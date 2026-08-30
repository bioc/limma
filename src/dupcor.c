/*
 * dupcor.c -- per-gene intra-block (duplicate) correlation by REML.
 *
 * C implementation of the per-gene work inside duplicateCorrelation(): it
 * reproduces statmod::mixedModel2Fit(y, X, Z, only.varcomp=TRUE) for each gene.
 * The consensus correlation (atanh transform + trimmed mean, and the rho
 * clamping) stays in R, in dups.R.
 *
 * Model for one gene:  y = X beta + Z u + e,  u ~ N(0, sb2 I), e ~ N(0, se2 I),
 * where X is the fixed-effect design and Z the block-indicator matrix. The
 * intra-block correlation reported per gene is
 *
 *     rho = sb2 / (sb2 + se2).
 *
 * REML estimation of (se2, sb2) without forming large matrices:
 *   1. QR of the (weighted) fixed-effect design X. Working in the residual
 *      space of dimension mq = nobs - rank removes the fixed effects (the
 *      "restricted" in REML): we obtain Q2^T y and A = Q2^T Z restricted to the
 *      last mq rows (dqrls returns these as the tail of its `effec` output).
 *   2. Diagonalise Cov(Q2^T y) = se2 I + sb2 A A^T without forming an mq x mq SVD.
 *      A thin QR A = Q R rotates the restricted response to g = Q^T (Q2^T y); in
 *      this basis the covariance se2 I + sb2 R R^T is block diagonal. The bottom
 *      mq-s coordinates of g (s = min(mq, nblocks)) are pure error, and the top
 *      s x nblocks block is diagonalised by a small SVD R = U_R D V^T (dgeqrf /
 *      dormqr + dgesvd on the s-row factor). The singular values d_i are those of
 *      A; only U_R (s x s) and d_i are needed, so V is skipped.
 *   3. In the rotated residual basis the squared components are independent with
 *           E[r_i^2] = se2 + sb2 * d_i^2,
 *      where r_top = U_R^T g_top (predictors d_i^2) and r_bottom = g_bottom (pure
 *      error, predictor 0). The mean is linear in (1, d_i^2). Because each r_i is
 *      ~ a scaled chi-square_1, Var(r_i^2) proportional to E[r_i^2]^2, so this is a gamma GLM
 *      with identity mean. Fitting it (OLS initialisation, then a damped gamma
 *      IRLS) gives  beta0 = se2,  beta1 = sb2,  and  rho = beta1 / (beta0 + beta1).
 *
 * Probe weights enter as sqrt(w) scaling of X and y (as in lm.c). Genes
 * with too few observations/blocks, a rank-deficient projection, or a failed
 * fit return NA_REAL (matching mixedModel2Fit's NA handling).
 */

#include <float.h>
#include <math.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "limma.h"

/*
 * Inputs:
 *   d, y: length-n vectors for the simple regression y_i ~ beta0 + beta1 d_i.
 * Outputs:
 *   beta0 and beta1 receive the intercept and slope; fitted receives length-n
 *   fitted means.
 * Returns:
 *   1 for a two-parameter fit; 0 when the predictor is effectively constant.
 * Notes:
 *   Used as the moment-based gamma-GLM start. On a constant predictor, beta1 is
 *   set to NA_REAL and fitted is the intercept-only mean.
 */
static int fitlm2(const double *d, const double *y, int n, double *beta0, double *beta1, double *fitted)
{
	int i;
	double sd = 0.0, sy = 0.0, sdd = 0.0, sdy = 0.0;
	double det;
	for(i = 0; i < n; i++)
	{
		sd += d[i];
		sy += y[i];
		sdd += d[i] * d[i];
		sdy += d[i] * y[i];
	}
	det = n * sdd - sd * sd;
	if(det <= DBL_EPSILON * fmax(1.0, n * sdd))
	{
		*beta0 = sy / n;
		*beta1 = NA_REAL;
		for(i = 0; i < n; i++)
			fitted[i] = *beta0;
		return 0;
	}
	*beta0 = (sy * sdd - sd * sdy) / det;
	*beta1 = (n * sdy - sd * sy) / det;
	for(i = 0; i < n; i++)
		fitted[i] = *beta0 + *beta1 * d[i];
	return 1;
}

/*
 * Inputs:
 *   y, mu: length-n non-negative responses and fitted means.
 * Returns:
 *   Gamma-family deviance D = 2 sum[(y-mu)/mu - log(y/mu)], or +Inf for
 *   non-finite or infeasible means.
 * Notes:
 *   Used only as the monotone objective for accepting damped gamma IRLS steps.
 */
static double gammadev(const double *y, const double *mu, int n)
{
	int i;
	double dev = 0.0;
	for(i = 0; i < n; i++)
	{
		if(mu[i] < 0.0 || !R_FINITE(mu[i]))
			return R_PosInf;
		if(y[i] < 1e-15 && mu[i] < 1e-15)
			continue;
		if(mu[i] <= 0.0)
			return R_PosInf;
		dev += (y[i] - mu[i]) / mu[i] - log(y[i] / mu[i]);
	}
	return 2.0 * dev;
}

/*
 * Inputs:
 *   d, y: length-n predictor d_i^2 and squared rotated residual components.
 *   beta0, beta1: starting values for se2 and sb2.
 * Outputs:
 *   beta0 and beta1 are updated in place; mu receives fitted gamma means.
 * Returns:
 *   1 on convergence or accepted early exit; 0 on numerical failure.
 * Notes:
 *   Fits identity-mean gamma IRLS with variance V(mu)=mu^2. A Levenberg
 *   damping parameter is increased until the gamma deviance does not increase.
 */
static int gammafit(const double *d, const double *y, int n, double *beta0, double *beta1, double *mu)
{
	int i, iter, lev;
	double dev, lambda = 0.0, maxinfo = 0.0, maxy = 0.0;

	for(i = 0; i < n; i++)
		if(y[i] > maxy)
			maxy = y[i];
	if(maxy == 0.0)
	{
		*beta0 = 0.0;
		*beta1 = 0.0;
		for(i = 0; i < n; i++)
			mu[i] = 0.0;
		return 1;
	}
	for(i = 0; i < n; i++)
		mu[i] = *beta0 + *beta1 * d[i];
	dev = gammadev(y, mu, n);

	for(iter = 1; iter <= 21; iter++)
	{
		double a = 0.0, b = 0.0, c = 0.0, score0 = 0.0, score1 = 0.0;
		double beta0_old = *beta0, beta1_old = *beta1, dev_old = dev;
		double dbeta0 = 0.0, dbeta1 = 0.0, maxmu = 0.0, maxv = 0.0;
		for(i = 0; i < n; i++)
		{
			double v = mu[i] * mu[i];
			if(mu[i] > maxmu)
				maxmu = mu[i];
			if(v > maxv)
				maxv = v;
		}
		/* For an identity-link gamma GLM, the derivative dmu/dbeta is the
		 * design vector (1,d_i). The Fisher scoring weight is therefore
		 * 1/V(mu_i) = 1/mu_i^2. The floor prevents near-zero fitted means from
		 * making a single component dominate the 2x2 system numerically. */
		double floor_v = maxv / 1000.0;
		for(i = 0; i < n; i++)
		{
			double v = mu[i] * mu[i];
			double invv, residual;
			if(v < floor_v)
				v = floor_v;
			if(v <= 0.0 || !R_FINITE(v))
				return 0;
			invv = 1.0 / v;
			residual = y[i] - mu[i];
			/* Accumulate the 2x2 Fisher information matrix and the score:
			 *   info = sum invv * [1,d_i]'[1,d_i]
			 *   score = sum (y_i-mu_i) * invv * [1,d_i].
			 * The damped solve below gives the proposed beta update. */
			a += invv;
			b += d[i] * invv;
			c += d[i] * d[i] * invv;
			score0 += residual * invv;
			score1 += d[i] * residual * invv;
		}
		maxinfo = fmax(a, c);
		if(iter == 1)
			lambda = fabs((a + c) / 2.0) / 2.0;

		/* Levenberg damping adds lambda to the information diagonal. If the
		 * proposed update increases the gamma deviance, lambda is doubled and the
		 * 2x2 system is solved again, shrinking the step toward steepest descent. */
		for(lev = 1; ; lev++)
		{
			double aa = a + lambda;
			double cc = c + lambda;
			double det = aa * cc - b * b;
			if(det <= 0.0 || !R_FINITE(det))
				return 0;
			dbeta0 = (cc * score0 - b * score1) / det;
			dbeta1 = (aa * score1 - b * score0) / det;
			*beta0 = beta0_old + dbeta0;
			*beta1 = beta1_old + dbeta1;
			for(i = 0; i < n; i++)
				mu[i] = *beta0 + *beta1 * d[i];
			dev = gammadev(y, mu, n);
			if(dev <= dev_old || dev / maxmu < 1e-15)
				break;
			if(lambda / maxinfo > 1e15)
			{
				*beta0 = beta0_old;
				*beta1 = beta1_old;
				return 1;
			}
			lambda *= 2.0;
		}
		if(lambda / maxinfo > 1e15)
			break;
		if(lev == 1)
			lambda /= 10.0;
		if(score0 * dbeta0 + score1 * dbeta1 < 1e-6 || dev / maxmu < 1e-15)
			break;
	}
	return R_FINITE(*beta0) && R_FINITE(*beta1);
}

/*
 * Inputs:
 *   m: ngenes x narrays expression matrix.
 *   design: narrays x nbeta fixed-effect design matrix.
 *   block: length-narrays one-based block labels.
 *   weights: optional ngenes x narrays prior-weight matrix, or NULL.
 *   gene: zero-based row of m; dws is reusable scratch for one thread.
 * Outputs:
 *   No caller-visible arrays are filled directly; all intermediate storage is
 *   in dws.
 * Returns:
 *   Per-gene REML correlation rho, or NA_REAL when the gene is not estimable or
 *   a numerical step fails.
 * Notes:
 *   Finite observations are compacted and weighted by sqrt(weight). Fixed
 *   effects are removed by dqrls, then Q2'Z is diagonalised by QR plus a small
 *   SVD before the gamma variance-component fit.
 */
double dupcorgene(const double *m, const double *design, const int *block, const double *weights, int ngenes, int narrays, int nbeta, int nblock_levels, int gene, dupws *dws)
{
	int array, i, j, nobs = 0, nblocks = 0, rank = 0, ny, mq, s, info, one = 1;
	int nonzero_d = 0;
	int lwork = dws->lwork;
	int ldvt = 1;
	double beta0, beta1, sum_s = 0.0, sum_ss = 0.0;
	double tol = 1e-7;
	const char jobu = 'A', jobvt = 'N', side = 'L', trans = 'T';

	const double *mptr = m + gene;
	const double *wptr = weights ? weights + gene : NULL;

	/* Scan observed arrays: relabel the blocks actually present to 0..nblocks-1
	 * (bmap), cache the response (mval) and prewhitening factor scale = sqrt(w). */
	for(i = 0; i < nblock_levels; i++)
		dws->bmap[i] = -1;
	for(array = 0; array < narrays; array++)
	{
		double value = *mptr;
		if(R_FINITE(value))
		{
			int original_block = block[array] - 1;
			if(dws->bmap[original_block] < 0)
				dws->bmap[original_block] = nblocks++;
			dws->mval[nobs] = value;
			dws->scale[nobs] = weights ? sqrt(*wptr) : 1.0;
			dws->obser[nobs++] = array;
		}
		mptr += ngenes;
		if(weights)
			wptr += ngenes;
	}
	/* need enough residual df and >1 block, and not a saturated block design */
	if(nobs <= nbeta + 2 || nblocks <= 1 || nblocks >= nobs - 1)
		return NA_REAL;

	/* weighted fixed-effect design X~ = sqrt(w) .* design[obs, ] */
	for(j = 0; j < nbeta; j++)
	{
		for(i = 0; i < nobs; i++)
		{
			array = dws->obser[i];
			dws->x[i + j * nobs] = design[array + j * narrays];
			if(!R_FINITE(dws->x[i + j * nobs]))
				return NA_REAL;
			dws->x[i + j * nobs] *= dws->scale[i];
		}
		dws->pivot[j] = j + 1;
	}
	/* rhs holds [ Z | y ]: nblocks block-indicator columns then the weighted
	 * response, all to be rotated by Q^T in one dqrls call (ny = nblocks+1) */
	ny = nblocks + 1;
	for(j = 0; j < nblocks; j++)
	{
		for(i = 0; i < nobs; i++)
		{
			array = dws->obser[i];
			dws->rhs[i + j * nobs] =
				dws->bmap[block[array] - 1] == j ? 1.0 : 0.0;
		}
	}
	for(i = 0; i < nobs; i++)
	{
		dws->rhs[i + nblocks * nobs] = dws->mval[i] * dws->scale[i];
	}

	/* dqrls (R LINPACK-based) solves real least-squares problems by
	 * rank-revealing QR decomposition with column pivoting. limma uses it here
	 * to remove fixed effects from the block indicators and response:
	 *   - dws->x (nobs x nbeta, lda=nobs) enters as X~ and is overwritten by
	 *     compact pivoted-QR storage;
	 *   - dws->rhs has ny=nblocks+1 right-hand sides, [Z | y], so the same Q^T
	 *     rotation is applied to all block-indicator columns and the response;
	 *   - coeff/resid/effec receive coefficient, residual and Q^T[Z|y] outputs;
	 *   - rank and the 1-based pivot vector describe the fixed-effect rank;
	 *   - qraux/qwork hold Householder and temporary QR data.
	 * The rows of effec below rank are the REML residual space: Q2^T Z and
	 * Q2^T y.
	 * Equivalent R operation: lm.fit(Xtilde, cbind(Z, y)) followed by extracting
	 * the residual-space effects below the fitted rank.
	 * Netlib references:
	 *   https://www.netlib.org/linpack/dqrdc.f
	 *   https://www.netlib.org/linpack/dqrsl.f */
	F77_CALL(dqrls)(dws->x, &nobs, &nbeta, dws->rhs, &ny, &tol, dws->coeff, dws->resid, dws->effec, &rank, dws->pivot, dws->qraux, dws->qwork);
	mq = nobs - rank;
	if(mq == 0)
		return NA_REAL;

	/* qtz = Q2^T Z (mq x nblocks): the block columns of effec below row `rank` */
	for(j = 0; j < nblocks; j++)
		for(i = 0; i < mq; i++)
			dws->qtz[i + j * mq] = dws->effec[rank + i + j * nobs];

	/* Diagonalise Cov(Q2^T y) = se2 I + sb2 (qtz)(qtz)^T without an mq x mq SVD.
	 * Thin-QR qtz = Q R, then rotate the restricted response g = Q^T (Q2^T y). In
	 * the Q basis Cov(g) = se2 I + sb2 R R^T is block diagonal: the bottom mq-s
	 * coordinates g[s..] are pure error (predictor 0), and the top s x nblocks
	 * block (s = min(mq, nblocks)) is diagonalised by the small SVD R = U_R D V^T.
	 * The singular values d_i are those of qtz; only U_R (s x s) is needed. */
	s = mq < nblocks ? mq : nblocks;

	/* Copy b = Q2^T y, the restricted response, before factorizing qtz. dgeqrf
	 * computes qtz = Q R using Householder reflectors, overwriting qtz with R in
	 * its upper triangle and reflectors below the diagonal while tau stores the
	 * reflector scalars. dormqr then applies Q^T to b in gvec without explicitly
	 * forming Q, giving g = Q^T b. Only the first s reflectors are needed because
	 * the effective rank of the thin QR block is at most s=min(mq,nblocks). */
	for(i = 0; i < mq; i++)
		dws->gvec[i] = dws->effec[rank + i + nblocks * nobs];

	/* DGEQRF (Netlib LAPACK) computes a QR factorization of a real M-by-N matrix
	 * A. limma uses it here to factor Q2^T Z before the smaller SVD:
	 *   - dws->qtz (mq x nblocks, lda=mq) enters as the residual-space block
	 *     design Q2^T Z;
	 *   - dws->qtz is overwritten: its upper trapezoid contains R, while entries
	 *     below the diagonal plus dws->tau store the Householder reflectors for Q;
	 *   - dws->swork/lwork provide LAPACK workspace;
	 *   - info<0 would indicate an illegal argument.
	 * The compact QR storage supplies R and the reflectors used to rotate Q2^T y.
	 * Equivalent R operation: qr(Q2tZ), retaining the compact QR representation.
	 * Netlib references:
	 *   https://www.netlib.org/lapack/explore-html/d0/da1/group__geqrf_gade26961283814bb4e62183d9133d8bf5.html */
	F77_CALL(dgeqrf)(&mq, &nblocks, dws->qtz, &mq, dws->tau, dws->swork, &lwork, &info);
	if(info != 0)
		return NA_REAL;
	/* DORMQR (Netlib LAPACK) multiplies a real matrix C by the orthogonal matrix
	 * Q from a QR factorization. limma uses it here to rotate the restricted
	 * response without forming Q explicitly:
	 *   - SIDE='L' applies Q from the left;
	 *   - TRANS='T' applies Q^T;
	 *   - A=dws->qtz and TAU=dws->tau hold the DGEQRF reflectors;
	 *   - C=dws->gvec is an mq x 1 matrix with ldc=mq;
	 *   - K=s uses the first s reflectors for the thin-Q action;
	 *   - dws->gvec is overwritten by Q^T (Q2^T y).
	 * The rotated response aligns with the triangular block whose SVD follows.
	 * Equivalent R operation: crossprod(qr.Q(qr(Q2tZ)), Q2ty).
	 * Netlib references:
	 *   https://www.netlib.org/lapack/explore-html/d7/d50/group__unmqr_ga768bd221f959be1b3d15bd177bb5c1b3.html */
	F77_CALL(dormqr)(&side, &trans, &mq, &one, &s, dws->qtz, &mq, dws->tau, dws->gvec, &mq, dws->swork, &lwork, &info FCONE FCONE);
	if(info != 0)
		return NA_REAL;

	/* rmat receives the top s rows of R from the QR factorization. Entries below
	 * the diagonal in qtz hold Householder reflector data, not R, so they are
	 * explicitly zeroed while copying.
	 *
	 * DGESVD (Netlib LAPACK) computes the singular value decomposition of a real
	 * M-by-N matrix A. limma uses it here to diagonalise the small triangular
	 * block from Q2^T Z:
	 *   - A=dws->rmat is s x nblocks with LDA=s;
	 *   - JOBU='A' stores all s left singular vectors in dws->u;
	 *   - JOBVT='N' skips right singular vectors, so VT is passed as NULL with a
	 *     dummy LDVT;
	 *   - dws->svals receives singular values in descending order;
	 *   - dws->rmat is overwritten by LAPACK work data.
	 * Only U and the singular values are needed because the REML gamma fit uses
	 * U^T g and d_i^2, not V.
	 * Equivalent R operation: svd(rmat, nu=s, nv=0).
	 * Netlib references:
	 *   https://www.netlib.org/lapack/explore-html/d1/d7f/group__gesvd_gac6bd5d4e645049e49bb70691180abf07.html */
	for(j = 0; j < nblocks; j++)
		for(i = 0; i < s; i++)
			dws->rmat[i + j * s] = (i <= j) ? dws->qtz[i + j * mq] : 0.0;
	F77_CALL(dgesvd)(&jobu, &jobvt, &s, &nblocks, dws->rmat, &s, dws->svals, dws->u, &s, NULL, &ldvt, dws->swork, &lwork, &info FCONE FCONE);
	if(info != 0)
		return NA_REAL;

	/* Rotated squared components: predictor d_i^2, response r_i^2 with
	 *   E[r_i^2] = se2 + sb2 * d_i^2.  Top block r = U_R^T g[0..s-1]; bottom block
	 * (pure error, predictor 0) is g[s..mq-1] directly. Store predictor in
	 * rhs[0..mq-1], response in rhs[mq..2mq-1]. */
	for(i = 0; i < s; i++)
	{
		double uqy = 0.0;
		for(j = 0; j < s; j++)
			uqy += dws->u[j + i * s] * dws->gvec[j];
		dws->rhs[i] = dws->svals[i] * dws->svals[i];
		dws->rhs[mq + i] = uqy * uqy;
	}
	for(i = s; i < mq; i++)
	{
		dws->rhs[i] = 0.0;
		dws->rhs[mq + i] = dws->gvec[i] * dws->gvec[i];
	}
	for(i = 0; i < mq; i++)
	{
		if(fabs(dws->rhs[i]) > 1e-15)
			nonzero_d++;
		sum_s += dws->rhs[i];
		sum_ss += dws->rhs[i] * dws->rhs[i];
	}

	/* moment-based start: OLS of the squared residuals on (1, d_i^2) */
	if(!fitlm2(dws->rhs, dws->rhs + mq, mq, &beta0, &beta1, dws->effec))
		return NA_REAL;

	/* refine by gamma GLM only when there is genuine spread in the predictor and
	 * the OLS fitted values are usable; otherwise keep the moment estimate */
	if(mq > 2 && nonzero_d > 1 && (sum_ss - sum_s * sum_s / mq) / (mq - 1) > 1e-15)
	{
		int nonnegative = 1;
		for(i = 0; i < mq; i++)
			if(dws->effec[i] < 0.0)
				nonnegative = 0;
		if(!nonnegative)
		{            /* infeasible OLS fit -> restart from mean-only */
			beta0 = 0.0;
			for(i = 0; i < mq; i++)
				beta0 += dws->rhs[mq + i];
			beta0 /= mq;
			beta1 = 0.0;
		}
		if(!gammafit(dws->rhs, dws->rhs + mq, mq, &beta0, &beta1, dws->effec))
			return NA_REAL;
	}

	/* rho = sb2 / (se2 + sb2) = beta1 / (beta0 + beta1) */
	if(ISNA(beta0) || ISNA(beta1))
		return NA_REAL;
	return beta1 / (beta0 + beta1);
}

static void allocdup(dupws *dws, int narrays, int nbeta, int nblock_levels)
{
	int min_dim = narrays < nblock_levels ? narrays : nblock_levels;
	int max_rhs = nblock_levels + 1;
	dws->x = R_Calloc((size_t) narrays * nbeta, double);
	dws->rhs = R_Calloc((size_t) narrays * max_rhs, double);
	dws->coeff = R_Calloc((size_t) nbeta * max_rhs, double);
	dws->resid = R_Calloc((size_t) narrays * max_rhs, double);
	dws->effec = R_Calloc((size_t) narrays * max_rhs, double);
	dws->qraux = R_Calloc(nbeta, double);
	dws->qwork = R_Calloc(2 * nbeta, double);
	dws->qtz = R_Calloc((size_t) narrays * nblock_levels, double);
	dws->tau = R_Calloc(min_dim, double);
	dws->rmat = R_Calloc((size_t) nblock_levels * nblock_levels, double);
	dws->svals = R_Calloc(min_dim, double);
	dws->u = R_Calloc((size_t) nblock_levels * nblock_levels, double);
	dws->gvec = R_Calloc(narrays, double);
	dws->mval = R_Calloc(narrays, double);
	dws->scale = R_Calloc(narrays, double);
	dws->obser = R_Calloc(narrays, int);
	dws->bmap = R_Calloc(nblock_levels, int);
	dws->pivot = R_Calloc(nbeta, int);
	dws->swork = NULL;
}

static void freedup(dupws *dws)
{
	R_Free(dws->x);
	R_Free(dws->rhs);
	R_Free(dws->coeff);
	R_Free(dws->resid);
	R_Free(dws->effec);
	R_Free(dws->qraux);
	R_Free(dws->qwork);
	R_Free(dws->qtz);
	R_Free(dws->tau);
	R_Free(dws->rmat);
	R_Free(dws->svals);
	R_Free(dws->u);
	R_Free(dws->gvec);
	R_Free(dws->swork);
	R_Free(dws->mval);
	R_Free(dws->scale);
	R_Free(dws->obser);
	R_Free(dws->bmap);
	R_Free(dws->pivot);
}

/*
 * Inputs:
 *   m: ngenes x narrays expression matrix.
 *   design: narrays x nbeta fixed-effect design matrix.
 *   block: length-narrays one-based block labels with nblock_levels levels.
 *   weights: optional ngenes x narrays prior-weight matrix, or NULL.
 *   nthreads: requested OpenMP thread count.
 * Outputs:
 *   rho: length-ngenes vector of per-gene correlations.
 * Returns:
 *   0 on success; 1 if LAPACK workspace query for the small SVD fails.
 * Notes:
 *   One dupws workspace is allocated per active thread. The consensus
 *   correlation and atanh trimming are handled in R, not here.
 */
int dupcor(const double *m, const double *design, const int *block, const double *weights, int ngenes, int narrays, int nbeta, int nblock_levels, int nthreads, double *rho)
{
	int gene, info, querylw = -1, lwork, one = 1, ldvt = 1, smax;
	double wq, wqmax;
	dupws *dws;
	const char jobu = 'A', jobvt = 'N', side = 'L', trans = 'T';

	nthreads = clampthreads(nthreads, ngenes);

	dws = R_Calloc(nthreads, dupws);
	for(gene = 0; gene < nthreads; gene++)
		allocdup(dws + gene, narrays, nbeta, nblock_levels);

	/* DGEQRF (Netlib LAPACK) computes a QR factorization, DORMQR (Netlib LAPACK)
	 * multiplies by a QR orthogonal factor, and DGESVD (Netlib LAPACK) computes
	 * a singular value decomposition. limma uses these calls only to query
	 * LAPACK workspace for the per-gene QR/SVD chain:
	 *   - LWORK=-1 requests query mode, so no factorization, multiplication or
	 *     SVD is performed;
	 *   - dws[0].qtz, dws[0].tau, dws[0].gvec, dws[0].rmat, dws[0].svals and
	 *     dws[0].u are dummy buffers with maximum narrays/nblock_levels sizes;
	 *   - WORK(1), here scalar wq, receives each routine's recommended workspace;
	 *   - info would report illegal dimensions if negative.
	 * The maximum queried wq becomes a safe swork length for every gene-specific
	 * mq and nblocks.
	 * Equivalent R operation: none; this asks LAPACK how much workspace to allocate.
	 * Netlib references:
	 *   https://www.netlib.org/lapack/explore-html/d0/da1/group__geqrf_gade26961283814bb4e62183d9133d8bf5.html
	 *   https://www.netlib.org/lapack/explore-html/d7/d50/group__unmqr_ga768bd221f959be1b3d15bd177bb5c1b3.html
	 *   https://www.netlib.org/lapack/explore-html/d1/d7f/group__gesvd_gac6bd5d4e645049e49bb70691180abf07.html */
	smax = narrays < nblock_levels ? narrays : nblock_levels;
	F77_CALL(dgeqrf)(&narrays, &nblock_levels, dws[0].qtz, &narrays, dws[0].tau, &wq, &querylw, &info);
	wqmax = wq;
	F77_CALL(dormqr)(&side, &trans, &narrays, &one, &smax, dws[0].qtz, &narrays, dws[0].tau, dws[0].gvec, &narrays, &wq, &querylw, &info FCONE FCONE);
	if(wq > wqmax)
		wqmax = wq;
	F77_CALL(dgesvd)(&jobu, &jobvt, &smax, &nblock_levels, dws[0].rmat, &smax, dws[0].svals, dws[0].u, &smax, NULL, &ldvt, &wq, &querylw, &info FCONE FCONE);
	if(info != 0)
	{
		for(gene = 0; gene < nthreads; gene++)
			freedup(dws + gene);
		R_Free(dws);
		return 1;
	}
	if(wq > wqmax)
		wqmax = wq;
	lwork = (int) wqmax;
	if(lwork < 1)
		lwork = 1;
	for(gene = 0; gene < nthreads; gene++)
	{
		dws[gene].lwork = lwork;
		dws[gene].swork = R_Calloc(lwork, double);
	}

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static) num_threads(nthreads)
	#endif
	for(gene = 0; gene < ngenes; gene++)
	{
		int thread = 0;
		#ifdef _OPENMP
		thread = omp_get_thread_num();
		#endif
		rho[gene] = dupcorgene(m, design, block, weights, ngenes, narrays, nbeta, nblock_levels, gene, dws + thread);
	}
	for(gene = 0; gene < nthreads; gene++)
		freedup(dws + gene);
	R_Free(dws);
	return 0;
}
