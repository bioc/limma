/*
 * awreml.c -- REML array-weight estimation with prior observation weights.
 *
 * C implementation of the Fisher-scoring repeat{} loop of .arrayWeightsPrWtsREML
 * (R/arrayWeightsPrWtsREML.R). Because each gene carries its own observation
 * weights w*weights[g,], every gene needs a separate weighted QR, so the inner
 * per-gene loop is the bottleneck; it is parallelised over genes with OpenMP.
 *
 * Model: array variances are modelled as log-linear in the variance design Z2
 * (columns sum to zero), w = exp(-Z2 gam). Each outer iteration accumulates, over
 * genes, the REML Fisher information info2 (ngam x ngam) and score-related vector
 * z (narrays), then takes a Fisher-scoring step gam += solve(info2, Z2' z).
 *
 * Per gene g, with Wg = w*weights[g,], sw = sqrt(Wg), Xw = sw.*design, yw = sw.*y[g,]:
 *   - QR of Xw gives Q1 (first p columns of Q); residual rw = (I - Q1 Q1') yw,
 *     s2 = ||rw||^2 / (narrays - p). The R score term w*weights*resid^2 equals
 *     rw^2 elementwise, so no original-scale residual is needed.
 *   - leverages h = diag(Q1 Q1'); Q2 = pairwise column products of Q1 (with the
 *     off-diagonal blocks scaled by sqrt(2)); with Z = [1 | Z2],
 *         info = Z' diag(1-2h) Z + (Q2' Z)'(Q2' Z),
 *     and info2_g = Schur complement of info removing the intercept.
 *   - info2 += info2_g; if s2 > 1e-15: z += rw^2/s2 - (1-h).
 * h, rw, s2 and info depend only on the projection P = Q1 Q1' (Q2 Q2' = P.*P), so
 * any stable QR reproduces the lm.wfit/qr.qy result of the R reference to rounding.
 *
 * Uses LAPACK dgeqrf/dorgqr (QR + form Q, no char args) and dgesv for the
 * ngam x ngam Fisher-scoring solve (matching R's solve()). All small dense
 * products are explicit loops.
 */

#include <math.h>
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "limma.h"

static void allocaw(awws *aws, int narrays, int p, int p2, int ngam, int lwork)
{
	int ng1 = ngam + 1;
	aws->xw = R_Calloc((size_t) narrays * p, double);
	aws->yw = R_Calloc(narrays, double);
	aws->rw = R_Calloc(narrays, double);
	aws->tau = R_Calloc(p, double);
	aws->cvec = R_Calloc(p, double);
	aws->q2 = R_Calloc((size_t) narrays * p2, double);
	aws->bmat = R_Calloc((size_t) p2 * ng1, double);
	aws->hvec = R_Calloc(narrays, double);
	aws->info = R_Calloc((size_t) ng1 * ng1, double);
	aws->pinfo2 = R_Calloc((size_t) ngam * ngam, double);
	aws->pz = R_Calloc(narrays, double);
	aws->lwork = lwork;
	aws->work = R_Calloc(lwork, double);
}

static void freeaw(awws *aws)
{
	R_Free(aws->xw);
	R_Free(aws->yw);
	R_Free(aws->rw);
	R_Free(aws->tau);
	R_Free(aws->cvec);
	R_Free(aws->q2);
	R_Free(aws->bmat);
	R_Free(aws->hvec);
	R_Free(aws->info);
	R_Free(aws->pinfo2);
	R_Free(aws->pz);
	R_Free(aws->work);
}

/* Accumulate gene g's contributions into aws->pinfo2 (ngam^2) and aws->pz (narrays). */
/*
 * Inputs:
 *   y: ngenes x narrays expression matrix.
 *   design: narrays x p mean-model design matrix.
 *   weights: ngenes x narrays prior observation weights.
 *   Z: narrays x (ngam+1) variance design with intercept in the first column.
 *   w: current length-narrays array weights; gene is the zero-based row.
 * Outputs:
 *   aws->pinfo2 accumulates this gene's Fisher information contribution.
 *   aws->pz accumulates this gene's score-related residual contribution.
 * Notes:
 *   All scratch is in aws and is private to one thread. The weighted QR forms
 *   Q1, leverages, residuals and Q1 column products for the REML update.
 */
void awremlgene(const double *y, const double *design, const double *weights, const double *Z, int ngenes, int narrays, int p, int p2, int ngam, const double *w, int gene, awws *aws)
{
	int i, j, k, a, c, r, col, info, ng1 = ngam + 1;
	int lwork = aws->lwork;
	double *q1 = aws->xw;          /* xw becomes Q1 after dorgqr */
	double rss, s2;

	/* weighted design Xw = sqrt(Wg).*design and response yw = sqrt(Wg).*y[g,] */
	const double *mptr = y + gene;
	const double *wptr = weights + gene;
	for(i = 0; i < narrays; i++)
	{
		double sw  = sqrt(w[i] * wptr[(size_t) i * ngenes]);
		aws->yw[i] = sw * mptr[(size_t) i * ngenes];
		for(j = 0; j < p; j++)
			aws->xw[i + (size_t) j * narrays] = sw * design[i + (size_t) j * narrays];
	}

	/* DGEQRF (Netlib LAPACK) computes a QR factorization of a real M-by-N matrix
	 * A, and DORGQR (Netlib LAPACK) generates the real orthonormal matrix Q from
	 * the Householder reflectors returned by DGEQRF. limma uses this pair here
	 * to build the fixed-effect basis for REML residual calculations:
	 *   - aws->xw (narrays x p, lda=narrays) enters as the weighted design
	 *     Xw = sqrt(Wg) * X and is overwritten by compact QR storage;
	 *   - after DGEQRF, the upper triangle of aws->xw contains R, while the lower
	 *     triangle plus aws->tau stores the Householder reflectors defining Q;
	 *   - aws->work/lwork provide LAPACK workspace, and info would report an
	 *     illegal argument if negative;
	 *   - DORGQR reads the reflectors from aws->xw/aws->tau and overwrites
	 *     aws->xw with the first p explicit columns of Q.
	 * After DORGQR, q1 points to Q1, the fitted fixed-effect space used below
	 * for Q1^T yw, Q1 Q1^T, REML residuals, leverages and score terms.
	 * Equivalent R operation: qr.Q(qr(Xw)).
	 * Netlib references:
	 *   https://www.netlib.org/lapack/explore-html/d0/da1/group__geqrf_gade26961283814bb4e62183d9133d8bf5.html
	 *   https://www.netlib.org/lapack/explore-html/d4/dfc/group__ungqr_ga9f4abcfe1543a5d6d90fcd4dd21f12f0.html */
	F77_CALL(dgeqrf)(&narrays, &p, aws->xw, &narrays, aws->tau, aws->work, &lwork, &info);
	F77_CALL(dorgqr)(&narrays, &p, &p, aws->xw, &narrays, aws->tau, aws->work, &lwork, &info);

	/* Project the weighted response onto the fitted mean-model space. cvec is
	 * Q1^T yw, the fitted vector is Q1 cvec, and rw is the REML residual in the
	 * orthogonal complement. The leverage h_i is row i of Q1 times itself, i.e.
	 * diag(Q1 Q1^T), and rss is ||rw||^2. */
	for(j = 0; j < p; j++)
	{
		double s = 0.0;
		const double *qj = q1 + (size_t) j * narrays;
		for(i = 0; i < narrays; i++)
			s += qj[i] * aws->yw[i];
		aws->cvec[j] = s;
	}
	rss = 0.0;
	for(i = 0; i < narrays; i++)
	{
		double fit = 0.0, hh = 0.0;
		for(j = 0; j < p; j++)
		{
			double qij = q1[i + (size_t) j * narrays];
			fit += qij * aws->cvec[j];
			hh  += qij * qij;
		}
		aws->rw[i] = aws->yw[i] - fit;
		rss += aws->rw[i] * aws->rw[i];
		aws->hvec[i] = hh;
	}
	s2 = rss / (narrays - p);

	/* Q2 stores the unique column-wise products of Q1. Diagonal products are
	 * q_a*q_a. Off-diagonal products are multiplied by sqrt(2) below so that
	 * inner products of Q2 columns reproduce the symmetric quadratic terms in
	 * the REML Fisher information without storing duplicate a,b and b,a pairs. */
	col = 0;
	for(k = 0; k < p; k++)
	{
		for(a = 0; a + k < p; a++)
		{
			double *q2c = aws->q2 + (size_t) col * narrays;
			const double *qa = q1 + (size_t) a * narrays;
			const double *qb = q1 + (size_t)(a + k) * narrays;
			for(i = 0; i < narrays; i++)
				q2c[i] = qa[i] * qb[i];
			col++;
		}
	}
	double sqrt2 = sqrt(2.0);
	for(col = p; col < p2; col++)
	{
		double *q2c = aws->q2 + (size_t) col * narrays;
		for(i = 0; i < narrays; i++)
			q2c[i] *= sqrt2;
	}

	/* bmat = Q2^T Z  (p2 x ng1), Z = [1 | Z2] */
	for(c = 0; c < ng1; c++)
	{
		const double *zc = Z + (size_t) c * narrays;
		for(r = 0; r < p2; r++)
		{
			const double *q2r = aws->q2 + (size_t) r * narrays;
			double s = 0.0;
			for(i = 0; i < narrays; i++)
				s += q2r[i] * zc[i];
			aws->bmat[r + (size_t) c * p2] = s;
		}
	}

	/* info = Z^T diag(1-2h) Z + bmat^T bmat   (ng1 x ng1) */
	for(c = 0; c < ng1; c++)
	{
		const double *zc = Z + (size_t) c * narrays;
		for(r = 0; r < ng1; r++)
		{
			const double *zr = Z + (size_t) r * narrays;
			double t1 = 0.0, t2 = 0.0;
			for(i = 0; i < narrays; i++)
				t1 += zr[i] * (1.0 - 2.0 * aws->hvec[i]) * zc[i];
			for(k = 0; k < p2; k++)
				t2 += aws->bmat[k + (size_t) r * p2] * aws->bmat[k + (size_t) c * p2];
			aws->info[r + (size_t) c * ng1] = t1 + t2;
		}
	}

	/* The first column of Z is the intercept for the variance model. The REML
	 * scoring step for the free variance parameters uses the Schur complement of
	 * the full information matrix after eliminating that intercept. Accumulate
	 * only the ngam x ngam block needed by the outer Fisher-scoring solve. */
	double i00 = aws->info[0];
	for(c = 0; c < ngam; c++)
		for(r = 0; r < ngam; r++)
		{
			double v = aws->info[(r + 1) + (size_t)(c + 1) * ng1] - aws->info[(r + 1)] * aws->info[(size_t)(c + 1) * ng1] / i00;
			aws->pinfo2[r + (size_t) c * ngam] += v;
		}

	/* The score vector contribution is rw_i^2/s2 - (1-h_i). In the original R
	 * expression this appears as w*weights*resid^2/s2, but rw is already the
	 * weighted residual, so no original-scale residual is needed here. */
	if(s2 > 1e-15)
		for(i = 0; i < narrays; i++)
			aws->pz[i] += aws->rw[i] * aws->rw[i] / s2 - (1.0 - aws->hvec[i]);
}

/*
 * Inputs:
 *   y: ngenes x narrays expression matrix.
 *   design: narrays x p mean-model design matrix.
 *   weights: ngenes x narrays prior observation weights.
 *   Z2: narrays x ngam variance design without intercept.
 *   prior_n, maxiter, tol, trace: REML prior and Fisher-scoring controls.
 *   nthreads: requested OpenMP thread count.
 * Outputs:
 *   w: length-narrays estimated array weights, filled in place.
 *   iter_out: number of Fisher-scoring iterations used.
 * Returns:
 *   0 for convergence, 1 for non-finite/failed solve, 2 for maxiter reached.
 * Notes:
 *   Per-gene contributions are accumulated in thread-local awws workspaces and
 *   reduced after the OpenMP loop before each Fisher-scoring solve.
 */
int awreml(const double *y, const double *design, const double *weights, const double *Z2, int ngenes, int narrays, int p, int ngam, double prior_n, int maxiter, double tol, int trace, int nthreads, double *w, int *iter_out)
{
	int i, c, r, t, gene, iter = 0, status = 0, info, one = 1, ng1 = ngam + 1;
	int p2 = p * (p + 1) / 2;
	int lwork = 1, querylw = -1;
	double wq, denom = ngenes + prior_n;
	double *Z, *Z2tZ2, *info2, *z, *dl, *Asolve, *sol, *gam;
	int *ipiv;
	awws *aws;

	nthreads = clampthreads(nthreads, ngenes);

	/* Z = [1 | Z2]  (narrays x ng1); Z2tZ2 = Z2^T Z2 (ngam x ngam, constant) */
	Z = R_Calloc((size_t) narrays * ng1, double);
	for(i = 0; i < narrays; i++)
		Z[i] = 1.0;
	for(c = 0; c < ngam; c++)
		for(i = 0; i < narrays; i++)
			Z[i + (size_t)(c + 1) * narrays] = Z2[i + (size_t) c * narrays];
	Z2tZ2 = R_Calloc((size_t) ngam * ngam, double);
	for(c = 0; c < ngam; c++)
		for(r = 0; r < ngam; r++)
		{
			double s = 0.0;
			for(i = 0; i < narrays; i++)
				s += Z2[i + (size_t) r * narrays] * Z2[i + (size_t) c * narrays];
			Z2tZ2[r + (size_t) c * ngam] = s;
		}

	info2 = R_Calloc((size_t) ngam * ngam, double);
	z = R_Calloc(narrays, double);
	dl = R_Calloc(ngam, double);
	Asolve = R_Calloc((size_t) ngam * ngam, double);
	sol = R_Calloc(ngam, double);
	gam = R_Calloc(ngam, double);
	ipiv = R_Calloc(ngam, int);

	for(i = 0; i < narrays; i++)
		w[i] = 1.0;

	/* DGEQRF (Netlib LAPACK) computes a QR factorization of a real M-by-N matrix
	 * A, and DORGQR (Netlib LAPACK) generates the real orthonormal matrix Q from
	 * DGEQRF reflectors. limma uses these calls only to query LAPACK workspace:
	 *   - LWORK=-1 requests query mode, so no QR factorization or Q construction
	 *     is performed;
	 *   - aq/tq are dummy A/TAU buffers with the maximum narrays x p dimensions
	 *     used by every gene;
	 *   - WORK(1), here scalar wq, receives each routine's recommended workspace;
	 *   - info would report illegal dimensions if negative.
	 * The larger DGEQRF/DORGQR query result becomes the per-thread aws->work length.
	 * Equivalent R operation: none; this asks LAPACK how much workspace to allocate.
	 * Netlib references:
	 *   https://www.netlib.org/lapack/explore-html/d0/da1/group__geqrf_gade26961283814bb4e62183d9133d8bf5.html
	 *   https://www.netlib.org/lapack/explore-html/d4/dfc/group__ungqr_ga9f4abcfe1543a5d6d90fcd4dd21f12f0.html */
	double *aq = R_Calloc((size_t) narrays * p, double);
	double *tq = R_Calloc(p, double);
	F77_CALL(dgeqrf)(&narrays, &p, aq, &narrays, tq, &wq, &querylw, &info);
	lwork = (int) wq;
	F77_CALL(dorgqr)(&narrays, &p, &p, aq, &narrays, tq, &wq, &querylw, &info);
	if((int) wq > lwork)
		lwork = (int) wq;
	if(lwork < 1)
		lwork = 1;
	R_Free(aq);
	R_Free(tq);

	aws = R_Calloc(nthreads, awws);
	for(t = 0; t < nthreads; t++)
		allocaw(aws + t, narrays, p, p2, ngam, lwork);

	for(iter = 1; iter <= maxiter; iter++)
	{
		double convcrit, num = 0.0;
		int solok = 1;

		/* priors: info2 starts at prior.n*Z2'Z2, z at prior.n*(w-1) */
		for(i = 0; i < ngam * ngam; i++)
			info2[i] = prior_n * Z2tZ2[i];
		for(i = 0; i < narrays; i++)
			z[i] = prior_n * (w[i] - 1.0);
		for(t = 0; t < nthreads; t++)
		{
			for(i = 0; i < ngam * ngam; i++)
				aws[t].pinfo2[i] = 0.0;
			for(i = 0; i < narrays; i++)
				aws[t].pz[i] = 0.0;
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
			awremlgene(y, design, weights, Z, ngenes, narrays, p, p2, ngam, w, gene, aws + thread);
		}

		/* reduce per-thread accumulators */
		for(t = 0; t < nthreads; t++)
		{
			for(i = 0; i < ngam * ngam; i++)
				info2[i] += aws[t].pinfo2[i];
			for(i = 0; i < narrays; i++)
				z[i] += aws[t].pz[i];
		}
		for(i = 0; i < ngam * ngam; i++)
			info2[i] /= denom;
		for(i = 0; i < narrays; i++)
			z[i] /= denom;

		/* DGESV (Netlib LAPACK) solves a real linear system A X = B by LU
		 * factorization with partial pivoting. limma uses it here for the dense
		 * Fisher-scoring update:
		 *   - A is Asolve, a copy of the ngam x ngam information matrix info2;
		 *   - NRHS=1 and B is sol, initially dl = Z2^T z;
		 *   - Asolve is overwritten by L and U factors of P*A;
		 *   - ipiv receives the row pivots defining P;
		 *   - sol is overwritten by the solution gamstep when info=0;
		 *   - info=i>0 means U(i,i) is exactly zero, so the step is singular.
		 * The solution is added to the variance-model coefficients gam.
		 * Equivalent R operation: solve(info2, dl).
		 * Netlib references:
		 *   https://www.netlib.org/lapack/explore-html/d8/da6/group__gesv_ga831ce6a40e7fd16295752d18aed2d541.html */
		for(c = 0; c < ngam; c++)
		{
			double s = 0.0;
			for(i = 0; i < narrays; i++)
				s += Z2[i + (size_t) c * narrays] * z[i];
			dl[c] = s;
			sol[c] = s;
		}
		for(i = 0; i < ngam * ngam; i++)
			Asolve[i] = info2[i];
		F77_CALL(dgesv)(&ngam, &one, Asolve, &ngam, ipiv, sol, &ngam, &info);
		if(info != 0)
			solok = 0;

		/* Update the variance-model coefficients and convert them back to array
		 * weights. Since log variance is Z2 gam, relative precision weights are
		 * exp(-Z2 gam). */
		for(c = 0; c < ngam; c++)
			gam[c] += sol[c];
		for(i = 0; i < narrays; i++)
		{
			double s = 0.0;
			for(c = 0; c < ngam; c++)
				s += Z2[i + (size_t) c * narrays] * gam[c];
			w[i] = exp(-s);
		}

		for(c = 0; c < ngam; c++)
			num += dl[c] * sol[c];
		convcrit = num / denom / ngam;

		if(trace)
		{
			double wmin = w[0], wmax = w[0];
			for(i = 1; i < narrays; i++)
			{
				if(w[i] < wmin)
					wmin = w[i];
				if(w[i] > wmax)
					wmax = w[i];
			}
			Rprintf("%d %g %g %g\n", iter, convcrit, wmin, wmax);
		}

		if(!solok || !R_FINITE(convcrit))
		{
			status = 1;
			break;
		}
		if(convcrit < tol)
		{
			status = 0;
			break;
		}
	}
	if(iter > maxiter)
	{
		iter = maxiter;
		status = 2;
	}

	for(t = 0; t < nthreads; t++)
		freeaw(aws + t);
	R_Free(aws);
	R_Free(Z);
	R_Free(Z2tZ2);
	R_Free(info2);
	R_Free(z);
	R_Free(dl);
	R_Free(Asolve);
	R_Free(sol);
	R_Free(gam);
	R_Free(ipiv);

	*iter_out = iter;
	return status;
}
