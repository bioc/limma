/*
 * gls.c -- gene-by-gene generalized least squares with a known
 * correlation structure. C implementation of the slow path of gls.series()
 * (the .Call("glsfit") branch), used when probe weights or missing
 * values are present; the weight-free fast path stays in R.
 *
 * Model for gene g:  y ~ N(X beta, sigma^2 V), where V is the correlation
 * matrix implied by the consensus correlation (block-diagonal for within-array
 * duplicates, or Z rho Z^T for general blocks; unit diagonal). GLS is solved by
 * prewhitening with the Cholesky factor of V:
 *     V = U^T U            (dpotrf, upper triangle)
 *     solve U^T y~ = y     (dtrsv, upper+transpose)   => y~ = U^{-T} y
 *     solve U^T X~ = X     (per design column)         => X~ = U^{-T} X
 *     then ordinary LS on (X~, y~) via qrfit().
 * This is the C analogue of R's backsolve(cholV, ., transpose=TRUE) + lm.fit().
 *
 * Probe weights enter through V: Cov(y_i) carries 1/w_i, so the (i,j) entry of
 * the working V is cormatrix[i,j] / (sqrt(w_i) sqrt(w_j)) (the diagonal is
 * 1/w_i). The per-observation sqrt(w) is cached during the scan as qr.scale.
 */

#include <math.h>
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "limma.h"

/*
 * Inputs:
 *   m: ngenes x narrays expression matrix in R column-major order.
 *   design: narrays x nbeta design matrix.
 *   cormatrix: narrays x narrays working correlation matrix.
 *   weights: optional ngenes x narrays prior-weight matrix, or NULL.
 *   gene: zero-based row of m to fit; gws is reusable scratch for one thread.
 * Outputs:
 *   beta, stdev_unscaled, sigma and df_residual are filled for this gene only.
 * Returns:
 *   0 on success or normal no-observation paths; otherwise the LAPACK dpotrf
 *   leading-minor index for a non-positive-definite working covariance.
 * Notes:
 *   Finite observations are compacted, the observed covariance submatrix is
 *   weighted when needed, and the model is fitted after Cholesky whitening.
 */
int glsgene(const double *m, const double *design, const double *cormatrix, const double *weights, int ngenes, int narrays, int nbeta, const double *contrasts, int ncont, int gene, double *beta, double *stdev_unscaled, double *sigma, int *df_residual, const double *cholc, glsws *gws)
{
	int array, i, j, nobs = 0, info = 0, all_zero = 1;
	int ncols = contrasts ? ncont : nbeta;
	const int inc = 1;
	const char upper = 'U', left = 'L', transpose = 'T', nonunit = 'N';
	const double one = 1.0;
	const double *mptr = m + gene;
	const double *wptr = weights ? weights + gene : NULL;
	double *bptr = beta + gene, *sptr = stdev_unscaled + gene;

	/* self-initialise this gene's output row (see lm.c for rationale);
	 * ncont columns under contrasts, otherwise nbeta */
	for(j = 0; j < ncols; j++, bptr += ngenes, sptr += ngenes)
	{
		*bptr = NA_REAL;
		*sptr = NA_REAL;
	}
	sigma[gene] = NA_REAL;
	df_residual[gene] = 0;

	/* keep finite observations; store y (unscaled here) and cache scale = sqrt(w) */
	for(array = 0; array < narrays; array++)
	{
		double value = *mptr;
		if(R_FINITE(value))
		{
			gws->qr.obser[nobs] = array;
			gws->qr.y[nobs] = value;
			gws->qr.scale[nobs] = weights ? sqrt(*wptr) : 1.0;
			nobs++;
		}
		mptr += ngenes;
		if(weights)
			wptr += ngenes;
	}
	if(nobs == 0)
		return 0;

	/* design rows for the observed arrays; all_zero flags a no-coefficient gene */
	for(j = 0; j < nbeta; j++)
	{
		for(i = 0; i < nobs; i++)
		{
			array = gws->qr.obser[i];
			gws->qr.x[i + j * nobs] = design[array + j * narrays];
			if(gws->qr.x[i + j * nobs] != 0.0)
				all_zero = 0;
		}
	}

	/* Complete gene (all arrays observed): reuse the one-time Cholesky factor of the
	 * full correlation matrix. Since chol(D^-1 C D^-1) = chol(C) D^-1 with
	 * D = diag(sqrt(w)), scaling the observed rows of y and X by sqrt(w) and whitening
	 * with the shared factor cholc reproduces the per-gene result without building V or
	 * calling dpotrf. Genes with missing values use a submatrix of C and fall through to
	 * the per-gene factorization below. cholc has leading dimension narrays. */
	if(cholc && nobs == narrays)
	{
		for(i = 0; i < nobs; i++)
			gws->qr.y[i] *= gws->qr.scale[i];
		/* DTRSV (Netlib BLAS) solves a real triangular matrix-vector system.
		 * limma uses it here to whiten the complete-gene response with the shared
		 * Cholesky factor:
		 *   - UPLO='U' uses the upper-triangular factor cholc;
		 *   - TRANS='T' solves cholc^T x = b;
		 *   - DIAG='N' treats cholc as non-unit triangular;
		 *   - N=nobs, LDA=narrays and INCX=1 describe the full-array layout;
		 *   - gws->qr.y is overwritten by y_tilde = cholc^{-T} y.
		 * The transformed response can be fitted by ordinary least squares.
		 * Equivalent R operation: backsolve(cholc, y, transpose=TRUE).
		 * Netlib references:
		 *   https://www.netlib.org/blas/dtrsv.f */
		F77_CALL(dtrsv)(&upper, &transpose, &nonunit, &nobs, cholc, &narrays, gws->qr.y, &inc FCONE FCONE FCONE);
		if(all_zero)
		{
			double sumsq = 0.0;
			for(i = 0; i < nobs; i++)
				sumsq += gws->qr.y[i] * gws->qr.y[i];
			df_residual[gene] = nobs;
			sigma[gene] = sqrt(sumsq / nobs);
			return 0;
		}
		for(j = 0; j < nbeta; j++)
			for(i = 0; i < nobs; i++)
				gws->qr.x[i + j * nobs] *= gws->qr.scale[i];
		/* DTRSM (Netlib BLAS) solves real triangular matrix equations with
		 * multiple right-hand sides. limma uses it here to whiten all complete-gene
		 * design columns with the shared Cholesky factor:
		 *   - SIDE='L' solves with cholc on the left;
		 *   - UPLO='U', TRANS='T' and DIAG='N' solve cholc^T X_tilde = X;
		 *   - M=nobs and N=nbeta describe the design block;
		 *   - ALPHA=1 leaves the right-hand side unscaled;
		 *   - gws->qr.x is overwritten by X_tilde = cholc^{-T} X.
		 * Ordinary least squares on (X_tilde, y_tilde) is the GLS fit.
		 * Equivalent R operation: backsolve(cholc, X, transpose=TRUE).
		 * Netlib references:
		 *   https://www.netlib.org/blas/dtrsm.f */
		F77_CALL(dtrsm)(&left, &upper, &transpose, &nonunit, &nobs, &nbeta, &one, cholc, &narrays, gws->qr.x, &nobs FCONE FCONE FCONE FCONE);
		return qrfit(nobs, nbeta, gene, ngenes, contrasts, ncont, beta, stdev_unscaled, sigma, df_residual, &gws->qr);
	}

	/* working covariance V on the observed arrays: cormatrix scaled by weights,
	 * V[i,j] = cormatrix[row_i, col_j] / (sqrt(w_i) sqrt(w_j)). Only the upper triangle
	 * is filled because dpotrf('U') reads just that. */
	for(j = 0; j < nobs; j++)
	{
		int col = gws->qr.obser[j];
		double col_scale = gws->qr.scale[j];
		for(i = 0; i <= j; i++)
		{
			int row = gws->qr.obser[i];
			double value = cormatrix[row + col * narrays];
			if(weights)
				value /= gws->qr.scale[i] * col_scale;
			gws->v[i + j * nobs] = value;
		}
	}

	/* DPOTRF (Netlib LAPACK) computes the Cholesky factorization of a real
	 * symmetric positive-definite matrix A. limma uses it here to factor the
	 * observed gene-specific covariance submatrix:
	 *   - UPLO='U' says the upper triangle of gws->v holds V on entry;
	 *   - N=nobs and LDA=nobs describe the observed covariance submatrix;
	 *   - gws->v is overwritten in its upper triangle by U where V = U^T U;
	 *   - INFO=0 means success; INFO=i>0 means the leading principal minor of
	 *     order i is not positive definite, which the R shim reports as an error.
	 * The factor U is reused below to whiten the response and design.
	 * Equivalent R operation: chol(V).
	 * Netlib references:
	 *   https://www.netlib.org/lapack/explore-html/d2/d09/group__potrf_ga84e90859b02139934b166e579dd211d4.html */
	F77_CALL(dpotrf)(&upper, &nobs, gws->v, &nobs, &info FCONE);
	if(info != 0)
		return info;

	/* DTRSV (Netlib BLAS) solves a real triangular matrix-vector system. limma
	 * uses it here to whiten the response with the gene-specific Cholesky factor:
	 *   - UPLO='U' uses the Cholesky factor U stored in gws->v;
	 *   - TRANS='T' solves U^T x = b;
	 *   - DIAG='N' treats U as non-unit triangular;
	 *   - N=nobs, LDA=nobs and INCX=1 describe the matrix/vector layout;
	 *   - gws->qr.y is overwritten by y_tilde = U^{-T} y.
	 * This is the response half of the GLS whitening transform.
	 * Equivalent R operation: backsolve(U, y, transpose=TRUE).
	 * Netlib references:
	 *   https://www.netlib.org/blas/dtrsv.f */
	F77_CALL(dtrsv)(&upper, &transpose, &nonunit, &nobs, gws->v, &nobs, gws->qr.y, &inc FCONE FCONE FCONE);

	/* If every observed design entry is zero, there are no estimable
	 * coefficients. The whitened response is still meaningful, so keep
	 * coefficients as NA, set df to nobs, and estimate sigma from ||y_tilde||. */
	if(all_zero)
	{
		double sumsq = 0.0;
		for(i = 0; i < nobs; i++)
			sumsq += gws->qr.y[i] * gws->qr.y[i];
		df_residual[gene] = nobs;
		sigma[gene] = sqrt(sumsq / nobs);
		return 0;
	}

	/* DTRSM (Netlib BLAS) solves real triangular matrix equations with multiple
	 * right-hand sides. limma uses it here to whiten all design columns at once:
	 *   - SIDE='L' solves with U on the left;
	 *   - UPLO='U', TRANS='T' and DIAG='N' solve U^T X_tilde = X;
	 *   - M=nobs and N=nbeta describe the observed design block;
	 *   - ALPHA=1 leaves the right-hand side unscaled;
	 *   - gws->qr.x is overwritten by X_tilde = U^{-T} X.
	 * Ordinary least squares on (X_tilde, y_tilde) is the GLS fit with covariance V.
	 * Equivalent R operation: backsolve(U, X, transpose=TRUE).
	 * Netlib references:
	 *   https://www.netlib.org/blas/dtrsm.f */
	F77_CALL(dtrsm)(&left, &upper, &transpose, &nonunit, &nobs, &nbeta, &one, gws->v, &nobs, gws->qr.x, &nobs FCONE FCONE FCONE FCONE);

	/* GLS reduces to ordinary LS on the prewhitened (X~, y~) */
	return qrfit(nobs, nbeta, gene, ngenes, contrasts, ncont, beta, stdev_unscaled, sigma, df_residual, &gws->qr);
}

static void allocgls(glsws *gws, int narrays, int nbeta)
{
	gws->qr.x = R_Calloc((size_t) narrays * nbeta, double);
	gws->qr.y = R_Calloc(narrays, double);
	gws->qr.coeff = R_Calloc(nbeta, double);
	gws->qr.resid = R_Calloc(narrays, double);
	gws->qr.effec = R_Calloc(narrays, double);
	gws->qr.qraux = R_Calloc(nbeta, double);
	gws->qr.work = R_Calloc(2 * nbeta, double);
	gws->qr.rinv = R_Calloc((size_t) nbeta * nbeta, double);
	gws->qr.cest = R_Calloc(nbeta, double);
	gws->qr.scale = R_Calloc(narrays, double);
	gws->qr.obser = R_Calloc(narrays, int);
	gws->qr.pivot = R_Calloc(nbeta, int);
	gws->qr.estpos = R_Calloc(nbeta, int);
	gws->v = R_Calloc((size_t) narrays * narrays, double);
}

static void freegls(glsws *gws)
{
	R_Free(gws->qr.x);
	R_Free(gws->qr.y);
	R_Free(gws->qr.coeff);
	R_Free(gws->qr.resid);
	R_Free(gws->qr.effec);
	R_Free(gws->qr.qraux);
	R_Free(gws->qr.work);
	R_Free(gws->qr.rinv);
	R_Free(gws->qr.cest);
	R_Free(gws->qr.scale);
	R_Free(gws->qr.obser);
	R_Free(gws->qr.pivot);
	R_Free(gws->qr.estpos);
	R_Free(gws->v);
}

/*
 * Inputs:
 *   m: ngenes x narrays expression matrix.
 *   design: narrays x nbeta design matrix.
 *   cormatrix: narrays x narrays correlation matrix shared by all genes.
 *   weights: optional ngenes x narrays prior-weight matrix, or NULL.
 *   nthreads: requested OpenMP thread count.
 * Outputs:
 *   beta and stdev_unscaled: ngenes x nbeta matrices.
 *   sigma and df_residual: length-ngenes vectors.
 * Returns:
 *   0 if all gene fits succeed; otherwise the first nonzero dpotrf status.
 * Notes:
 *   One glsws scratch workspace is allocated per active thread. Any covariance
 *   failure is collected during the parallel loop and reported by the R shim.
 */
int gls(const double *m, const double *design, const double *cormatrix, const double *weights, int ngenes, int narrays, int nbeta, const double *contrasts, int ncont, int nthreads, double *beta, double *stdev_unscaled, double *sigma, int *df_residual)
{
	int gene, failed = 0, cinfo = 0;
	int *status;
	double *cholc;
	const double *shared;
	const char upper = 'U';
	glsws *gws;

	nthreads = clampthreads(nthreads, ngenes);

	gws = R_Calloc(nthreads, glsws);
	status = R_Calloc(ngenes, int);
	for(gene = 0; gene < nthreads; gene++)
		allocgls(gws + gene, narrays, nbeta);

	/* Factor the full correlation matrix once; complete genes reuse it instead of a
	 * per-gene dpotrf (see glsgene). If cormatrix is not positive definite the factor
	 * is unavailable (shared = NULL) and every gene takes the per-gene path, which
	 * reports the failure for any complete gene exactly as before. */
	cholc = R_Calloc((size_t) narrays * narrays, double);
	int ci, cj;
	for(cj = 0; cj < narrays; cj++)
		for(ci = 0; ci <= cj; ci++)
			cholc[ci + cj * narrays] = cormatrix[ci + cj * narrays];
	/* DPOTRF (Netlib LAPACK) computes the Cholesky factorization of a real
	 * symmetric positive-definite matrix A. limma uses it here to pre-factor
	 * the full correlation matrix for complete genes:
	 *   - UPLO='U' says the upper triangle of cholc holds cormatrix on entry;
	 *   - N=narrays and LDA=narrays describe the full matrix;
	 *   - cholc is overwritten in its upper triangle by U where C = U^T U;
	 *   - cinfo=0 means the shared factor is available; cinfo>0 leaves the
	 *     complete genes to the per-gene path, which reports the same failure.
	 * The shared U avoids repeating the same Cholesky factorization for every
	 * complete gene.
	 * Equivalent R operation: chol(cormatrix).
	 * Netlib references:
	 *   https://www.netlib.org/lapack/explore-html/d2/d09/group__potrf_ga84e90859b02139934b166e579dd211d4.html */
	F77_CALL(dpotrf)(&upper, &narrays, cholc, &narrays, &cinfo FCONE);
	shared = cinfo == 0 ? cholc : NULL;

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static) num_threads(nthreads)
	#endif
	for(gene = 0; gene < ngenes; gene++)
	{
		int thread = 0;
		#ifdef _OPENMP
		thread = omp_get_thread_num();
		#endif
		status[gene] = glsgene(m, design, cormatrix, weights, ngenes, narrays, nbeta, contrasts, ncont, gene, beta, stdev_unscaled, sigma, df_residual, shared, gws + thread);
	}
	for(gene = 0; gene < ngenes; gene++)
	{
		if(status[gene])
		{
			failed = status[gene];
			break;
		}
	}
	for(gene = 0; gene < nthreads; gene++)
		freegls(gws + gene);
	R_Free(gws);
	R_Free(status);
	R_Free(cholc);
	return failed;
}
