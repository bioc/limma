/*
 * lm.c -- gene-by-gene linear model fit with probe weights and/or
 * missing values. This is the C implementation of the slow path of lm.series()
 * (the .Call("lmfit") branch): R takes it only when probe-wise weights
 * or non-finite values make the QR differ between genes. The weight-free,
 * fully-observed "fast path" (one shared lm.fit sweep) stays in R.
 *
 * Model for gene g:  y_g = X beta_g + e,  e ~ N(0, sigma^2 W^{-1}),
 * with W = diag(probe weights). Multiplying each observation by sqrt(w)
 * (prewhitening) turns weighted LS into ordinary LS:
 *     y~ = sqrt(w) .* y,   X~ = sqrt(w) .* X,   then solve via qrfit().
 * Only finite observations are kept, so each gene has its own design submatrix.
 */

#include <math.h>
#include <R.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "limma.h"

/*
 * Inputs:
 *   m: ngenes x narrays expression matrix in R column-major order.
 *   design: narrays x nbeta design matrix in column-major order.
 *   weights: optional ngenes x narrays prior-weight matrix, or NULL.
 *   gene: zero-based row of m to fit; lws is reusable scratch for one thread.
 * Outputs:
 *   beta, stdev_unscaled, sigma and df_residual are filled for this gene only.
 * Returns:
 *   Status from qrfit(); currently 0 for all normal early-return and fit paths.
 * Notes:
 *   Non-finite expression values are dropped. Weighted fits prewhiten the
 *   response and selected design rows by sqrt(weight).
 */
int lmgene(const double *m, const double *design, const double *weights, int ngenes, int narrays, int nbeta, const double *contrasts, int ncont, int gene, double *beta, double *stdev_unscaled, double *sigma, int *df_residual, lmws *lws)
{
	int array, i, j, nobs = 0;
	int ncols = contrasts ? ncont : nbeta;
	/* column-major walks for this gene's row of the data/weights and outputs */
	const double *mptr = m + gene;
	const double *wptr = weights ? weights + gene : NULL;
	double *bptr = beta + gene, *sptr = stdev_unscaled + gene;

	/* Self-initialise this gene's output row to NA/0 so every early-return path
	 * (no observations, rank 0) and every dropped column/contrast is well defined.
	 * This reproduces the prefilled NA of the R reference loop. The output has
	 * ncont columns under contrasts, otherwise nbeta. */
	for(j = 0; j < ncols; j++, bptr += ngenes, sptr += ngenes)
	{
		*bptr = NA_REAL;
		*sptr = NA_REAL;
	}
	sigma[gene] = NA_REAL;
	df_residual[gene] = 0;

	/* Collect finite observations and their prewhitening factor scale = sqrt(w);
	 * y~ = scale * y. scale is cached so the design rows below reuse it. */
	for(array = 0; array < narrays; array++)
	{
		double value = *mptr;
		if(R_FINITE(value))
		{
			double scale = weights ? sqrt(*wptr) : 1.0;
			lws->obser[nobs] = array;
			lws->scale[nobs] = scale;
			lws->y[nobs] = value * scale;
			nobs++;
		}
		mptr += ngenes;
		if(weights)
			wptr += ngenes;
	}

	if(nobs == 0)
		return 0;                       /* nothing observed; leave NA/0 row */

	/* Build the prewhitened design X~[obs, ] = scale .* design[obs, ] (column-major) */
	for(j = 0; j < nbeta; j++)
	{
		for(i = 0; i < nobs; i++)
		{
			array = lws->obser[i];
			lws->x[i + j * nobs] = design[array + j * narrays] * lws->scale[i];
		}
	}

	/* Ordinary LS on (X~, y~) -> beta, stdev.unscaled, sigma, df for this gene.
	 * When contrasts is non-NULL, qrfit returns contrast-space estimates instead. */
	return qrfit(nobs, nbeta, gene, ngenes, contrasts, ncont, beta, stdev_unscaled, sigma, df_residual, lws);
}

static void alloclm(lmws *lws, int narrays, int nbeta)
{
	lws->x = R_Calloc((size_t) narrays * nbeta, double);
	lws->y = R_Calloc(narrays, double);
	lws->coeff = R_Calloc(nbeta, double);
	lws->resid = R_Calloc(narrays, double);
	lws->effec = R_Calloc(narrays, double);
	lws->qraux = R_Calloc(nbeta, double);
	lws->work = R_Calloc(2 * nbeta, double);
	lws->rinv = R_Calloc((size_t) nbeta * nbeta, double);
	lws->cest = R_Calloc(nbeta, double);
	lws->scale = R_Calloc(narrays, double);
	lws->obser = R_Calloc(narrays, int);
	lws->pivot = R_Calloc(nbeta, int);
	lws->estpos = R_Calloc(nbeta, int);
}

static void freelm(lmws *lws)
{
	R_Free(lws->x);
	R_Free(lws->y);
	R_Free(lws->coeff);
	R_Free(lws->resid);
	R_Free(lws->effec);
	R_Free(lws->qraux);
	R_Free(lws->work);
	R_Free(lws->rinv);
	R_Free(lws->cest);
	R_Free(lws->scale);
	R_Free(lws->obser);
	R_Free(lws->pivot);
	R_Free(lws->estpos);
}

/*
 * Inputs:
 *   m: ngenes x narrays expression matrix.
 *   design: narrays x nbeta design matrix.
 *   weights: optional ngenes x narrays prior-weight matrix, or NULL.
 *   nthreads: requested OpenMP thread count.
 * Outputs:
 *   beta and stdev_unscaled: ngenes x nbeta matrices.
 *   sigma and df_residual: length-ngenes vectors.
 * Notes:
 *   Genes are independent. One lmws scratch workspace is allocated per active
 *   thread and reused across genes; builds without OpenMP run serially.
 */
void lm(const double *m, const double *design, const double *weights, int ngenes, int narrays, int nbeta, const double *contrasts, int ncont, int nthreads, double *beta, double *stdev_unscaled, double *sigma, int *df_residual)
{
	int gene;
	lmws *lws;

	nthreads = clampthreads(nthreads, ngenes);

	/* one reusable scratch workspace per thread (allocated once, not per gene) */
	lws = R_Calloc(nthreads, lmws);
	for(gene = 0; gene < nthreads; gene++)
		alloclm(lws + gene, narrays, nbeta);

	/* genes are independent; each thread fits its share using its own workspace */
	#ifdef _OPENMP
	#pragma omp parallel for schedule(static) num_threads(nthreads)
	#endif
	for(gene = 0; gene < ngenes; gene++)
	{
		int thread = 0;
		#ifdef _OPENMP
		thread = omp_get_thread_num();
		#endif
		lmgene(m, design, weights, ngenes, narrays, nbeta, contrasts, ncont, gene, beta, stdev_unscaled, sigma, df_residual, lws + thread);
	}
	for(gene = 0; gene < nthreads; gene++)
		freelm(lws + gene);
	R_Free(lws);
}
