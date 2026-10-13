/*
 * gls.c -- gene-by-gene generalized least squares with a known
 * correlation structure. C implementation of the slow path of gls.series()
 * (the .Call("glsfit") branch), used when probe weights or missing
 * values are present; the weight-free fast path stays in R.
 *
 * Model for gene g:  y ~ N(X beta, sigma^2 V) with V = W^{-1/2} R W^{-1/2},
 * where W holds the probe weights and R is block diagonal with exchangeable
 * blocks: R_b = (1 - rho) I + rho J (unit diagonal, correlation rho within a
 * block, zero between blocks). Within-array duplicates are blocks of size ndups.
 *
 * GLS is solved by prewhitening with T = R^{-1/2} W^{1/2}, which satisfies
 * T^T T = V^{-1}, followed by ordinary LS on (T X, T y) via qrfit(). The
 * symmetric inverse square root of an exchangeable block has the closed form
 *     R_b^{-1/2} z = (z - c_b S_b 1) / s,
 *     s = sqrt(1 - rho), lambda_b = 1 + (n_b - 1) rho,
 *     c_b = rho / (lambda_b + s sqrt(lambda_b)),
 * where n_b is the number of observed arrays in block b and S_b = 1^T z is the
 * block sum. The observed arrays of a block are again an exchangeable block,
 * so missing values only change n_b. Each gene therefore costs O(nobs nbeta)
 * for the whitening instead of a Cholesky factorization of the observed V.
 * This is the C analogue of R's .invSqrtBlockCor() + lm.fit().
 */

#include <math.h>
#include <R.h>
#include <R_ext/RS.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "limma.h"

/*
 * Inputs:
 *   m: ngenes x narrays expression matrix in R column-major order.
 *   design: narrays x nbeta design matrix.
 *   block: length-narrays one-based block labels in 1..nblocks.
 *   correlation: common within-block correlation rho; the caller guarantees
 *     rho < 1 and 1 + (n_b - 1) rho > 0 for every full block.
 *   weights: optional ngenes x narrays prior-weight matrix, or NULL.
 *   gene: zero-based row of m to fit; gws is reusable scratch for one thread.
 * Outputs:
 *   beta, stdev_unscaled, sigma and df_residual are filled for this gene only.
 * Returns:
 *   Status from qrfit(); 0 on success or normal no-observation paths.
 * Notes:
 *   Finite observations are compacted, scaled by sqrt(w), whitened block by
 *   block with the closed-form R^{-1/2}, and fitted by ordinary least squares.
 */
int glsgene(const double *m, const double *design, const int *block, int nblocks, double correlation, const double *weights, int ngenes, int narrays, int nbeta, const double *contrasts, int ncont, int gene, double *beta, double *stdev_unscaled, double *sigma, int *df_residual, glsws *gws)
{
	int array, b, i, j, nobs = 0, all_zero = 1;
	int ncols = contrasts ? ncont : nbeta;
	const double s = sqrt(1.0 - correlation);
	const double *mptr = m + gene;
	const double *wptr = weights ? weights + gene : NULL;
	double *bptr = beta + gene, *sptr = stdev_unscaled + gene;
	double *xptr, *sumptr;

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

	/* W^{1/2} scaling of y and X, accumulating observed block counts and the
	 * block sums of the scaled [y | X] columns (column 0 is y, column j+1 is X_j) */
	for(b = 0; b < nblocks; b++)
		gws->bcount[b] = 0;
	for(j = 0; j < (nbeta + 1) * nblocks; j++)
		gws->bsum[j] = 0.0;
	for(i = 0; i < nobs; i++)
	{
		double w = gws->qr.scale[i];
		b = block[gws->qr.obser[i]] - 1;
		gws->bcount[b]++;
		gws->qr.y[i] *= w;
		gws->bsum[b] += gws->qr.y[i];
		for(j = 0, xptr = gws->qr.x + i, sumptr = gws->bsum + nblocks + b; j < nbeta; j++, xptr += nobs, sumptr += nblocks)
		{
			*xptr *= w;
			*sumptr += *xptr;
		}
	}

	/* closed-form coefficient c_b for each block with observed arrays */
	for(b = 0; b < nblocks; b++)
	{
		if(gws->bcount[b] > 0)
		{
			double lambda = 1.0 + (gws->bcount[b] - 1) * correlation;
			gws->bcoef[b] = correlation / (lambda + s * sqrt(lambda));
		}
	}

	/* whiten the response: y_tilde = R^{-1/2} (W^{1/2} y) */
	for(i = 0; i < nobs; i++)
	{
		b = block[gws->qr.obser[i]] - 1;
		gws->qr.y[i] = (gws->qr.y[i] - gws->bcoef[b] * gws->bsum[b]) / s;
	}

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

	/* whiten the design columns the same way: X_tilde = R^{-1/2} (W^{1/2} X) */
	for(i = 0; i < nobs; i++)
	{
		double c;
		b = block[gws->qr.obser[i]] - 1;
		c = gws->bcoef[b];
		for(j = 0, xptr = gws->qr.x + i, sumptr = gws->bsum + nblocks + b; j < nbeta; j++, xptr += nobs, sumptr += nblocks)
			*xptr = (*xptr - c * *sumptr) / s;
	}

	/* GLS reduces to ordinary LS on the prewhitened (X~, y~) */
	return qrfit(nobs, nbeta, gene, ngenes, contrasts, ncont, beta, stdev_unscaled, sigma, df_residual, &gws->qr);
}

static void allocgls(glsws *gws, int narrays, int nbeta, int nblocks)
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
	gws->bsum = R_Calloc((size_t) nblocks * (nbeta + 1), double);
	gws->bcoef = R_Calloc(nblocks, double);
	gws->bcount = R_Calloc(nblocks, int);
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
	R_Free(gws->bsum);
	R_Free(gws->bcoef);
	R_Free(gws->bcount);
}

/*
 * Inputs:
 *   m: ngenes x narrays expression matrix.
 *   design: narrays x nbeta design matrix.
 *   block: length-narrays one-based block labels in 1..nblocks.
 *   correlation: common within-block correlation rho.
 *   weights: optional ngenes x narrays prior-weight matrix, or NULL.
 *   nthreads: requested OpenMP thread count.
 * Outputs:
 *   beta and stdev_unscaled: ngenes x nbeta matrices.
 *   sigma and df_residual: length-ngenes vectors.
 * Returns:
 *   0 if all gene fits succeed; 1 if the correlation matrix is not positive
 *   definite (rho >= 1 or 1 + (n_b - 1) rho <= 0 for some full block), in which
 *   case no gene is fitted; otherwise the first nonzero gene status.
 * Notes:
 *   Observed sub-blocks are never larger than full blocks, so checking the full
 *   block sizes once guarantees lambda_b > 0 for every gene. One glsws scratch
 *   workspace is allocated per active thread.
 */
int gls(const double *m, const double *design, const int *block, int nblocks, double correlation, const double *weights, int ngenes, int narrays, int nbeta, const double *contrasts, int ncont, int nthreads, double *beta, double *stdev_unscaled, double *sigma, int *df_residual)
{
	int array, b, gene, failed = 0;
	int *status, *bsize;
	glsws *gws;

	if(correlation >= 1.0)
		return 1;
	bsize = R_Calloc(nblocks, int);
	for(array = 0; array < narrays; array++)
		bsize[block[array] - 1]++;
	for(b = 0; b < nblocks; b++)
	{
		if(1.0 + (bsize[b] - 1) * correlation <= 0.0)
			failed = 1;
	}
	R_Free(bsize);
	if(failed)
		return failed;

	nthreads = clampthreads(nthreads, ngenes);

	gws = R_Calloc(nthreads, glsws);
	status = R_Calloc(ngenes, int);
	for(gene = 0; gene < nthreads; gene++)
		allocgls(gws + gene, narrays, nbeta, nblocks);

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static) num_threads(nthreads)
	#endif
	for(gene = 0; gene < ngenes; gene++)
	{
		int thread = 0;
		#ifdef _OPENMP
		thread = omp_get_thread_num();
		#endif
		status[gene] = glsgene(m, design, block, nblocks, correlation, weights, ngenes, narrays, nbeta, contrasts, ncont, gene, beta, stdev_unscaled, sigma, df_residual, gws + thread);
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
	return failed;
}
