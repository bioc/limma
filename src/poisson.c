/*
 * poisson.c -- per-gene Poisson GLM (log link) with a per-sample offset.
 *
 * For each gene (row of X) fit  y_j ~ Poisson(mu_j), log(mu_j) = offset_j +
 * (design beta)_j  by IRLS and return the fitted means mu (ngenes x narrays).
 * voomLmFit uses this in place of edgeR's glmFit(dispersion=0) to flag the
 * structural-zero fitted values among count rows that contain many zeros.
 *
 * Start: the IRLS begins from the warm-start linear predictor `start` (the
 * non-offset part X beta, ngenes x narrays) via eta = offset + start; voomLmFit
 * passes the lmFit fit rescaled to the log-count scale. The converged Poisson MLE
 * is unique, so the start only affects the iteration count.
 *
 * For log-link Poisson the working weight is a_j = mu_j and the working response
 * is u_j = (design beta)_j + (y_j - mu_j)/mu_j, so each step is a weighted
 * least-squares solve of (sqrt(a).design, sqrt(a).u). mu is floored at MU_FLOOR in
 * the weight and working response so a structural-zero group (mu -> 0) cannot
 * produce 0/0; the reported fitted value mu = exp(eta) may still be ~0.
 *
 * Both the start and every WLS step use R's LINPACK dqrls (R_ext/Applic.h) -- the
 * rank-revealing routine lm.fit uses, resident in libR. Genes are independent, so
 * the loop is parallelised over genes with OpenMP.
 */

#include <math.h>
#include <R.h>
#include <R_ext/Applic.h>   /* dqrls */
#include <R_ext/RS.h>       /* F77_CALL */
#ifdef _OPENMP
#include <omp.h>
#endif
#include "limma.h"

#define MU_FLOOR 1e-10
#define POIS_MAXIT 25
#define POIS_TOL 1e-8

static void allocpois(poisws *pws, int narrays, int p)
{
	pws->y = R_Calloc(narrays, double);
	pws->uw = R_Calloc(narrays, double);
	pws->eta = R_Calloc(narrays, double);
	pws->mu = R_Calloc(narrays, double);
	pws->xw = R_Calloc((size_t) narrays * p, double);
	pws->coef = R_Calloc(p, double);
	pws->beta = R_Calloc(p, double);
	pws->resid = R_Calloc(narrays, double);
	pws->effec = R_Calloc(narrays, double);
	pws->qraux = R_Calloc(p, double);
	pws->work = R_Calloc(2 * p, double);
	pws->pivot = R_Calloc(p, int);
}

static void freepois(poisws *pws)
{
	R_Free(pws->y);
	R_Free(pws->uw);
	R_Free(pws->eta);
	R_Free(pws->mu);
	R_Free(pws->xw);
	R_Free(pws->coef);
	R_Free(pws->beta);
	R_Free(pws->resid);
	R_Free(pws->effec);
	R_Free(pws->qraux);
	R_Free(pws->work);
	R_Free(pws->pivot);
}

/*
 * Inputs:
 *   y, mu: length-n count and fitted-mean vectors.
 * Returns:
 *   Poisson deviance without the conventional factor of 2.
 * Notes:
 *   Used only as the IRLS convergence statistic, so only relative changes
 *   matter. The y == 0 contribution is mu.
 */
static double poisdev(const double *y, const double *mu, int n)
{
	int i;
	double s = 0.0;
	for(i = 0; i < n; i++)
		s += (y[i] > 0.0) ? (y[i] * log(y[i] / mu[i]) - (y[i] - mu[i])) : mu[i];
	return s;
}

/*
 * Inputs:
 *   pws->xw: narrays x p weighted design matrix, overwritten by dqrls.
 *   pws->uw: length-narrays weighted working response.
 * Outputs:
 *   pws->beta receives coefficients in original column order.
 * Notes:
 *   Aliased coefficients are set to 0. This is an internal IRLS WLS solve, not
 *   a user-visible coefficient fit.
 */
static void wlssolve(poisws *pws, int narrays, int p)
{
	int rank = 0, ny = 1, j;
	double tol = 1e-7;
	for(j = 0; j < p; j++)
		pws->pivot[j] = j + 1;

	/* dqrls (R LINPACK-based) solves real least-squares problems by
	 * rank-revealing QR decomposition with column pivoting. limma uses it here
	 * to solve one Poisson IRLS weighted least-squares subproblem:
	 *   - xw (narrays x p, lda=narrays) enters as sqrt(mu) * design and is
	 *     overwritten by compact QR storage;
	 *   - uw enters as sqrt(mu) times the working response, with ny=1;
	 *   - tol=1e-7 controls rank detection;
	 *   - coef receives estimable coefficients in pivot order;
	 *   - resid and effec receive residual and Q^T uw working outputs;
	 *   - rank and the 1-based pivot vector identify estimable columns;
	 *   - qraux/work hold Householder and temporary QR data.
	 * Aliased columns are set to zero when beta is unpivoted below.
	 * Equivalent R operation: lm.fit(xw, uw) for the current IRLS working
	 * response and weighted design.
	 * Netlib references:
	 *   https://www.netlib.org/linpack/dqrdc.f
	 *   https://www.netlib.org/linpack/dqrsl.f */
	F77_CALL(dqrls)(pws->xw, &narrays, &p, pws->uw, &ny, &tol, pws->coef, pws->resid, pws->effec, &rank, pws->pivot, pws->qraux, pws->work);
	for(j = 0; j < p; j++)
		pws->beta[j] = 0.0;
	for(j = 0; j < rank; j++)
		pws->beta[pws->pivot[j] - 1] = pws->coef[j];
}

/*
 * Inputs:
 *   design: narrays x p design matrix.
 *   offset: length-narrays offset vector.
 *   pws->beta: length-p current coefficients.
 * Outputs:
 *   pws->eta and pws->mu receive eta = offset + design beta and exp(eta).
 * Returns:
 *   Nothing.
 * Notes:
 *   Called after each WLS solve to refresh the Poisson mean.
 */
static void etamu(const double *design, const double *offset, int narrays, int p, poisws *pws)
{
	int i, j;
	for(i = 0; i < narrays; i++)
	{
		double e = offset[i];
		for(j = 0; j < p; j++)
			e += design[i + (size_t) j * narrays] * pws->beta[j];
		pws->eta[i] = e;
		pws->mu[i] = exp(e);
	}
}

/*
 * Inputs:
 *   X: ngenes x narrays count matrix.
 *   design: narrays x p design matrix.
 *   offset: length-narrays offset vector.
 *   start: ngenes x narrays warm-start non-offset linear predictor.
 *   gene: zero-based row of X; pws is reusable scratch for one thread.
 * Outputs:
 *   out: ngenes x narrays fitted Poisson means, filled for this gene only.
 * Returns:
 *   Nothing.
 * Notes:
 *   Runs up to POIS_MAXIT IRLS iterations with dqrls WLS steps. The working
 *   mean is floored for numerical stability, but the stored output is exp(eta).
 */
void poisgene(const double *X, const double *design, const double *offset, const double *start, int ngenes, int narrays, int p, int gene, double *out, poisws *pws)
{
	int i, j, it;
	double dev, devnew;
	const double *yptr = X + gene;
	const double *sptr = start + gene;

	/* Counts are read from this gene's row. The warm start stores only the
	 * non-offset linear predictor, so eta starts at offset + start before
	 * exponentiation to the current Poisson mean. */
	for(i = 0; i < narrays; i++)
	{
		pws->y[i] = yptr[(size_t) i * ngenes];
		pws->eta[i] = offset[i] + sptr[(size_t) i * ngenes];
		pws->mu[i] = exp(pws->eta[i]);
	}
	dev = poisdev(pws->y, pws->mu, narrays);

	for(it = 0; it < POIS_MAXIT; it++)
	{
		/* working WLS problem: sqrt(mu).design, sqrt(mu).[ (eta-offset) + (y-mu)/mu ];
		 * mu floored in weight and working response to stay finite as mu -> 0 */
		for(i = 0; i < narrays; i++)
		{
			double mui = pws->mu[i];
			double muf = mui < MU_FLOOR ? MU_FLOOR : mui;
			double sa = sqrt(muf);
			pws->uw[i] = sa * ((pws->eta[i] - offset[i]) + (pws->y[i] - mui) / muf);
			for(j = 0; j < p; j++)
				pws->xw[i + (size_t) j * narrays] = sa * design[i + (size_t) j * narrays];
		}
		wlssolve(pws, narrays, p);
		etamu(design, offset, narrays, p, pws);
		devnew = poisdev(pws->y, pws->mu, narrays);
		if(fabs(devnew - dev) / (fabs(devnew) + 0.1) < POIS_TOL)
		{
			dev = devnew;
			break;
		}
		dev = devnew;
	}

	for(i = 0; i < narrays; i++)
		out[gene + (size_t) i * ngenes] = pws->mu[i];
}

/*
 * Inputs:
 *   X: ngenes x narrays count matrix.
 *   design: narrays x p design matrix.
 *   offset: length-narrays offset vector.
 *   start: ngenes x narrays warm-start non-offset linear predictor.
 *   nthreads: requested OpenMP thread count.
 * Outputs:
 *   out: ngenes x narrays fitted Poisson means.
 * Notes:
 *   Genes are independent. One poisws scratch workspace is allocated per active
 *   thread and reused across genes.
 */
void pois(const double *X, const double *design, const double *offset, const double *start, int ngenes, int narrays, int p, int nthreads, double *out)
{
	int gene, t;
	poisws *pws;

	nthreads = clampthreads(nthreads, ngenes);

	pws = R_Calloc(nthreads, poisws);
	for(t = 0; t < nthreads; t++)
		allocpois(pws + t, narrays, p);

	#ifdef _OPENMP
	#pragma omp parallel for schedule(static) num_threads(nthreads)
	#endif
	for(gene = 0; gene < ngenes; gene++)
	{
		int thread = 0;
		#ifdef _OPENMP
		thread = omp_get_thread_num();
		#endif
		poisgene(X, design, offset, start, ngenes, narrays, p, gene, out, pws + thread);
	}

	for(t = 0; t < nthreads; t++)
		freepois(pws + t);
	R_Free(pws);
}
