/*
 * threads.c -- shared OpenMP thread-count resolver for the genewise kernels.
 *
 * lm, gls, dupcor, poisson and awreml each call clampthreads() to turn a
 * requested thread count into the count actually used. This replaces the guard
 * that was inlined identically in every kernel and makes the C side robust to
 * any requested value on its own: the R wrappers no longer pre-validate
 * nthreads (the former .checkNThreads).
 */

#ifdef _OPENMP
#include <omp.h>
#endif
#include "limma.h"

/*
 * Inputs:
 *   nthreads: requested thread count (any int; not pre-validated by R).
 *   nwork: number of independent work items (genes) the caller will loop over.
 * Returns:
 *   The thread count to use: at least 1, capped at the available OpenMP threads,
 *   dropped to 1 when the request exceeds the work, and forced to 1 without
 *   OpenMP or inside an existing parallel region (nested parallelism stays
 *   serial for reproducibility).
 */
int clampthreads(int nthreads, int nwork)
{
	if(nthreads < 1)
		nthreads = 1;
	if(nwork > 0 && nthreads > nwork)
		nthreads = 1;
	#ifdef _OPENMP
	int hi = omp_get_max_threads();
	if(hi >= 1 && nthreads > hi)
		nthreads = hi;
	if(omp_in_parallel())
		nthreads = 1;
	#else
	nthreads = 1;
	#endif
	return nthreads;
}
