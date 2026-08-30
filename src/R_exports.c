/*
 * R_exports.c -- the .Call SEXP shims for limma's C kernels. Each function
 * unpacks its R arguments, calls the matching kernel (declared in limma.h) and
 * packs the result; all numerical work lives in the kernel .c files. The
 * routines are registered in init.c.
 */

#include <R.h>
#include <Rinternals.h>
#include "limma.h"
#include "R_exports.h"

/*
 * Inputs:
 *   m: numeric ngenes x narrays matrix.
 *   design: numeric narrays x nbeta design matrix.
 *   weight: optional numeric ngenes x narrays matrix, or R_NilValue.
 *   contrasts: optional numeric nbeta x ncont contrast matrix, or R_NilValue.
 *   nthr: requested thread count.
 * Outputs:
 *   Returns a named list with coefficients, stdev.unscaled, sigma and
 *   df.residual. With contrasts the coefficient/stdev matrices have ncont
 *   columns in contrast space; otherwise nbeta columns.
 * Returns:
 *   SEXP list allocated for R.
 * Notes:
 *   Inputs are coerced to REALSXP without modifying caller-owned objects, then
 *   passed to lm(), which fills newly allocated result vectors.
 */
SEXP lmfit(SEXP m, SEXP design, SEXP weight, SEXP contrasts, SEXP nthr)
{
	SEXP output, beta, stdev, sigma, df, names;
	SEXP mdim, ddim;
	int ngenes, narrays, nbeta, nthreads, ncont = 0, ncols;
	const double *cptr = NULL;
	int nprotect = 8;

	m = PROTECT(coerceVector(m, REALSXP));
	design = PROTECT(coerceVector(design, REALSXP));
	if(weight != R_NilValue)
	{
		weight = PROTECT(coerceVector(weight, REALSXP));
		nprotect++;
	}
	if(contrasts != R_NilValue)
	{
		contrasts = PROTECT(coerceVector(contrasts, REALSXP));
		nprotect++;
		ncont = INTEGER(getAttrib(contrasts, R_DimSymbol))[1];
		cptr = REAL(contrasts);
	}
	mdim = getAttrib(m, R_DimSymbol);
	ddim = getAttrib(design, R_DimSymbol);

	ngenes = INTEGER(mdim)[0];
	narrays = INTEGER(mdim)[1];
	nbeta = INTEGER(ddim)[1];
	nthreads = asInteger(nthr);
	ncols = (contrasts != R_NilValue) ? ncont : nbeta;

	output = PROTECT(allocVector(VECSXP, 4));
	beta = PROTECT(allocMatrix(REALSXP, ngenes, ncols));
	stdev = PROTECT(allocMatrix(REALSXP, ngenes, ncols));
	sigma = PROTECT(allocVector(REALSXP, ngenes));
	df = PROTECT(allocVector(INTSXP, ngenes));
	names = PROTECT(allocVector(STRSXP, 4));

	lm(REAL(m), REAL(design), weight == R_NilValue ? NULL : REAL(weight), ngenes, narrays, nbeta, cptr, ncont, nthreads, REAL(beta), REAL(stdev), REAL(sigma), INTEGER(df));

	SET_VECTOR_ELT(output, 0, beta);
	SET_VECTOR_ELT(output, 1, stdev);
	SET_VECTOR_ELT(output, 2, sigma);
	SET_VECTOR_ELT(output, 3, df);
	SET_STRING_ELT(names, 0, mkChar("coefficients"));
	SET_STRING_ELT(names, 1, mkChar("stdev.unscaled"));
	SET_STRING_ELT(names, 2, mkChar("sigma"));
	SET_STRING_ELT(names, 3, mkChar("df.residual"));
	setAttrib(output, R_NamesSymbol, names);

	UNPROTECT(nprotect);
	return output;
}

/*
 * Inputs:
 *   m: numeric ngenes x narrays matrix.
 *   design: numeric narrays x nbeta design matrix.
 *   cormat: numeric narrays x narrays correlation matrix.
 *   weight: optional numeric ngenes x narrays matrix, or R_NilValue.
 *   contrasts: optional numeric nbeta x ncont contrast matrix, or R_NilValue.
 *   nthr: requested thread count.
 * Outputs:
 *   Returns a named list with coefficients, stdev.unscaled, sigma and
 *   df.residual. With contrasts the coefficient/stdev matrices have ncont
 *   columns in contrast space; otherwise nbeta columns.
 * Returns:
 *   SEXP list allocated for R.
 * Notes:
 *   Calls gls(); a nonzero Cholesky status is converted to an R error.
 */
SEXP glsfit(SEXP m, SEXP design, SEXP cormat, SEXP weight, SEXP contrasts, SEXP nthr)
{
	SEXP output, beta, stdev, sigma, df, names;
	SEXP mdim, ddim;
	int ngenes, narrays, nbeta, nthreads, ncont = 0, ncols, nprotect = 9, failed;
	const double *cptr = NULL;

	m = PROTECT(coerceVector(m, REALSXP));
	design = PROTECT(coerceVector(design, REALSXP));
	cormat = PROTECT(coerceVector(cormat, REALSXP));
	if(weight != R_NilValue)
	{
		weight = PROTECT(coerceVector(weight, REALSXP));
		nprotect++;
	}
	if(contrasts != R_NilValue)
	{
		contrasts = PROTECT(coerceVector(contrasts, REALSXP));
		nprotect++;
		ncont = INTEGER(getAttrib(contrasts, R_DimSymbol))[1];
		cptr = REAL(contrasts);
	}
	mdim = getAttrib(m, R_DimSymbol);
	ddim = getAttrib(design, R_DimSymbol);

	ngenes = INTEGER(mdim)[0];
	narrays = INTEGER(mdim)[1];
	nbeta = INTEGER(ddim)[1];
	nthreads = asInteger(nthr);
	ncols = (contrasts != R_NilValue) ? ncont : nbeta;

	output = PROTECT(allocVector(VECSXP, 4));
	beta = PROTECT(allocMatrix(REALSXP, ngenes, ncols));
	stdev = PROTECT(allocMatrix(REALSXP, ngenes, ncols));
	sigma = PROTECT(allocVector(REALSXP, ngenes));
	df = PROTECT(allocVector(INTSXP, ngenes));
	names = PROTECT(allocVector(STRSXP, 4));

	failed = gls(REAL(m), REAL(design), REAL(cormat), weight == R_NilValue ? NULL : REAL(weight), ngenes, narrays, nbeta, cptr, ncont, nthreads, REAL(beta), REAL(stdev), REAL(sigma), INTEGER(df));
	if(failed)
		error("the leading minor of order %d is not positive", failed);

	SET_VECTOR_ELT(output, 0, beta);
	SET_VECTOR_ELT(output, 1, stdev);
	SET_VECTOR_ELT(output, 2, sigma);
	SET_VECTOR_ELT(output, 3, df);
	SET_STRING_ELT(names, 0, mkChar("coefficients"));
	SET_STRING_ELT(names, 1, mkChar("stdev.unscaled"));
	SET_STRING_ELT(names, 2, mkChar("sigma"));
	SET_STRING_ELT(names, 3, mkChar("df.residual"));
	setAttrib(output, R_NamesSymbol, names);

	UNPROTECT(nprotect);
	return output;
}

/*
 * Inputs:
 *   m: numeric ngenes x narrays matrix.
 *   design: numeric narrays x nbeta fixed-effect design matrix.
 *   block: integer length-narrays one-based block labels.
 *   nlev: number of block levels; weight: optional prior weights.
 *   nthr: requested thread count.
 * Outputs:
 *   Returns a length-ngenes numeric vector of per-gene correlations.
 * Returns:
 *   SEXP REALSXP vector allocated for R.
 * Notes:
 *   Consensus trimming and clipping are performed in R. A workspace-query
 *   failure in dupcor() is converted to an R error.
 */
SEXP dupcorfit(SEXP m, SEXP design, SEXP block, SEXP nlev, SEXP weight, SEXP nthr)
{
	SEXP rho, mdim, ddim;
	int ngenes, narrays, nbeta, nblocks, nthreads, nprotect = 3, info;

	m = PROTECT(coerceVector(m, REALSXP));
	design = PROTECT(coerceVector(design, REALSXP));
	if(weight != R_NilValue)
	{
		weight = PROTECT(coerceVector(weight, REALSXP));
		nprotect++;
	}
	mdim = getAttrib(m, R_DimSymbol);
	ddim = getAttrib(design, R_DimSymbol);
	ngenes = INTEGER(mdim)[0];
	narrays = INTEGER(mdim)[1];
	nbeta = INTEGER(ddim)[1];
	nblocks = asInteger(nlev);
	nthreads = asInteger(nthr);

	rho = PROTECT(allocVector(REALSXP, ngenes));
	if(ngenes == 0)
	{
		UNPROTECT(nprotect);
		return rho;
	}

	info = dupcor(REAL(m), REAL(design), INTEGER(block), weight == R_NilValue ? NULL : REAL(weight), ngenes, narrays, nbeta, nblocks, nthreads, REAL(rho));
	if(info)
		error("unable to determine SVD workspace");

	UNPROTECT(nprotect);
	return rho;
}

/*
 * Inputs:
 *   y: numeric ngenes x narrays expression matrix.
 *   design: numeric narrays x p mean-model design matrix.
 *   weight: numeric ngenes x narrays prior-weight matrix.
 *   vd: numeric narrays x ngam variance-design matrix.
 *   priorn, maxit, tol, trace, nthr: REML controls.
 * Outputs:
 *   Returns a named list with w, iter and status.
 * Returns:
 *   SEXP list allocated for R.
 * Notes:
 *   Calls awreml(), which fills the array-weight vector and convergence status.
 */
SEXP awremlfit(SEXP y, SEXP design, SEXP weight, SEXP vd, SEXP priorn, SEXP maxit, SEXP tol, SEXP trace, SEXP nthr)
{
	SEXP out, w, iter, status, names;
	SEXP ydim, ddim, vdim;
	int ngenes, narrays, p, ngam, it = 0, st, nprotect = 9;

	y = PROTECT(coerceVector(y, REALSXP));
	design = PROTECT(coerceVector(design, REALSXP));
	weight = PROTECT(coerceVector(weight, REALSXP));
	vd = PROTECT(coerceVector(vd, REALSXP));

	ydim = getAttrib(y, R_DimSymbol);
	ddim = getAttrib(design, R_DimSymbol);
	vdim = getAttrib(vd, R_DimSymbol);
	ngenes = INTEGER(ydim)[0];
	narrays = INTEGER(ydim)[1];
	p = INTEGER(ddim)[1];
	ngam = INTEGER(vdim)[1];

	w = PROTECT(allocVector(REALSXP, narrays));

	st = awreml(REAL(y), REAL(design), REAL(weight), REAL(vd), ngenes, narrays, p, ngam, asReal(priorn), asInteger(maxit), asReal(tol), asLogical(trace), asInteger(nthr), REAL(w), &it);

	iter = PROTECT(ScalarInteger(it));
	status = PROTECT(ScalarInteger(st));
	out = PROTECT(allocVector(VECSXP, 3));
	names = PROTECT(allocVector(STRSXP, 3));
	SET_VECTOR_ELT(out, 0, w);
	SET_VECTOR_ELT(out, 1, iter);
	SET_VECTOR_ELT(out, 2, status);
	SET_STRING_ELT(names, 0, mkChar("w"));
	SET_STRING_ELT(names, 1, mkChar("iter"));
	SET_STRING_ELT(names, 2, mkChar("status"));
	setAttrib(out, R_NamesSymbol, names);

	UNPROTECT(nprotect);
	return out;
}

/*
 * Inputs:
 *   y: numeric ngenes x narrays count matrix.
 *   design: numeric narrays x p design matrix.
 *   offset: numeric length-narrays offset vector.
 *   start: numeric ngenes x narrays warm-start linear predictor.
 *   nthr: requested thread count.
 * Outputs:
 *   Returns an ngenes x narrays numeric matrix of fitted Poisson means.
 * Returns:
 *   SEXP matrix allocated for R.
 * Notes:
 *   The numerical IRLS work is delegated to pois().
 */
SEXP poisfit(SEXP y, SEXP design, SEXP offset, SEXP start, SEXP nthr)
{
	SEXP out, ydim, ddim;
	int ngenes, narrays, p;

	y = PROTECT(coerceVector(y, REALSXP));
	design = PROTECT(coerceVector(design, REALSXP));
	offset = PROTECT(coerceVector(offset, REALSXP));
	start = PROTECT(coerceVector(start, REALSXP));

	ydim = getAttrib(y, R_DimSymbol);
	ddim = getAttrib(design, R_DimSymbol);
	ngenes = INTEGER(ydim)[0];
	narrays = INTEGER(ydim)[1];
	p = INTEGER(ddim)[1];

	out = PROTECT(allocMatrix(REALSXP, ngenes, narrays));
	pois(REAL(y), REAL(design), REAL(offset), REAL(start), ngenes, narrays, p, asInteger(nthr), REAL(out));

	UNPROTECT(5);
	return out;
}

/*
 * Inputs:
 *   covariate, response, weight: equal-length numeric vectors prepared by R.
 *   span, iter, delta: smoothing controls.
 * Outputs:
 *   Returns a two-component list: fitted values and final robustness weights.
 * Returns:
 *   SEXP list allocated for R.
 * Notes:
 *   weightedLowess() validates sorting, lengths and types before calling this
 *   shim; lowess() performs the numerical smoothing.
 */
SEXP weighted_lowess(SEXP covariate, SEXP response, SEXP weight, SEXP span, SEXP iter, SEXP delta)
{
	int npts = length(covariate);
	SEXP out = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(out, 0, allocVector(REALSXP, npts));
	SET_VECTOR_ELT(out, 1, allocVector(REALSXP, npts));

	lowess(REAL(covariate), REAL(response), REAL(weight), npts, asReal(span), asInteger(iter), asReal(delta), REAL(VECTOR_ELT(out, 0)), REAL(VECTOR_ELT(out, 1)));

	UNPROTECT(1);
	return out;
}
