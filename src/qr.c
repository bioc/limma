/*
 * qr.c -- pivoted (weighted) least-squares solve for a single gene.
 *
 * Shared kernel used by both lm.c and gls.c. The caller has
 * already placed the (optionally prewhitened) response in lws->y and the
 * design in lws->x; this routine solves the ordinary least-squares problem
 *
 *     beta_hat = argmin_beta || y - X beta ||^2
 *
 * by a rank-revealing pivoted QR factorisation X = Q R, using the LINPACK
 * routine dqrls -- the very routine base R's lm.fit()/lm.wfit() dispatch to, so
 * the numerical results agree to rounding error. From X = Q R we have
 *
 *     beta_hat        = R^{-1} Q^T y
 *     Var(beta_hat)   = sigma^2 (X^T X)^{-1} = sigma^2 R^{-1} R^{-T}
 *     stdev.unscaled  = sqrt(diag((X^T X)^{-1})) = sqrt(rowSumSq(R^{-1}))
 *     sigma^2         = RSS / (nobs - rank),  RSS = sum_{i>=rank} (Q^T y)_i^2
 *
 * dqrls returns the rotated response effec = Q^T y, so the residual sum of
 * squares is just the tail of effec (no explicit residual vector needed). This
 * reproduces, per gene, the lm.fit()/lm.wfit() call inside the lm.series and
 * gls.series reference loops, and chol2inv(qr$qr) for the unscaled covariance.
 *
 * Rank deficiency: dqrls calls dqrdc2, which keeps the linearly-independent
 * columns in their original order and pushes collinear columns to the end.
 * Hence pivot[0..rank-1] is ascending and lists the retained columns; the
 * dropped columns are left at their caller-initialised NA (mirroring lm.fit's
 * aliased-coefficient NAs).
 */

#include <math.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/RS.h>
#include "limma.h"

/*
 * Inputs:
 *   r: QR storage whose leading rank x rank upper triangle is R.
 *   ldr: leading dimension of r; rank: numerical rank to invert.
 * Outputs:
 *   rinv: compact rank x rank column-major inverse of R.
 * Notes:
 *   Uses column-wise back-substitution to solve R rinv = I. The row sums of
 *   squares of rinv give diag((X'X)^-1) for estimable coefficients.
 */
static void invup(const double *r, int ldr, int rank, double *rinv)
{
	int i, j, k;
	for(j = 0; j < rank; j++)
	{
		for(i = 0; i < rank; i++)
			rinv[i + j * rank] = 0.0;
		/* back-substitute up column j: rinv[i,j] = (delta_ij - sum_{k>i} R[i,k] rinv[k,j]) / R[i,i] */
		for(i = j; i >= 0; i--)
		{
			double value = (i == j) ? 1.0 : 0.0;
			for(k = i + 1; k <= j; k++)
				value -= r[i + k * ldr] * rinv[k + j * rank];
			rinv[i + j * rank] = value / r[i + i * ldr];
		}
	}
}

/*
 * Inputs:
 *   nobs: number of compacted observations in lws->x and lws->y.
 *   nbeta: number of design columns; gene/ngenes identify the output row.
 *   contrasts: optional nbeta x ncont contrast matrix (column-major), or NULL.
 *   ncont: number of contrast columns (ignored when contrasts is NULL).
 *   lws->x: nobs x nbeta design matrix, overwritten by dqrls QR storage.
 *   lws->y: length-nobs response vector.
 * Outputs:
 *   beta and stdev_unscaled are filled for this gene. Without contrasts these
 *   hold the nbeta estimable coefficients and sqrt(diag((X'X)^-1)). With
 *   contrasts they hold the ncont contrast values C'beta and the exact contrast
 *   standard deviations sqrt(diag(C'(X'X)^-1 C)) for this gene.
 *   sigma and df_residual are filled for this gene when residual df is positive.
 * Returns:
 *   0 for all paths.
 * Notes:
 *   Caller pre-initialises output rows to NA/0. Columns/contrasts that are not
 *   estimable for this gene remain NA. The routine uses R's dqrls tolerance
 *   1e-7, matching lm.fit/lm.wfit.
 */
int qrfit(int nobs, int nbeta, int gene, int ngenes, const double *contrasts, int ncont, double *beta, double *stdev_unscaled, double *sigma, int *df_residual, lmws *lws)
{
	int i, j, next = 0, rank = 0, ny = 1;
	double tol = 1e-7;
	/* Outputs are column-major matrices with genes as rows. These pointers walk
	 * the current gene across coefficient (or contrast) columns by jumping ngenes
	 * entries. */
	double *bptr = beta + gene;
	double *sptr = stdev_unscaled + gene;

	if(nbeta == 0)
		return 0;

	/* dqrls (R LINPACK-based) solves real least-squares problems by
	 * rank-revealing QR decomposition with column pivoting. limma uses it here
	 * to match the numerical path used by lm.fit() for one gene:
	 *   - x (nobs x nbeta, lda=nobs) enters as the design matrix and is
	 *     overwritten by compact pivoted-QR storage;
	 *   - y enters as the response, with ny=1 right-hand side;
	 *   - tol controls numerical rank detection;
	 *   - coeff receives coefficients in pivot order;
	 *   - resid and effec receive residual and Q^T y working outputs;
	 *   - rank receives the detected rank;
	 *   - pivot is 1-based and maps pivot-order columns back to original columns;
	 *   - qraux/work store Householder and temporary QR data.
	 * The overwritten x and effec outputs provide the triangular R factor and
	 * Q^T y used below for stdev.unscaled, sigma and residual df.
	 * Equivalent R operation: lm.fit(x, y), or qr.coef(qr(x), y) together with
	 * qr.qty(qr(x), y) for the effects.
	 * Netlib references:
	 *   https://www.netlib.org/linpack/dqrdc.f
	 *   https://www.netlib.org/linpack/dqrsl.f */
	for(j = 0; j < nbeta; j++)
		lws->pivot[j] = j + 1;
	F77_CALL(dqrls)(lws->x, &nobs, &nbeta, lws->y, &ny, &tol, lws->coeff, lws->resid, lws->effec, &rank, lws->pivot, lws->qraux, lws->work);
	if(rank == 0)
		return 0;

	/* Nothing is estimable when rank is zero, so the caller-initialised NA/0 row
	 * is left unchanged. Otherwise invert the estimable triangular factor R from
	 * dqrls' overwritten design matrix. */
	invup(lws->x, nobs, rank, lws->rinv);
	if(contrasts == NULL)
	{
		/* Walk all nbeta original columns; a column is retained iff it is the next
		 * pivot entry (pivot[0..rank-1] is ascending). Retained column `next` maps to
		 * coeff[next] and to column `next` of R^{-1}; dropped columns stay NA. */
		for(j = 0; j < nbeta; j++, bptr += ngenes, sptr += ngenes)
		{
			if(next < rank && lws->pivot[next] - 1 == j)
			{
				double value = 0.0;
				/* row-sum-of-squares of R^{-1} = diagonal of (X^T X)^{-1} */
				for(i = next; i < rank; i++)
				{
					double entry = lws->rinv[next + i * rank];
					value += entry * entry;
				}
				*bptr = lws->coeff[next];
				*sptr = sqrt(value);
				next++;
			}
		}
	}
	else
	{
		/* Contrast space. estpos[c] gives the pivot position of original column c
		 * (0..rank-1) or -1 when that column was dropped as non-estimable for this
		 * gene. cest holds the contrast restricted to the estimable columns, in
		 * pivot order. For contrast k: coef = sum cest*coeff and the unscaled
		 * variance is c'(X'X)^-1 c = ||R^{-T} cest||^2 = sum_c (sum_a cest[a]
		 * R^{-1}[a,c])^2, the exact per-gene contrast s.d. A contrast that loads on
		 * any column dropped for this gene is not estimable and stays NA. */
		for(j = 0; j < nbeta; j++)
			lws->estpos[j] = -1;
		for(i = 0; i < rank; i++)
			lws->estpos[lws->pivot[i] - 1] = i;
		for(j = 0; j < ncont; j++, bptr += ngenes, sptr += ngenes)
		{
			const double *ck = contrasts + (size_t) j * nbeta;
			int estimable = 1, a, c;
			double coefval = 0.0, var = 0.0;
			for(i = 0; i < nbeta; i++)
				if(ck[i] != 0.0 && lws->estpos[i] < 0)
				{
					estimable = 0;
					break;
				}
			if(!estimable)
				continue;                     /* leave NA/NA for this contrast */
			for(i = 0; i < rank; i++)
			{
				lws->cest[i] = ck[lws->pivot[i] - 1];
				coefval += lws->cest[i] * lws->coeff[i];
			}
			for(c = 0; c < rank; c++)
			{
				double w = 0.0;
				for(a = 0; a <= c; a++)       /* R^{-1} is upper triangular */
					w += lws->cest[a] * lws->rinv[a + c * rank];
				var += w * w;
			}
			*bptr = coefval;
			*sptr = sqrt(var);
		}
	}

	/* The residual degrees of freedom are nobs-rank. dqrls has already rotated
	 * the response to Q^T y in effec; the fitted part is in the first rank
	 * entries, so RSS is the squared tail beyond rank. */
	df_residual[gene] = nobs - rank;
	if(df_residual[gene] > 0)
	{
		double rss = 0.0;
		for(i = rank; i < nobs; i++)
			rss += lws->effec[i] * lws->effec[i];
		sigma[gene] = sqrt(rss / df_residual[gene]);
	}
	return 0;
}
