/*
 * weighted_lowess.c -- C version of the local weighted regression (lowess)
 * trend fitting algorithm, based on the Fortran lowess.f from
 * http://www.netlib.org/go written by Cleveland. Consideration of non-equal
 * prior weights is added to the span calculations and linear regression. These
 * weights are intended to have the equivalent effect of frequency weights (at
 * least, in the integer case; extended by analogy to all non-negative values).
 *
 * The kernel fills two caller-provided npts-length outputs: the fitted trend
 * and the final robustness weights. All scratch lives in one lowessws workspace
 * (allocated once, freed once). The SEXP shim is in R_exports.c.
 */

#include <math.h>
#include <R.h>
#include "limma.h"

#define THRESHOLD 0.0000001

static void alloclowess(lowessws *ws, int npts)
{
	ws->seed = R_Calloc(npts, int);
	ws->fstart = R_Calloc(npts, int);
	ws->fend = R_Calloc(npts, int);
	ws->fdist = R_Calloc(npts, double);
	ws->work = R_Calloc(npts, double);
	ws->ror = R_Calloc(npts, int);
}

static void freelowess(lowessws *ws)
{
	R_Free(ws->seed);
	R_Free(ws->fstart);
	R_Free(ws->fend);
	R_Free(ws->fdist);
	R_Free(ws->work);
	R_Free(ws->ror);
}

/*
 * Inputs:
 *   xptr: sorted length-npts covariate vector.
 *   delta: minimum covariate spacing between seed points.
 * Outputs:
 *   ws->seed receives the selected seed indices.
 * Returns:
 *   Number of seed points written to ws->seed.
 * Notes:
 *   The first and last observations are always seeds.
 */
static int find_seeds(lowessws *ws, const double *xptr, int npts, double delta)
{
	int pt, last_pt = 0;
	int total = 1;
	ws->seed[0] = 0;
	for(pt = 1; pt < npts - 1; ++pt)
	{
		if(xptr[pt] - xptr[last_pt] > delta)
		{
			ws->seed[total] = pt;
			last_pt = pt;
			++total;
		}
	}
	ws->seed[total] = npts - 1;
	++total;
	return total;
}

/*
 * Inputs:
 *   ws->seed: num seed indices.
 *   xptr, wptr: sorted covariates and non-negative prior weights.
 *   spanweight: target total prior weight in each local span.
 * Outputs:
 *   ws->fstart/ws->fend receive inclusive span limits for each seed.
 *   ws->fdist receives the maximum covariate distance for each local span.
 * Returns:
 *   Nothing.
 * Notes:
 *   Ties at span boundaries are included. The implementation recomputes spans
 *   directly instead of using Cleveland's update recurrence because weighted
 *   updates are less numerically stable in double precision.
 *
 */
static void find_limits(lowessws *ws, int num, const double *xptr, const double *wptr, int npts, double spanweight)
{
	int curx;
	for(curx = 0; curx < num; ++curx)
	{
		const int curpt = ws->seed[curx];
		int left = curpt, right = curpt;
		double curw = wptr[curpt];
		int ende = (curpt == npts - 1), ends = (curpt == 0);
		double mdist = 0, ldist, rdist;

		while(curw < spanweight && (!ende || !ends))
		{
			if(ende)
			{
				/* Can only extend backwards. */
				--left;
				curw += wptr[left];
				if(left == 0)
					ends = 1;
				ldist = xptr[curpt] - xptr[left];
				if(mdist < ldist)
					mdist = ldist;
			}
			else if(ends)
			{
				/* Can only extend forwards. */
				++right;
				curw += wptr[right];
				if(right == npts - 1)
					ende = 1;
				rdist = xptr[right] - xptr[curpt];
				if(mdist < rdist)
					mdist = rdist;
			}
			else
			{
				/* Can do either; extending by the one that minimizes the curpt mdist. */
				ldist = xptr[curpt] - xptr[left - 1];
				rdist = xptr[right + 1] - xptr[curpt];
				if(ldist < rdist)
				{
					--left;
					curw += wptr[left];
					if(left == 0)
						ends = 1;
					if(mdist < ldist)
						mdist = ldist;
				}
				else
				{
					++right;
					curw += wptr[right];
					if(right == npts - 1)
						ende = 1;
					if(mdist < rdist)
						mdist = rdist;
				}
			}
		}

		/* Extending to ties. */
		while(left > 0 && xptr[left] == xptr[left - 1])
			--left;
		while(right < npts - 1 && xptr[right] == xptr[right + 1])
			++right;

		/* Recording */
		ws->fstart[curx] = left;
		ws->fend[curx] = right;
		ws->fdist[curx] = mdist;
	}
}

/*
 * Inputs:
 *   xptr, yptr, wptr, rwptr: covariate, response, prior and robust weights.
 *   curpt: point to fit; left/right: inclusive local span; dist: span distance.
 * Outputs:
 *   work[left..right] receives combined local regression weights.
 * Returns:
 *   Fitted value at curpt.
 * Notes:
 *   Uses tricube distance weights times prior and robustness weights. If the
 *   local x variance is too small, returns the weighted local mean.
 */
static double lowess_fit(const double *xptr, const double *yptr, const double *wptr, const double *rwptr, int curpt, int left, int right, double dist, double *work)
{
	double ymean = 0, allweight = 0;
	int pt;
	if(dist < THRESHOLD)
	{
		for(pt = left; pt <= right; ++pt)
		{
			work[pt] = wptr[pt] * rwptr[pt];
			ymean += yptr[pt] * work[pt];
			allweight += work[pt];
		}
		ymean /= allweight;
		return ymean;
	}
	double xmean = 0;
	for(pt = left; pt <= right; ++pt)
	{
		work[pt] = pow(1 - pow(fabs(xptr[curpt] - xptr[pt]) / dist, 3.0), 3.0) * wptr[pt] * rwptr[pt];
		xmean += work[pt] * xptr[pt];
		ymean += work[pt] * yptr[pt];
		allweight += work[pt];
	}
	xmean /= allweight;
	ymean /= allweight;

	double var = 0, covar = 0, temp;
	for(pt = left; pt <= right; ++pt)
	{
		temp = xptr[pt] - xmean;
		var += temp * temp * work[pt];
		covar += temp * (yptr[pt] - ymean) * work[pt];
	}
	if(var < THRESHOLD)
		return ymean;

	const double slope = covar / var;
	const double intercept = ymean - slope * xmean;
	return slope * xptr[curpt] + intercept;
}

/*
 * Inputs:
 *   x, y, w: sorted covariate, response and non-negative prior-weight vectors.
 *   span: target fraction of total prior weight per local fit.
 *   niter: number of robustness iterations; delta: seed spacing threshold.
 * Outputs:
 *   fitted: length-npts fitted trend values.
 *   robust: length-npts final robustness weights.
 * Returns:
 *   Nothing.
 * Notes:
 *   Allocates one lowessws workspace internally. Seed fits are linearly
 *   interpolated between seed points; robustness weights use a weighted MAD.
 */
void lowess(const double *x, const double *y, const double *w, int npts, double span, int niter, double delta, double *fitted, double *robust)
{
	lowessws ws;
	int pt, it, nseeds;
	double totalweight = 0;

	alloclowess(&ws, npts);

	/* The lowess span is interpreted as a fraction of total prior weight, not a
	 * raw count of observations. Each local window is expanded until its prior
	 * weights sum to at least spanweight. */
	for(pt = 0; pt < npts; ++pt)
		totalweight += w[pt];
	double spanweight = totalweight * span;
	const double subrange = (x[npts - 1] - x[0]) / npts;

	/* delta selects a subset of seed covariate points where full local
	 * regressions are evaluated. find_limits records the inclusive window for
	 * each seed and the maximum covariate distance used by the tricube kernel. */
	nseeds = find_seeds(&ws, x, npts, delta);
	find_limits(&ws, nseeds, x, w, npts, spanweight);

	for(pt = 0; pt < npts; ++pt)
		robust[pt] = 1;

	/* Each robustness iteration refits the seed points using the current robust
	 * weights, interpolates fitted values between seeds, then recomputes robust
	 * weights from absolute residuals. */
	for(it = 0; it < niter; ++it)
	{
		int cur_seed, last_pt = 0, subpt;
		double current;

		/* Compute fitted values at seed points. Non-seed points between adjacent
		 * seeds are filled by linear interpolation in x, following the original
		 * lowess algorithm's delta shortcut. */
		fitted[0] = lowess_fit(x, y, w, robust, 0, ws.fstart[0], ws.fend[0], ws.fdist[0], ws.work);
		for(cur_seed = 1; cur_seed < nseeds; ++cur_seed)
		{
			pt = ws.seed[cur_seed];
			fitted[pt] = lowess_fit(x, y, w, robust, pt, ws.fstart[cur_seed], ws.fend[cur_seed], ws.fdist[cur_seed], ws.work);

			if(pt - last_pt > 1)
			{
				/* Some protection is provided against infinite slopes. This shouldn't be
				 * a problem for non-zero delta; the only concern is at the final point
				 * where the covariate distance may be zero. Besides, if delta is not
				 * positive, pt-last_pt could never be 1 so we'd never reach this point.
				 */
				current = x[pt] - x[last_pt];
				if(current > THRESHOLD * subrange)
				{
					const double slope = (fitted[pt] - fitted[last_pt]) / current;
					const double intercept = fitted[pt] - slope * x[pt];
					for(subpt = last_pt + 1; subpt < pt; ++subpt)
						fitted[subpt] = slope * x[subpt] + intercept;
				}
				else
				{
					const double endave = 0.5 * (fitted[pt] + fitted[last_pt]);
					for(subpt = last_pt + 1; subpt < pt; ++subpt)
						fitted[subpt] = endave;
				}
			}
			last_pt = pt;
		}

		/* Compute a weighted median absolute residual. rsort_with_index sorts
		 * residual magnitudes in ws.work and carries their original indices in
		 * ws.ror, so prior weights can be accumulated in sorted-residual order. */
		double resid_scale = 0;
		for(pt = 0; pt < npts; ++pt)
		{
			ws.work[pt] = fabs(y[pt] - fitted[pt]);
			resid_scale += ws.work[pt];
			ws.ror[pt] = pt;
		}
		resid_scale /= npts;
		rsort_with_index(ws.work, ws.ror, npts);

		current = 0;
		double cmad = 0;
		const double halfweight = totalweight / 2;
		for(pt = 0; pt < npts; ++pt)
		{
			current += w[ws.ror[pt]];
			if(current == halfweight)
			{
				/* Exact half-weight matches are rare; average the adjacent
				 * residual magnitudes, then multiply by the bisquare constant. */
				cmad = 3 * (ws.work[pt] + ws.work[pt + 1]);
				break;
			}
			else if(current > halfweight)
			{
				cmad = 6 * ws.work[pt];
				break;
			}
		}

		/* If it's too small, then robustness weighting will have no further effect.
		 * Any points with large residuals would already be pretty lowly weighted.
		 * This is based on a similar step in lowess.c in the core R code.
		 */
		if(cmad <= THRESHOLD * resid_scale)
			break;

		/* Tukey bisquare robustness weights. Points with residual magnitude at or
		 * above cmad receive zero robust weight in the next iteration. */
		for(pt = 0; pt < npts; ++pt)
		{
			if(ws.work[pt] < cmad)
				robust[ws.ror[pt]] = pow(1 - pow(ws.work[pt] / cmad, 2.0), 2.0);
			else
				robust[ws.ror[pt]] = 0;
		}
	}

	freelowess(&ws);
}
