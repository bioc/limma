weightedLowess <- function(x, y, weights=rep(1, length(y)), delta=NULL, npts=200, span=0.3, iterations=4L) 
# This function clusters points by average linkage and fits a lowess curve of degree 1 to the cluster
# midpoints. Fitted values are computed by linear interpolation of the fitted coefficients (i.e. 
# quadratic interpolation between points. Several iterations of robustification are performed
# using the fitted residuals.
#
# Created by Aaron Lun 13 Jan 2014.
# Last modified 16 Feb 2020.
{
	x <- as.double(x)
	y <- as.double(y)
	weights <- as.double(weights)
	if(!identical(length(y),length(x))) stop("x and y should have same length")
	if(!identical(length(y),length(weights))) stop("weights should have same length as x and y")
	iterations <- as.integer(iterations)

#	Choosing an appropriate 'delta' for approximation. We assume that the covariates
#	have some cluster structure, where clusters are defined by partitioning on the 
#	'numclusters' largest distances between points. We want 'npts' evenly spaced points 
#	across the cluster range (i.e. the total range minus the partitioned distances).
#	Each cluster must also have at least one point; so, we compute the spacing (and
#	thus delta) as the ratio of the cluster range to the remaining number of points
#	(after each additional cluster beyond the required first one has eaten 1 point). 
	o <- order(x)
	x <- x[o]
	if (is.null(delta)) {
		npts <- as.integer(npts+0.5)
		if (npts < 1L) { stop("number of points should be a positive integer") }
		if (npts >= length(x)) { 
			delta <- 0
		} else {
			dx <- sort(diff(x))
			cumrange <- cumsum(dx)
			numclusters <- 1:npts-1L
			delta <- min(cumrange[length(dx)-numclusters]/(npts-numclusters))
		}
	}
	delta <- as.double(delta)
	
#	Running the smoothing procedure with specified values.
	out <- .Call("weighted_lowess", x, y[o], weights[o], span, iterations, delta, PACKAGE="limma")
	names(out) <- c("fitted", "weights")
	out$fitted[o] <- out$fitted
	out$residuals <- y-out$fitted
	out$weights[o] <- out$weights
	out$delta <- delta

	out
}

