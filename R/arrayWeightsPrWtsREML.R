.arrayWeightsPrWtsREML <- function(y, design=NULL, weights, var.design, prior.n=10, maxiter=50L, tol=1e-6, trace=FALSE, nthreads=1L)
#	Estimate array weights by REML allowing for prior observation weights
#	Probes with missing or infinite values are removed.
#	Gordon Smyth
#	Created 12 Feb 2019 from .arrayWeightsREML.
#	C implement and parallel support by Lizhong Chen
#	Last revised 18 June 2026.
{
#	y should be a numeric matrix
	y <- as.matrix(y)

#	Columns of var.design should sum to zero, and intercept column should be omitted.
	Z2 <- as.matrix(var.design)

	fit <- .Call("awremlfit",y,design,weights,Z2,as.double(prior.n),as.integer(maxiter),as.double(tol),as.logical(trace),nthreads,PACKAGE="limma")

	if(fit$status==1L) warning("convergence tolerance not achievable, stopping prematurely")
	if(fit$status==2L) warning("iteration limit reached")

	fit$w
}
