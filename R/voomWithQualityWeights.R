voomWithQualityWeights <- function(counts, design=NULL, lib.size=NULL, normalize.method="none", plot=FALSE, span=0.5, var.design=NULL, var.group=NULL, method="genebygene", maxiter=50, tol=1e-5, trace=FALSE, col=NULL, ...)
#	Combine voom weights with sample-specific weights estimated by arrayWeights() function for RNA-seq data
#	Matt Ritchie, Cynthia Liu, Gordon Smyth
#	Created 22 Sept 2014. Last modified 7 June 2019.
{
#	Setup side-by-side plots showing (1) the voom trend and (2) the array weights
	if(plot) {
		oldpar <- par(mfrow=c(1,2))
		on.exit(par(oldpar))
	}

#	Voom without array weights
	v <- voom(counts, design=design, lib.size=lib.size, normalize.method=normalize.method, plot=FALSE, span=span, ...)

#	Estimate array weights on top of voom weights
	aw <- arrayWeights(v, design=design, method=method, maxiter=maxiter, tol=tol, var.design=var.design, var.group=var.group)

#	Update voom weights now using the array weights, plotting trend if requested
	v <- voom(counts, design=design, weights=aw, lib.size=lib.size, normalize.method=normalize.method, plot=plot, span=span, ...)

#	Update array weights
	aw <- arrayWeights(v, design=design, method=method, maxiter=maxiter, tol=tol, trace=trace, var.design=var.design, var.group=var.group)

#	Incorporate the array weights into the voom weights
	v$weights <- t(aw * t(v$weights))
	v$targets$sample.weights <- aw

#	Plot array weights
	if(plot) {
		barplot(aw, names=1:length(aw), main="Sample-specific weights", ylab="Weight", xlab="Sample", col=col)
		abline(h=1, col=2, lty=2)
	}

	v
}

