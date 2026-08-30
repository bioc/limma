normalizeBetweenSamples <- function(object, method="quantile", cyclic.method="fast", ...)
#	Normalize between samples.
#	Same as normalizeBetweenArrays() but for EList objects.
#	Gordon Smyth
#	Created 30 May 2026.  Last revised 30 May 2026.
{
#	Check for data.frame input
	if(is.data.frame(object)) {
		object <- as.matrix(object)
		if(mode(object) != "numeric") stop("'object' is a data.frame and not all columns are numeric")
	}

#	Check method
	choices <- c("none","scale","quantile","cyclicloess")
	method <- match.arg(method,choices)

#	Methods for matrices
	if(is(object,"matrix")) {
		return(switch(method,
			none = object,
			scale = normalizeMedianValues(object),
			quantile = normalizeQuantiles(object, ...),
			cyclicloess = normalizeCyclicLoess(object,method=cyclic.method, ...)
		))
	}

#	EList objects
	if(is(object,"EList")) {
		object$E <- Recall(object$E,method=method,cyclic.method=cyclic.method,...)
		return(object)
	} else {
		stop("class of object not supported")
	}
}
