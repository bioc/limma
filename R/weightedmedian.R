weighted.median <- function (x, w, na.rm = FALSE)
#	Weighted median
#	Gordon Smyth
#	Created 30 June 2005. Last revised 9 June 2020.
{
#	Check weights have correct length
	if (missing(w)) {
		w <- rep_len(1, length.out=length(x))
	} else {
		if(!identical(length(w),length(x))) stop("'x' and 'w' must have the same length")
	}

#	Remove NAs if necessary
	if(na.rm) {
		i <- !is.na(x)
		x <- x[i]
		w <- w[i]
	}

#	Check content of weights
	r <- range(w)
	if(anyNA(r)) stop("NA weights not allowed")
	if(r[1]<0) stop("Negative weights not allowed")
	if(r[2]==0) {
		warning("All weights are zero")
		return(NA_real_)
	}

#	Remove zero weights
	if(r[1]==0) {
		i <- which(w==0)
		x <- x[-i]
		w <- w[-i]
	}

#	Return the median of the discrete distribution with weights as probabilities
	o <- order(x)
	x <- x[o]
	w <- w[o]
	p <- cumsum(w)/sum(w)
	n <- sum(p<0.5)
	if(p[n+1L] > 0.5)
		x[n+1L]
	else
		(x[n+1L]+x[n+2L])/2
}
