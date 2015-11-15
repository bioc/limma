fitFDist <- function(x,df1,covariate=NULL)
#	Moment estimation of the parameters of a scaled F-distribution.
#	The numerator degrees of freedom are given, the denominator is to be estimated.
#	Gordon Smyth and Belinda Phipson
#	8 Sept 2002.  Last revised 27 Oct 2012.
{
#	Check covariate
	if(!is.null(covariate)) {
		if(length(covariate) != length(x)) stop("covariate and x must be of same length")
		if(any(is.na(covariate))) stop("NA covariate values not allowed")
		isfin <- is.finite(covariate)
		if(!all(isfin)) {
			if(!any(isfin))
				covariate <- sign(covariate)
			else {
				r <- range(covariate[isfin])
				covariate[covariate == -Inf] <- r[1]-1
				covariate[covariate == Inf] <- r[2]+1
			}
		}
		splinedf <- min(4,length(unique(covariate)))
		if(splinedf < 2) covariate <- NULL
	}
#	Remove missing or infinite values and zero degrees of freedom
	ok <- is.finite(x) & is.finite(df1) & (x > -1e-15) & (df1 > 1e-15)
	notallok <- !all(ok)
	if(notallok) {
		x <- x[ok]
		df1 <- df1[ok]
		if(!is.null(covariate)) {
			covariate2 <- covariate[!ok]
			covariate <- covariate[ok]
		}
	}
	n <- length(x)
	if(n==0) return(list(scale=NA,df2=NA))

#	Avoid exactly zero values
	x <- pmax(x,0)
	m <- median(x)
	if(m==0) {
		warning("More than half of residual variances are exactly zero: eBayes unreliable")
		m <- 1
	} else {
		if(any(x==0)) warning("Zero sample variances detected, have been offset",call.=FALSE)
	}
	x <- pmax(x, 1e-5 * m)

#	Better to work on with log(F)
	z <- log(x)
	e <- z-digamma(df1/2)+log(df1/2)

	if(is.null(covariate)) {
		emean <- mean(e)
		evar <- sum((e-emean)^2)/(n-1)
	} else {
		if(!requireNamespace("splines",quietly=TRUE)) stop("splines package required but is not available")
		design <- try(splines::ns(covariate,df=splinedf,intercept=TRUE),silent=TRUE)
		if(is(design,"try-error")) stop("Problem with covariate")
		fit <- lm.fit(design,e)
		if(notallok) {
			design2 <- predict(design,newx=covariate2)
			emean <- rep.int(0,n+length(covariate2))
			emean[ok] <- fit$fitted
			emean[!ok] <- design2 %*% fit$coefficients
		} else {
			emean <- fit$fitted
		}
		evar <- mean(fit$residuals[-(1:fit$rank)]^2)
	}
	evar <- evar - mean(trigamma(df1/2))
	if(evar > 0) {
		df2 <- 2*trigammaInverse(evar)
		s20 <- exp(emean+digamma(df2/2)-log(df2/2))
	} else {
		df2 <- Inf
		s20 <- exp(emean)
	}
	list(scale=s20,df2=df2)
}


trigammaInverse <- function(x) {
#	Solve trigamma(y) = x for y
#	Gordon Smyth
#	8 Sept 2002.  Last revised 12 March 2004.

#	Non-numeric or zero length input
	if(!is.numeric(x)) stop("Non-numeric argument to mathematical function")
	if(length(x)==0) return(numeric(0))

#	Treat out-of-range values as special cases
	omit <- is.na(x)
	if(any(omit)) {
		y <- x
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 0)
	if(any(omit)) {
		y <- x
		y[omit] <- NaN
		warning("NaNs produced")
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x > 1e7)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/sqrt(x[omit])
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 1e-6)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/x[omit]
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}

#	Newton's method
#	1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
#	so iteration to solve 1/x = 1/trigamma is monotonically convergent
	y <- 0.5+1/x
	iter <- 0
	repeat {
		iter <- iter+1
		tri <- trigamma(y)
		dif <- tri*(1-tri/x)/psigamma(y,deriv=2)
		y <- y+dif
		if(max(-dif/y) < 1e-8) break
		if(iter > 50) {
			warning("Iteration limit exceeded")
			break
		}
	}
	y
}
