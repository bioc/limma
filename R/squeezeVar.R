#	EMPIRICAL BAYES SQUEEZING OF VARIANCES

squeezeVar <- function(var, df, covariate=NULL, robust=FALSE, winsor.tail.p=c(0.05,0.1))
#	Empirical Bayes posterior variances
#	Gordon Smyth
#	2 March 2004.  Last modified 2 Dec 2013.
{
	n <- length(var)
	if(n == 0) stop("var is empty")
	if(n == 1) return(list(var.post=var,var.prior=var,df.prior=0))
	if(length(df)==1) { 
		df <- rep.int(df,n)
	} else {
		if(length(df) != n) stop("lengths differ")
	}

#	Estimate prior var and df
	if(robust)
		fit <- fitFDistRobustly(var, df1=df, covariate=covariate, winsor.tail.p=winsor.tail.p)
	else
		fit <- fitFDist(var, df1=df, covariate=covariate)

#	Prior var will be vector if robust=TRUE, otherwise scalar
 	var.prior <- fit$scale

#	Prior df will be vector if covariate is non-NULL, otherwise scalar
	df.prior <- fit$df2.shrunk
	if(is.null(df.prior)) df.prior <- fit$df2

#	Check estimated prior df
	if(is.null(df.prior) || any(is.na(df.prior))) stop("Could not estimate prior df")

#	Squeeze the posterior variances
	df.total <- df + df.prior
	var[df==0] <- 0 # guard against missing or infinite values
	Infdf <- df.prior==Inf
	if(any(Infdf)) {
		var.post <- rep(var.prior,length.out=n)
		i <- which(!Infdf)
		if(length(i)) {
			if(is.null(covariate))
				s02 <- var.prior
			else
				s02 <- var.prior[i]
			var.post[i] <- (df[i]*var[i] + df.prior[i]*s02) / df.total[i]
		}
	} else {
		var.post <- (df*var + df.prior*var.prior) / df.total
	}

	list(df.prior=df.prior,var.prior=var.prior,var.post=var.post)
}
