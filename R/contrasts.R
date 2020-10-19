#  CONTRASTS

contrasts.fit <- function(fit,contrasts=NULL,coefficients=NULL)
#	Convert coefficients and std deviations in fit object to reflect contrasts of interest
#	Note: does not completely take probe-wise weights into account
#	because this would require refitting the linear model for each probe
#	Gordon Smyth
#	Created 13 Oct 2002.  Last modified 19 Oct 2020.
{
#	Check number of arguments
	if(identical(is.null(contrasts),is.null(coefficients))) stop("Must specify exactly one of contrasts or coefficients")

#	If coefficients are input, just subset
	if(!is.null(coefficients)) return(fit[,coefficients])

#	Check for valid fit object
	if(is.null(fit$coefficients)) stop("fit must contain coefficients component")
	if(is.null(fit$stdev.unscaled)) stop("fit must contain stdev.unscaled component")

#	Remove test statistics in case eBayes() has previously been run on the fit object
	fit$t <- NULL
	fit$p.value <- NULL
	fit$lods <- NULL
	fit$F <- NULL
	fit$F.p.value <- NULL

#	Number of coefficients in fit
	ncoef <- ncol(fit$coefficients)

#	Check contrasts.
	if(anyNA(contrasts)) stop("NAs not allowed in contrasts")
	contrasts <- as.matrix(contrasts)
	if(!identical(nrow(contrasts),ncoef)) stop("Number of rows of contrast matrix must match number of coefficients in fit")
	rn <- rownames(contrasts)
	cn <- colnames(fit$coefficients)
	if(!is.null(rn) && !is.null(cn) && !identical(rn,cn)) warning("row names of contrasts don't match col names of coefficients")
	fit$contrasts <- contrasts

#	Special case of contrast matrix with 0 columns
	if(!ncol(contrasts)) return(fit[,0])

#	Correlation matrix of estimable coefficients
#	Test whether design was orthogonal
	if(is.null(fit$cov.coefficients)) {
		warning("cov.coefficients not found in fit - assuming coefficients are orthogonal",call.=FALSE)
		var.coef <- colMeans(fit$stdev.unscaled^2)
		fit$cov.coefficients <- diag(var.coef,nrow=ncoef)
		cormatrix <- diag(nrow=ncoef)
		orthog <- TRUE
	} else {
		cormatrix <- cov2cor(fit$cov.coefficients)
		if(length(cormatrix) < 2) {
			orthog <- TRUE
		} else {
			orthog <- sum(abs(cormatrix[lower.tri(cormatrix)])) < 1e-12
		}
	}

#	If design matrix was singular, reduce to estimable coefficients
	r <- nrow(cormatrix)
	if(r < ncoef) {
		if(is.null(fit$pivot)) stop("cor.coef not full rank but pivot column not found in fit")
		est <- fit$pivot[1:r]
		if(any(contrasts[-est,]!=0)) stop("trying to take contrast of non-estimable coefficient")
		contrasts <- contrasts[est,,drop=FALSE]
		fit$coefficients <- fit$coefficients[,est,drop=FALSE]
		fit$stdev.unscaled <- fit$stdev.unscaled[,est,drop=FALSE]
		ncoef <- r
	}

#	Remove coefficients that don't appear in any contrast
#	(Not necessary but can make function faster)
	ContrastsAllZero <- which(rowSums(abs(contrasts))==0)
	if(length(ContrastsAllZero)) {
		contrasts <- contrasts[-ContrastsAllZero,,drop=FALSE]
		fit$coefficients <- fit$coefficients[,-ContrastsAllZero,drop=FALSE]
		fit$stdev.unscaled <- fit$stdev.unscaled[,-ContrastsAllZero,drop=FALSE]
		fit$cov.coefficients <- fit$cov.coefficients[-ContrastsAllZero,-ContrastsAllZero,drop=FALSE]
		cormatrix <- cormatrix[-ContrastsAllZero,-ContrastsAllZero,drop=FALSE]
		ncoef <- ncol(fit$coefficients)
	}

#	Replace NA coefficients with large (but finite) standard deviations
#	to allow zero contrast entries to clobber NA coefficients.
	NACoef <- anyNA(fit$coefficients)
	if(NACoef) {
		i <- is.na(fit$coefficients)
		fit$coefficients[i] <- 0
		fit$stdev.unscaled[i] <- 1e30
	}

#	New coefficients
	fit$coefficients <- fit$coefficients %*% contrasts

#	Test whether design was orthogonal
	if(length(cormatrix) < 2) {
		orthog <- TRUE
	} else {
		orthog <- all(abs(cormatrix[lower.tri(cormatrix)]) < 1e-14)
	}

#	New correlation matrix
	R <- chol(fit$cov.coefficients)
	fit$cov.coefficients <- crossprod(R %*% contrasts)
#	fit$pivot <- NULL

#	New standard deviations
	if(orthog)
		fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
	else {
		R <- chol(cormatrix)
		ngenes <- NROW(fit$stdev.unscaled)
		ncont <- NCOL(contrasts)
		U <- matrix(1,ngenes,ncont,dimnames=list(rownames(fit$stdev.unscaled),colnames(contrasts)))
		o <- array(1,c(1,ncoef))
		for (i in 1:ngenes) {
			RUC <- R %*% .vecmat(fit$stdev.unscaled[i,],contrasts)
			U[i,] <- sqrt(o %*% RUC^2)
		}
		fit$stdev.unscaled <- U
	}

#	Replace NAs if necessary
	if(NACoef) {
		i <- (fit$stdev.unscaled > 1e20)
		fit$coefficients[i] <- NA
		fit$stdev.unscaled[i] <- NA
	}

	fit
}
