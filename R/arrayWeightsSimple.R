arrayWeightsQuick <- function(y, fit)
#	Compute approximate array quality weights
#	Gordon Smyth
#	25 Oct 2004.  Last revised 28 Oct 2004.
{
	if(!is.null(fit$weights)) warning("spot quality weights found but not taken into account")
	res <- as.matrix(y)- fit$coef %*% t(fit$design)
	h <- hat(fit$design, intercept=FALSE)
	mures2 <- fit$sigma^2 %*% array(1-h,c(1,length(h)))
	1/colMeans(res*res/mures2,na.rm=TRUE)
}

arrayWeightsSimple <- function(object,design=NULL,var.group=NULL,maxiter=100L,tol=1e-6,maxratio=100L,trace=FALSE)
#	Estimate array weights by REML.
#	Assumes no prior weights and any probes with missing or infinite values are removed.
#	Uses an exact Fisher scoring algorithm similar to statmod::remlscor.
#	Gordon Smyth
#	Created 13 Dec 2005. Last revised 31 Jan 2019.
{
	M <- as.matrix(object)
	narrays <- ncol(M)
	if(narrays < 3L) stop("too few arrays")

#	Remove rows with missing or infinite values
	allfin <- is.finite(rowSums(M))
	if(!all(allfin)) {
		nrowna <- sum(!allfin)
		M <- M[allfin,,drop=FALSE]
		message(as.integer(nrowna)," rows with missing or inf values removed")
	}
	ngenes <- nrow(M)
	if(ngenes < narrays) stop("too few probes")

#	Default design is just an intercept
	if(is.null(design)) design <- matrix(1,narrays,1)
	p <- ncol(design)

#	Setup genewise variance design matrices, with and without intercept
	if(is.null(var.group)) {
		Z1 <- contr.sum(narrays)
		Z <- cbind(1,Z1)
	} else {
		var.group <- droplevels(as.factor(var.group))
		if(length(var.group) != narrays) stop("var.group has wrong length")
		contrasts(var.group) <- contr.sum(levels(var.group))
		Z <- model.matrix(~var.group)
		Z1 <- Z[,-1,drop=FALSE]
	}

#	Ratio of genes to arrays, used for convergence criterion
#	The sqrt is used to reduce iterations with very large datasets
	m <- sqrt(ngenes)/narrays

#	Starting values
	gam <- rep_len(0,ncol(Z1))
	w <- rep_len(1,narrays)
	if(trace) cat("iter convcrit range(w)\n")

	iter <- 0L
	p2 <- p * (p+1L) %/% 2L
	Q2 <- array(0,c(narrays,p2))
	repeat {
		iter <- iter+1L

#		Fit weighted linear models and extract residual variances
		fitm <- lm.wfit(design, t(M), w)
		Effects <- fitm$effects[(fitm$rank+1):narrays,,drop=FALSE]
		s2 <- colMeans(Effects^2)

#		Fisher information matrix for variance parameters (including intercept)
		Q <- qr.qy(fitm$qr,diag(1,nrow=narrays,ncol=p))
		j0 <- 0
		for (k in 0:(p-1)) {
			Q2[,(j0+1):(j0+p-k)] <- Q[,1:(p-k)]*Q[,(k + 1):p]
			j0 <- j0 + p - k
		}
		if(p > 1) Q2[,(p+1):p2] <- sqrt(2)*Q2[,(p+1):p2]
		h <- rowSums(Q2[,1:p,drop=FALSE])
		info <- crossprod(Z,(1-2*h)*Z) + crossprod(crossprod(Q2,Z))

#		Fisher information excluding intercept (i.e., for gam)
		info1 <- info[-1,-1,drop=FALSE] - (info[-1,1,drop=FALSE]/info[1,1]) %*% info[1,-1,drop=FALSE]

#		Score vector (log-lik derivative) for gam
		fitvald <- matrix(1/w,narrays,1)%*%s2
		dl1 <- crossprod(Z1, rowMeans(fitm$residuals^2/fitvald - (1-h)) )

#		Fisher scoring
		gamstep <- solve(info1,dl1)
		gam <- gam + gamstep

#		Update array weights
		w <- drop(exp(Z1 %*% (-gam)))

#		Test for convergence
		convcrit <- m*crossprod(dl1,gamstep)
		if(trace) cat(iter,convcrit,range(w),"\n")
		if(is.na(convcrit)) {
			warning("convergence tolerance not achievable, stopping prematurely")
			break
		}
		if(convcrit < tol) break

#		Check for other exceedences
		if(max(w)/min(w) > maxratio) {
			warning("maximum divergence of weights reached, stopping before convergence")
			break
		}
		if(iter==maxiter) {
			warning("iteration limit reached")
			break
		}
	}
	w
}
