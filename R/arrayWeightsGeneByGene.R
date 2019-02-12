.arrayWeightsGeneByGene <- function(E, design=NULL, weights=NULL, var.design=NULL, prior.n=10, trace=FALSE)
#	Estimate array variances via gene-by-gene update algorithm
#	Created by Matt Ritchie 7 Feb 2005.
#	Gordon Smyth simplified argument checking to use getEAWP, 9 Mar 2008.
#	Cynthia Liu added var.design argument so that variance model can be modified by user, 22 Sep 2014
#	Last modified 19 Oct 2016.
{
	ngenes <- nrow(E)
	narrays <- ncol(E)
	if(is.null(design)) design <- matrix(1,narrays,1)
	nparams <- ncol(design)

#	Columns of var.design should sum to zero
	if(is.null(var.design)) {
		Z2 <- contr.sum(narrays)
	} else {
		Z2 <- var.design
	}

#	Intialise array gammas to zero (with prior weight of prior.n genes having leverage=0)
	ngam <- ncol(Z2)
	arraygammas <- rep_len(0, ngam)
	Zinfo <- prior.n*crossprod(Z2)
	if(trace) cat("gene convcrit range(w)\n")

#	Step progressive algorithm once through all genes
	for(i in 1:ngenes) {
		if(!all(is.finite(arraygammas))) stop("convergence problem at gene ", i, ": array weights not estimable")
	 	aw <- drop(exp(Z2 %*% (-arraygammas)))
		if(is.null(weights)) {
			w <- aw
		} else {
			w <- aw*weights[i,]
		}
		y <- E[i,]
		obs <- is.finite(y)
		if (sum(obs) > 1) {
			if(sum(obs) == narrays)	{
				X <- design
			} else {
#				remove missing/infinite values
				X <- design[obs, , drop = FALSE]
				y <- y[obs]
				vary <- vary[obs]
				Z2obs <- Z2[obs,,drop=FALSE]
			}
			out <- lm.wfit(X, y, w[obs])
			d <- rep(0, narrays)
			d[obs] <- w[obs]*out$residuals^2
			s2 <- sum(d[obs])/out$df.residual
			h <- rep_len(1, narrays)
			h[obs] <- hat(out$qr)
			Agam <- crossprod(Z2, (1-h)*Z2)
			Agam.del <- crossprod(t(rep(h[narrays], length(arraygammas))-h[1:(length(narrays)-1)]))
			Agene.gam <- (Agam - 1/out$df.residual*Agam.del) # 1/(narrays-nparams)
			if(is.finite(sum(Agene.gam)) && sum(obs) == narrays) {
				Zinfo <- Zinfo + Agene.gam
				R <- chol(Zinfo)
				Zinfoinv <- chol2inv(R)
				zd <- d/s2 - 1 + h
				Zzd <- crossprod(Z2, zd)
				gammas.iter <- Zinfoinv%*%Zzd
				arraygammas <- arraygammas + gammas.iter
			}
			if(is.finite(sum(Agene.gam)) && sum(obs) < narrays && sum(obs) > 2) { 
				Zinfo <- Zinfo + Agene.gam
				A1 <- (diag(1, sum(obs))-1/sum(obs)*matrix(1, sum(obs), sum(obs)))%*%Z2obs
				A1 <- A1[-sum(obs),] # remove last row
				R <- chol(Zinfo)
				Zinfoinv <- chol2inv(R)
				zd <- d/s2 - 1 + h
				Zzd <- A1%*%crossprod(Z2, zd)
				Zinfoinv.A1 <- A1%*%Zinfoinv%*%t(A1)
				alphas.old <- A1%*%arraygammas
				alphas.iter <- Zinfoinv.A1%*%Zzd
				alphas.new <- alphas.old + alphas.iter
				Us <- rbind(diag(1, sum(obs)-1), -1)
				Usalphas <- Us%*%(alphas.new-alphas.old)
				Usgammas <- Z2%*%arraygammas
				Usgammas[obs] <- Usgammas[obs] + Usalphas
				arraygammas <- Usgammas[1:(narrays-1)]
			}

			if(trace && (i%%1000L==1L)) {
				convcrit <- crossprod(Zzd, gammas.iter) / narrays
				cat(i,convcrit,range(1/vary),"\n")
			}
		}
	}

	drop(exp(Z2 %*% (-arraygammas)))
}
