library(limma)

lm.series.loop <- function(M,design,ndups=1,spacing=1,weights=NULL)
{
	M <- as.matrix(M)
	if(!is.null(weights)) {
		weights <- asMatrixWeights(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep_len(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}
	ngenes <- nrow(M)
	nbeta <- ncol(design)
	beta <- stdev.unscaled <- matrix(NA_real_,ngenes,nbeta)
	sigma <- rep_len(NA_real_,ngenes)
	df.residual <- rep_len(0,ngenes)
	for (i in 1:ngenes) {
		y <- as.vector(M[i,])
		obs <- is.finite(y)
		if(sum(obs) > 0) {
			X <- design[obs,,drop=FALSE]
			y <- y[obs]
			if(is.null(weights))
				out <- lm.fit(X,y)
			else
				out <- lm.wfit(X,y,as.vector(weights[i,obs]))
			est <- !is.na(out$coefficients)
			beta[i,] <- out$coefficients
			if(out$rank > 0) {
				stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
				df.residual[i] <- out$df.residual
				if(df.residual[i] > 0) sigma[i] <- sqrt(mean(out$effects[-(1:out$rank)]^2))
			}
		}
	}
	list(coefficients=beta,stdev.unscaled=stdev.unscaled,sigma=sigma,df.residual=df.residual)
}

check.lm.series <- function(M,design,ndups=1,spacing=1,weights=NULL)
{
	actual <- lm.series(M,design,ndups=ndups,spacing=spacing,weights=weights,nthreads=1L)
	parallel <- lm.series(M,design,ndups=ndups,spacing=spacing,weights=weights,nthreads=2L)
	expected <- lm.series.loop(M,design,ndups=ndups,spacing=spacing,weights=weights)
	for(n in names(expected)) {
		stopifnot(isTRUE(all.equal(actual[[n]],expected[[n]],tolerance=1e-12,check.attributes=FALSE)))
		stopifnot(isTRUE(all.equal(parallel[[n]],actual[[n]],tolerance=1e-12,check.attributes=FALSE)))
	}
}

set.seed(1)
M <- matrix(rnorm(8*6),8,6)
design <- cbind(Int=1,Group=c(0,0,0,1,1,1))
M[c(1,10)] <- NA
M[19] <- Inf
check.lm.series(M,design)

weights <- matrix(runif(8*6),8,6)
weights[c(1,8,15,22)] <- c(0,-1,NA,Inf)
check.lm.series(M,design,weights=weights)
check.lm.series(M,cbind(design,design[,2]),weights=weights)

M.rank <- M
M.rank[1,4:6] <- NA
check.lm.series(M.rank,design)
check.lm.series(rbind(M,M),design,ndups=2,weights=rbind(weights,weights))
check.lm.series(matrix(c(NA,1),2,1),matrix(1,1,1))
check.lm.series(matrix(c(1,NA),1,2),matrix(1,2,1))
check.lm.series(matrix(NA_real_,1,2),matrix(1,2,1))

# Rank-0 gene: an intercept-free design whose observed rows are all zero for one
# gene leaves that gene NA (df=0) and continues, rather than aborting the fit.
M.deg <- matrix(rnorm(2*5),2,5)
M.deg[1,3:5] <- NA
design.deg <- matrix(c(0,0,1,1,1),5,1)
check.lm.series(M.deg,design.deg)

M <- matrix(rnorm(8*6),8,6)
fit.fast <- lm.series(M,design)
fit.slow <- lm.series(M,design,weights=matrix(1,8,6))
for(n in c("coefficients","stdev.unscaled","sigma","df.residual"))
	stopifnot(isTRUE(all.equal(fit.fast[[n]],fit.slow[[n]],tolerance=1e-12,check.attributes=FALSE)))

fit.one <- lmFit(M,design,weights=matrix(1,8,6),nthreads=1L)
fit.two <- lmFit(M,design,weights=matrix(1,8,6),nthreads=2L)
stopifnot(isTRUE(all.equal(fit.one$coefficients,fit.two$coefficients,tolerance=1e-12)))

# invalid nthreads is clamped to a valid count in C (no error), so results match nthreads=1L
# the NA forces the genewise C path where the clamp happens
M.ck <- M
M.ck[1] <- NA
fit.ref <- lm.series(M.ck,design,nthreads=1L)
for(nthreads in list(0L,-1L,NA_integer_,1.5,1:2)) {
	fit.clamp <- lm.series(M.ck,design,nthreads=nthreads)
	for(n in c("coefficients","stdev.unscaled","sigma","df.residual"))
		stopifnot(isTRUE(all.equal(fit.clamp[[n]],fit.ref[[n]],tolerance=1e-12,check.attributes=FALSE)))
}
