library(limma)

gls.series.loop <- function(M,design,cormatrix,weights=NULL)
{
	ngenes <- nrow(M)
	nbeta <- ncol(design)
	beta <- stdev.unscaled <- matrix(NA_real_,ngenes,nbeta)
	sigma <- rep_len(NA_real_,ngenes)
	df.residual <- rep_len(0,ngenes)
	for (i in 1:ngenes) {
		y <- drop(M[i,])
		o <- is.finite(y)
		y <- y[o]
		n <- length(y)
		if(n > 0) {
			X <- design[o,,drop=FALSE]
			V <- cormatrix[o,o,drop=FALSE]
			if(!is.null(weights)) {
				wrs <- 1/sqrt(drop(weights[i,o]))
				V <- wrs * t(wrs * t(V))
			}
			cholV <- chol(V)
			y <- backsolve(cholV,y,transpose=TRUE)
			if(all(X==0)) {
				df.residual[i] <- n
				sigma[i] <- sqrt(array(1/n,c(1,n)) %*% y^2)
			} else {
				X <- backsolve(cholV,X,transpose=TRUE)
				out <- lm.fit(X,y)
				est <- !is.na(out$coefficients)
				beta[i,] <- out$coefficients
				stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
				df.residual[i] <- out$df.residual
				if(df.residual[i] > 0)
					sigma[i] <- sqrt(array(1/out$df.residual,c(1,n)) %*% out$residuals^2)
			}
		}
	}
	list(coefficients=beta,stdev.unscaled=stdev.unscaled,sigma=sigma,df.residual=df.residual)
}

check.gls.series <- function(M,design,cormatrix,weights=NULL)
{
	actual <- .Call("glsfit",M,design,cormatrix,weights,NULL,1L,PACKAGE="limma")
	parallel <- .Call("glsfit",M,design,cormatrix,weights,NULL,2L,PACKAGE="limma")
	expected <- gls.series.loop(M,design,cormatrix,weights)
	for(n in names(expected)) {
		stopifnot(isTRUE(all.equal(actual[[n]],expected[[n]],tolerance=1e-11,check.attributes=FALSE)))
		stopifnot(isTRUE(all.equal(parallel[[n]],actual[[n]],tolerance=1e-11,check.attributes=FALSE)))
	}
}

set.seed(2)
M <- matrix(rnorm(8*6),8,6)
M[c(1,10)] <- NA
M[19] <- Inf
design <- cbind(Int=1,Group=c(0,0,0,1,1,1))
block <- c(1,1,2,2,3,3)
cormatrix <- outer(block,block,"==") * 0.3
diag(cormatrix) <- 1
check.gls.series(M,design,cormatrix)
fit <- gls.series(M,design,block=block,correlation=0.3)
expected <- gls.series.loop(M,design,cormatrix)
for(n in names(expected))
	stopifnot(isTRUE(all.equal(fit[[n]],expected[[n]],tolerance=1e-11,check.attributes=FALSE)))

weights <- matrix(runif(8*6,0.2,2),8,6)
check.gls.series(M,design,cormatrix,weights)
check.gls.series(M,cbind(design,design[,2]),cormatrix,weights)
check.gls.series(M,matrix(0,6,2),cormatrix,weights)
check.gls.series(matrix(c(1,NA),1,2),matrix(1,2,1),diag(2))
check.gls.series(matrix(NA_real_,1,2),matrix(1,2,1),diag(2))

M.rank <- M
M.rank[1,4:6] <- NA
check.gls.series(M.rank,design,cormatrix,weights)

weights.raw <- weights
weights.raw[c(1,8,15)] <- c(0,-1,NA)
weights.raw[is.na(weights.raw)] <- 0
M.weights <- M
M.weights[weights.raw < 1e-15] <- NA
weights.raw[weights.raw < 1e-15] <- NA
check.gls.series(M.weights,design,cormatrix,weights.raw)

M.dups <- matrix(rnorm(8*3),8,3)
M.dups[c(2,13)] <- NA
design.dups <- cbind(Int=1,Group=c(0,1,1))
M.dups <- unwrapdups(M.dups,ndups=2)
design.dups <- design.dups %x% rep_len(1,2)
cormatrix.dups <- diag(rep_len(0.4,3),nrow=3,ncol=3) %x% array(1,c(2,2))
diag(cormatrix.dups) <- 1
check.gls.series(M.dups,design.dups,cormatrix.dups)

M.complete <- matrix(rnorm(8*6),8,6)
fit.fast <- gls.series(M.complete,design,block=block,correlation=0.3)
fit.slow <- gls.series(M.complete,design,block=block,correlation=0.3,weights=matrix(1,8,6))
for(n in c("coefficients","stdev.unscaled","sigma","df.residual"))
	stopifnot(isTRUE(all.equal(fit.fast[[n]],fit.slow[[n]],tolerance=1e-11,check.attributes=FALSE)))

bad.cormatrix <- matrix(c(1,2,2,1),2,2)
bad.one <- try(.Call("glsfit",matrix(1,4,2),matrix(1,2,1),bad.cormatrix,NULL,NULL,1L,PACKAGE="limma"),silent=TRUE)
bad.two <- try(.Call("glsfit",matrix(1,4,2),matrix(1,2,1),bad.cormatrix,NULL,NULL,2L,PACKAGE="limma"),silent=TRUE)
stopifnot(inherits(bad.one,"try-error"), identical(as.character(bad.one),as.character(bad.two)))

fit.one <- lmFit(M,design,block=block,correlation=0.3,weights=weights,nthreads=1L)
fit.two <- lmFit(M,design,block=block,correlation=0.3,weights=weights,nthreads=2L)
stopifnot(isTRUE(all.equal(fit.one$coefficients,fit.two$coefficients,tolerance=1e-11)))

fit.one <- suppressWarnings(gls.series(M,design,block=block,weights=weights,nthreads=1L))
fit.two <- suppressWarnings(gls.series(M,design,block=block,weights=weights,nthreads=2L))
stopifnot(isTRUE(all.equal(fit.one$coefficients,fit.two$coefficients,tolerance=1e-11)))

# invalid nthreads is clamped to a valid count in C (no error), so results match nthreads=1L
fit.ref <- gls.series(M,design,block=block,correlation=0.3,nthreads=1L)
for(nthreads in list(0L,-1L,NA_integer_,1.5,1:2)) {
	fit.clamp <- gls.series(M,design,block=block,correlation=0.3,nthreads=nthreads)
	stopifnot(isTRUE(all.equal(fit.clamp$coefficients,fit.ref$coefficients,tolerance=1e-11,check.attributes=FALSE)))
}
