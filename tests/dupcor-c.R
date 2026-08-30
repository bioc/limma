library(limma)
library(statmod)

dupcor.loop <- function(M,design,Array,weights=NULL)
{
	ngenes <- nrow(M)
	nbeta <- ncol(design)
	rho <- rep_len(NA_real_,ngenes)
	nafun <- function(e) NA
	for (i in seq_len(ngenes)) {
		y <- drop(M[i,])
		o <- is.finite(y)
		A <- factor(Array[o])
		nobs <- sum(o)
		nblocks <- length(levels(A))
		if(nobs>(nbeta+2L) && nblocks>1L && nblocks<(nobs-1L)) {
			y <- y[o]
			X <- design[o,,drop=FALSE]
			Z <- model.matrix(~0+A)
			if(!is.null(weights)) {
				w <- drop(weights[i,])[o]
				s <- tryCatch(suppressWarnings(mixedModel2Fit(y,X,Z,w,only.varcomp=TRUE,maxit=20)$varcomp),error=nafun)
			} else
				s <- tryCatch(suppressWarnings(mixedModel2Fit(y,X,Z,only.varcomp=TRUE,maxit=20)$varcomp),error=nafun)
			if(!is.na(s[1])) rho[i] <- s[2]/sum(s)
		}
	}
	rho
}

check.dupcor <- function(M,design,Array,weights=NULL)
{
	Array <- as.integer(factor(Array))
	nblocks <- max(Array)
	actual <- .Call("dupcorfit",M,design,Array,nblocks,weights,1L,PACKAGE="limma")
	parallel <- .Call("dupcorfit",M,design,Array,nblocks,weights,2L,PACKAGE="limma")
	expected <- dupcor.loop(M,design,Array,weights)
	stopifnot(identical(is.na(actual),is.na(expected)))
	stopifnot(identical(is.na(parallel),is.na(actual)))
	stopifnot(isTRUE(all.equal(actual,expected,tolerance=1e-8,check.attributes=FALSE)))
	stopifnot(isTRUE(all.equal(parallel,actual,tolerance=1e-8,check.attributes=FALSE)))
}

set.seed(3)
M <- matrix(rnorm(40*12),40,12)
design <- model.matrix(~factor(rep(1:3,4)))
block <- rep(1:4,each=3)
check.dupcor(M,design,block)
rho <- dupcor.loop(M,design,block)
rho[rho < -0.49] <- -0.49
rho[rho > 0.99] <- 0.99
expected <- list(
	consensus.correlation=tanh(mean(atanh(rho),trim=0.15,na.rm=TRUE)),
	atanh.correlations=atanh(rho)
)
actual <- duplicateCorrelation(M,design,block=block)
stopifnot(isTRUE(all.equal(actual$consensus.correlation,expected$consensus.correlation,tolerance=1e-8)))
stopifnot(isTRUE(all.equal(actual$atanh.correlations,expected$atanh.correlations,tolerance=1e-8)))

weights <- matrix(runif(length(M),0.1,2),nrow(M))
check.dupcor(M,design,block,weights)

M[sample(length(M),40)] <- NA
M[1,] <- NA
M[2,-c(1,2,4)] <- NA
check.dupcor(M,design,block,weights)
check.dupcor(M,cbind(design,design[,2]),block,weights)

M <- matrix(rnorm(80*6),80,6)
design <- cbind(Int=1,Group=c(0,0,0,1,1,1))
M <- unwrapdups(M,ndups=2)
design <- design %x% rep_len(1,2)
Array <- rep(1:6,each=2)
check.dupcor(M,design,Array)

M.zero <- matrix(0,3,6)
check.dupcor(M.zero,matrix(1,6,1),rep(1:3,each=2))

# invalid nthreads is clamped to a valid count in C (no error), so results match nthreads=1L
M.ck <- matrix(rnorm(12*6),12,6)
block.ck <- rep(1:3,each=2)
cor.ref <- duplicateCorrelation(M.ck,block=block.ck,nthreads=1L)
for(nthreads in list(0L,-1L,NA_integer_,1.5,1:2)) {
	cor.clamp <- duplicateCorrelation(M.ck,block=block.ck,nthreads=nthreads)
	stopifnot(isTRUE(all.equal(cor.clamp$consensus.correlation,cor.ref$consensus.correlation,tolerance=1e-11)))
}
