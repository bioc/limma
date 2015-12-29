#  Plot regularized linear discriminant functions

plotRLDF <- function(y,design=NULL,z=NULL,labels.y=NULL,labels.z=NULL,col.y="black",col.z="black",df.prior=5,show.dimensions=c(1,2),ndim=max(show.dimensions),nprobes=100,plot=TRUE,...)
#	Regularized linear discriminant function
#	Di Wu and Gordon Smyth
#	29 June 2009.  Last revised 10 Dec 2015.
{
#	Check y
	y <- as.matrix(y)
	g <- nrow(y)
	n <- ncol(y)

#	Check labels.y
	if(is.null(labels.y))
		labels.y <- colnames(y)
	else {
		if(length(labels.y) != n) stop("length(labels.y) doesn't agree with ncol(y)")
		labels.y <- as.character(labels.y)
	}
	if(is.null(labels.y)) labels.y <- as.character(1:n)

#	Check design
	if(is.null(design)) {
		if(is.null(labels.y)) stop("groups not specified")
		if(!anyDuplicated(labels.y)) stop("design not specified and all labels.y are different")
		f <- as.factor(labels.y)
		design <- model.matrix(~f)
	} else {
		design <- as.matrix(design)
		if(nrow(design) != n) stop("nrow(design) doesn't match ncol(y)")
	}

#	Check show.dimensions
	show.dimensions <- as.integer(show.dimensions)[1:2]
	if(show.dimensions[1]==show.dimensions[2]) stop("show.dimensions must specify two different columns")
	if(any(show.dimensions<1) || any(show.dimensions)>n) stop("show.dimensions must be a column number of y")

#	Check ndim
	if(any(show.dimensions)>ndim) ndim <- max(show.dimensions)

#	Check nprobes
	if(nprobes<1) stop("'nprobes' must be at least 1")

#	Project onto between and within spaces
#	Discard first column as intercept
	qrd <- qr(design)
	p <- qrd$rank
	if(p==n) stop("No residual degrees of freedom")
	U <- qr.qty(qrd, t(y))
	UB <- U[2:p,,drop=FALSE]
	UW <- U[(p+1):n,,drop=FALSE]

#	Prior variance
	s <- colMeans(UW*UW)
	s0 <- median(s)

#	Select probes by moderated F
	if(g>nprobes) {
		modF <- colMeans(UB*UB)/(s+df.prior*s0)
		o <- order(modF,decreasing=TRUE)
		top <- o[1:nprobes]
		y <- y[top,,drop=FALSE]
		if(!is.null(z)) z <- z[top,,drop=FALSE]
		UB <- UB[,top,drop=FALSE]
		UW <- UW[,top,drop=FALSE]
		g <- nprobes
	} else {
		top <- 1:nprobes
	}

#	Within group SS
	W <- crossprod(UW)

#	Regularized within group SS
	Wreg <- W+diag(df.prior*s0,g,g)

#	Ratio of between to within SS
	WintoB <- backsolve(chol(Wreg),t(UB),transpose=TRUE)

#	Linear discriminant gene weights
	d1 <- show.dimensions[1]
	d2 <- show.dimensions[2]
	SVD <- svd(WintoB,nu=ndim,nv=0)
	metagenes <- SVD$u

	rank <- min(dim(WintoB))
	if(ndim > rank) message("Note: only ",rank," discriminant function(s) have any predictive power")

#	LDF for training set
	d.y <- t(y) %*% metagenes
	d1.y <- d.y[,d1]
	d2.y <- d.y[,d2]

#	LDF for classified set
	if(!is.null(z)) {
		z <- as.matrix(z)
		d.z <- t(z) %*% metagenes
		d1.z <- d.z[,d1]
		d2.z <- d.z[,d2]
		if(is.null(labels.z)) labels.z <- colnames(z)
		if(is.null(labels.z)) labels.z <- letters[1:ncol(z)]
		if(length(labels.z) != ncol(z)) stop("length(labels.z) doesn't agree with ncol(z)")
		labels.z <- as.character(labels.z)
	} else {
		d1.z <- d2.z <- numeric(0)
		labels.z <- character(0)
	}

#	Make plot
	if(plot) {
		plot(c(d1.y,d1.z),c(d2.y,d2.z),type="n", xlab=paste("Discriminant Function",d1), ylab=paste("Discriminant Function",d2))
		text(d1.y,d2.y,labels=labels.y,col=col.y,...)
		if(!is.null(z)) text(d1.z,d2.z,labels=labels.z,col=col.z,...)
	}

#	Output
	out <- list(training=d.y,top=top,metagenes=metagenes,singular.values=SVD$d,rank=rank)
	if(!is.null(z)) out$predicting <- d.z
	invisible(out)
}

