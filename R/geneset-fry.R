fry <- function(y,...) UseMethod("fry")

fry.default <- function(y,index=NULL,design=NULL,contrast=ncol(design),weights=NULL,sort="directional",...)
#	Quick version of roast gene set test assuming equal variances between genes
#	The up and down p-values are equivalent to those from roast with nrot=Inf
#	in the special case of prior.df=Inf.
#	Gordon Smyth and Goknur Giner
#	Created 30 January 2015.  Last modified 1 December 2016
{
#	Issue warning if extra arguments found
	dots <- names(list(...))
	if(length(dots)) warning("Extra arguments disregarded: ",sQuote(dots))

#	Extract components from y
	y <- getEAWP(y)
	G <- nrow(y$exprs)
	n <- ncol(y$exprs)

#	Check index
	if(is.null(index)) index <- list(set1=1L:G)
	if(!is.list(index)) index <- list(set1=index)
	nsets <- length(index)
	if(nsets==0) stop("index is empty")

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design))
		design <- matrix(1,n,1)
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
	}
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")
	ne <- nonEstimable(design)
	if(!is.null(ne)) cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")

	p <- ncol(design)
	df.residual <- n-p

#	Check weights
	if(is.null(weights)) weights <- y$weights

#	Reduce to numeric expression matrix
	y <- y$exprs

#	Check weights
	if(!is.null(weights)) {
		if(any(weights<=0)) stop("weights must be positive")
		if(length(weights)==n) {
			sw <- sqrt(weights)
			y <- t(t(y)*sw)
			design <- design*sw
			weights <- NULL
		}
	}
	if(!is.null(weights)) {
		if(length(weights)==G) weights <- matrix(weights,G,n)
		weights <- as.matrix(weights)
		if(any( dim(weights) != dim(y) )) stop("weights not conformal with y")
	}

#	Reform design matrix so that contrast of interest is last column
	if(is.character(contrast)) {
		contrast <- which(contrast==colnames(design))
		if(length(contrast)==0) stop("coef ",contrast," not found")
	}
	if(length(contrast)==1) {
		j <- c((1:p)[-contrast], contrast)
		if(contrast<p) design <- design[,j]
	} else {
		QR <- qr(contrast)
		design <- t(qr.qty(QR,t(design)))
		if(sign(QR$qr[1,1]<0)) design[,1] <- -design[,1]
		design <- design[,c(2:p,1)]
	}

#	Compute effects matrix
	if(is.null(weights)) {
		QR <- qr(design)
		if(QR$rank<p) stop("design matrix is not of full rank")
		effects <- qr.qty(QR,t(y))
		unscaledt <- effects[p,]
		if(QR$qr[p,p]<0) unscaledt <- -unscaledt
	} else {
		effects <- matrix(0,n,G)
		unscaledt <- rep(0,n)
		sw <- sqrt(weights)
		yw <- y*sw
		for (g in 1:G) {
			xw <- design*sw[g,]
			QR <- qr(xw)
			if(QR$rank<p) stop("weighted design matrix not of full rank for gene ",g)
			effects[,g] <- qr.qty(QR,yw[g,])
			unscaledt[g] <- effects[p,g]
			if(QR$qr[p,p]<0) unscaledt[g] <- -unscaledt[g]
		}
	}

#	Standardized residuals
	U <- t(effects[-(1:(p-1)),,drop=FALSE])

#	Global statistics
	nsets <- length(index)
	NGenes <- rep.int(0L,nsets)
	Direction <- rep.int("",nsets)
	PValue.Mixed <- PValue <- rep.int(0,nsets)
	for (i in 1:nsets) {
		iset <- index[[i]]
		USet <- U[iset,,drop=FALSE]
		NGenes[i] <- nrow(USet)
		MeanUSet <- colMeans(USet)
		t.stat <- MeanUSet[1L] / sqrt(mean(MeanUSet[-1L]^2L))
		if(t.stat>0) Direction[i] <- "Up" else Direction[i] <- "Down"
		PValue[i] <- 2*pt(-abs(t.stat),df=df.residual)

		if(NGenes[i]==1) {
			PValue.Mixed[i] <- PValue[i]
		} else {
			SVD <- svd(USet,nu=0)
			A <- SVD$d^2
			d1 <- length(A)
			d <- d1-1L
			beta.mean <- 1/d1
			beta.var <- d/d1/d1/(d1/2+1)
			Fobs <- (sum(USet[,1]^2)-A[d1]) / (A[1]-A[d1])
			Frb.mean <- (sum(A) * beta.mean - A[d1]) / (A[1]-A[d1])
			COV <- matrix(-beta.var/d,d1,d1)
			diag(COV) <- beta.var
			Frb.var <- (A %*% COV %*% A ) / (A[1]-A[d1])^2
			alphaplusbeta <- Frb.mean*(1-Frb.mean)/Frb.var-1
			alpha <- alphaplusbeta*Frb.mean
			beta <- alphaplusbeta-alpha
			PValue.Mixed[i] <- pbeta(Fobs,shape1=alpha,shape2=beta,lower.tail=FALSE)
		}
	}

#	Add FDR
	if(nsets>1) {
		FDR <- p.adjust(PValue,method="BH")
		FDR.Mixed <- p.adjust(PValue.Mixed,method="BH")
		tab <- data.frame(NGenes=NGenes,Direction=Direction,PValue=PValue,FDR=FDR,PValue.Mixed=PValue.Mixed,FDR.Mixed=FDR.Mixed)
	} else {
		tab <- data.frame(NGenes=NGenes,Direction=Direction,PValue=PValue,PValue.Mixed=PValue.Mixed)
	}
	rownames(tab) <- names(index)

#	Sort results
	if(is.logical(sort)) if(sort) sort <- "directional" else sort <- "none"
	sort <- match.arg(sort,c("directional","mixed","none"))
	if(sort=="none") return(tab)
	if(sort=="directional") {
		o <- order(tab$PValue,-tab$NGenes,tab$PValue.Mixed)
	} else {
		o <- order(tab$PValue.Mixed,-tab$NGenes,tab$PValue)
	}
	tab[o,,drop=FALSE]
}
