voom <- function(
	counts, design=NULL,
	lib.size=NULL, offset=NULL, offset.prior=NULL,
	normalize.method="none",
	block=NULL, correlation=NULL, weights=NULL,
	span=0.5, adaptive.span=TRUE, plot=FALSE, save.plot=FALSE
)
#	Linear modelling of count data with mean-variance modelling at the observation level.
#	Creates an EList object with precision weights ready for entry to lmFit() etc in the limma pipeline.
#	Gordon Smyth and Charity Law
#	Created 22 Jun 2011.
#	Support for offset.prior matrix added 11 Apr 2026.
#	Deprecation message added 8 Oct 2026.
#	Support for offset matrix added 1 Nov 2026.
{
#	Function is deprecated in favor of voomLmFit()
	message("voomLmFit is now recommended over voom")

#	Create output list
	out <- list()

#	Extract counts from known data objects
	if(is(counts,"DGEList")) {
		out$genes <- counts$genes
		out$targets <- counts$samples
		if(is.null(design) && diff(range(as.numeric(counts$samples$group)))>0)
			design <- model.matrix(~group,data=counts$samples)
#		Normalized library sizes
		if(is.null(lib.size)) {
			nf <- counts$samples$norm.factors
			if(is.null(nf)) nf <- 1
			lib.size <- counts$samples$lib.size * nf
		}
		if(is.null(offset)) offset <- counts[["offset"]]
		if(is.null(offset.prior)) offset.prior <- counts[["offset.prior"]]
		counts <- counts$counts
	} else if(is(counts,"SummarizedExperiment")) {
		if(!requireNamespace("SummarizedExperiment",quietly=TRUE))
			stop("SummarizedExperiment package required but is not installed (or can't be loaded)")
		se <- counts
		an <- SummarizedExperiment::assayNames(se)
		if("counts" %in% an) {
			counts <- SummarizedExperiment::assay(se,"counts")
		} else {
			counts <- SummarizedExperiment::assay(se,1L)
		}
		counts <- as.matrix(counts)
		rd <- SummarizedExperiment::rowData(se)
		if(ncol(rd)) out$genes <- as.data.frame(rd)
		cd <- SummarizedExperiment::colData(se)
		if(ncol(cd)) out$targets <- as.data.frame(cd)
		if(is.null(design) && ("group" %in% colnames(cd)) && diff(range(as.numeric(as.factor(cd$group))))>0)
			design <- model.matrix(~group,data=as.data.frame(cd))
	} else if(is(counts,"eSet")) {
		if(!requireNamespace("Biobase",quietly=TRUE))
			stop("Biobase package required but is not installed (or can't be loaded)")
		if(length(Biobase::fData(counts))) out$genes <- Biobase::fData(counts)
		if(length(Biobase::pData(counts))) out$targets <- Biobase::pData(counts)
		if("counts" %in% names(Biobase::assayData(counts))) {
			counts <- get("counts",Biobase::assayData(counts))
		} else {
			counts <- Biobase::exprs(counts)
			message("Extracting exprs(counts). These should be raw counts.")
		}
	} else {
		counts <- as.matrix(counts)
	}

#	Check counts
	ngenes <- nrow(counts)
	if(ngenes < 2L) stop("Need at least two genes to fit a mean-variance trend")
	m <- min(counts)
	if(is.na(m)) stop("NA counts not allowed")
	if(m < 0) stop("Negative counts not allowed")
	nsamples <- ncol(counts)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,ncol(counts),1)
		rownames(design) <- colnames(counts)
		colnames(design) <- "GrandMean"
	}

#	Library sizes default to column sums
	if(is.null(lib.size)) lib.size <- colSums(counts)

#	Combine library sizes with offsets.
#	offset takes precedence over lib.size and offset.prior.
#	offset.prior adds to log(lib.size) while offset replaces them.
	if(is.null(offset)) {
		if(is.null(offset.prior)) {
			lib.size.matrix <- matrix(lib.size,ngenes,nsamples,byrow=TRUE)
		} else {
			message("Using offset.prior matrix")
			if(!identical(dim(counts),dim(offset.prior))) stop("counts and offset.prior must have equal dimensions.")
			lib.size.matrix <- exp(matrix(log(lib.size),ngenes,nsamples,byrow=TRUE) + offset.prior)
		}
	} else {
		message("Using offset matrix")
#		Offset can be a matrix or a sample-wise vector.
		if(is.matrix(offset)) {
			if(!identical(dim(counts),dim(offset))) stop("counts and offset must have equal dimensions")
			lib.size.matrix <- exp(offset)
		} else {
			if(!identical(length(offset),nsamples)) stop("if offset is a vector, its length must be the number of samples")
			lib.size.matrix <- matrix(exp(offset),ngenes,nsamples,byrow=TRUE)
		}
	}

#	Choose span based on the number of genes
	if(adaptive.span) span <- chooseLowessSpan(ngenes, small.n=50, min.span=0.3, power=1/3)

#	Fit linear model to log2-counts-per-million
	y <- log2((counts+0.5)/(lib.size.matrix+1)*1e6)
	y <- normalizeBetweenArrays(y,method=normalize.method)
	fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights)
	if(is.null(fit$Amean)) fit$Amean <- rowMeans(y,na.rm=TRUE)

#	If no replication found, set all weight to 1
	NWithReps <- sum(fit$df.residual > 0L)
	if(NWithReps < 2L) {
		if(NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
		if(NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
		out$E <- y
		out$weights <- y
		out$weights[] <- 1
		out$design <- design
		if(is.null(out$targets))
			out$targets <- data.frame(lib.size=lib.size)
		else
			out$targets$lib.size <- lib.size
		out$other$offset.prior <- offset.prior
		return(new("EList",out))
	}

#	Fit lowess trend to sqrt-standard-deviations by log-count-size
	sx <- fit$Amean+mean(log2(lib.size+1))-log2(1e6)
	sy <- sqrt(fit$sigma)
	allzero <- rowSums(counts)==0
	if(any(allzero)) {
		sx <- sx[!allzero]
		sy <- sy[!allzero]
	}
	l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title("voom: Mean-variance trend")
		lines(l,col="red")
	}

#	Make interpolating rule
#	Special treatment of zero counts is now removed;
#	instead zero counts get same variance as smallest gene average.
#	l$x <- c(0.5^0.25, l$x)
#	l$x <- c(log2(0.5), l$x)
#	var0 <- var(log2(0.5*1e6/(lib.size+0.5)))^0.25
#	var0 <- max(var0,1e-6)
#	l$y <- c(var0, l$y)
	f <- approxfun(l, rule=2, ties=list("ordered",mean))

#	Find individual quarter-root fitted counts
	if(fit$rank < ncol(design)) {
		j <- fit$pivot[1:fit$rank]
		fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
	} else {
		fitted.values <- fit$coefficients %*% t(fit$design)
	}
	fitted.cpm <- 2^fitted.values
	fitted.count <- 1e-6 * fitted.cpm * (lib.size.matrix+1)
	fitted.logcount <- log2(fitted.count)

#	Apply trend to individual observations
	w <- 1/f(fitted.logcount)^4
	dim(w) <- dim(fitted.logcount)

#	Output
	out$E <- y
	out$weights <- w
	out$design <- design
	if(is.null(out$targets)) {
		out$targets <- data.frame(lib.size=lib.size)
	} else {
		out$targets$lib.size <- lib.size
	}
	if(adaptive.span) out$span <- span
	if(save.plot) {
		out$voom.xy <- list(x=sx,y=sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		out$voom.line <- l
	}
	out$other$offset <- offset
	out$other$offset.prior <- offset.prior

	new("EList",out)
}
