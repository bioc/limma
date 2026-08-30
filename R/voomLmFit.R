voomLmFit <- function(
	counts, design=NULL, contrasts=NULL, block=NULL, prior.weights=NULL,
	sample.weights=FALSE, var.design=NULL, var.group=NULL, prior.n=10,
	lib.size=NULL, offset=NULL, offset.prior=NULL, normalize.method="none",
	span=0.5, adaptive.span=TRUE, plot=FALSE, save.plot=FALSE, keep.EList=TRUE,
	nthreads=1L
)
#	Implements the voom+lmFit+contrasts.fit pipeline for counts while taking
#	into account fitted values and residuals that are exactly zero and
#	contribute no information about the residual variances.
#	Creates an MArrayLM object for entry to eBayes() etc in the limma pipeline.
#	Columns will be contrasts if the contrasts argument is set.
#
#	This version is self-contained within limma, although the
#	SummarizedExperiment and Biobase packages are assumed to be available if
#	counts is of SummarizedExperiment class.
#
#	Gordon Smyth. Self-contained limma version by Lizhong Chen, including C code
#	and contrasts argument.
#	Created 21 Jan 2020.  Last modified 16 Aug 2026.
{
	Block <- !is.null(block)
	PriorWeights <- !is.null(prior.weights)
	SampleWeights <- sample.weights || !is.null(var.design) || !is.null(var.group)

#	Can't specify prior weights and ask for sample weights to be estimated as well
	if(PriorWeights && SampleWeights) stop("Can't specify prior.weights and estimate sample weights")

#	Create output object
	out <- list()

#	Extract counts from known data objects
	if(is(counts,"DGEList")) {
		out$genes <- counts$genes
		out$targets <- counts$samples
		if(is.null(design) && diff(range(as.numeric(counts$samples$group)))>0) design <- model.matrix(~group,data=counts$samples)
#		Effective (normalized) library sizes, as getNormLibSizes() returns for a DGEList
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
		if(!is.null(an) && "counts" %in% an)
			counts <- SummarizedExperiment::assay(se,"counts")
		else
			counts <- SummarizedExperiment::assay(se,1L)
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
	G <- nrow(counts)
	if(G < 2L) stop("Need at least two genes to fit a mean-variance trend")
	m <- min(counts)
	if(is.na(m)) stop("NA counts not allowed")
	if(m < 0) stop("Negative counts not allowed")
	n <- ncol(counts)

#	Check design
	if(is.null(design)) {
		design <- matrix(1,n,1)
		rownames(design) <- colnames(counts)
		colnames(design) <- "GrandMean"
	}

#	Check contrasts. If given, the final fit is returned in contrast space, with
#	exact per-gene stdev.unscaled because voom weights are probe-specific.
	if(!is.null(contrasts)) {
		contrasts <- as.matrix(contrasts)
		if(!is.numeric(contrasts)) stop("contrasts must be a numeric matrix")
		if(anyNA(contrasts)) stop("NAs not allowed in contrasts")
		if(!identical(nrow(contrasts),ncol(design))) stop("Number of rows of contrasts must match number of columns of design")
	}

#	Check library sizes
	if(is.null(lib.size)) lib.size <- colSums(counts)

#	Offset matrix takes precedence over lib.size and offset.prior. Otherwise, lib.sizes default to column sums.
#	offset.prior adds to log(lib.size) while offset replaces them.
	if(is.null(offset)) {
		if(is.null(offset.prior)) {
			lib.size.matrix <- matrix(lib.size,G,n,byrow=TRUE)
		} else {
			message("Using offset.prior matrix")
			if(!identical(dim(counts),dim(offset.prior))) stop("counts and offset.prior must have equal dimensions.")
			lib.size.matrix <- exp(matrix(log(lib.size),G,n,byrow=TRUE) + offset.prior)
		}
	} else {
		message("Using offset matrix")
#		Offset can be a matrix or a row vector.
		if(is.matrix(offset)) {
			if(!identical(dim(counts),dim(offset))) stop("counts and offset must have equal dimensions")
		} else {
			if(!identical(length(offset),n)) stop("if offset is a vector, its length must be the number of samples")
		}
		lib.size.matrix <- exp(offset)
	}

#	Expand prior.weights if necessary
	if(!is.null(prior.weights)) prior.weights <- asMatrixWeights(prior.weights,dim(counts))

#	Choose span based on the number of genes
	if(adaptive.span) span <- chooseLowessSpan(nrow(counts), small.n=50, min.span=0.3, power=1/3)

#	log2-counts-per-million
	y <- log2((counts+0.5)/(lib.size.matrix+1)*1e6)

#	Microarray-style normalization
	y <- normalizeBetweenArrays(y,method=normalize.method)

#	Fit linear model
	fit <- lmFit(y,design,weights=prior.weights,nthreads=nthreads)

#	Find largest leverage value of design matrix
	if(is.null(fit$qr))
		h <- hat(design,intercept=FALSE)
	else
		h <- hat(fit$qr)
	MinGroupSize <- 1/max(h)

#	Identify fitted values that are exactly zero and should not contribute to the genewise variances
#	Note that a single zero is never a problem
	eps <- 1e-4
	RowHasZero <- which(rowSums(counts < eps) > (max(2,MinGroupSize)-eps))
	AnyZeroRows <- as.logical(length(RowHasZero))
	if(AnyZeroRows) {
		countsZero <- counts[RowHasZero,,drop=FALSE]
#		Poisson fit (log link, offset=log(lib.size), dispersion=0) to flag structural-zero
#		fitted values, warm-started from the lmFit coefficients rescaled to the log-count
#		scale: log(mu) = log(2)*log2CPM + log(lib.size) - log(1e6)
		startZero <- fit$coefficients[RowHasZero,,drop=FALSE] %*% t(fit$design)
		startZero <- log(2)*startZero - log(1e6)
		PoissonFitted <- .Call("poisfit",countsZero,design,log(lib.size),startZero,nthreads,PACKAGE="limma")
		IsZero <- (PoissonFitted < eps & countsZero < eps)
		RowHasExactZero <- which(rowSums(IsZero) > eps)
#		If any exact zero fits, then rerun the linear model for those rows with NAs
		if(length(RowHasExactZero)) {
			RowHasZero <- RowHasZero[RowHasExactZero]
			IsZero <- IsZero[RowHasExactZero,,drop=FALSE]
			yNAshort <- y[RowHasZero,,drop=FALSE]
			yNAshort[IsZero] <- NA
			fitNA <- suppressWarnings(lmFit(yNAshort,design,weights=prior.weights[RowHasZero,,drop=FALSE],nthreads=nthreads))
			fit$df.residual[RowHasZero] <- fitNA$df.residual
			fit$sigma[RowHasZero] <- fitNA$sigma
#			If blocking or sample weights are present, then we will later on need a full length copy of y with NAs inserted
			if(Block || SampleWeights) {
				yNAfull <- y
				yNAfull[RowHasZero,] <- yNAshort
			}
		} else {
			AnyZeroRows <- FALSE
		}
	}

#	If no replication found, assume all weights are 1 and return fit already computed
	HasRep <- (fit$df.residual > 0L)
	NWithReps <- sum(HasRep)
	if(NWithReps < 2L) {
		if(NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
		if(NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
		fit$genes <- out$genes
		if(!is.null(contrasts)) fit <- contrasts.fit(fit,contrasts)
		return(fit)
	}

#	Fit lowess trend to sqrt-standard-deviations by log-count-size
	Amean <- Amean2 <- rowMeans(y)
	if(AnyZeroRows) Amean2[RowHasZero] <- rowMeans(yNAshort,na.rm=TRUE)
	sx <- Amean2[HasRep]+mean(log2(lib.size+1))-log2(1e6)
	sy <- sqrt(fit$sigma[HasRep])
	if(AnyZeroRows)
		l <- weightedLowess(sx,sy,span=span,weights=fit$df.residual[HasRep],output.style="lowess")
	else
		l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title("voom: Mean-variance trend")
		lty <- ifelse(Block || SampleWeights,2,1)
		lines(l,col="red",lty=lty)
	}

#	Make interpolating rule
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

#	Apply trend to individual observations to get voom weights
	w <- 1/f(fitted.logcount)^4
	dim(w) <- dim(fitted.logcount)

#	Add voom weights to prior weights
	if(PriorWeights) {
		weights <- w * prior.weights
		attr(weights,"arrayweights") <- NULL
	} else
		weights <- w

#	Estimate sample weights?
	if(SampleWeights) {
		if(AnyZeroRows) {
			sw <- arrayWeights(yNAfull,design,weights=weights,var.design=var.design,var.group=var.group,prior.n=prior.n,method="genebygene",nthreads=nthreads)
		} else {
			sw <- arrayWeights(y,design,weights=weights,var.design=var.design,var.group=var.group,prior.n=prior.n,method="reml",nthreads=nthreads)
		}
		message("First sample weights (min/max) ", paste(format(range(sw)),collapse="/") )
		if(Block) weights <- t(sw * t(weights))
	}

#	Estimate correlation?
	if(Block) {
		if(AnyZeroRows) {
			dc <- duplicateCorrelation(yNAfull,design,block=block,weights=weights,nthreads=nthreads)
		} else {
			dc <- duplicateCorrelation(y,design,block=block,weights=weights,nthreads=nthreads)
		}
		correlation <- dc$consensus.correlation
		if(is.na(correlation)) {
			warning("Intra-block correlation not estimable, setting to zero.", call. = FALSE)
			correlation <- 0
		}
		if(identical(correlation,0)) {
			Block <- FALSE
			if(plot) lines(l,col="red",lty=1)
		} else {
			message("First intra-block correlation  ",format(correlation))
		}
	} else {
		correlation <- NULL
	}

#	Second iteration to refine intra-block correlation or sample weights
	if(Block || SampleWeights) {
#		Rerun voom weights with new correlation and sample weights
		if(SampleWeights)
			weights <- asMatrixWeights(sw,dim(y))
		else
			weights <- prior.weights
		fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights,nthreads=nthreads)
		if(AnyZeroRows) {
			fitNA <- suppressWarnings(lmFit(yNAshort,design,block=block,correlation=correlation,weights=weights[RowHasZero,,drop=FALSE],nthreads=nthreads))
			fit$df.residual[RowHasZero] <- fitNA$df.residual
			fit$sigma[RowHasZero] <- fitNA$sigma
		}
		sy <- sqrt(fit$sigma[HasRep])
		if(AnyZeroRows)
			l <- weightedLowess(sx,sy,span=span,weights=fit$df.residual[HasRep],output.style="lowess")
		else
			l <- lowess(sx,sy,f=span)
		if(plot) {
			lines(l,col="red")
			legend("topright",lty=c(2,1),col="red",legend=c("First","Final"))
		}
		f <- approxfun(l, rule=2, ties=list("ordered",mean))
		if(fit$rank < ncol(design)) {
			j <- fit$pivot[1:fit$rank]
			fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
		} else {
			fitted.values <- fit$coefficients %*% t(fit$design)
		}
		fitted.cpm <- 2^fitted.values
		fitted.count <- 1e-6 * fitted.cpm * (lib.size.matrix+1)
		fitted.logcount <- log2(fitted.count)
		w <- 1/f(fitted.logcount)^4
		dim(w) <- dim(fitted.logcount)
		if(PriorWeights) {
			weights <- w * prior.weights
			attr(weights,"arrayweights") <- NULL
		} else
			weights <- w
		if(SampleWeights) {
			if(AnyZeroRows) {
				sw <- arrayWeights(yNAfull,design,weights=weights,var.design=var.design,var.group=var.group,prior.n=prior.n,method="genebygene",nthreads=nthreads)
			} else {
				sw <- arrayWeights(y,design,weights=weights,var.design=var.design,var.group=var.group,prior.n=prior.n,method="reml",nthreads=nthreads)
			}
			message("Final sample weights (min/max) ", paste(format(range(sw)),collapse="/") )
			weights <- t(sw * t(weights))
		}
		if(Block) {
			if(AnyZeroRows) {
				dc <- suppressWarnings(duplicateCorrelation(yNAfull,design,block=block,weights=weights,nthreads=nthreads))
			} else {
				dc <- suppressWarnings(duplicateCorrelation(y,design,block=block,weights=weights,nthreads=nthreads))
			}
			correlation <- dc$consensus.correlation
			if(is.na(correlation)) {
				warning("Intra-block correlation not estimable, setting to zero.")
				correlation <- 0
			}
			message("Final intra-block correlation  ",format(correlation))
		}
	}

#	Final linear model fit with voom weights. When contrasts are supplied they are
#	applied here, so the per-gene contrast stdev.unscaled is exact under voom weights.
	fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights,nthreads=nthreads,contrasts=contrasts)
	if(is.null(fit$Amean)) fit$Amean <- Amean
	if(AnyZeroRows) {
		fitNA <- suppressWarnings(lmFit(yNAshort,design,block=block,correlation=correlation,weights=weights[RowHasZero,,drop=FALSE],nthreads=nthreads))
		fit$df.residual[RowHasZero] <- fitNA$df.residual
		fit$sigma[RowHasZero] <- fitNA$sigma
	}

#	Output
	fit$genes <- out$genes
	fit$targets <- out$targets
	if(is.null(fit$targets)) {
		fit$targets <- data.frame(lib.size=lib.size)
		row.names(fit$targets) <- colnames(y)
	}
	if(SampleWeights) fit$targets$sample.weight <- sw
	if(save.plot) {
		fit$voom.xy <- list(x=sx,y=sy,xlab="log2( count size + 0.5 )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		fit$voom.line <- l
	}
	if(keep.EList) {
		fit$EList <- new("EList",list(E=y,weights=weights,genes=out$genes))
	}
	fit
}
