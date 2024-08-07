voomaLmFit <- function(
	y, design=NULL, prior.weights=NULL,
	block=NULL,
	sample.weights=FALSE, var.design=NULL, var.group=NULL, prior.n=10,
	predictor=NULL, span=NULL, legacy.span=FALSE, plot=FALSE, save.plot=FALSE,
	keep.EList=TRUE
)
#	vooma+lmFit with iteration of sample weights and intrablock correlation.
#	Creates an MArrayLM object for entry to eBayes() etc in the limma pipeline.
#	Mengbo Li and Gordon Smyth
#	Created 24 Nov 2023. Last modifed 15 Jul 2024.
{
	Block <- !is.null(block)
	PriorWeights <- !is.null(prior.weights)
	SampleWeights <- sample.weights || !is.null(var.design) || !is.null(var.group)
	
#	Can't specify prior weights and ask for sample weights to be estimated as well
	if(PriorWeights && SampleWeights) stop("Can't specify prior.weights and estimate sample weights")
	
#	Check y
	if(!is(y,"EList")) y <- new("EList",list(E=as.matrix(y)))
	narrays <- ncol(y)
	ngenes <- nrow(y)
	
#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design)) {
		design <- matrix(1,narrays,1)
		rownames(design) <- colnames(y)
		colnames(design) <- "GrandMean"
	}
	
#	Expand prior.weights if necessary
	if(!is.null(prior.weights)) prior.weights <- asMatrixWeights(prior.weights,dim(y))
	
#	Fit linear model
	fit <- lmFit(y, design, weights = prior.weights)
	
#	Find individual fitted values
	if(fit$rank < ncol(design)) {
		j <- fit$pivot[1:fit$rank]
		fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
	} else {
		fitted.values <- fit$coefficients %*% t(fit$design)
	}
	
#	Prepare to fit lowess trend
	sx <- rowMeans(y$E)
	sy <- sqrt(fit$sigma)
	mu <- fitted.values
	
#	Optionally combine ave log intensity with precision predictor
	if(!is.null(predictor)) {
		sxc <- rowMeans(predictor)
		vartrend <- lm.fit(cbind(1,sx,sxc),sy)
		beta <- coef(vartrend)
		sx <- vartrend$fitted.values
		mu <- beta[1] + beta[2]*mu + beta[3]*predictor
		xlab <- "Combined predictor"
		main.title <- "vooma variance trend"
	} else {
		xlab <- "Average log2 expression"
		main.title <- "vooma mean-variance trend"
	}

#	Choose span based on the number of genes
	if(is.null(span))
		if(legacy.span)
			span <- chooseLowessSpan(ngenes,small.n=10,min.span=0.3,power=0.5) 
		else
			span <- chooseLowessSpan(ngenes,small.n=50,min.span=0.3,power=1/3)
	
#	Fit lowess trend
	l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab=xlab,ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title(main.title)
		lty <- ifelse(Block || SampleWeights,2,1)
		lines(l,col="red",lty=lty)
	}
	
#	Make interpolating rule
	f <- approxfun(l, rule=2, ties=list("ordered",mean))
	
#	Apply trend to individual observations to get vooma weights
	w <- 1/f(mu)^4
	dim(w) <- dim(y)
	colnames(w) <- colnames(y)
	rownames(w) <- rownames(y)
	
#	Add voom weights to prior weights
	if(PriorWeights) {
		weights <- w * prior.weights
		attr(weights,"arrayweights") <- NULL
	} else {
		weights <- w
	}
	
#	Estimate sample weights?
	if(SampleWeights) {
		sw <- arrayWeights(y,design,weights=weights,var.design=var.design,var.group=var.group,prior.n=prior.n)
		message("First sample weights (min/max) ", paste(format(range(sw)),collapse="/") )
		if(Block) weights <- t(sw * t(weights))
	}
	
#	Estimate correlation?
	if(Block) {
		dc <- suppressWarnings(duplicateCorrelation(y,design,block=block,weights=weights))
		correlation <- dc$consensus.correlation
		if(is.na(correlation)) correlation <- 0
		message("First intra-block correlation  ",format(correlation))
	} else {
		correlation <- NULL
	}
	
#	Second iteration to refine intra-block correlation or sample weights
	if(Block || SampleWeights) {
#		Rerun voom weights with new correlation and sample weights
		if(SampleWeights) {
			weights <- asMatrixWeights(sw,dim(y))
		} else {
			weights <- prior.weights
		}
		fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights)
#		Find individual fitted values
		if(fit$rank < ncol(design)) {
			j <- fit$pivot[1:fit$rank]
			fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
		} else {
			fitted.values <- fit$coefficients %*% t(fit$design)
		}
#		Prepare to fit NEW lowess trend
		sx <- rowMeans(y$E)
		sy <- sqrt(fit$sigma)
		mu <- fitted.values
#		Optionally combine ave log intensity with precision predictor for the NEW trend
		if(!is.null(predictor)) {
			vartrend <- lm.fit(cbind(1,sx,sxc),sy)
			beta <- coef(vartrend)
			sx <- vartrend$fitted.values
			mu <- beta[1] + beta[2]*mu + beta[3]*predictor
		}
		l <- lowess(sx,sy,f=span)
		if(plot) {
			lines(l,col="red")
			legend("topleft",lty=c(2,1),col="red",legend=c("First","Final"))
		}
#		Make NEW interpolating rule
		f <- approxfun(l, rule=2, ties=list("ordered",mean))
#		Apply NEW trend to individual observations to get vooma weights
		w <- 1/f(mu)^4
		dim(w) <- dim(y)
		colnames(w) <- colnames(y)
		rownames(w) <- rownames(y)
#		Add voom weights to prior weights
		if(PriorWeights) {
			weights <- w * prior.weights
			attr(weights,"arrayweights") <- NULL
		} else {
			weights <- w
		}
		if(SampleWeights) {
			sw <- arrayWeights(y,design,weights=weights,var.design=var.design,var.group=var.group,prior.n=prior.n)
			message("Final sample weights (min/max) ", paste(format(range(sw)),collapse="/") )
			weights <- t(sw * t(weights))
		}
		if(Block) {
			dc <- suppressWarnings(duplicateCorrelation(y,design,block=block,weights=weights))
			correlation <- dc$consensus.correlation
			if(is.na(correlation)) correlation <- 0
			message("Final intra-block correlation  ",format(correlation))
		}
	}
	
#	Final linear model fit with voom weights
	fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights)
	
#	Output
	if(!is.null(y$genes)) {
		fit$genes <- y$genes
	}
	if(SampleWeights) {
		fit$targets <- y$targets
		if(is.null(fit$targets)) {
			fit$targets <- data.frame(sample.weight=sw)
			row.names(fit$targets) <- colnames(y)
		} else {
			fit$targets$sample.weight <- sw
		}
	}
	fit$span <- span
	if(save.plot) {
		fit$voom.xy <- list(x=sx,y=sy,xlab=xlab,ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		fit$voom.line <- l
		fit$voom.line$col <- "red"
	}
	if(keep.EList) {
		y$weights <- weights
		fit$EList <- y
	}
	fit
}
