vooma <- function(y,design=NULL,block=NULL,correlation,predictor=NULL,span=NULL,legacy.span=FALSE,plot=FALSE,save.plot=FALSE)
#	Linear modelling of continuous log-expression data with mean-variance modeling at the observational level.
#	Analogous to voom() but for non-count data.
#	y must not contain NAs.
#	Gordon Smyth, Charity Law, Mengbo Li.
#	Created 31 July 2012.  Last modified 13 Feb 2024.
{
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

#	Fit linear model
	if(is.null(block)) {
		fit <- lm.fit(design,t(y$E))
		mu <- fit$fitted.values
	} else {
		block <- as.vector(block)
		if(length(block)!=narrays) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,narrays,nblocks)==matrix(ub,narrays,nblocks,byrow=TRUE)
		cormatrix <- Z%*%(correlation*t(Z))
		diag(cormatrix) <- 1
		cholV <- chol(cormatrix)
		z <- backsolve(cholV,t(y$E),transpose=TRUE)
		X <- backsolve(cholV,design,transpose=TRUE)
		fit <- lm.fit(X,z)
		mu <- crossprod(cholV,fit$fitted.values)
	}
	mu <- t(mu)
	s2 <- colMeans(fit$effects[-(1:fit$rank),,drop=FALSE]^2)

#	Prepare to predict sqrt-standard-deviations by ave log intensity
	sx <- rowMeans(y$E)
	sy <- sqrt(sqrt(s2))

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
		xlab <- "Average log-expression"
		main.title <- "vooma mean-variance trend"
	}

#	Choose span based on the number of genes
	if(is.null(span))
		if(legacy.span)
			if(ngenes<=10) span <- 1 else span <- 0.3+0.7*(10/ngenes)^0.5
		else
			if(ngenes<=50) span <- 1 else span <- 0.3+0.7*(50/ngenes)^0.4

#	Fit lowess trend
	l <- lowess(sx,sy,f=span)
	if(plot) {
		plot(sx,sy,xlab=xlab,ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title(main.title)
		lines(l,col="red")
	}

#	Make interpolating rule
	f <- approxfun(l, rule=2, ties=list("ordered",mean))

#	Apply trend to individual observations
	w <- 1/f(mu)^4
	dim(w) <- dim(y)
	colnames(w) <- colnames(y)
	rownames(w) <- rownames(y)

#	Output
	y$weights <- w
	y$design <- design
	y$span <- span
	if(save.plot) {
		fit$voom.xy <- list(x=sx,y=sy,xlab=xlab,ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		fit$voom.line <- l
		fit$voom.line$col <- "red"
	}
	y
}

voomaByGroup <- function(y,group,design=NULL,block=NULL,correlation,
  span=NULL,legacy.span=FALSE,plot=FALSE,col=NULL,lwd=1,
  pch=16,cex=0.3,alpha=0.5,legend="topright")
#	Vooma by group
#	Linear modelling of microarray data with mean-variance modeling
#	at the observational level by fitting group-specific trends.
#	Creates an EList object for entry to lmFit() etc in the limma pipeline.
#	Charity Law and Gordon Smyth
#	Created 13 Feb 2013.  Modified 13 Sept 2024.
{
#	Check y
	if(!is(y,"EList")) y <- new("EList",list(E=as.matrix(y)))
	ngenes <- nrow(y)
	narrays <- ncol(y)

#	Check group
	group <- as.factor(group)
	intgroup <- as.integer(group)
	levgroup <- levels(group)
	ngroups <- length(levgroup)

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design)) design <- model.matrix(~group)

#	Check color
	if(is.null(col))
		if(ngroups==2L)
			col <- c("red","blue")
		else
			col <- 1L+1L:ngroups
	col <- rep_len(col,ngroups)

	w <- y$E
	sx <- sy <- matrix(0, nrow=nrow(y), ncol=nlevels(group))
	colnames(sx) <- levgroup
	rownames(sx) <- rownames(y)
	for (lev in 1L:ngroups) {
		i <- intgroup==lev
		yi <- y[,i]
		designi <- design[i,,drop=FALSE]
		voomi <- vooma(y=yi,design=designi,correlation=correlation, block=block[i], plot=FALSE, span=span, legacy.span=legacy.span)
		w[,i] <- voomi$weights
		sx[,lev] <- voomi$meanvar.trend$x
		sy[,lev] <- voomi$meanvar.trend$y	
	}
	span <- voomi$span

# 	Voom plot	
	if(plot) {
		RGB <- col2rgb(col)/255
		plot(sx,sy,xlab="Average log2 expression",ylab="Sqrt( standard deviation )",main="voom: Mean-variance trend",type="n")
		for (lev in 1:nlevels(group)) {
			coli.transparent <- rgb(RGB[1,lev],RGB[2,lev],RGB[3,lev],alpha=alpha)
			points(sx[,lev],sy[,lev],pch=pch,cex=cex,col=coli.transparent)
			l <- lowess(sx[,lev],sy[,lev],f=span)
			lines(l,col=col[lev],lwd=lwd)
		}
		if(is.character(legend)) legend(legend, levels(group), col=col, lty=1, lwd=lwd)
	}

# 	Output	
	y$meanvar.trend <- list(x=sx,y=sy)
	y$weights <- w
	y$design <- design
	y$span <- span
	y
}
