#  VENN DIAGRAM COUNTS AND PLOTS

vennCounts <- function(classification,include="both") {
#	Venn diagram counts
#	Gordon Smyth
#	4 July 2003.  Last modified 15 September 2003.

	if(is(classification,"classification")) classification <- classification$classification
	include <- match.arg(include,c("both","up","down"))
	classification <- switch(include,
		both = as.matrix(abs(classification) > 0.5),
		up = as.matrix(classification > 0.5),
		down = as.matrix(classification < -0.5)
	)
	ngenes <- nrow(classification)
	ncontrasts <- ncol(classification)
	names <- colnames(classification)
	if(is.null(names)) names <- paste("Group",1:ncontrasts)
	noutcomes <- 2^ncontrasts
	outcomes <- matrix(0,noutcomes,ncontrasts)
	colnames(outcomes) <- names
	for (j in 1:ncontrasts)
		outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
	counts <- rep(0,noutcomes)
	one <- rep(1,ngenes)
	for (i in 1:noutcomes)
		counts[i] <- sum(apply(classification==(one%*%outcomes[i,,drop=FALSE]),1,all),na.rm=TRUE)
	structure(cbind(outcomes,Counts=counts),class="VennCounts")
}

vennDiagram <- function(object,include="both",names,...) {
#	Plot Venn diagram
#	Gordon Smyth and James Wettenhall
#	4 July 2003.  Last modified 4 September 2003.

	if(class(object) != "VennCounts") object <- vennCounts(object,include=include)
	nsets <- ncol(object)-1
	if(nsets > 3) stop("Can't plot Venn diagram for more than 3 sets")
	if(missing(names)) names <- colnames(object)[1:nsets]
	counts <- object[,"Counts"]
	theta <- 2*pi*(1:360)/360
	xcentres <- list(0,c(-1,1),c(-1,1,0))[[nsets]]
	ycentres <- list(0,c(0,0),c(1/sqrt(3),1/sqrt(3),-2/sqrt(3)))[[nsets]]
	r <- c(1.5,1.5,1.5)[nsets]
	xtext <- list(-1.2,c(-1.2,1.2),c(-1.2,1.2,0))[[nsets]]
	ytext <- list(1.8,c(1.8,1.8),c(2.4,2.4,-3))[[nsets]]
	cex <- 1.5
	plot(x=0,y=0,type="n",xlim=c(-4,4),ylim=c(-4,4),xlab="",ylab="",axes=FALSE,...)
	for(circle in 1:nsets) {
		lines(xcentres[circle]+r*cos(theta),ycentres[circle]+r*sin(theta))
		text(xtext[circle],ytext[circle],names[circle],cex=cex)
	}
	switch(nsets,
		{
			rect(-3,-2.5,3,2.5)
			text(2.3,-2.1,counts[1],cex=cex)
			text(0,0,counts[2],cex=cex)
		}, {
			rect(-3,-2.5,3,2.5)
		    text(2.3,-2.1,counts[1],cex=cex)
			text(1.5,0.1,counts[2],cex=cex)
			text(-1.5,0.1,counts[3],cex=cex)
			text(0,0.1,counts[4],cex=cex)
		}, {
			rect(-3,-3.5,3,3.3)
			text(2.5,-3,counts[1],cex=cex)
			text(0,-1.7,counts[2],cex=cex)
			text(1.5,1,counts[3],cex=cex)
			text(.75,-.35,counts[4],cex=cex)
			text(-1.5,1,counts[5],cex=cex)
			text(-.75,-.35,counts[6],cex=cex)
			text(0,.9,counts[7],cex=cex)
			text(0,0,counts[8],cex=cex)
		}
	)
	invisible()
}
