contrastAsCoef <- function(design, contrast=NULL, first=TRUE)
#	Reform a design matrix so that one or more contrasts become simple coefficients
#	Gordon Smyth
#	31 August 2013. Last modified 1 Oct 2025.
{
	design <- as.matrix(design)
	if(is.null(contrast)) return(design)
	contrast <- as.matrix(contrast)
	if(ncol(design) != nrow(contrast)) stop("Length of contrast doesn't match ncol(design)")
	QRc <- qr(contrast)
	ncontrasts <- QRc$rank
	if(identical(ncontrasts,0L)) stop("contrast is all zero")
	coef <- 1:ncontrasts
	designT <- qr.qty(QRc,t(design))
	designT[coef,] <- backsolve(QRc$qr,designT,k=ncontrasts)
	design <- t(designT)
	colnames(design) <- paste0("Q",1:ncol(design))
	cn <- colnames(contrast)
	if(is.null(cn)) cn <- paste0("C",QRc$pivot[coef])
	colnames(design)[coef] <- cn
	if(!first) {
		design <- cbind(design[,-coef,drop=FALSE],design[,coef,drop=FALSE])
		coef <- rev( ncol(design)-coef+1 )
	}
	list(design=design,coef=coef,qr=QRc)
}
