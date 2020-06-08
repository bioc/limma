topTableF <- function(fit,number=10,genelist=fit$genes,adjust.method="BH",sort.by="F",p.value=1,lfc=0)
#	Summary table of top genes by F-statistic
#	Gordon Smyth
#	27 August 2006. Last modified 24 June 2014.
{
#	Deprecated message added 6 June 2020
	message("topTableF is obsolete and will be removed in a future version of limma. Please considering using topTable instead.")

#	Check fit
	if(is.null(fit$coefficients)) stop("Coefficients not found in fit")
	M <- as.matrix(fit$coefficients)
	rn <- rownames(M)
	if(is.null(colnames(M))) colnames(M) <- paste("Coef",1:ncol(M),sep="")
	Amean <- fit$Amean
	Fstat <- fit$F
	Fp <- fit$F.p.value
	if(is.null(Fstat)) stop("F-statistics not found in fit")

#	Ensure genelist is a data.frame
	if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ProbeID=genelist,stringsAsFactors=FALSE)

#	Check rownames
	if(is.null(rn))
		rn <- 1:nrow(M)
	else
		if(anyDuplicated(rn)) {
			if(is.null(genelist))
				genelist <- data.frame(ID=rn,stringsAsFactors=FALSE)
			else
				if("ID" %in% names(genelist))
					genelist$ID0 <- rn
				else
					genelist$ID <- rn
			rn <- 1:nrow(M)
		}

#	Check sort.by
	sort.by <- match.arg(sort.by,c("F","none"))

#	Apply multiple testing adjustment
	adj.P.Value <- p.adjust(Fp,method=adjust.method)

#	Thin out fit by lfc and p.value thresholds
	if(lfc > 0 || p.value < 1) {
		if(lfc>0)
			big <- rowSums(abs(M)>lfc,na.rm=TRUE)>0
		else
			big <- TRUE
		if(p.value<1) {
			sig <- adj.P.Value <= p.value
			sig[is.na(sig)] <- FALSE
		} else
			sig <- TRUE
		keep <- big & sig
		if(!all(keep)) {
			M <- M[keep,,drop=FALSE]
			rn <- rn[keep]
			Amean <- Amean[keep]
			Fstat <- Fstat[keep]
			Fp <- Fp[keep]
			genelist <- genelist[keep,,drop=FALSE]
			adj.P.Value <- adj.P.Value[keep]
		}
	}

#	Enough rows left?
	if(nrow(M) < number) number <- nrow(M)
	if(number < 1) return(data.frame())

#	Find rows of top genes
	if(sort.by=="F")
		o <- order(Fp,decreasing=FALSE)[1:number]
	else
		o <- 1:number

#	Assemble data.frame
	if(is.null(genelist))
		tab <- data.frame(M[o,,drop=FALSE])
	else
		tab <- data.frame(genelist[o,,drop=FALSE],M[o,,drop=FALSE])
	tab$AveExpr <- Amean[o]
	tab <- data.frame(tab,F=Fstat[o],P.Value=Fp[o],adj.P.Val=adj.P.Value[o])
	rownames(tab) <- rn[o]
	tab
}
