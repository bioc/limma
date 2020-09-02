#  OUTPUT

write.fit <- function(fit, results=NULL, file, digits=NULL, adjust="none", method="separate", F.adjust="none", quote=FALSE, sep="\t", row.names=TRUE, ...)
#	Write an MArrayLM fit to a file
#	Gordon Smyth
#	14 Nov 2003.  Last modified 2 Sep 2020.
{
	if(!is(fit, "MArrayLM")) stop("fit should be an MArrayLM object")
	if(!is.null(results) && !is(results,"TestResults")) stop("results should be a TestResults object")
	if(is.null(fit$t) || is.null(fit$p.value)) fit <- eBayes(fit)
	method <- match.arg(method, c("separate","global"))

	p.value <- as.matrix(fit$p.value)
	if(adjust=="none") {
		p.value.adj <- NULL
	} else {
		p.value.adj <- p.value
		if(method=="separate") for (j in 1:ncol(p.value)) p.value.adj[,j] <- p.adjust(p.value[,j],method=adjust)
		if(method=="global") p.value.adj[] <- p.adjust(p.value,method=adjust)
	}
	if(F.adjust=="none" || is.null(fit$F.p.value))
		F.p.value.adj <- NULL
	else
		F.p.value.adj <- p.adjust(fit$F.p.value,method=F.adjust)

#	Optionally, round results for easy reading
	if(is.null(digits)) {
		rn <- function(x,digits=digits) x
	} else {
		rn <- function(x,digits=digits)
				if(is.null(x))
					NULL
				else {
					if(is.matrix(x) && ncol(x)==1) x <- x[,1]
					round(x,digits=digits)
				}
	}

#	Prepare output data.frame
	tab <- list()
	tab$A <- rn(fit$Amean,digits=digits-1)
	tab$Coef <- rn(fit$coef,digits=digits)
	tab$t <- rn(fit$t,digits=digits-1)
	tab$p.value <- rn(p.value,digits=digits+2)
	tab$p.value.adj <- rn(p.value.adj,digits=digits+3)
	tab$F <- rn(fit$F,digits=digits-1)
	tab$F.p.value <- rn(fit$F.p.value,digits=digits+2)
	tab$F.p.value.adj <- rn(F.p.value.adj,digits=digits+3)
	tab$Res <- unclass(results)
	tab$Genes <- fit$genes
	tab <- data.frame(tab,check.names=FALSE)

#	Unlike write.table, the row.names argument must be a logical value
	if(is.character(row.names)) {
		warning("attempt to set new row.names ignored. row.names argument should be TRUE or FALSE.")
		row.names <- TRUE
	} else {
		if(!is.logical(row.names)) stop("row.names should be logical value")
	}

#	If row.names=TRUE but fit doesn't contain row.names, then override argument and issue warning
	if(is.null(row.names(fit))) {
		if(row.names) {
			warning("fit doesn't contain row.names")
			row.names <- FALSE
		}
	} else {
		row.names(tab) <- row.names(fit)
	}

#	This function treats col.names similarly to write.csv for col.names
#	and ensures a blank column name for the row.names column if present.
    Call <- match.call(expand.dots = TRUE)
    if (!is.null(Call[["col.names"]])) warning("attempt to set 'col.names' ignored")
    if (!is.null(Call[["qmethod"]])) warning("attempt to set 'qmethod' ignored")
	if(row.names) col.names <- NA else col.names=TRUE

	write.table(tab,file=file,quote=quote,row.names=row.names,col.names=col.names,sep=sep,qmethod="double",...)
}

