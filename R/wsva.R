wsva <- function(y, design, n.sv=1L, weight.by.sd=FALSE, ...)
#	Weighted surrogate variable analysis
#	Yifang Hu and Gordon Smyth
#	Created 26 Nov 2015.  Last modified 17 Aug 2016.
{
	ngenes <- nrow(y)
	narrays <- ncol(y)
	p <- ncol(design)
	for(i in 1L:n.sv) {
		Effects <- .lmEffects(y, design, ...)[,-1L]
		if(weight.by.sd) {
			s <- sqrt(rowMeans(Effects^2))
			Effects <- s * Effects
		}
		u <- drop(svd(Effects,nu=1L,nv=0L)$u)
		if(weight.by.sd) u <- u*s
		sv <- colSums(u*y)
		design <- cbind(design, sv)
	}
	j <- 1L:n.sv
	colnames(design)[p+j] <- paste0("SV",j)
	design
}
