cumOverlap <- function(ol1, ol2)
#	Cumulative overlap analysis
#	Di Wu and Gordon Smyth
#	Createdin 2007. Last modified 16 April 2010
#	Based on hypergeometric distribution. Test whether two ordered lists of genes are significantly overlapped. 
{
#	Check for duplicates
	if(anyDuplicated(ol1)) stop("Duplicate IDs found in ol1")
	if(anyDuplicated(ol2)) stop("Duplicate IDs found in ol2")

#	Reduce to IDs found in both lists
	m1 <- match(ol1,ol2)
	redo <- FALSE
	if(anyNA(m1)) {
		ol1 <- ol1[!is.na(m1)]
		redo <- TRUE
	}
	m2 <- match(ol2,ol1)
	if(anyNA(m2)) {
		ol2 <- ol2[!is.na(m2)]
		redo <- TRUE
	}

#	Match ol1 to ol2
	if(redo) m1 <- match(ol1,ol2)

#	Count overlaps
	ngenes <- length(ol1)
	i <- 1L:ngenes
	inoverlap <- m1 <= i
	noverlap <- cumsum(inoverlap)

#	Hypergeometric p-valules
	p <- phyper(noverlap-0.5,m=i,n=ngenes-i,k=i,lower.tail=FALSE)

#	Bonferroni
	p.b <- p*i
	nmin <- which.min(p.b)
	p.b <- pmin(p.b,1)

#	Which are ids contribute to overlap
	idoverlap <- ol1[(1:nmin)[inoverlap[1:nmin]]]

	list(n.min=nmin,p.min=p.b[nmin],n.overlap=noverlap,id.overlap=idoverlap,p.value=p,adj.p.value=p.b)
}
