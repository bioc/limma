cameraPR <- function(statistic,...) UseMethod("cameraPR")

cameraPR.default <- function(statistic,index,use.ranks=FALSE,inter.gene.cor=0.01,sort=TRUE,directional=TRUE,...)
#	Competitive gene set test allowing for correlation between genes.
#	This version uses pre-ranked genes (i.e., with a provided vector of test statistics)
#	and uses pre-specified inter-gene correlations for the sets.
#	Gordon Smyth
#	Created 18 Apr 2017. Last modified 4 Oct 2026.
{
#	Issue warning if extra arguments found
	dots <- names(list(...))
	if(length(dots)) warning("Extra arguments disregarded: ",sQuote(dots))

#	Check statistic
	if(is.list(statistic)) stop("statistic should be a numeric vector")
	storage.mode(statistic) <- "numeric"
	if(anyNA(statistic)) stop("NA values for statistic not allowed")
	G <- length(statistic)
	ID <- names(statistic)
	if(G<3) stop("Too few genes in dataset: need at least 3")

#	Check index
	if(!is.list(index) || is.data.frame(index)) index <- list(set1=index)
	nsets <- length(index)

#	Check inter.gene.cor
	if(anyNA(inter.gene.cor)) stop("NA inter.gene.cor not allowed")
	if(max(abs(inter.gene.cor) >= 1)) stop("`inter.gene.cor` must be strictly between -1 and 1")
	if(length(inter.gene.cor) > 1L) {
		if(length(inter.gene.cor) != nsets) stop("Length of `inter.gene.cor` doesn't match number of sets")
		fixed.cor <- FALSE
	} else {
		fixed.cor <- TRUE
		inter.gene.cor <- rep_len(inter.gene.cor,nsets)
	}

#	Check directional
	if(!directional) {
		if(min(statistic) < 0) {
			statistic <- abs(statistic)
			warning("Converting `statistic` to absolute values for non-directional tests",call.=FALSE)
		}
		if(!use.ranks) {
			use.ranks <- TRUE
			warning("Setting `use.ranks=TRUE` for non-directional tests",call.=FALSE)
		}
	}

#	Set df
	if(use.ranks)
		df.camera <- Inf
	else
		df.camera <- G-2L

#	Global statistics
	meanStat <- mean(statistic)
	varStat <- var(statistic)

	NGenes <- Down <- Up <- rep_len(0,nsets)
	for (i in 1:nsets) {
		iset <- index[[i]]
		if(is.data.frame(iset)) {
			if(ncol(iset) < 2) stop("index ",i," is a data.frame but doesn't contain gene weights")
			if(!is.numeric(iset[,2])) stop("index ",i," is a data.frame but second column isn't numeric")
			iw <- iset[,2]
			iset <- iset[,1]
			if(is.factor(iset)) iset <- as.character(iset)
			if(is.character(iset)) {
				if(is.null(ID)) stop("names(statistic) is NULL, so no gene names to match to")
				if(anyDuplicated(iset)) stop("Duplicate gene ids in set ",i)
				iset <- match(iset,ID)
				if(anyNA(iset)) {
					o <- which(!is.na(iset))
					iset <- iset[o]
					iw <- iw[o]
				}
			}
			if(max(iset) > G) stop("set ",i," contains indices too large")
			if(min(iset) < 1) stop("set ",i," contains negative indices")
			NGenes[i] <- length(iset)
			if(length(iset)) {
				x <- rep_len(0,G)
				x[iset] <- iw
				Z <- .invSqrtCorTimes(cbind(1,x,statistic),sign(x),inter.gene.cor[i])
				fit <- lm.fit(Z[,1:2],Z[,3])
				s2 <- mean(fit$effects[3:G]^2)
				t.stat <- sign(fit$qr$qr[2,2])*fit$effects[2]/sqrt(s2)
				Down[i] <- pt(t.stat,df=G-2L)
				Up[i] <- pt(t.stat,df=G-2L,lower.tail=FALSE)
			} else {
				Down[i] <- 1
				Up[i] <- 1
			}
		} else {
			if(is.character(iset)) iset <- which(ID %in% iset)
			StatInSet <- statistic[iset]
			m <- length(StatInSet)
			NGenes[i] <- m
			if(use.ranks) {
				p.value <- rankSumTestWithCorrelation(iset,statistics=statistic,correlation=inter.gene.cor[i],df=df.camera)
				Down[i] <- p.value[1]
				Up[i] <- p.value[2]
			} else {	
				vif <- 1+(m-1)*inter.gene.cor[i]
				m2 <- G-m
				meanStatInSet <- mean(StatInSet)
				delta <- G/m2*(meanStatInSet-meanStat)
				varStatPooled <- ( (G-1L)*varStat - delta^2*m*m2/G ) / (G-2L)
				two.sample.t <- delta / sqrt( varStatPooled * (vif/m + 1/m2) )
				Down[i] <- pt(two.sample.t,df=df.camera)
				Up[i] <- pt(two.sample.t,df=df.camera,lower.tail=FALSE)
			}
		}
	}
	TwoSided <- 2*pmin(Down,Up)

#	Assemble into data.frame
	if(directional) {
		D <- (Down < Up)
		Direction <- rep_len("Up",nsets)
		Direction[D] <- "Down"
		if(fixed.cor)
			tab <- data.frame(NGenes=NGenes,Direction=Direction,PValue=TwoSided,stringsAsFactors=FALSE)
		else
			tab <- data.frame(NGenes=NGenes,Correlation=inter.gene.cor,Direction=Direction,PValue=TwoSided,stringsAsFactors=FALSE)
	} else {
		if(fixed.cor)
			tab <- data.frame(NGenes=NGenes,PValue=Up,stringsAsFactors=FALSE)
		else
			tab <- data.frame(NGenes=NGenes,Correlation=inter.gene.cor,PValue=Up,stringsAsFactors=FALSE)
	}
	rownames(tab) <- names(index)

#	Add FDR
	if(nsets>1L) tab$FDR <- p.adjust(tab$PValue,method="BH")

#	Sort by p-value
	if(sort && nsets>1L) {
		o <- order(tab$PValue)
		tab <- tab[o,]
	}

	tab
}

.invSqrtCorTimes <- function(B, u, rho)
#	Compute V^(-1/2) %*% B where
#	V = (1-rho)*I + rho*u u^T + rho*diag(1-u^2),
#	u has elements 0, 1 or -1, and V^(1/2) is the symmetric square root.
#	Algorithm is an adaption of the Sherman-Morrison formula.
#	Gordon Smyth and Claude.
#	Created 4 Oct 2026. Last modified 5 Oct 2026.
{
	B <- as.matrix(B)
	if(length(u) != nrow(B)) stop("length(u) must equal nrow(B)")
	i <- which(u != 0)
	k <- length(i)
	if(k < 2L) return(B)
	if(rho >= 1) stop("rho must be less than 1")
	if(rho <= -1/(k-1)) stop("V is not positive definite")
	s <- u[i]
	sa <- sqrt(1 - rho)
	sd <- sqrt(1 + (k-1)*rho)
	gamma <- -rho / (sa * sd * (sa + sd))
	Bs <- B[i, , drop=FALSE]
	stB <- crossprod(s, Bs)
	B[i, ] <- Bs / sa + outer(gamma * s, stB[1, ])
	B
}
