#  SCORE.R

zscore <- function(q, distribution=NULL, ...) 
#  Z-score equivalents for deviates from specified distribution
#  Gordon Smyth
#  13 June 2012
{
	z <- q
	n <- length(q)
	pdist <- get(paste("p",as.character(distribution),sep=""))
	pupper <- pdist(q,...,lower.tail=FALSE,log.p=TRUE)
	plower <- pdist(q,...,lower.tail=TRUE,log.p=TRUE)
	up <- pupper<plower
	if(any(up)) z[up] <- qnorm(pupper[up],lower.tail=FALSE,log.p=TRUE)
	if(any(!up)) z[!up] <- qnorm(plower[!up],lower.tail=TRUE,log.p=TRUE)
	z
}

zscoreGamma <- function(q, shape, rate = 1, scale = 1/rate) 
#  Z-score equivalents for gamma deviates
#  Gordon Smyth
#  1 October 2003
{
	z <- q
	n <- length(q)
	shape <- rep(shape,length.out=n)
	scale <- rep(scale,length.out=n)
	up <- (q > shape*scale)
	if(any(up)) z[up] <- qnorm(pgamma(q[up],shape=shape[up],scale=scale[up],lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE)
	if(any(!up)) z[!up] <- qnorm(pgamma(q[!up],shape=shape[!up],scale=scale[!up],lower.tail=TRUE,log.p=TRUE),lower.tail=TRUE,log.p=TRUE)
	z
}

zscoreT <- function(x, df, approx=FALSE)
#  Z-score equivalents for t distribution deviates
#  Gordon Smyth
#  Created 24 August 2003. Last modified 15 July 2019.
{
	if(length(x)==0L) return(x)

	if(length(df)==1L) {
		if(is.na(df)) {
			x[] <- NA
			return(x)
		}
		if(approx) {
			if(df > 1e300) return(x)
			if(df > 1e5 || df < 1) {
				z <- .zscoreTWallace(x=x,df=df)
			} else {
				z <- .zscoreTHill(x=x,df=df)
			}
		} else {
			z <- x
			pos <- (x>0) & is.finite(x)
			neg <- (x<0) & is.finite(x)
			z[pos] <- qnorm(pt(x[pos],df=df,lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE) 
			z[neg] <- qnorm(pt(x[neg],df=df,lower.tail=TRUE,log.p=TRUE),lower.tail=TRUE,log.p=TRUE)
		}
		return(z)
	}

	if(length(df) != length(x)) stop("length of df doesn't match length of x")
	if(approx) {
		if(anyNA(df)) {
			i <- is.na(df)
			x[i] <- NA
			df[i] <- Inf
		}
		z <- x
		VBig <- (df > 1e300)
		Extreme <- ((df > 1e5) & !VBig) | (df < 1)
		Mid <- !(Extreme | VBig) 
		z[Extreme] <- .zscoreTWallace(x=x[Extreme],df=df[Extreme])
		z[Mid] <- .zscoreTHill(x=x[Mid],df=df[Mid])
	} else {
		z <- x
		pos <- (x>0) & is.finite(x)
		neg <- (x<0) & is.finite(x)
		z[pos] <- qnorm(pt(x[pos],df=df[pos],lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE) 
		z[neg] <- qnorm(pt(x[neg],df=df[neg],lower.tail=TRUE,log.p=TRUE),lower.tail=TRUE,log.p=TRUE)
	}
	z
}

.zscoreTWallace <- function(x, df)
#  Z-score equivalents for t distribution deviates using an upper bound from Wallace (1959).
#  Wallace, D. L. (1959). Bounds on normal approximations to Student's and the chi-square distributions. The Annals of Mathematical Statistics, 30(4), 1121-1130.
#  Gordon Smyth
#  Created 16 July 2019.
{
	z <- sqrt(df*log1p(x/df*x))
	z[x<0] <- -z[x<0]
	z
}

.zscoreTHill <- function(x, df)
#  Z-score equivalents for t distribution deviates using Hill's 1970 approximation:
#  Hill, G. W. (1970). Algorithm 396: Student's t-quantiles. Communications of the ACM, 13(10), 619-620.
#  The approx requires df > 0.5 and is best for large df.
#  Gordon Smyth
#  Created 3 June 2014. Last modified 15 July 2019.
{
	A <- df-0.5
	B <- 48*A*A
	z <- A*log1p(x/df*x)
	z <- (((((-0.4*z-3.3)*z-24)*z-85.5)/(0.8*z*z+100+B)+z+3)/B+1)*sqrt(z)
	z[x<0] <- -z[x<0]
	z
}

tZscore <- function(x, df)
#  t-statistics equivalents for z-scores deviates
#  Gordon Smyth
#  1 June 2004
{
	z <- x
	df <- rep(df,length.out=length(x))
	pos <- x>0
	if(any(pos)) z[pos] <- qt(pnorm(x[pos],lower.tail=FALSE,log.p=TRUE),df=df[pos],lower.tail=FALSE,log.p=TRUE) 
	if(any(!pos)) z[!pos] <- qt(pnorm(x[!pos],lower.tail=TRUE,log.p=TRUE),df=df[!pos],lower.tail=TRUE,log.p=TRUE)
	z
}
