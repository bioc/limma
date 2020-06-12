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
#  1 October 2003. Last modified 9 June 2020.
{
	z <- q
	n <- length(q)
	shape <- rep_len(shape,length.out=n)
	scale <- rep_len(scale,length.out=n)
	up <- (q > shape*scale)
	if(any(up)) z[up] <- qnorm(pgamma(q[up],shape=shape[up],scale=scale[up],lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE)
	if(any(!up)) z[!up] <- qnorm(pgamma(q[!up],shape=shape[!up],scale=scale[!up],lower.tail=TRUE,log.p=TRUE),lower.tail=TRUE,log.p=TRUE)
	z
}

zscoreT <- function(x, df, approx=FALSE, method="bailey")
#  Z-score equivalents of t distribution deviates
#  Gordon Smyth
#  Created 24 August 2003 with exact quantile method only.
#  Hill method added 3 June 2014.
#  Bailey and Wallace methods added 21 July 2019.
{
	if(approx) {
		df <- pmin(df,1e100)
		method <- match.arg(method,c("bailey","hill","wallace"))
		if(method=="bailey") return(.zscoreTBailey(x=x,df=df))
		if(method=="hill") return(.zscoreTHill(x=x,df=df))
		if(method=="wallace") return(.zscoreTWallace(x=x,df=df))
	} else {
		return(.zscoreTQuantile(x=x,df=df))
	}
}

.zscoreTQuantile <- function(x, df)
#  Z-score equivalents of t distribution deviates using an approximatiuon from Wallace (1959).
#  Wallace, D. L. (1959). Bounds on normal approximations to Student's and the chi-square distributions. The Annals of Mathematical Statistics, 30(4), 1121-1130.
#  Gordon Smyth
#  Created 21 July 2019 by modifying zscoreT code written 24 August 2003.
{
	qnorm(pt(abs(x),df=df,lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE) * sign(x)
}

.zscoreTWallace <- function(x, df)
#  Z-score equivalents of t distribution deviates using an approximatiuon from Wallace (1959).
#  Wallace, D. L. (1959). Bounds on normal approximations to Student's and the chi-square distributions. The Annals of Mathematical Statistics, 30(4), 1121-1130.
#  Gordon Smyth
#  Created 16 July 2019.
{
	((df+0.125)/(df+0.375)) * sqrt(df*log1p(x/df*x)) * sign(x)
}

.zscoreTBailey <- function(x, df)
#  Z-score equivalents of t distribution deviates using an approximatiuon from Wallace (1959).
#  Wallace, D. L. (1959). Bounds on normal approximations to Student's and the chi-square distributions. The Annals of Mathematical Statistics, 30(4), 1121-1130.
#  Gordon Smyth
#  Created 16 July 2019.
{
	((df+0.125)/(df+1.125)) * sqrt((df+19/12)*log1p(x/(df+1/12)*x)) * sign(x)
}

.zscoreTHill <- function(x, df)
#  Z-score equivalents for t distribution deviates using Hill's 1970 approximation:
#  Hill, G. W. (1970). Algorithm 396: Student's t-quantiles. Communications of the ACM, 13(10), 619-620.
#  The approx requires df > 0.5 and gives good accuracy for df >= 2.
#  Gordon Smyth
#  Created 3 June 2014. Last modified 21 July 2019.
{
	A <- df-0.5
	B <- 48*A*A
	z <- A*log1p(x/df*x)
	z <- (((((-0.4*z-3.3)*z-24)*z-85.5)/(0.8*z*z+100+B)+z+3)/B+1)*sqrt(z)
	z * sign(x)
}

tZscore <- function(z, df)
#  t-statistic equivalents of z-score deviates
#  Gordon Smyth
#  Created 1 June 2004. Modified 21 July 2019.
{
	qt(pnorm(abs(z),lower.tail=FALSE,log.p=TRUE),df=df,lower.tail=FALSE,log.p=TRUE) * sign(z)
}
