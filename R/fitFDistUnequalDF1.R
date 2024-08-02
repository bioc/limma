fitFDistUnequalDF1 <- function(x,df1,covariate=NULL,robust=FALSE,prior.weights=NULL)
# Robust estimation of the parameters of a scaled F-distribution given df1.
# This version gives special attention to the possibility that df1 may vary
# substantially between observations.
# Gordon Smyth and Lizhong Chen
# Created 18 Jul 2024. Last modified 31 Jul 2024.
{
  n <- length(x)

# Check df1. Can be unit vector.
  if(all(length(df1) != c(1L,n))) stop("x and df1 are different lengths")
  if(anyNA(df1)) stop("NA df1 values")

# Check covariate
  if(!is.null(covariate)) {
    if(!identical(length(covariate),n)) stop("x and covariate are different lengths")
    if(anyNA(covariate)) stop("prior.weights contain NA values")
  }

# Check prior.weights
  if(!is.null(prior.weights)) {
    if(!identical(length(prior.weights),n)) stop("x and covariate are different lengths")
    m <- min(prior.weights)
    if(is.na(m)) stop("prior.weights contain NA values")
    if(m < 0) stop("prior.weights are negative")
  }

# Check for NAs in x
  if(anyNA(x)) {
    i <- is.na(x)
    if(is.null(prior.weights))
      prior.weights <- as.numeric(!i)
    else
      prior.weights[i] <- 0
    x[i] <- 0
   }

# Treat small df1 values as un-informative
  if(min(df1) < 0.01) {
    i <- (df1 < 0.01)
    if(is.null(prior.weights))
      prior.weights <- as.numeric(!i)
    else
      prior.weights[i] <- 0
    df1[i] <- 1
  }

  PriorWeights <- !is.null(prior.weights)

# Check there are some informative x values
  i <- (x>0)
  if(PriorWeights) i[prior.weights==0] <- FALSE
  n.informative <- sum(i)
  if(n.informative<2) return(list(scale=NA_real_,df2=NA_real_))
  if(n.informative==2) {
    covariate <- NULL
    robust <- FALSE
    prior.weights <- NULL
  }

# Avoid exactly zero x values for moment estimation
  m <- median(x[i])
  xpos <- pmax(x, 1e-12 * m)

# Work on with log(F)
  z <- log(xpos)

# Average log(F) adjusted for d1
  d1 <- df1/2
  e <- z+logmdigamma(d1)
  w <- 1/trigamma(d1)
  if(length(w) < n) w <- rep_len(w,n)
  if(PriorWeights) w <- w*prior.weights
  if(is.null(covariate)) {
    emean <- sum(w*e)/sum(w)
  } else {
    span <- chooseLowessSpan(n,small.n=500)
    emean <- loessFit(e, covariate, weights=w/quantile(w,probs=0.75), min.weight=1e-8, max.weight=1e2, span=span, iterations=1)$fitted
  }

# Log-likelihood function
  d1x <- d1*xpos
  if(PriorWeights) {
    minusTwiceLogLik <- function(par) {
      d2 <- par/(1-par)
      d2s20 <- d2*exp(emean-logmdigamma(d2))
      -2*sum(prior.weights*(-(d1+d2)*log1p(d1x/d2s20)-d1*log(d2s20)+lgamma(d1+d2)-lgamma(d2)))
    }
  } else {
    minusTwiceLogLik <- function(par) {
      d2 <- par/(1-par)
      d2s20 <- d2*exp(emean-logmdigamma(d2))
      -2*sum(-(d1+d2)*log1p(d1x/d2s20)-d1*log(d2s20)+lgamma(d1+d2)-lgamma(d2))
    }
  }

# Optimization
  out <- optimize(minusTwiceLogLik, c(1/2,0.9998))
  d2 <- out$minimum/(1-out$minimum)
  s20 <- exp(emean-logmdigamma(d2))

# Finish here if robust=FALSE
  if(!robust) return(list(scale=s20,df2=2*d2))

# Use FDR to identify two-sided outliers
# FDR is approximate probability of not being an outlier
  df2 <- 2*d2
  FStat <- x/s20
  RightP <- pf(FStat,df1=df1,df2=df2,lower.tail=FALSE)
  LeftP <- 1-RightP
  if(min(LeftP) < 0.001) {
    i <- (LeftP < 0.001)
    LeftP[i] <- pf(FStat[i],df1=.subsetUnitVector(df1,i),df2=df2,lower.tail=TRUE)
  }
  TwoSidedP <- 2*pmin(LeftP,RightP)
  FDR <- p.adjust(TwoSidedP, method="BH")
  FDR[FDR > 0.3] <- 1

# If no outliers, return non-robust estimates
  if(identical(min(FDR),1)) return(list(scale=s20,df2=df2))

# Refit F-distribution with FDR as prior weights
  outpw <- Recall(x=x,df1=df1,covariate=covariate,robust=FALSE,prior.weights=FDR)
  s20 <- outpw$scale
  df2 <- outpw$df2

# Use qqplot-type method to identify right outliers
  r <- rank(FStat)
  UniformP <- (n-r+0.5)/n
  ProbNotOutlier <- pmin(RightP/UniformP,1)

# If no right outliers, return robust estimates without df2 shrinkage
  if(identical(min(ProbNotOutlier),1)) return(outpw)

# Posterior df for right outliers
# Find df2.outlier to make maxFStat the median of the distribution
# Exploit fact that LogRightP is nearly linear with positive 2nd deriv as a function of df2
  i <- which.min(RightP)
  minRightP <- RightP[i]
  if(identical(minRightP,0)) {
    df2.outlier <- 0
    df2.shrunk <- ProbNotOutlier*df2
  } else {
    df2.outlier <- log(0.5)/log(minRightP)*df2
#   Iterate for accuracy
    NewLogRightP <- pf(FStat[i],df1=.subsetUnitVector(df1,i),df2=df2.outlier,lower.tail=FALSE,log.p=TRUE)
    df2.outlier <- log(0.5)/NewLogRightP*df2.outlier
    df2.shrunk <- ProbNotOutlier*df2+(1-ProbNotOutlier)*df2.outlier
  }

# Force df2.shrunk to be monotonic in RightP
  o <- order(RightP)
  df2.ordered <- df2.shrunk[o]
  m <- cumsum(df2.ordered)
  m <- m/(1:n)
  imin <- which.min(m)
  df2.ordered[1:imin] <- m[imin]
  df2.shrunk[o] <- cummax(df2.ordered)

  list(scale=s20,df2=df2,df2.outlier=df2.outlier,df2.shrunk=df2.shrunk)
}

.subsetUnitVector <- function(x,i)
# Subset vector treating unit vector as vector of equal values
{
  if(identical(length(x),1L)) x else x[i]
}
