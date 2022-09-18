goanaTrend <- function(index.de, covariate, n.prior=10, plot=FALSE, xlab="Covariate Rank", ylab="Probability gene is DE", main="DE status vs covariate")
# Estimate probability of DE given a covariate
# Gordon Smyth
# Created 18 Sep 2022
{
  if(anyNA(index.de)) stop("index.de should not contain missing values")
  ngenes <- length(covariate)
  if(identical(ngenes,0L)) return(numeric(0))
  p <- isDE <- rep_len(0,ngenes)
  isDE[index.de] <- 1
  nDE <- sum(isDE)
  p.mean <- max(nDE/ngenes, 1e-5)
  if(anyNA(covariate)) {
    isna <- is.na(covariate)
    p[isna] <- p.mean
    names(p) <- names(covariate)
    covariate <- covariate[!isna]
    index.de <- as.logical(isDE[!isna])
    p[!isna] <- Recall(index.de, covariate, n.prior=n.prior, plot=plot, xlab=xlab, ylab=ylab, main=main)
    return(p)
  }
  o <- order(covariate)
  isDE.o <- isDE[o]
  span <- approx(x=c(20,200),y=c(1,0.5),xout=sum(isDE),rule=2,ties=list("ordered",mean))$y
  p.o <- tricubeMovingAverage(isDE.o,span=span)
  p.o <- (p.mean*n.prior + p.o*nDE)/(n.prior + nDE)
  p[o] <- p.o
  if(plot) {
    ylim <- c(0, ceiling(max(p)*ngenes)/ngenes )
    x <- 1:ngenes
    plot(x,p.o,type="l",ylim=ylim,xlab=xlab,ylab=ylab,main=main,lwd=2)
    xseg <- x[as.logical(isDE.o)]
    yseg0 <- rep_len(0,nDE)
    yseg1 <- rep_len(ylim[2]/9,nDE)
    if(nDE < 100) lwd <- 2 else lwd <- 1
    segments(xseg,yseg0,xseg,yseg1,lwd=lwd)
  }
  names(p) <- names(covariate)
  invisible(p)
}
