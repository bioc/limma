\name{ebayes}
\alias{ebayes}
\title{Empirical Bayes Statistics for Differential Expression (Deprecated)}
\description{Given a microarray linear model fit, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a common value.
This function is now deprecated.
Please use \code{\link{eBayes}} instead.
}
\usage{
ebayes(fit, proportion=0.01, stdev.coef.lim=c(0.1,4),
       trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
}
\arguments{
  \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit}.
  For \code{ebayes} only, \code{fit} can alternatively be an unclassed list produced by \code{lm.series}, \code{gls.series} or \code{mrlm} containing components \code{coefficients}, \code{stdev.unscaled}, \code{sigma} and \code{df.residual}.}
  \item{proportion}{numeric value between 0 and 1, assumed proportion of genes which are differentially expressed}
  \item{stdev.coef.lim}{numeric vector of length 2, assumed lower and upper limits for the standard deviation of log2-fold-changes for differentially expressed genes}
  \item{trend}{logical, should an intensity-trend be allowed for the prior variance? Default is that the prior variance is constant.}
  \item{robust}{logical, should the estimation of \code{df.prior} and \code{var.prior} be robustified against outlier sample variances?}
  \item{winsor.tail.p}{numeric vector of length 1 or 2, giving left and right tail proportions of \code{x} to Winsorize. Used only when \code{robust=TRUE}.}
}
\value{
A list with components:
  \item{t}{numeric vector or matrix of moderated t-statistics}
  \item{p.value}{numeric vector of p-values corresponding to the t-statistics}
  \item{s2.prior}{estimated prior value for \code{sigma^2}. A vector if \code{covariate} is non-\code{NULL}, otherwise a scalar.}
  \item{df.prior}{degrees of freedom associated with \code{s2.prior}}
  \item{df.total}{numeric vector of total degrees of freedom associated with t-statistics and p-values. Equal to \code{df.prior+df.residual} or \code{sum(df.residual)}, whichever is smaller.}
  \item{s2.post}{numeric vector giving the posterior values for \code{sigma^2}}
  \item{lods}{numeric vector or matrix giving the log-odds of differential expression (natural log scale).}
  \item{var.prior}{estimated prior value for the variance of the log2-fold-change for differentially expressed gene}
}
\details{
This function is now deprecated.
Please use \code{\link{eBayes}} instead.
}
\seealso{
\code{\link{eBayes}}
}
\author{Gordon Smyth}

\references{
Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\url{http://projecteuclid.org/euclid.aoas/1469199900}

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology} 3, Article 3.
\url{http://www.statsci.org/smyth/pubs/ebayes.pdf}
}
