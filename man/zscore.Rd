\name{zscore}
\alias{zscore}
\alias{zscoreGamma}
\alias{zscoreHyper}

\title{Transform non-normal random deviates to standard normal}

\description{
Compute z-score equivalents of non-normal random deviates.
}

\usage{
zscore(q, distribution, \dots)
zscoreGamma(q, shape, rate = 1, scale = 1/rate) 
zscoreHyper(q, m, n, k) 
}

\arguments{
  \item{q}{numeric vector or matrix giving deviates of a random variable}
  \item{distribution}{character name of probabability distribution for which a cumulative distribution function exists}
  \item{\ldots}{other arguments specify distributional parameters and are passed to the cumulative distribution function}
  \item{shape}{gamma shape parameter (>0)}
  \item{rate}{gamma rate parameter (>0)}
  \item{scale}{gamma scale parameter (>0)}
  \item{m}{as for \code{\link{qhyper}}}
  \item{n}{as for \code{\link{qhyper}}}
  \item{k}{as for \code{\link{qhyper}}}
}

\value{
Numeric vector or matrix of equivalent deviates from the standard normal distribution.
}

\details{
These functions compute the standard normal deviates which have the same quantiles as the given values in the specified distribution.
For example, if \code{z <- zscoreGamma(x, shape, rate)} then \code{pnorm(z)} equals \code{pgamma(x, shape, rate)}.

\code{zscore} works for any distribution for which a cumulative distribution function (like \code{pnorm}) exists in R.
The argument \code{distribution} is the name of the cumulative distribution function with the \code{"p"} removed.

\code{zscoreGamma} and \code{zscoreHyper} are specific functions for the gamma and hypergeometric distributions respectively.

The transformation to z-scores is done by converting to log tail probabilities, and then using \code{qnorm}.
For numerical accuracy, the left or right tail is used, depending on which is likely to be smaller.
}

\author{Gordon Smyth}

\seealso{
\code{\link{zscoreT}}.

\code{zscoreNBinom} in the edgeR package.

\code{\link{qnorm}} in the stats package.
}

\examples{
# These are all equivalent
zscore(c(1,2.5), dist="gamma", shape=0.5, scale=2)
zscore(c(1,2.5), dist="chisq", df=1)
zscoreGamma(c(1,2.5), shape=0.5, scale=2)
}

\keyword{distributions}
