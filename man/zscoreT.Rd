\name{zscoreT}
\alias{zscoreT}
\alias{tZscore}

\title{Transform t-statistics to standard normal}

\description{
Compute z-score equivalents of t-distributed random deviates.
}

\usage{
zscoreT(x, df, approx=FALSE, method = "bailey")
tZscore(z, df)
}

\arguments{
  \item{x}{numeric vector or matrix of values from a t-distribution.}
  \item{df}{degrees of freedom (>0) of the t-distribution.}
  \item{approx}{logical. If \code{TRUE} then a fast approximation is used otherwise exact z-scores are computed.}
  \item{method}{
    character string specifying transformation to be used when \code{approx=TRUE},
    options being \code{"bailey"}, \code{"hill"} or \code{"wallace"}.
  }
  \item{z}{numeric vector or matrix of values from the standard normal distribution.}
}

\value{
Numeric vector or matrix of z-scores or t-distribution deviates.
}

\details{
\code{zscoreT} transforms t-distributed values to standard normal.
Each value is converted to the equivalent quantile of the normal distribution so that
if \code{z <- zscoreT(x, df=df)} then \code{pnorm(z)} equals \code{pt(x, df=df)}.

\code{tZscore} is the inverse of \code{zscoreT} and computes t-distribution equivalents of standard normal deviates.

If \code{approx=FALSE}, the transformation is done by converting to log tail probabilities using \code{pt} or \code{pnorm} and then converting back to quantiles using \code{qnorm} or \code{qt}.
For numerical accuracy, the smaller of the two tail probabilities is used for each deviate.

If \code{approx=TRUE}, then an approximate closed-form transformation is used to convert t-statistics to z-scores directly without computing tail probabilities.
The \code{method} argument provides a choice of three transformations.
\code{method="bailey"} is equation (5) of Bailey (1980) or equation (7) of Brophy (1987).
\code{method="hill"} is from Hill (1970) as given by equation (5) of Brophy (1987).
\code{method="wallace"} is from Wallace (1959) as given by equation equation (2) of Brophy (1987).
Bailey's transformation is a modification of Wallace's approximation.
The Hill approximation is generally the most accurate for df > 2 but is poor for df < 1.
Bailey's approximation is faster than Hill's and gives acceptable two-figure accuracy throughout.
Bailey's approximation also works for some extreme values, with very large \code{x} or \code{df}, for which Hill's approximation fails due to overflow.
}

\note{
The default approximation used when \code{approx=TRUE} was changed from Hill to Bailey in limma version 3.41.13.
}

\author{Gordon Smyth}

\references{
Bailey, B. J. R. (1980).
Accurate normalizing transformations of a Student's t variate.
\emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} 29(3), 304--306.

Hill, GW (1970). Algorithm 395: Student's t-distribution.
\emph{Communications of the ACM} 13, 617--620.

Brophy, AL (1987).
Efficient estimation of probabilities in the t distribution.
\emph{Behavior Research Methods} 19, 462--466.

Wallace, D. L. (1959).
Bounds on normal approximations to Student's and the chi-square distributions.
\emph{The Annals of Mathematical Statistics}, 30(4), 1121--1130.
}

\seealso{
\code{\link{zscore}}.

\code{zscoreNBinom} in the edgeR package.
}

\examples{
zscoreT(4, df=3)
zscoreT(4, df=3, approx=TRUE)
zscoreT(4, df=Inf)
tZscore(2.2, df=3)
}

\keyword{distributions}
