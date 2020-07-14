\name{chooseLowessSpan}
\alias{chooseLowessSpan}
\title{Choose Span for Local-Weighted Regression Smoothing}

\description{
Choose an optimal span, depending on the number of points, for lowess smoothing of variance trends.
}

\usage{
chooseLowessSpan(n=1000, small.n=25, min.span=0.2, power=1/3)
}

\arguments{
  \item{n}{the number of points the lowess curve will be applied to.}
  \item{small.n}{the span will be set to 1 for any \code{n} less than or equal to size value.}
  \item{min.span}{the minimum span for large \code{n}.}
  \item{power}{numeric power between 0 and 1 how fast the chosen span decreases with \code{n>}.}
}

\details{
The span is the proportion of points used for each of the local regressions.
When there a few points, a large span should be used to ensure a smooth curve.
When there are a large number of points, smaller spans can be used because each span window still contains good coverage.
By default, the chosen span decreases as the cube-root of the number of points, a rule that is motivated by analogous rules to choose the number of bins for a histogram (Scott, 1979; Freedman & Diaconis, 1981; Hyndman, 1995).

The span returned is essentially
\code{min.span + (1-min.span) * (small.n/n)^power}.
The span is set to 1 for any \code{n} less than \code{small.n}.

The function is tuned for smoothing of mean-variance trends, for which the trend is usually monotonic, so preference is given to moderately large spans.
Even for the large datasets, the span is always greater than \code{min.span}.

This function is used to create the default span for \code{vooma}, \code{eBayes}, \code{squeezeVar} and \code{fitFDistRobustly}.
}

\value{
A numeric vector of length 1 containing the span value.
}

\author{Gordon Smyth}

\references{
Freedman, D. and Diaconis, P. (1981). On the histogram as a density estimator: L_2 theory.
\emph{Zeitschrift fur Wahrscheinlichkeitstheorie und verwandte Gebiete}, 57, 453-476.

Hyndman, R. J. (1995). The problem with Sturges' rule for constructing histograms.
\url{http://robjhyndman.com/papers/sturges.pdf}.

Scott, D. W. (1979). On optimal and data-based histograms.
\emph{Biometrika}, 66, 605-610.
}

\seealso{
\code{\link{loessFit}}, \code{\link{weightedLowess}}, \code{\link{lowess}}, \code{\link{loess}}.

\code{\link{vooma}}, \code{\link{eBayes}}, \code{\link{squeezeVar}}, \code{\link{fitFDistRobustly}}.
}

\examples{
chooseLowessSpan(100)
chooseLowessSpan(1e6)

n <- 10:5000
span <- chooseLowessSpan(n)
plot(n,span,type="l")
}

\concept{locally-weighted regression smoothing}
