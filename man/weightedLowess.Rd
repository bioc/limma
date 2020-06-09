\name{weightedLowess}
\alias{weightedLowess}

\title{LOWESS Smoother with Prior Weights}

\description{
This function generalizes the original LOWESS smoother (locally-weighted regression) to incorporate prior weights while preserving the original algorithm design and efficiency as closely as possible.
}

\usage{
weightedLowess(x, y, weights = NULL,
               delta = NULL, npts = 200, span = 0.3, iterations = 4,
               output.style = "loess")
}

\arguments{
  \item{x}{a numeric vector of values for the covariate or x-axis coordinates.}
  \item{y}{a numeric vector of response values or y-axis coordinates, of same length as \code{x}.}
  \item{weights}{
    a numeric vector containing non-negative prior weights, of same length as \code{x}.
    Defaults to a constant vector.
  }
  \item{delta}{
    a numeric scalar specifying the maximum distance between successive anchor x-values where a local regression will be computed.
    Roughly corresponds to \code{diff(range(x))/npts} if the x-values are equally spaced.
    Setting \code{delta=0} forces every distinct x-value to be an anchor point.
    If \code{NULL} then a suitable \code{delta} value will be computed from \code{npts}.
  }
  \item{npts}{
    an integer scalar specifying the approximate number of anchor x-values at which local regressions will be computed.
    Ignored if \code{delta} is not \code{NULL}.
  }
  \item{span}{
    a numeric scalar between 0 and 1 specifying the width of the smoothing window as a proportion of the total weight.
  }
  \item{iterations}{
    an integer scalar specifying the number of iterations.
    \code{iterations=1} corresponds to local least squares regression without robustifying weights.
    Each additional iteration incorporates robustifying weights.
  }
  \item{output.style}{
    character string indicating whether the output should be in the style of \code{"loess"} or of \code{"lowess"}.
  }
}

\value{
If \code{output.style="loess"}, then a list with the following components:
  \item{fitted}{numeric vector of smoothed y-values (in the same order as the input vectors).}
  \item{residuals}{numeric vector or residuals.}
  \item{weights}{numeric vector of robustifying weights used in the most recent iteration.}
  \item{delta}{the delta used, either the input value or the value derived from \code{npts}.}

If \code{output.style="lowess"}, then a list with the following components:
  \item{x}{numeric vector of x-values in ascending order.}
  \item{y}{numeric vector or smoothed y-values.}
  \item{delta}{the delta used, either the input value or the value derived from \code{npts}.}
}

\details{
This function extends the LOWESS algorithm of Cleveland (1979, 1981) to handle non-negative prior weights.

The LOWESS method consists of computing a series of local linear regressions, with each local regression restricted to a window of x-values.
Smoothness is achieved by using overlapping windows and by gradually down-weighting points in each regression according to their distance from the anchor point of the window (tri-cube weighting).

To conserve running time and memory, locally-weighted regressions are computed at only a limited number of anchor x-values, either \code{npts} or the number of distinct x-values, whichever is smaller.
Anchor points are defined exactly as in the original LOWESS algorithm.
Any x-value within distance \code{delta} of an anchor point is considered adjacent to it.
The first anchor point is \code{min(x)}.
With the x-values sorted in ascending order, successive anchor points are defined as follows.
The next anchor point is the smallest x-value not adjacent to any previous anchor points.
The last anchor point is \code{max(x)}.

For each anchor point, a weighted linear regression is performed for a window of neighboring points.
The neighboring points consist of the smallest set of closest neighbors such as the sum of weights is greater than or equal to \code{span} times the total weight of all points.
Each local regression produces a fitted value for that anchor point.
Fitted values for other x-values are then obtained by linear interpolation between anchor points.

For the first iteration, the local linear regressions use weights equal to prior weights times the tri-cube distance weights.
Subsequent iterations multiple these weights by robustifying weights.
Points with residuals greater than 6 times the median absolute residual are assigned weights of zero and otherwise Tukey's biweight function is applied to the residuals to obtain the robust weights.
More iterations produce greater robustness.

In summary, the prior weights are used in two ways.
First, the prior weights are used during the span calculations such that the points included in the window for each local regression must account for the specified proportion of the total sum of weights.
Second, the weights used for the local regressions are the product of the prior weights, tri-cube local weights and biweight robustifying weights.
Hence a point with prior weight equal to an integer n has the same influence as n points with unit weight and the same x and y-values.

See also \code{\link{loessFit}}, which is is essentially a wrapper function for \code{lowess} and \code{\link{weightedLowess}} with added error checking.

\strong{Relationship to lowess and loess}

The stats package provides two functions \code{lowess} and \code{loess}.
\code{lowess} implements the original LOWESS algorithm of Cleveland (1979, 1981) designed for scatterplot smoothing with  single x-variable while \code{loess} implements the more complex algorithm by Cleveland et al (1988, 1992) designed to fit multivariate surfaces.
The \code{loess} algorithm is more general than \code{lowess} in a number of ways, notably because it allows prior weights and up to four numeric predictors.
On the other hand, \code{loess} is necessarily slower and uses more memory than \code{lowess}.
Furthermore, it has less accurate interpolation than \code{lowess} because it uses a cruder algorithm to choose the anchor points whereby anchor points are equi-spaced in terms of numbers of points rather than in terms of x-value spacing.
\code{lowess} and \code{loess} also have different defaults and input parameters.
See Smyth (2003) for a detailed discussion.

Another difference between \code{lowess} and \code{loess} is that \code{lowess} returns the x and y coordinates of the fitted curve, with x in ascending order, whereas \code{loess} returns fitted values and residuals in the original data order.

The purpose of the current function is to incorporate prior weights but keep the algorithmic advantages of the original \code{lowess} code for scatterplot smoothing.
The current function therefore generalizes the \code{span} and interpolation concepts of \code{lowess} differently to \code{loess}.

When \code{output.style="loess"}, \code{weightedLowess} outputs results in original order similar to \code{loessFit} and \code{loess}.
When \code{output.style="lowess"}, \code{weightedLowess} outputs results in sorted order the same as \code{lowess}.

The \code{span} argument corresponds to the \code{f} argument of \code{lowess} and the \code{span} argument of \code{loess}.
The \code{delta} argument is the same as the \code{delta} argument of \code{lowess}.
The \code{npts} argument is new and amounts to a more convenient way to specify \code{delta}.
The \code{iterations} argument is the same as the corresponding argument of \code{loess} and is equivalent to \code{iter+1} where \code{iter} is the \code{lowess} argument.
}

\author{C code and R function by Aaron Lun.}

\seealso{
\code{\link{lowess}},
\code{\link{loess}},
\code{\link{loessFit}},
\code{\link{tricubeMovingAverage}}.
}

\examples{
y <- rt(100,df=4)
x <- runif(100)
w <- runif(100)
l <- weightedLowess(x, y, w, span=0.7, output.style="lowess")
plot(x, y, cex=w)
lines(l, col = "red")
}

\references{
Cleveland, W.S. (1979).
Robust Locally Weighted Regression and Smoothing Scatterplots.
\emph{Journal of the American Statistical Association} 74(368), 829-836.

Cleveland, W.S. (1981).
LOWESS: A program for smoothing scatterplots by robust locally weighted regression.
\emph{The American Statistician} 35(1), 54.

Cleveland, W.S., and Devlin, S.J. (1988).
Locally-weighted regression: an approach to regression analysis by local fitting.
\emph{Journal of the American Statistical Association} 83(403), 596-610.

Cleveland, W.S., Grosse, E., and Shyu, W.M. (1992).
Local regression models.
Chapter 8 In:
\emph{Statistical Models in S} edited by J.M. Chambers and T.J. Hastie, Chapman & Hall/CRC, Boca Raton.

Smyth, G.K. 2003.
\emph{lowess vs. loess}.
Answer on the Bioconductor Support forum \url{https://support.bioconductor.org/p/2323/}.
}

\concept{locally-weighted regression smoothing}
