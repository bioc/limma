\name{vooma}
\alias{vooma}
\alias{voomaByGroup}

\title{Convert Variance Trend to Observation Precision Weights for Microarray-style Data}

\description{
Estimate the variance trend for microarray data and use it to compute appropriate observational-level weights.
The variance trend optionally depends on a second predictor as well as on average log-expression.
}

\usage{
vooma(y, design = NULL, block = NULL, correlation,
      predictor = NULL, span = NULL, legacy.span = FALSE, plot = FALSE, save.plot = FALSE)

voomaByGroup(y, group, design = NULL, block = NULL, correlation,
             span = NULL, legacy.span = FALSE, plot = FALSE, col = NULL, lwd = 1,
             pch = 16, cex = 0.3, alpha = 0.5, legend = "topright")
}

\arguments{
  \item{y}{a numeric \code{matrix}, \code{EList} object, or any object containing log-expression data that can be coerced to a matrix. Rows correpond to genes and columns to samples.}
  \item{design}{design matrix with rows corresponding to samples and columns to coefficients to be estimated.  Defaults to the unit vector meaning that samples are treated as replicates.}
  \item{block}{vector or factor specifying a blocking variable on the arrays. Has length equal to the number of arrays.}
  \item{correlation}{intra-block correlation}
  \item{predictor}{
    precision predictor.
    Numeric matrix of the same dimensions as \code{y} that predicts the precision of each log-expression value.
    Is used as a second covariate together with the log-intensities to predict the variances and produce the final precision weights.
  }
  \item{span}{
    width of the smoothing window, as a proportion of the data set.
    Defaults to a value that depends the number of genes (\code{nrow(y)}).
    Equal to 1 if the number of genes is less than or equal to 50, then decreases slowly to 0.3 if the number of genes is very large.
  }
  \item{legacy.span}{logical. If \code{TRUE}, then the original default setting will be used for \code{span}, which is slightly smaller than the new default.}
  \item{plot}{logical. If \code{TRUE}, a plot of the mean-variance trend is displayed.}
  \item{save.plot}{
    logical, should the coordinates and line of the plot be saved in the output?
  }
  \item{group}{categorical vector or factor giving group membership of columns of \code{y}.}
  \item{col}{vector of colors for plotting group trends}
  \item{lwd}{line width for plotting group trends}
  \item{pch}{plotting character. Default is integer code 16, which gives a solid circle. If a vector, then should be of length \code{nrow(y)}.}
  \item{cex}{numeric vector of plot symbol expansions.  If a vector, then should be of length equal to number of groups.}
  \item{alpha}{transparency of points, on scale from \code{0} for fully transparent to \code{1} for fully opaque.}
  \item{legend}{character string giving position to place legend.}
}

\details{
\code{vooma} is an acronym for "mean-variance modelling at the observational level for arrays".
It is analogous to \code{voom} but for continuous log-expression values rather than for sequencing counts.

\code{vooma} estimates the mean-variance relationship in the data and uses it to compute appropriate precision weights for each observation.
The mean-variance trend is estimated from gene-level data but is extrapolated back to individual observations to obtain a precision weight (inverse variance) for each observation.
The weights can then used by other functions such as \code{lmFit} to adjust for heteroscedasticity.

If \code{span=NULL}, then an optimal span value is estimated depending on \code{nrow(y)}.
The span is chosen by \code{chooseLowessSpan} with \code{n=nrow(y)}, \code{small.n=50}, \code{min.span=0.3} and \code{power=1.3}.
If \code{legacy.span=TRUE}, then the \code{chooseLowessSpan} arguments are reset to \code{small.n=10}, \code{min.span=0.3} and \code{power=0.5} to match the settings used by \code{vooma} in limma version 3.59.1 and earlier.

The variance trend can be modeled using a second optional \code{predictor} as well as in terms of log-expression.
If \code{predictor} is not \code{NULL}, then the variance trend is modeled as a function of both the mean log-expression and the \code{predictor} using a multiple linear regression with the two predictors.
In this case, the \code{predictor} is assumed to be some prior predictor of the precision or standard deviation of each log-expression value.
Any \code{predictor} that is correlated with the precision of each observation should give good results.

\code{voomaByGroup} estimates precision weights separately for different groups of samples.
In other words, it allows for different mean-variance curves in different groups.
\code{voomaByGroup} has a quite simple implementation and will not work work with complex design matrices.
It essentially assumes that the \code{design} matrix defines a oneway layout and it works best when the oneway layout refers to the same groups as defined by the `group` argument.
}

\value{
An EList object with the following components:
\item{E}{numeric matrix of log-expression values. Equal to \code{y} for \code{vooma} or log2-counts-per-million for \code{voomaByGroup}.}
\item{weights}{numeric matrix of observation precision weights.}
\item{design}{numeric matrix of experimental design.}
\item{genes}{data-frame of gene annotation, only if \code{counts} was a \code{DGEList} object.}
\item{voom.xy}{if \code{save.plot}, list containing x and y coordinates for points in mean-variance plot}
\item{voom.line}{if \code{save.plot}, list containing coordinates of lowess line in the mean-variance plot}
}

\author{Charity Law, Gordon Smyth and Mengbo Li. Mengbo Li contributed the functionality associated with the \code{predictor} argument.}

\references{
Law CW (2013).
\emph{Precision weights for gene expression analysis}.
PhD Thesis. University of Melbourne, Australia.
\url{http://hdl.handle.net/11343/38150}
}

\seealso{
\code{\link{voom}}, \code{\link{arrayWeights}}
}

\examples{
group <- gl(2,4)
design <- model.matrix(~group)
y <- matrix(rnorm(500*8),500,8)
u <- matrix(runif(length(y)),500,8)
yu <- y*u
v <- vooma(yu,design,plot=TRUE,predictor=u)
}
