\name{arrayWeights}
\alias{arrayWeights}
\title{Array Quality Weights}

\description{
Estimate relative quality weights for each array or group in a multi-array
experiment.
}

\usage{
arrayWeights(object, design = NULL, weights = NULL,
     var.design = NULL, var.group = NULL, prior.n = 10,
     method = "auto", maxiter = 50L, tol = 1e-5, trace = FALSE,
     nthreads = 1L)
}

\arguments{
  \item{object}{any matrix-like object containing log-expression values or log-ratio expression values, for example an \code{EList} or \code{ExpressionSet} object.
    See \code{help("getEAWP")} for a list of possible classes.}
  \item{design}{the design matrix of the microarray experiment, with rows
           corresponding to arrays and columns to coefficients to be
           estimated.  Defaults to the unit vector meaning that the
           arrays are treated as replicates.}
  \item{weights}{numeric matrix containing prior weights for each expression value.}
  \item{var.design}{design matrix for the variance model. Defaults to the sample-specific model 
           whereby each sample has a distinct variance.}
  \item{var.group}{vector or factor indicating groups to have different array weights. This is another way to specify \code{var.design} for groupwise variance models.}
  \item{prior.n}{prior number of genes. Larger values squeeze the array weights more strongly towards equality.}
  \item{method}{character string specifying the estimating algorithm to be used. Choices
          are \code{"genebygene"}, \code{"reml"} or \code{"auto"}.}
  \item{maxiter}{maximum number of iterations allowed when \code{method="reml"}.}
  \item{tol}{convergence tolerance when \code{method="reml"}.}
  \item{trace}{logical. If \code{TRUE} then progress information is printed
          at each iteration of the \code{"reml"} algorithm or at every 1000th gene for the 
          \code{"genebygene"} algorithm.}
  \item{nthreads}{positive integer giving the number of OpenMP threads used by the REML algorithm. Defaults to 1.}
}

\details{
The relative reliability of each array is estimated by measuring how
well the expression values for that array follow the linear model.
Arrays that tend to have larger residuals are assigned lower weights.

The basic method is described by Ritchie et al (2006) and the extension to custom variance models by Liu et al (2015).
A weighted linear model is fitted to the expression values for 
each gene.
The variance model is fitted to the squared residuals from the linear model fit and is updated either by full REML 
scoring iterations (\code{method="reml"}) or using an efficient gene-by-gene update algorithm (\code{method="genebygene"}).
The final estimates of these array variances are converted to weights.
The gene-by-gene algorithm is described by Ritchie et al (2006) while the REML algorithm is an adaptation of the algorithm of Smyth (2002).

For stability, the array weights are squeezed slightly towards equality.
This is done by adding a prior likelihood corresponding to unit array weights equivalent to \code{prior.n} genes.
The gene-by-gene algorithm is started from the prior genes while the REML algorithm adds the prior to the log-likelihood derivatives.

By default, \code{arrayWeights} chooses between the REML and gene-by-gene algorithms automatically (\code{method="auto"}).
REML is chosen if there are no prior weights or missing values and otherwise the gene-by-gene algorithm is used.

The input \code{object} is interpreted as for \code{lmFit} and \code{getEAWP}.
In particular, the arguments \code{design} and \code{weights} will be extracted from the data 
\code{object} if available and do not normally need to be set explicitly in 
the call; if any of these are set in the call then they will over-ride 
the slots or components in the data \code{object}.
}

\value{
	A numeric vector of array weights, which multiply to 1.
 }

\references{
Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat ML, Smyth GK, Ritchie ME (2015).
Why weight? Combining voom with estimates of sample quality improves power in RNA-seq analyses.
\emph{Nucleic Acids Research} 43, e97.
\doi{https://doi.org/10.1093/nar/gkv412}.

Ritchie ME, Diyagama D, Neilson J, van Laar R, Dobrovic A, Holloway A, Smyth GK (2006).
Empirical array quality weights in the analysis of microarray data.
\emph{BMC Bioinformatics} 7, 261.
\doi{10.1186/1471-2105-7-261}.

Smyth GK (2002).
An efficient algorithm for REML in heteroscedastic regression.
\emph{Journal of Computational and Graphical Statistics} 11, 836-847.
\url{https://gksmyth.github.io/pubs/remlalgo.pdf}
}

\seealso{
\code{\link{arrayWeightsQuick}}, \code{\link{voomWithQualityWeights}}

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\examples{
ngenes <- 1000
narrays <- 6
y <- matrix(rnorm(ngenes*narrays), ngenes, narrays)
var.group <- c(1,1,1,2,2,2)
y[,var.group==1] <- 2*y[,var.group==1]
arrayWeights(y, var.group=var.group)
}

\author{Matthew Ritchie and Gordon Smyth}
\keyword{regression}
\keyword{models}
