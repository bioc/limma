\name{dupcor}
\alias{duplicateCorrelation}
\title{Correlation Between Duplicates or Within Blocks}
\description{
Estimate the intra-block correlation given a block structure for the arrays or samples.
}
\usage{
duplicateCorrelation(object, design=NULL, ndups=2, spacing=1, block=NULL,
                     trim=0.15, weights=NULL)
}
\arguments{
  \item{object}{A matrix-like data object containing log-ratios or log-expression values for a series of samples, with rows corresponding to genes and columns to samples.
  Any type of data object that can be processed by \code{\link{getEAWP}} is acceptable.}
  \item{design}{the design matrix of the microarray experiment, with rows corresponding to arrays and columns to comparisons to be estimated. The number of rows must match the number of columns of \code{object}. Defaults to the unit vector meaning that the arrays are treated as replicates.} 
  \item{ndups}{a positive integer giving the number of times each gene is printed on an array. \code{nrow(object)} must be divisible by \code{ndups}.
  Ignored if \code{block} is specified.}
  \item{spacing}{the spacing between the rows of \code{object} corresponding to duplicate spots, \code{spacing=1} for consecutive spots}
  \item{block}{vector or factor specifying a blocking variable}
  \item{trim}{the fraction of observations to be trimmed from each end of \code{tanh(all.correlations)} when computing the trimmed mean.}
  \item{weights}{an optional numeric matrix of the same dimension as \code{object} containing weights for each spot. If smaller than \code{object} then it will be filled out to the same size.}
}
\value{
  A list with components
  \item{consensus.correlation}{the average estimated inter-duplicate correlation. The average is the trimmed mean of the individual correlations on the atanh-transformed scale.}
  \item{cor}{same as \code{consensus.correlation}, for compatibility with earlier versions of the software}
  \item{atanh.correlations}{numeric vector of length \code{nrow(object)/ndups} giving the individual genewise atanh-transformed correlations.}
}
\details{
When \code{block=NULL}, this function estimates the correlation between duplicate spots (regularly spaced within-array replicate spots).
If \code{block} is not null, this function estimates the correlation between repeated observations on the blocking variable.
Typically the blocks are biological replicates and repeated observations on the same block may be correlated.
In either case, the correlation is estimated by fitting a mixed linear model by REML individually for each gene.
The function also returns a consensus correlation, which is a robust average of the individual correlations, intended for input to functions such as \code{lmFit}, \code{gls.series} or \code{voom}.

It is not possible to estimate correlations between duplicate spots and with sample blocks simultaneously.
If \code{block} is not null, then the function will set \code{ndups=1}, which is equivalent to ignoring duplicate spots.

For this function to return statistically useful results, there must be at least two more arrays than the number of coefficients to be estimated, i.e., two more than the column rank of \code{design}.

The function may take long time to execute as it fits a mixed linear model for each gene using an iterative algorithm.

If present, \code{ndups} and \code{spacing} will be extracted from \code{object$printer$ndups} and \code{object$printer$spacing}.
}
\seealso{
These functions use \code{\link[statmod:mixedmodel]{mixedModel2Fit}} from the statmod package.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth}
\references{
Smyth, G. K., Michaud, J., and Scott, H. (2005). The use of within-array replicate spots for assessing differential expression in microarray experiments. \emph{Bioinformatics} 21(9), 2067-2075.
[\url{http://bioinformatics.oxfordjournals.org/content/21/9/2067}]
[Preprint with corrections: \url{http://www.statsci.org/smyth/pubs/dupcor.pdf}]
}

\examples{
# Simulate a paired experiment with incomplete blocks
Block <- c(1,1,2,2,3,3,4,4,5,6,7,8)
Treat <- factor(c(1,2,1,2,1,2,1,2,1,2,1,2))
design <- model.matrix(~Treat)

ngenes <- 50
nsamples <- 12
y <- matrix(rnorm(ngenes*nsamples),ngenes,nsamples)
rownames(y) <- paste0("Gene",1:ngenes)

# Estimate the within-block correlation
dupcor <- duplicateCorrelation(y,design,block=Block)
dupcor$consensus.correlation

# Estimate the treatment effect using both complete and incomplete blocks
fit <- lmFit(y,design,block=Block,correlation=dupcor$consensus)
fit <- eBayes(fit)
topTable(fit,coef=2)
}

\keyword{models}
\keyword{regression}
