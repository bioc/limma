\name{removeBatchEffect}
\alias{removeBatchEffect}
\title{Remove Batch Effect}
\description{
Remove batch effects from expression data.
}
\usage{
removeBatchEffect(x, batch=NULL, batch2=NULL, covariates=NULL,
                  design=matrix(1,ncol(x),1), \dots)
}
\arguments{
  \item{x}{numeric matrix, or any data object that can be processed by \code{\link{getEAWP}} containing log-expression values for a series of samples.
  Rows correspond to probes and columns to samples.}
  \item{batch}{factor or vector indicating batches.}
  \item{batch2}{factor or vector indicating a second series of batches.}
  \item{covariates}{matrix or vector of numeric covariates to be adjusted for.}
  \item{design}{design matrix relating to treatment conditions to be preserved, usually the design matrix with all experimental factors other than the batch effects.}
  \item{\dots}{other arguments are passed to \code{\link{lmFit}}.}
}
\value{
A numeric matrix of log-expression values with batch and covariate effects removed.
}
\details{
This function is useful for removing unwanted batch effects, associated with hybridization time or other technical variables, ready for plotting or unsupervised analyses such as PCA, MDS or heatmaps.
The design matrix is used to describe comparisons between the samples, for example treatment effects, that should not be removed.
The function (in effect) fits a linear model to the data, including both batches and regular treatments, then removes the component due to the batch effects.

In most applications, only the first \code{batch} argument will be needed.
This case covers the situation where the data has been collected in a series of separate batches.

The \code{batch2} argument is used when there is a second series of batch effects, independent of the first series.
For example, \code{batch} might correspond to time of data collection while \code{batch2} might correspond to operator or some other change in operating characteristics.
If \code{batch2} is included, then the effects of \code{batch} and \code{batch2} are assumed to be additive.

The \code{covariates} argument allows correction for one or more continuous numeric effects, similar to the analysis of covariance method in statistics.
If \code{covariates} contains more than one column, then the columns are assumed to have additive effects.
Setting \code{covariates} to be a design matrix constructed from batch effects and technical effects allows very general batch effects to be accounted for.

The data object \code{x} can be of any class for which \code{lmFit} works.
If \code{x} contains weights, then these will be used in estimating the batch effects.
}

\note{
This function is not intended to be used prior to linear modelling.
For linear modelling, it is better to include the batch factors in the linear model.
}

\seealso{
\link{05.Normalization}
}
\author{Gordon Smyth and Carolyn de Graaf}

\examples{
y <- matrix(rnorm(10*9),10,9)
y[,1:3] <- y[,1:3] + 5
batch <- c("A","A","A","B","B","B","C","C","C")
y2 <- removeBatchEffect(y, batch)
par(mfrow=c(1,2))
boxplot(as.data.frame(y),main="Original")
boxplot(as.data.frame(y2),main="Batch corrected")
}
