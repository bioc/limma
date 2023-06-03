\name{removeBatchEffect}
\alias{removeBatchEffect}
\title{Remove Batch Effect}
\description{
Remove batch effects from expression data.
}
\usage{
removeBatchEffect(x, batch = NULL, batch2 = NULL, covariates = NULL,
                  design = matrix(1,ncol(x),1), group = NULL, \dots)
}
\arguments{
  \item{x}{numeric matrix, or any data object that can be processed by \code{\link{getEAWP}} containing log-expression values for a series of samples.
  Rows correspond to probes and columns to samples.}
  \item{batch}{factor or vector indicating batches.}
  \item{batch2}{factor or vector indicating a second series of batches.}
  \item{covariates}{matrix or vector of numeric covariates to be adjusted for.}
  \item{design}{design matrix relating to experimental conditions to be preserved, usually the design matrix with all experimental factors other than the batch effects. Ignored if \code{group} is not \code{NULL}.}
  \item{group}{factor defining the experimental conditions to be preserved. An alternative way to specify the design matrix \code{design}.}
  \item{\dots}{other arguments are passed to \code{\link{lmFit}}.}
}
\value{
A numeric matrix of log-expression values with batch and covariate effects removed.
}
\details{
This function is useful for removing unwanted batch effects, associated with hybridization time or other technical variables, ready for plotting or unsupervised analyses such as PCA, MDS or heatmaps.
The design matrix or group factor is used to define comparisons between the samples, for example treatment effects, that should not be removed.
The function fits a linear model to the data, including both batches and regular treatments, then removes the component due to the batch effects.

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
This function is intended for plotting and data exploration purposes.
This function is not intended to be used to prepare data for linear modeling by \code{lmFit}.
For linear modeling, it is better to include the batch factors in the linear model so fhat \code{lmFit} can correctly assess the standard errors of the linear model parameters.
}

\seealso{
\link{05.Normalization}
}
\author{Gordon Smyth and Carolyn de Graaf}

\examples{
ngenes <- 10
nsamples <- 8
y <- matrix(rnorm(ngenes*nsamples),ngenes,nsamples)
group <- factor(c("A","A","A","A","B","B","B","B"))
batch <- factor(c(1,1,2,2,1,1,2,2))
colnames(y) <- paste(group,batch,sep=".")
y[,batch==2] <- y[,batch==2] + 5
y[,group=="B"] <- y[,group=="B"] + 1
y.corrected <- removeBatchEffect(y, batch=batch, group=group)
oldpar <- par(mfrow=c(1,2))
plotMDS(y,main="Original")
plotMDS(y.corrected,main="Batch corrected")
par(oldpar)
devAskNewPage(FALSE)
}
