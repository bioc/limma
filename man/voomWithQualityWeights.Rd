\name{voomWithQualityWeights}
\alias{voomWithQualityWeights}
\title{Combining observational-level with sample-specific quality weights for RNA-seq analysis}
\description{
Combine voom observational-level weights with sample-specific quality weights in a designed experiment.
}
\usage{
voomWithQualityWeights(counts, design = NULL, lib.size = NULL, normalize.method = "none",
             plot = FALSE, span = 0.5, var.design = NULL, var.group = NULL,
             method = "genebygene", maxiter = 50, tol = 1e-5, trace = FALSE, col = NULL, ...) 
}
\arguments{
\item{counts}{
  a numeric \code{matrix} containing raw counts, or an \code{ExpressionSet} containing raw counts, or a \code{DGEList} object.
}
\item{design}{
  design matrix with rows corresponding to samples and columns to coefficients to be estimated.  
  Defaults to the unit vector meaning that samples are treated as replicates.
}
\item{lib.size}{
  numeric vector containing total library sizes for each sample.
  If \code{NULL} and \code{counts} is a \code{DGEList} then, the normalized library sizes are taken from \code{counts}.
  Otherwise library sizes are calculated from the columnwise counts totals.
}
\item{normalize.method}{
  normalization method to be applied to the logCPM values.
  Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.
}
\item{plot}{
  logical, should a plot of the mean-variance trend and sample-specific weights be displayed?
}
\item{span}{
  width of the lowess smoothing window as a proportion.
}
\item{var.design}{
  design matrix for the variance model.
  Defaults to the sample-specific model whereby each sample has a distinct quality weight.
}
\item{var.group}{
  vector or factor indicating groups to have different quality weights.
  This is another way to specify \code{var.design} for groupwise variance models.
}
\item{method}{
  character string specifying the method used to estimate the quality weights.
  Choices are \code{"genebygene"} or \code{"reml"}.
}
\item{maxiter}{
  maximum number of iterations allowed for quality weight estimation when \code{method = "reml"}.
}
\item{tol}{
  convergence tolerance for quality weight estimation when \code{method = "reml"}.
}
\item{trace}{
  logical.
  If \code{TRUE} then diagnostic information is output at each iteration of the \code{"reml"} algorithm, or at every 1000th iteration of the \code{"genebygene"} algorithm.
}
\item{col}{
  colours to use in the barplot of sample-specific weights if \code{plot=TRUE}).
  If \code{NULL}, then bars are plotted in grey.
}
\item{\dots}{
  other arguments are passed to \code{voom} and hence to \code{lmFit}.
}
}

\details{
This function is an alternative to \code{voom} and, like \code{voom}, is intended to process RNA-seq data prior to linear modeling in limma.
It combines observational-level weights from \code{voom} with sample-specific weights estimated using the \code{arrayWeights} function.
The method is described by Liu et al (2015).
}

\value{
An \code{\link[limma:EList]{EList}} object similar to that from \code{\link{voom}},
with an extra column \code{sample.weights} containing the vector of sample quality factors added to the \code{targets} data.frame.
The \code{weights} component combines the sample weights and the usual voom precision weights.
}

\author{Matthew Ritchie, Cynthia Liu, Gordon Smyth}

\references{
Liu, R., Holik, A. Z., Su, S., Jansz, N., Chen, K., Leong, H. S., Blewitt, M. E., Asselin-Labat, M.-L., Smyth, G. K., Ritchie, M. E. (2015).
Why weight? Combining voom with estimates of sample quality improves power in RNA-seq analyses.
\emph{Nucleic Acids Research} 43, e97.
\url{http://nar.oxfordjournals.org/content/43/15/e97}
}

\seealso{
\code{\link{voom}}, \code{\link{arrayWeights}}

See also \code{voomLmFit} in the edgeR package.

A summary of limma functions for RNA-seq analysis is given in \link{11.RNAseq}.
}

\keyword{regression}
\keyword{rna-seq}
