\name{voomWithQualityWeights}
\alias{voomWithQualityWeights}
\title{Combining observational-level with sample-specific quality weights for RNA-seq analysis}
\description{
Combine voom observational-level weights with sample-specific quality weights in a designed experiment.
}
\usage{
voomWithQualityWeights(counts, design=NULL, lib.size=NULL, normalize.method="none",
             plot=FALSE, span=0.5, var.design=NULL, method="genebygene", maxiter=50,
             tol=1e-10, trace=FALSE, replace.weights=TRUE, col=NULL, ...) 
}
\arguments{
 \item{counts}{a numeric \code{matrix} containing raw counts, or an \code{ExpressionSet} containing raw counts, or a \code{DGEList} object.}
 \item{design}{design matrix with rows corresponding to samples and columns to coefficients to be estimated.  
               Defaults to the unit vector meaning that samples are treated as replicates.}
 \item{lib.size}{numeric vector containing total library sizes for each sample.
 If \code{NULL} and \code{counts} is a \code{DGEList} then, the normalized library sizes are taken from \code{counts}.
 Otherwise library sizes are calculated from the columnwise counts totals.}
 \item{normalize.method}{normalization method to be applied to the logCPM values.
 Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.}
 \item{plot}{\code{logical}, should a plot of the mean-variance trend and sample-specific weights be displayed?}
 \item{span}{width of the lowess smoothing window as a proportion.}
 \item{var.design}{design matrix for the variance model. Defaults to the sample-specific model (i.e. each sample has a distinct variance) when \code{NULL}.}
 \item{method}{character string specifying the estimating algorithm to be used. Choices
          are \code{"genebygene"} and \code{"reml"}.}
 \item{maxiter}{maximum number of iterations allowed.}
 \item{tol}{convergence tolerance.}
 \item{trace}{logical variable. If true then output diagnostic information
          at each iteration of the '"reml"' algorithm, or at every 1000th iteration of the 
          '"genebygene"' algorithm.}
 \item{replace.weights}{logical variable. If TRUE then the weights in the voom object will be replaced with 
       the combined voom and sample-specific weights and the \code{\link[limma:EList]{EList}} object from voom is returned.
       If FALSE, then a matrix of combined weights is returned.}
 \item{col}{colours to use in the barplot of sample-specific weights (only used if \code{plot=TRUE}). If \code{NULL}, bars are plotted in grey.}
 \item{\dots}{other arguments are passed to \code{lmFit}.}
}
\details{
This function is intended to process RNA-Seq data prior to linear modelling in limma.

It combines observational-level weights from \code{voom} with sample-specific weights estimated using the \code{arrayWeights} function.
}
\value{
	Either an \code{\link[limma:EList]{EList}} object with the following components:
\item{E}{numeric matrix of normalized expression values on the log2 scale}
\item{weights}{numeric matrix of inverse variance weights}
\item{design}{design matrix}
\item{lib.size}{numeric vector of total normalized library sizes}
\item{genes}{dataframe of gene annotation extracted from \code{counts}}
       or a matrix of combined voom and sample-specific weights with same dimension as \code{counts}.
}

\author{Matthew Ritchie and Cynthia Liu}

\references{
Liu, R., Holik, A. Z., Su, S., Jansz, N., Chen, K., Leong, H. S., Blewitt, M. E., Asselin-Labat, M.-L., Smyth, G. K., Ritchie, M. E. (2014).
Why weight? Combining voom with estimates of sample quality improves power in RNA-seq analyses.
(in preparation)

Ritchie, M. E., Diyagama, D., Neilson, van Laar, R., J., Dobrovic, A., Holloway, A., and Smyth, G. K. (2006).
Empirical array quality weights in the analysis of microarray data.
\emph{BMC Bioinformatics} 7, 261.
\url{http://www.biomedcentral.com/1471-2105/7/261}

Law, C. W., Chen, Y., Shi, W., Smyth, G. K. (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
\url{http://genomebiology.com/2014/15/2/R29}
}

\seealso{
\code{\link{voom}}, \code{\link{arrayWeights}}

An overview of linear model functions in limma is given by \link{06.LinearModels}.
A \code{voomWithQualityWeights} case study is given in the User's Guide.
}

\keyword{regression}
\keyword{rna-seq}