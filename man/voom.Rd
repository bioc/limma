\name{voom}
\alias{voom}
\title{Transform RNA-Seq Counts for Linear Modeling With Precision Weights}
\description{
Transform count data to log2 counts-per-million (logCPM), estimate the mean-variance relationship and use it to compute observation-level precision weights.
The logCPM and associated precision weights are then ready for linear modeling.
}

\usage{
voom(counts, design = NULL, lib.size = NULL, normalize.method = "none",
     block = NULL, correlation = NULL, weights = NULL,
     span = 0.5, adaptive.span = FALSE, plot = FALSE, save.plot = FALSE)
}

\arguments{
  \item{counts}{
    a numeric \code{matrix} containing raw counts, or an \code{ExpressionSet} containing raw counts, or a \code{DGEList} object.
    Counts must be non-negative and NAs are not permitted.
  }
  \item{design}{
    design matrix with rows corresponding to samples and columns to coefficients to be estimated.
    Defaults to \code{model.matrix(~group, data=counts$samples)} if \code{counts} is a DGEList, otherwise defaults to the unit vector meaning that all samples are treated as replicates.
  }
  \item{lib.size}{
    numeric vector containing the library sizes for each sample.
    Defaults to the columnwise count totals if \code{counts} is a matrix or
    to \code{normLibSizes(counts)} if \code{counts} is a \code{DGEList}.
  }
  \item{normalize.method}{
    the microarray-style normalization method to be applied to the logCPM values.
    Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.
  }
  \item{block}{
    vector or factor specifying a blocking variable on the samples.
    Has length equal to the number of samples (\code{ncol(counts)}).
  }
  \item{correlation}{
    the intrablock correlation.
    Normally a single numeric value between -1 and 1, but a vector of genewise correlations is also allowed.
  }
  \item{weights}{
    prior weights.
    Can be a numeric matrix of individual weights of same dimensions as the \code{counts},
    or a numeric vector of sample weights with length equal to \code{ncol(counts)},
    or a numeric vector of gene weights with length equal to \code{nrow(counts)}.
  }
  \item{span}{
    width of the smoothing window used for the lowess mean-variance trend.
    Expressed as a proportion between 0 and 1.
  }
  \item{adaptive.span}{
    logical.
    If \code{TRUE}, then an optimal value for \code{span} will be chosen depending on the number of genes.
  }
  \item{plot}{
    logical, should a plot of the mean-variance trend be displayed?
  }
  \item{save.plot}{
    logical, should the coordinates and line of the plot be saved in the output?
  }
}

\details{
This function processes sequence count data from technologies such as RNA-seq or ChIP-seq to make it ready for linear modeling in limma.

\code{voom} is an acronym for "mean-variance modeling at the observational level".
The idea is to estimate the mean-variance relationship in the data, then use this to compute an appropriate precision weight for each observation.
Count data always show marked mean-variance relationships.
Raw counts show increasing variance with increasing count size, while log-counts typically show a decreasing mean-variance trend.
This function estimates the mean-variance trend for log-counts, then assigns a weight to each observation based on its predicted variance.
The weights are then used in the linear modeling process to adjust for heteroscedasticity. 
The mean-variance trend is estimated from gene-level data but is extrapolated back to individual observations to obtain a precision weight (inverse variance) for each observation.

\code{voom} performs the following specific calculations.
First, the counts are converted to logCPM values, adding 0.5 to all the counts to avoid taking the logarithm of zero.
The logCPM calculation uses normalized library sizes if \code{counts} is a DGEList or simply the column sums if \code{counts} is a matrix.
A microarray-style normalization method can also be optionally applied to the matrix of logCPM values.
The \code{lmFit} function is used to fit row-wise linear models.
The \code{lowess} function is then used to fit a trend to the square-root residual standard deviations as a function of an average log-count measure.
The trend line is then used to predict the variance of each logCPM value as a function of its fitted value on the count scale, and the inverse variances become the estimated precision weights.

The optional arguments \code{block}, \code{correlation} and \code{weights} are passed to \code{\link{lmFit}} in the above calling sequence, so they influence the row-wise standard deviations to which the mean-variance trend is fitted.
The arguments \code{block} and \code{correlation} have the same meaning as for \code{\link{lmFit}}.
Most users will not need to specify the \code{weights} argument but, if it is included, then the output \code{weights} are taken to modify the input prior weights in a multiplicative fashion.

For good results, the \code{counts} matrix should be filtered to remove rows with very low counts before running voom().
The \code{filterByExpr} function in the edgeR package can be used for that purpose.

If \code{counts} is a \code{DGEList} object from the edgeR package, then voom will use the normalization factors found in the object when computing the logCPM values.
In other words, the logCPM values are computed from the effective library sizes rather than the raw library sizes.
If the \code{DGEList} object has been scale-normalized in edgeR, then it is usual to leave \code{normalize.method="none"} in voom, i.e., the logCPM values should not usually be re-normalized in the \code{voom} call.

The \code{voom} method is similar in purpose to the limma-trend method, which uses \code{\link{eBayes}} or \code{\link{treat}} with \code{trend=TRUE}.
The voom method incorporates the mean-variance trend into the precision weights, whereas limma-trend incorporates the trend into the empirical Bayes moderation.
The voom method takes into account the sequencing depths (library sizes) of the individual columns of \code{counts} and applies the mean-variance trend on an individual observation basis.
limma-trend, on the other hand, assumes that the library sizes are not wildly different and applies the mean-variance trend on a genewise basis.
As noted by Law et al (2014), voom should be more powerful than limma-trend if the library sizes are very different but, otherwise, the two methods should give similar results.

Note that \code{edgeR::voomLmFit} is a further-developed version \code{voom} with more functionality and convenience.
\code{voomLmFit} is now recommended over \code{voom}, particularly if an intrablock correlation needs to be estimated or if the counts are sparse with a high proportion of zeros.
}

\note{
\code{voom} is designed to accept counts.
Usually these will be sequence read counts, but counts of species abundance or other biological quantities might also be appropriate.
Estimated counts are also acceptable provided that the column sums are representative of the total library size (total number of reads) for that sample.
\code{voom} can analyze scaled counts provided that the column sums remain proportional to the total library sizes.
\code{voom} is designed to take account of sample-specific library sizes and hence \code{voom} should not be used to analyze quantities that have been normalized for library size such as RPKM, transcripts per million (TPM) or counts per million (CPM).
Such quantities prevent \code{voom} from inferring the correct library sizes and hence the correct precision with which each value was measured.
}

\value{
An \code{\link[limma:EList]{EList}} object with the following components:
\item{E}{numeric matrix of normalized expression values on the log2 scale}
\item{weights}{numeric matrix of inverse variance weights}
\item{design}{design matrix}
\item{lib.size}{numeric vector of total normalized library sizes}
\item{genes}{data-frame of gene annotation extracted from \code{counts}}
\item{voom.xy}{if \code{save.plot}, list containing x and y coordinates for points in mean-variance plot}
\item{voom.line}{if \code{save.plot}, list containing coordinates of loess line in the mean-variance plot}
 }

\author{Charity Law and Gordon Smyth}

\references{
Law CW (2013).
\emph{Precision weights for gene expression analysis}.
PhD Thesis. University of Melbourne, Australia.
\url{http://hdl.handle.net/11343/38150}

Law CW, Chen Y, Shi W, Smyth GK (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
\doi{10.1186/gb-2014-15-2-r29}.
See also the Preprint Version at \url{https://gksmyth.github.io/pubs/VoomPreprint.pdf} incorporating some notational corrections.

Law CW, Alhamdoosh M, Su S, Smyth GK, Ritchie ME (2016).
RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR.
\emph{F1000Research} 5, 1408.
\url{https://f1000research.com/articles/5-1408}

Law CW, Alhamdoosh M, Su S, Dong X, Tian L, Smyth GK, Ritchie ME (2018).
RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR.
\emph{Bioconductor Workflow Package}.
\url{https://www.bioconductor.org/packages/RNAseq123/}
}

\seealso{
\code{\link{lmFit}} and \code{\link{eBayes}} are downstream of \code{voom}.

\code{\link{voomWithQualityWeights}} and \code{edgeR::voomLmFit} are further developed versions of \code{voom} with extra capabilities.
Either can be used as a replacement for \code{voom}.
\code{voomWithQualityWeights} estimates samples weights.
\code{voomLmFit} estimates both sample weights and intrablock correlation and also improves variance estimation for sparse data.

\code{\link{vooma}} is analogous to \code{voom} but for continuous log-expression data instead of RNA-seq counts.

A summary of limma functions for RNA-seq analysis is given in \link{11.RNAseq}.
}

\examples{
\dontrun{
keep <- filterByExpr(counts, design)
v <- voom(counts[keep,], design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust=TRUE)}
}
\keyword{rna-seq}
