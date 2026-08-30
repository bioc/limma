\name{voomLmFit}
\alias{voomLmFit}
\title{Voom Variance Modelling with Sample Weights, Block Correlation and Sparse Data}
\description{
Transform count data to log2-counts per million (logCPM), estimate voom precision weights and fit limma linear models while allowing for loss of residual degrees of freedom due to exact zeros.
}

\usage{
voomLmFit(counts, design = NULL, contrasts = NULL, block = NULL, prior.weights = NULL,
     sample.weights = FALSE, var.design = NULL, var.group = NULL, prior.n = 10,
     lib.size = NULL, offset = NULL, offset.prior = NULL, normalize.method = "none",
     span = 0.5, adaptive.span = TRUE, plot = FALSE, save.plot = FALSE,
     keep.EList = TRUE, nthreads = 1L)
}
\arguments{
  \item{counts}{
    a numeric \code{matrix} containing raw counts, or a \code{DGEList} object, or a \code{SummarizedExperiment} object containing raw counts.
    Counts must be non-negative. Fractional counts are permitted but NAs are not.
  }
  \item{design}{
    design matrix with rows corresponding to samples and columns to coefficients to be estimated.
    Defaults to the unit vector meaning that samples are treated as replicates.
  }
  \item{contrasts}{
    optional numeric matrix specifying coefficient contrasts.
    Rows correspond to coefficients (columns of \code{design}) and each column is a contrast.
    If specified, then the output object will contain the contrasts instead of the original coefficients.
	May be a vector if there is only one contrast.
	NAs are not allowed.
  }
  \item{block}{
    vector or factor specifying a blocking variable on the samples.
    Has length equal to \code{ncol(counts)}.
    Samples within each block are assumed to be correlated.
  }
  \item{prior.weights}{
    prior weights.
    Can be a numeric matrix of individual weights of same dimensions as the \code{counts},
    or a numeric vector of sample weights with length equal to \code{ncol(counts)},
    or a numeric vector of gene weights with length equal to \code{nrow(counts)}.
  }
  \item{sample.weights}{
    logical value, if \code{TRUE} then empirical sample quality weights will be estimated.
  }
  \item{var.design}{
    optional design matrix for the sample weights.
    Defaults to the sample-specific model whereby each sample has a distinct variance.
  }
  \item{var.group}{
    optional vector or factor indicating groups to have different array weights.
    This is another way to specify \code{var.design} for groupwise sample weights.
  }
  \item{prior.n}{
    prior number of genes for squeezing the weights towards equality.
    Larger values squeeze the sample weights more strongly towards equality.
  }
  \item{lib.size}{
    numeric vector containing total library sizes for each sample.
    Defaults to the normalized (effective) library sizes in \code{counts} if \code{counts} is a \code{DGEList} or to the columnwise count totals if \code{counts} is a matrix.
  }
  \item{offset}{
    an edgeR-style offset matrix of same dimensions as \code{counts} specifying observation-specific effective library sizes.
    Provides a way to input observation-wise normalization schemes.
    If not NULL, then the library sizes will be set to \code{exp(offset)}.
   }
  \item{offset.prior}{
    an alternative way to specify the offset matrix.
    \code{offset.prior} specifies effective library sizes relative to the raw library sizes instead of replacing the raw library sizes.
    The rows of \code{offset.prior} should sum to zero on average.
   }
  \item{normalize.method}{
    the microarray-style normalization method to be applied to the logCPM values (if any).
    Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.
    Any normalization factors found in \code{counts} will still be used even if \code{normalize.method="none"}.
  }
  \item{span}{
    width of the smoothing window used for the lowess mean-variance trend, expressed as a proportion between 0 and 1.
    Ignored if \code{adaptive.span=TRUE}.
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
  \item{keep.EList}{
    logical. If \code{TRUE}, then the normalized log2-CPM values and voom weights will be saved in the component \code{EList} of the output object.
  }
  \item{nthreads}{
    positive integer giving the number of OpenMP threads used by native-code fitting steps. Defaults to 1.
  }
}

\details{
"voom" is an acronym for "mean-variance modeling at the observational level".
This function processes count data from technologies such as RNA-seq or ChIP-seq to make it ready for linear modeling in limma.
The idea is to estimate the mean-variance relationship in the data, then use this to compute an appropriate precision weight for each observation.
Count data always show marked mean-variance relationships.
Raw counts show increasing variance with increasing count size, while log-counts typically show a decreasing mean-variance trend.
This function estimates the mean-variance trend for log-counts, then assigns a weight to each observation based on its predicted variance.
The weights are then used in the linear modeling process to adjust for heteroscedasticity. 
The mean-variance trend is estimated from gene-level data but is extrapolated back to individual observations to obtain a precision weight (inverse variance) for each observation.

The original voom method was described by Law et al (2014) and was implemented in the \code{\link{voom}} function.
Later, voom was combined with empirical sample quality weights by Liu et al (2015) and implemented in the \code{\link{voomWithQualityWeights}} function.
\code{voomLmFit} supersedes those earlier functions with additional features and ease of use (Baldoni et al, 2025).
It is particularly designed to give improved variance modeling for genes with large numbers of zero counts.
It allows for loss of residual degrees of freedom due to exact zero fitted values using the method of Lun & Smyth (2017).
The loss of residual degrees of freedom occurs when all the counts in a group are zero or when there are blocking factors that can fit zero counts exactly, meaning that not all the counts contribute effectively to variance estimation.
The refinement implemented in \code{voomLmFit} ensures that the residual variances remain approximately unbiased, even for low count genes with lots of zeros, whereas the residual variances could be underestimated by the original \code{voom} function for such genes.

This function is analogous to calling \code{voom} followed by \code{arrayWeights}, \code{duplicateCorrelation}, \code{lmFit} and \code{contrasts.fit}, except for the modified \code{df.residual} and \code{sigma} values.
This function returns \code{df.residual} values that are less than or equal to those from \code{voom} and \code{lmFit}, and \code{sigma} values that are greater than or equal to those from \code{lmFit}.
\code{voomLmFit} is more robust to zero counts than calling \code{voom} and \code{lmFit} separately and provides more rigorous error rate control.

\code{voomLmFit} replaces the original limma functions \code{voom}, \code{arrayWeights}, \code{duplicateCorrelation}, \code{lmFit} and \code{contrasts.fit} for the purpose of modeling of count data.
Compared to the original \code{voom} function, \code{voomLmFit} adds the following functionality:
  \itemize{
    \item Sample weights are estimated if \code{sample.weights=TRUE} or if \code{var.design} or \code{var.group} are non-NULL.
    \item An intra-block correlation is estimated if \code{block} is non-NULL.
    \item \code{residual.df} will be decreased for genes where all the counts are zero for particular experimental conditions.
    As a result, different genes may have different residual df, whereas the residual df from \code{voom} were all equal.
    \item The output is a \code{MArrayLM} linear model object rather than an \code{EList} data object.
    \item If \code{contrasts} is non-NULL, then the output object will contain contrasts rather than coefficients, superseding the need for downstream use of \code{contrasts.fit}.
    Importantly, the contrast standard errors from \code{voomLmFit} are exact whereas those from \code{contrasts.fit} are often approximate.
  }

\code{voomLmFit} performs the following specific calculations.
First, the counts are converted to logCPM values, adding 0.5 to all the counts to avoid taking the logarithm of zero.
The logCPM calculation uses normalized library sizes if \code{counts} is a DGEList or simply the column sums if \code{counts} is a matrix.
A microarray-style normalization method can also be optionally applied to the matrix of logCPM values.
The \code{lmFit} function is used to fit row-wise linear models.
The \code{lowess} function is then used to fit a trend to the square-root residual standard deviations as a function of an average log-count measure.
The trend line is then used to predict the variance of each logCPM value as a function of its fitted value on the count scale, and the inverse variances become the estimated precision weights.

The smoothness of the trend line is determined by \code{span}.
If \code{adaptive.span} is \code{TRUE}, then an optimal value for \code{span} is estimated by \code{chooseLowessSpan(nrow(counts), small.n=50, min.span=0.3, power=1/3)}.

If \code{sample.weights=TRUE}, then empirical sample quality weights are estimated using the \code{arrayWeights} function.
If \code{var.group} is specified, then a unique quality weight will be estimated for each sample group rather than for every individual sample.
Alternatively, if \code{var.design} is a numeric vector or design matrix, then the sample weights will be modeled as a log-linear function of that argument.
This functionality allows the sample weights to be modeled, for example, as a function of prior numeric quality scores for the samples.
The \code{lmFit} and \code{arrayWeights} fits will be iterated twice to ensure approximate convergence.

If \code{block} is non-NULL, then \code{duplicateCorrelation} will be called to estimate the intra-block correlation.
The \code{lmFit} and \code{duplicateCorrelation} fits will be iterated twice to ensure approximate convergence.

The \code{prior.weights} argument can be used to specify \emph{a priori} estimates of observational precision, if they exist, in which case the voom weights will modify these prior weights.
Most users will not need to use this argument.
If it is used, then the prior weights should arise from external annotation rather than being estimated from the data.

For good results, the \code{counts} matrix is usually filtered to remove rows with very low counts before running \code{voomLmFit}, and the \code{filterByExpr} function in the edgeR package can be used for that purpose.
However, \code{voomLmFit} is relatively robust to low counts and requires much less filtering for reliable results than did the original \code{voom} function.

If \code{counts} is a \code{DGEList} object from the edgeR package, then voomLmFit will use the normalization factors found in the object when computing the logCPM values.
In other words, the logCPM values are computed from the effective library sizes rather than the raw library sizes.
If the \code{DGEList} object contains library normalization factors, then it is usual to leave \code{normalize.method="none"} in voomLmFit, i.e., the logCPM values should not usually be re-normalized in the \code{voomLmFit} call.

If the \code{DGEList} object contains an offset matrix, then observation-specific effective library sizes will be set to \code{exp(offset)}.
If an \code{offset.prior} matrix is found, then \code{offset.prior} will be added to the log library sizes to make the offset matrix.

The \code{voom} method is similar in purpose to the limma-trend method, which uses \code{\link{eBayes}} or \code{\link{treat}} with \code{trend=TRUE}.
The voom method incorporates the mean-variance trend into the precision weights, whereas limma-trend incorporates the trend into the empirical Bayes moderation.
The voom method takes into account the sequencing depths (library sizes) of the individual columns of \code{counts} and applies the mean-variance trend on an individual observation basis.
limma-trend, on the other hand, assumes that the library sizes are not wildly different and applies the mean-variance trend on a genewise basis.
As noted by Law et al (2014), voom should be more powerful than limma-trend if the library sizes are very different but, otherwise, the two methods should give similar results.

\code{voomLmFit} is usually followed by running \code{\link{eBayes}} or \code{\link{treat}} on the fitted model object.
}

\value{
An MArrayLM object containing linear model fits for each row of data.
The object includes a \code{targets} data.frame component containing sample annotation and library sizes.
Columns of \code{targets} include \code{lib.size} and \code{sample.weight} (if \code{sample.weights=TRUE}).

If \code{block} was set, then the \code{correlation} component stores the consensus intra-block correlation.

If \code{save.plot=TRUE} then the output object will include components \code{voom.xy} and \code{voom.line}.
\code{voom.xy} contains the x and y coordinates of the points in the voom mean-variance plot in the same format as produced by \code{\link{xy.coords}} and \code{voom.line} contains the estimated trend curve.

If \code{keep.EList=TRUE} then the output object includes component \code{EList}, which is an \code{EList} object in the same format as produced by \code{voom} containing the voom log2-CPM values and the voom weights.
}

\author{
Original version by Gordon Smyth as a function in the edgeR package.
Self-contained limma version by Lizhong Chen, including the \code{contrasts} and \code{nthreads} arguments and C code backends for \code{lmFit}, \code{arrayWeights} and \code{duplicateCorrelation}.
}

\references{
Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{Nucleic Acids Research} 53(22), gkaf1305.
\doi{10.1093/nar/gkaf1305}.

Law CW, Chen Y, Shi W, Smyth GK (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
\doi{10.1186/gb-2014-15-2-r29}.
See also the Preprint Version at \url{https://gksmyth.github.io/pubs/VoomPreprint.pdf} incorporating some notational corrections.

Lun ATL, Smyth GK (2017).
No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data.
\emph{Statistical Applications in Genetics and Molecular Biology} 16(2), 83-93.
\doi{10.1515/sagmb-2017-0010}.

Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat ML, Smyth GK, Ritchie ME (2015).
Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses.
\emph{Nucleic Acids Research} 43, e97.
\doi{10.1093/nar/gkv412}.
}

\seealso{
\code{\link{voom}},
\code{\link{lmFit}},
\code{\link{contrasts.fit}},
\code{\link{voomWithQualityWeights}},
\code{\link{duplicateCorrelation}},
\code{\link{arrayWeights}},
\code{\link{MArrayLM-class}},
\code{\link{chooseLowessSpan}}.
}

\concept{Model fit}
