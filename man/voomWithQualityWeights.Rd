\name{voomWithQualityWeights}
\alias{voomWithQualityWeights}
\title{Transform RNA-Seq Counts for Linear Modeling With Precision and Sample Weights (Deprecated)}
\description{
Combine voom observational-level precision weights with sample-specific quality weights in a designed experiment.
}
\usage{
voomWithQualityWeights(counts, design = NULL, lib.size = NULL, normalize.method = "none",
             plot = FALSE, span = 0.5, adaptive.span = TRUE,
             var.design = NULL, var.group = NULL,
             method = "genebygene", maxiter = 50, tol = 1e-5, trace = FALSE,
             col = NULL, \dots)
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
    numeric vector containing total library sizes for each sample.
    If \code{NULL} and \code{counts} is a \code{DGEList}, then the normalized library sizes are taken from \code{counts}.
    Otherwise library sizes are calculated from the columnwise count totals.
  }
  \item{normalize.method}{
    the microarray-style normalization method to be applied to the logCPM values.
    Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.
  }
  \item{plot}{
    logical, should a plot of the mean-variance trend and sample-specific weights be displayed?
  }
  \item{span}{
    width of the smoothing window used for the lowess mean-variance trend.
    Expressed as a proportion between 0 and 1.
  }
  \item{adaptive.span}{
    logical.
    If \code{TRUE}, then an optimal value for \code{span} will be chosen depending on the number of genes, and the value set by the \code{span} argument will be ignored.
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
    colors to use in the barplot of sample-specific weights if \code{plot=TRUE}.
    If \code{NULL}, then bars are plotted in grey.
  }
  \item{\dots}{
    other arguments are passed to \code{voom} and hence to \code{lmFit}.
  }
}

\details{
This was the original function implementing the "voom" method with sample weights, although we now recommend \code{\link{voomLmFit}} instead.
\code{\link{voomLmFit}} is a further developed version of \code{voomWithQualityWeights} with extra capabilities and ease of use.
\code{voomLmFit} automates the computation of sample weights, intrablock correlation and contrasts, and also improves variance estimation for sparse data.

Like \code{voom}, \code{voomWithQualityWeights} is intended to process RNA-seq data prior to linear modeling in limma.
It combines observational-level weights from \code{voom} with sample-specific weights estimated using the \code{arrayWeights} function.
The method is described by Liu et al (2015).
It consists of a call to \code{voom}, followed by \code{arrayWeights}, and then another iteration of both functions.
}

\value{
An \code{\link[limma:EList]{EList}} object similar to that from \code{\link{voom}},
with an extra column, \code{sample.weights}, added to the \code{targets} data.frame, containing the vector of sample quality factors.
The \code{weights} component combines the sample weights and the usual voom precision weights.
}

\author{Matthew Ritchie, Cynthia Liu, Gordon Smyth}

\references{
Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat ML, Smyth GK, Ritchie ME (2015).
Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses.
\emph{Nucleic Acids Research} 43(15), e97.
\doi{10.1093/nar/gkv412}.
}

\seealso{
\code{\link{voom}}, \code{\link{arrayWeights}}, \code{\link{voomLmFit}}.

A summary of limma functions for RNA-seq analysis is given in \link{11.RNAseq}.
}

\keyword{regression}
\keyword{rna-seq}
