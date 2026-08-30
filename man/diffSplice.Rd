\name{diffSplice}
\alias{diffSplice}
\alias{diffSplice.MArrayLM}
\title{Differential Transcript or Exon Usage}
\description{Given a limma fit object with multiple features per gene, test for differential usage of those features within genes between experimental conditions.
The features may be any set of isoform-identifying features such as transcripts, exons or exon-junctions.}
\usage{
\method{diffSplice}{MArrayLM}(fit, geneid, exonid = NULL,
     robust = FALSE, legacy = FALSE, verbose = TRUE, \dots)
}
\arguments{
  \item{fit}{an \code{MArrayLM} fitted model object produced by \code{voomLmFit}, \code{lmFit} or \code{contrasts.fit}. Rows should correspond to transcripts, for a DTU analysis, or to exons and exon-exon junctions for a DEU analysis.}
  \item{geneid}{gene identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the gene identifiers. Rows with the same ID are assumed to belong to the same gene.}
  \item{exonid}{feature or exon identifiers. Either a vector of length \code{nrow(fit)} or the name of the column of \code{fit$genes} containing the feature identifiers.}
  \item{robust}{logical, should the estimation of the empirical Bayes prior parameters be robustified against outlier sample variances?}
  \item{legacy}{logical. If \code{FALSE} then the new empirical Bayes hyperparameter estimation (introduced in limma 3.61.8) will be used, if \code{TRUE} the earlier hyperparameter estimation will be used. The new method is particularly appropriate when the residual degrees of freedom are not all equal, which is likely to be the case for \code{diffSplice}.}
  \item{verbose}{logical, if \code{TRUE} some diagnostic information about the number of genes and exons is output.}
  \item{\dots}{other arguments are not currently used.}
}
\value{
An object of class \code{MArrayLM} containing both exon level and gene level tests.
Results are sorted by geneid and by exonid within gene.
  \item{coefficients}{numeric matrix of coefficients of same dimensions as \code{fit}. Each coefficient is the difference between the log-fold-change for that exon versus the average log-fold-change for all other exons for the same gene.}
  \item{t}{numeric matrix of moderated t-statistics, of same dimensions as \code{fit}.}
  \item{p.value}{numeric vector of p-values corresponding to the t-statistics}
  \item{genes}{data.frame of exon annotation}
  \item{genecolname}{character string giving the name of the column of \code{genes} containing gene IDs}
  \item{gene.F}{numeric matrix of moderated F-statistics, one row for each gene.}
  \item{gene.F.p.value}{numeric matrix of p-values corresponding to \code{gene.F}}
  \item{gene.simes.p.value}{numeric matrix of Simes adjusted p-values, one row for each gene.}
  \item{gene.bonferroni.p.value}{numeric matrix of Bonferroni adjusted p-values, one row for each gene.}
  \item{gene.genes}{data.frame of gene annotation.}
}

\details{
Given a limma fit object with multiple features per gene, test for differential usage of those features within genes between experimental conditions.
The features may be any set of isoform-identifying features such as transcripts, exons or exon-junctions.
The features can be transcripts for a differential transcript usage (DTU) analysis (Baldoni et al 2025), or can be a combination of exons and exon-exon junctions for a differential exon usage (DEU) analysis.
Another possibility is for the features to be transcript-overlap equivalence classes (Cmero et al 2019).

Testing for differential usage is equivalent to testing whether the log-fold-changes in the \code{fit} differ between features for the same gene.
Two different tests are provided.
The first is a F-test for differences between the log-fold-changes for each gene.
This is equivalent to testing for interaction between the features for that gene and the coefficient of the linear model.
The other is a series of t-tests in which each transcript is compared to the weighted average of all other features for the same gene.
The feature-level t-tests are converted into genewise tests by adjusting the p-values for the same gene by Simes' method.
As an alternative, the feature-level t-tests are also converted into genewise tests by adjusting the smallest p-value for each gene by Bonferroni's method.

This function can be used on transcript level RNA-seq counts from RSEM, kallisto, or Salmon, after pre-processing by the edgeR functions \code{catchSalmon()} (or \code{catchKallisto} or \code{catchRSEM}) and \code{voomLmFit()}, as described by Baldoni et al (2025).
It can also be used on equivalence-class counts from Salmon or kallisto, after pre-processing by voomLmFit(), as described by Cmero et al (2019).
It can also be used on exon-level read counts or on data from an exon microarray.
}

\note{
This function is not designed for situations with a very high level of multi-counting of RNA-seq reads that overlap two or more exons.
In particular, it is not designed for use with "chopped" exons, where overlapping exons belonging to different transcripts of the same gene are chopped up into unique sub-exons, because artificial exons of this sort lead to high levels of multi-counting.
}

\seealso{
\code{\link{topSplice}} and \code{\link{plotSplice}} are downstream functions that operate on the output from \code{diffSplice}.

Also see \code{diffSplice.DGEGLM} in the edgeR package, which has comparable functionality but for edgeR fit objects.

A summary of functions available in LIMMA for RNA-seq analysis is given in \link{11.RNAseq}.
}

\author{Gordon Smyth and Charity Law}

\references{
Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{Nucleic Acids Research} 53(22), gkaf1305.
\doi{10.1093/nar/gkaf1305}.

Cmero M, Davidson NM, Oshlack A (2019).
Using equivalence class counts for fast and accurate testing of differential transcript usage.
\emph{F1000Research} 8, 265.
\doi{10.12688/f1000research.18276.2}.
}

\examples{
\dontrun{
fit <- voomLmFit(dge, design)
ex <- diffSplice(fit, geneid="GeneID")
topSplice(ex)
plotSplice(ex, xlab="Transcript")
}
}

\concept{RNA-seq}
\concept{Differential transcript usage}
\concept{Differential usage}
