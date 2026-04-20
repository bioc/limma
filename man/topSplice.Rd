\title{Top Table of Differential Transcript Usage}
\name{topSplice}
\alias{topSplice}
\description{
Top table of results from a \code{diffSplice} differential transcript usage analysis.
}
\usage{
topSplice(fit, coef = ncol(fit), test = "F", number = 10, FDR=1, sort.by = "p", treat.lfc = 0)
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}.}
  \item{coef}{the coefficient (column) of fit for which differential usage is assessed.}
  \item{test}{character string specifying which statistical test to apply.
    Possible values are \code{"simes"}, \code{"F"} or \code{"t"}.
    \code{"F"} gives F-tests for each gene.
    \code{"t"} gives t-tests for each transcript or exon.
    \code{"simes"} gives genewise p-values derived from the t-tests after Simes adjustment for each gene.
  }
  \item{number}{integer, maximum number of rows to output.}
  \item{FDR}{numeric, only show exons or genes with false discovery rate less than this cutoff.}
  \item{sort.by}{character string specifying which column to sort results by.
    Possible values for \code{"p"}, \code{"logFC"}, \code{"NExons"} or \code{"none"}.
    \code{"logFC"} is only available if \code{test="t"} and \code{"NExons"} is only available if \code{test="simes"} or \code{test="F"}.
  }
  \item{treat.lfc}{
    a treat-style log2-fold-change threshold to apply to significant testing when \code{test="t"}.
    The null hypothesis will be that the true log2-fold-change is no larger than this theshold.
  }
}

\details{
This function ranks results from a differential transcript usage analysis conducted by \code{diffSplice}.
The output object from \code{diffSplice} includes results at the gene level and results at the isoform level.
The isoform-level results can be transcripts, if a DTU analysis was done, or can be exons and exon-junctions, if a DEU analysis was done.

Choosing \code{test="F"} or \code{test="simes"} gives gene-level results, where the null hypothesis is that all transcripts or exons are used equally.
Choosing \code{test="t"} gives transcript or exon level results.

The F-statistic option chosen by \code{test="F"} conducts an ANOVA-like F-tests of differential usage for each gene.
Simes method chosen by \code{test="simes"} processes the isoform-level p-values for each gene to give an overall call of differential usage for that gene.
The Simes p-values is likely to be more powerful than F-tests when only one transcript of the gene is differentially expressed.
The F-tests are likely to be powerful than Simes when there multiple transcripts of the same gene differentially expressed in different directions.
}

\value{A data.frame with annotation columns found in \code{fit} plus the following columns
  \item{logFC}{log2-fold change of transcript vs other transcripts for the same gene (if \code{test="t"})}
  \item{t}{moderated t-statistic (if \code{test="t"})}
  \item{F}{moderated F-statistic (if \code{test="F"})}
  \item{P.Value}{p-value}
  \item{FDR}{false discovery rate bound, computed by the Benjamini & Hochberg method.}
}

\author{Gordon Smyth}

\references{
Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{bioRxiv}
\doi{10.1101/2025.04.07.647659}.
}

\seealso{
\code{\link{diffSplice}}, \code{\link{plotSplice}}.
See also the \code{diffSplice} method for \code{DGEGLM} objects in the edgeR package.

A summary of functions available in LIMMA for RNA-seq analysis is given in \link{11.RNAseq}.
}

\examples{# See diffSplice}

\keyword{rna-seq}
\concept{differential usage}
