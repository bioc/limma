\title{Top table of differentially spliced genes or exons}
\name{topSplice}
\alias{topSplice}
\description{
Top table ranking the most differentially spliced genes or exons.
}
\usage{
topSplice(fit, coef = ncol(fit), test = "F", number = 10, FDR=1, sort.by = "p")
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}.}
  \item{coef}{the coefficient (column) of fit for which differentially splicing is assessed.}
  \item{test}{character string specifying which statistical test to apply.
    Possible values are \code{"simes"}, \code{"F"} or \code{"t"}.
    \code{"F"} gives F-tests for each gene.
    \code{"t"} gives t-tests for each exon.
    \code{"simes"} gives genewise p-values derived from the t-tests after Simes adjustment for each gene.
  }
  \item{number}{integer, maximum number of rows to output.}
  \item{FDR}{numeric, only show exons or genes with false discovery rate less than this cutoff.}
  \item{sort.by}{character string specifying which column to sort results by.
    Possible values for \code{"p"}, \code{"logFC"}, \code{"NExons"} or \code{"none"}.
    \code{"logFC"} is only available if \code{test="t"} and \code{"NExons"} is only available if \code{test="simes"} or \code{test="F"}.
  }
}

\details{
Ranks genes or exons by evidence of differential exon-usage arising from differential splicing.
Choosing \code{test="F"} or \code{test="simes"} gives gene-level tests, where the null hypothesis is that all exons are used equally.
Choosing \code{test="t"} gives exon-level tests of usage for each exon vs all other exons for the same gene.

The F-statistic option chosen by \code{test="F"} conducts an ANOVA-like F-tests of differential usage for each gene.
Simes method chosen by \code{test="simes"} processes the exon-level p-values for each gene to give an overall call of differential splicing for that gene.
The F-tests are likely to be powerful for genes in which several exons are differentially splices.
The Simes p-values is likely to be more powerful when only a minority of the exons for a gene are differentially spliced.

The F-test option gives the best error rate control.
Simes method and the exon-level t-tests are both slightly liberal in terms of error rate control.
}

\value{A data.frame with any annotation columns found in \code{fit} plus the following columns
  \item{logFC}{log2-fold change of exon vs other exons for the same gene (if \code{level="exon"})}
  \item{t}{moderated t-statistic (if \code{level="exon"})}
  \item{F}{moderated F-statistic (if \code{level="gene"})}
  \item{P.Value}{p-value}
  \item{FDR}{false discovery rate}
}

\author{Gordon Smyth}
\seealso{
\code{\link{diffSplice}}, \code{\link{plotSplice}}

A summary of functions available in LIMMA for RNA-seq analysis is given in \link{11.RNAseq}.
}
\examples{# See diffSplice}

\keyword{rna-seq}
