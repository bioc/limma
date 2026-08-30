\title{Differential usage plot}
\name{plotSplice}
\alias{plotSplice}
\description{
Plot relative log-fold changes by exons (or transcripts) for the specified gene and highlight differentially used exons (or transcripts).
}
\usage{
plotSplice(fit, coef = ncol(fit), geneid = NULL, genecolname = NULL,
           rank = 1L, FDR = 0.05, xlab = "Exon")
}
\arguments{
  \item{fit}{\code{MArrayLM} fit object produced by \code{diffSplice}.}
  \item{coef}{the coefficient (column) of fit for which differential splicing is assessed.}
  \item{geneid}{character string, ID of the gene to plot.}
  \item{genecolname}{column name of \code{fit$genes} containing gene IDs. Defaults to \code{fit$genecolname}.}
  \item{rank}{integer, if \code{geneid=NULL} then this ranked gene will be plotted.}
  \item{FDR}{numeric, highlight exons as red dots with false discovery rate less than this cutoff. The FDR of the individual exon is calculated based on the exon-level t-statistics test for differences between each exon and all other exons for the same gene.}
  \item{xlab}{x-axis lab for the plot. Typical choices are \code{"Exon"} for DEU analyses or \code{"Transcript"} for DTU analyses.}
}

\details{
Plot relative log2-fold-changes by exon (or transcript) for the specified gene.
The relative logFC is the difference between that exon (or transcript)'s logFC and the overall logFC for the gene, as computed by \code{diffSplice}.
Individual exons (or transcripts) that are significantly differentially used are highlighted as red dots.
The size of the red dots increases with the statistical significance of the exon (or transcript).
}

\value{A plot is created on the current graphics device.}

\author{Gordon Smyth and Yifang Hu}

\references{
Baldoni PL, Chen L, Li M, Chen Y, Smyth GK (2025).
Dividing out quantification uncertainty enables assessment of differential transcript usage with limma and edgeR.
\emph{bioRxiv}
\doi{10.1101/2025.04.07.647659}.

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\doi{10.1093/nar/gkv007}.
}

\seealso{
\code{\link{diffSplice}}, \code{\link{topSplice}}

A summary of functions available in LIMMA for RNA-seq analysis is given in \link{11.RNAseq}.
}

\examples{# See diffSplice}

\concept{differential usage}
\keyword{hplot}
\keyword{rna-seq}
