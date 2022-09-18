\name{goanaTrend}
\alias{goanaTrend}
\title{Estimate DE Trend for Gene Ontology or KEGG Pathway Analysis}

\description{
Given a list of differentially expressed (DE) genes and a covariate, estimate the probability of a gene being called significant as a function of the covariate.
This function is typically used to estimate the gene length or gene abundance bias for a pathway analysis.
}

\usage{
goanaTrend(index.de, covariate, n.prior = 10, plot = FALSE,
           xlab = "Covariate Rank", ylab = "Probability gene is DE",
           main="DE status vs covariate")
}

\arguments{
  \item{index.de}{an index vector specifying which genes are significantly DE.
  Can be a vector of integer indices, or a logical vector of length \code{nrow(covariate)}, or any vector such as \code{covariate[index]} selects the DE genes.}

  \item{covariate}{numeric vector, length equal to the number of genes in the analysis.
  Usually equal to gene length or average log-expression but can be any meaningful genewise covariate.}

  \item{n.prior}{prior number of genes using for moderating the trend towards constancy, for stability when the number of DE genes is small.}

  \item{plot}{if \code{TRUE}, plot the estimated tend.}
  \item{xlab}{label for x-axis of plot.}
  \item{ylab}{label for y-axis of plot.}
  \item{main}{main title for the plot.}
}

\details{
\code{goanaTrend} is called by \code{goana} and \code{kegga} when the \code{trend} argument is used
to correct for unwanted trends in the differential expression analysis associated with gene length, gene abundance or any other covariate (Young et al, 2010).

This function is analogous to the \code{nullp} function of the goseq package but the
trend is estimated using \code{\link{tricubeMovingAverage}} instead of by monotonic regression.
While \code{\link{tricubeMovingAverage}} does not enforce strict monotonicity, it has the advantage of numerical stability and statistical robustness when there are only a small number of DE genes.

This function also moderates the estimated trend slightly towards constancy to provide more stability.
The degree of moderation is determined by the \code{n.prior} argument relative to the number of DE genes.
}

\value{
Numeric vector of same length as \code{covariate} giving estimated probabilities.
}

\references{
  Young MD, Wakefield MJ, Smyth GK, Oshlack A (2010).
  Gene ontology analysis for RNA-seq: accounting for selection bias.
  \emph{Genome Biology} 11, R14.
  \doi{10.1186/gb-2010-11-2-r14}
}

\seealso{
\code{\link{goana}}, \code{\link{kegga}}

See \link{10.GeneSetTests} for a description of other functions used for gene set testing.
}

\author{Gordon Smyth and Yifang Hu}

\examples{
x <- runif(100)
i <- 1:10
goanaTrend(i, x, plot=TRUE)
}

\concept{gene set tests}
