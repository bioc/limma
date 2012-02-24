\name{camera}
\alias{camera}
\alias{interGeneCorrelation}
\title{Competitive Gene Set Test Accounting for Inter-gene Correlation}
\description{
Test whether a set of genes is highly ranked relative to other genes in terms of differential expression, accounting for inter-gene correlation.
}
\usage{
camera(indices, y, design, contrast=ncol(design), statistic="modt", use.ranks=FALSE, trend.var=FALSE)
interGeneCorrelation(y, design)
}
\arguments{
  \item{indices}{an index vector or a list of index vectors.  Can be any vector such that \code{y[indices,]} selects the rows corresponding to the test set.}
  \item{y}{numeric matrix giving log-expression or log-ratio values for a series of microarrays, or any object that can coerced to a matrix including \code{ExpressionSet}, \code{MAList}, \code{EList} or \code{PLMSet} objects.
  Rows correspond to probes and columns to samples.}
  \item{design}{design matrix.}
  \item{contrast}{contrast of the linear model coefficients for which the test is required. Can be an integer specifying a column of \code{design}, or else a numeric vector of same length as the number of columns of \code{design}.}
  \item{statistic}{character string indicating the genewise statistic to use, either \code{"modt"} or \code{"logFC"}.}
  \item{use.ranks}{should a rank sum test (\code{TRUE}) or a parametric t-test (\code{FALSE} be done?}
  \item{trend.var}{logical, should an empirical Bayes trend be estimated?  See \code{\link{eBayes}} for details.}
} 

\value{
A numeric matrix with a row for each set and the following columns:
\item{NGenes}{number of genes in set}
\item{Correlation}{inter-gene correlation}
\item{Down}{left-tail p-value}
\item{Up}{right-tail p-value}
\item{TwoSided}{two-tailed p-value}
}

\details{
This function implements the camera gene set test proposed by Wu and Smyth (2012).
\code{camera} performs a \emph{competitive} test in the sense defined by Goeman and Buhlmann (2007).
It tests whether the genes in the set are highly ranked in terms of differential expression relative to genes not in the set.
It has similar aims to \code{geneSetTest} but accounts for inter-gene correlation.
See \code{\link{roast}} for an analogous \emph{self-contained} gene set test.

The function can be used for any microarray experiment which can be represented by a linear model.
The design matrix for the experiment is specified as for the \code{\link{lmFit}} function, and the contrast of interest is specified as for the \code{\link{contrasts.fit}} function.
This allows users to focus on differential expression for any coefficient or contrast in a linear model by giving the vector of test statistic values.

\code{camera} estimates p-values after adjusting the variance of test statistics by an estimated variance inflation factor.
The inflation factor depends on estimated genewise correlation and the number of genes in the gene set. 
}

\seealso{
\code{\link{rankSumTestWithCorrelation}},
\code{\link{geneSetTest}},
\code{\link{roast}},
\code{\link{romer}}.
}

\author{Di Wu and Gordon Smyth}

\references{
Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
Submitted.

Goeman, JJ, and Buhlmann, P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 
}

\examples{
y <- matrix(rnorm(1000*6),1000,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))

# First set of 20 genes are genuinely differentially expressed
iset1 <- 1:20
y[iset1,4:6] <- y[iset1,4:6]+1

# Second set of 20 genes are not DE
iset2 <- 21:40
 
camera(iset1, y, design)
camera(iset2, y, design)

camera(list(set1=iset1,set2=iset2), y, design)
}
\keyword{htest}
