\name{camera}
\alias{camera}
\alias{camera.default}
\alias{interGeneCorrelation}
\alias{cameraPR}
\alias{cameraPR.default}
\title{Competitive Gene Set Test Accounting for Inter-gene Correlation}
\description{
Test whether a set of genes is highly ranked relative to other genes in terms of differential expression, accounting for inter-gene correlation.
}
\usage{
\method{camera}{default}(y, index, design, contrast = ncol(design), weights = NULL,
       use.ranks = FALSE, allow.neg.cor=FALSE, inter.gene.cor=0.01, trend.var = FALSE,
       sort = TRUE, directional = TRUE, \dots)
\method{cameraPR}{default}(statistic, index, use.ranks = FALSE, inter.gene.cor=0.01,
       sort = TRUE, directional = TRUE, \dots)
interGeneCorrelation(y, design)
}
\arguments{
  \item{y}{
    a numeric matrix of log-expression values or log-ratios of expression values, or any data object containing such a matrix.
    Rows correspond to probes and columns to samples.
    Any type of object that can be processed by \code{\link{getEAWP}} is acceptable.
    \code{NA} or infinite values are not allowed.
  }
  \item{statistic}{
    a numeric vector of genewise statistics.
    If \code{index} contains gene IDs, then \code{statistic} should be a named vector with the gene IDs as names.
  }
  \item{index}{
    a list of gene sets.
    Each gene set is either an index vector \code{i} or a data.frame where the first column is the index vector and the second column gives the positive or negative change magnitude associated with each gene in the gene set.
    The index vector can contain gene IDs or row numbers, provided that \code{y[i,]} of \code{statistic[i]} will subset to rows corresponding to genes in the gene set.
  }
  \item{design}{design matrix.}
  \item{contrast}{contrast of the linear model coefficients for which the test is required. Can be an integer specifying a column of \code{design}, or else a numeric vector of same length as the number of columns of \code{design}.}
  \item{weights}{numeric matrix of precision weights. Can be a matrix of the same size as \code{y}, or a numeric vector of array weights with length equal to \code{ncol(y)}, or a numeric vector of gene weights with length equal to \code{nrow(y)}.}
  \item{use.ranks}{do a rank-based test (\code{TRUE}) or a parametric test (\code{FALSE})?}
  \item{allow.neg.cor}{should reduced variance inflation factors be allowed for negative correlations?}
  \item{inter.gene.cor}{numeric, optional preset value for the inter-gene correlation within tested sets.  If \code{NA} or \code{NULL}, then an inter-gene correlation will be estimated for each tested set.}
  \item{trend.var}{logical, should an empirical Bayes trend be estimated?  See \code{\link{eBayes}} for details.}
  \item{sort}{logical, should the results be sorted by p-value?}
  \item{directional}{are the gene sets directional? If \code{TRUE}, will test for genes changing in the same direction within each set. If \code{FALSE}, will test for large effects without regard to direction in each set.}
  \item{\dots}{other arguments are not currently used}
}

\details{
\code{camera} and \code{interGeneCorrelation} implement methods proposed by Wu and Smyth (2012).
\code{camera} performs a \emph{competitive} test in the sense defined by Goeman and Buhlmann (2007).
It tests whether the genes in the set are highly ranked in terms of differential expression relative to genes not in the set.
It has similar aims to \code{geneSetTest} but accounts for inter-gene correlation.
See \code{\link{roast}} for an analogous \emph{self-contained} gene set test.

The function can be used for any microarray experiment which can be represented by a linear model.
The design matrix for the experiment is specified as for the \code{\link{lmFit}} function, and the contrast of interest is specified as for the \code{\link{contrasts.fit}} function.
This allows users to focus on differential expression for any coefficient or contrast in a linear model by giving the vector of test statistic values.

The argument \code{index} specifies one or more gene sets.
The gene sets can be simple lists of genes, as indicated by a vector of gene IDs.
For such gene sets, all the genes in the set are expected to change in the same direction and by the same amount when the molecular pathway is active.
Alternatively, the gene sets can be provided as expression signatures, whereby an expected directional change is provided for each gene in the gene set.
For example, the directional changes could be the log-fold-changes observed in previous experiments when a particular biological pathway was perturbed.
Expression signature gene sets are provided as data.frames, with the gene indices in the first column and the directional changes in the second column.

One convenient way to construct the \code{index} argument is by using \code{\link{ids2indices}}.

\code{camera} estimates p-values after adjusting the variance of test statistics by an estimated variance inflation factor.
The inflation factor depends on estimated genewise correlation and the number of genes in the gene set.

By default, \code{camera} uses \code{interGeneCorrelation} to estimate the mean pair-wise correlation within each set of genes.
\code{camera} can alternatively be used with a preset correlation specified by \code{inter.gene.cor} that is shared by all sets.
This usually works best with a small value, say \code{inter.gene.cor=0.01}.

If \code{inter.gene.cor=NA}, then \code{camera} will estimate the inter-gene correlation for each set.
In this mode, \code{camera} gives rigorous error rate control for all sample sizes and all gene sets.
However, in this mode, highly co-regulated gene sets that are biologically interpretable may not always be ranked at the top of the list.

With the default value \code{inter.gene.cor=0.01}, \code{camera} will rank biologically interpretable sets more highly.
This gives a useful compromise between strict error rate control and interpretable gene set rankings.

\code{cameraPR} is a "pre-ranked" version of \code{camera} where the genes are pre-ranked according to a pre-computed statistic.

If \code{directional=TRUE}, then the gene sets are assumed to contain genes changing in the same direction and the functions will do two-sided directional tests for each gene set.
If \code{directional=FALSE}, then the gene sets are assumed to be non-directional, so that significant genes in each set could be changing in different directions.
In this case, the functions will test for large changes in each set without regard to direction of change.
If \code{directional=FALSE}, then \code{use.ranks} should be \code{TRUE} and \code{statistic} should be non-negative.
}

\note{The default settings for \code{inter.gene.cor} and \code{allow.neg.cor} were changed to the current values in limma 3.29.6.
Previously, the default was to estimate an inter-gene correlation for each set.
To reproduce the earlier default, use \code{allow.neg.cor=TRUE} and \code{inter.gene.cor=NA}.}

\value{
\code{camera} and \code{cameraPR} return a data.frame with a row for each set and the following columns:
\item{NGenes}{number of genes in set.}
\item{Correlation}{inter-gene correlation (only included if the \code{inter.gene.cor} was not preset).}
\item{Direction}{direction of change (\code{"Up"} or \code{"Down"}).}
\item{PValue}{two-tailed p-value.}
\item{FDR}{Benjamini and Hochberg FDR adjusted p-value.}

\code{interGeneCorrelation} returns a list with components:
\item{vif}{variance inflation factor.}
\item{correlation}{inter-gene correlation.}
}

\author{Di Wu and Gordon Smyth}

\references{
Wu D, Smyth GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation.
\emph{Nucleic Acids Research} 40, e133.
\doi{10.1093/nar/gks461}.

Goeman JJ, Buhlmann P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 
}

\seealso{
\code{\link{getEAWP}}

\code{\link{rankSumTestWithCorrelation}},
\code{\link{geneSetTest}},
\code{\link{roast}},
\code{\link{fry}},
\code{\link{romer}},
\code{\link{ids2indices}}.

There is a topic page on \link{10.GeneSetTests}.
}

\examples{
y <- matrix(rnorm(1000*6),1000,6)
rownames(y) <- paste0("G",1:1000)
colnames(y) <- paste0("S",1:6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))

# First set of 20 genes are genuinely up-regulated
index1 <- 1:20
y[index1,4:6] <- y[index1,4:6] + 1

# Second set of 20 genes are genuinely down-regulated
index2 <- 21:40
y[index2,4:6] <- y[index1,4:6] - 1

# Other genes are not DE
index3 <- 41:60
 
# Directional gene set containing both up and down genes
index4 <- data.frame(
  ID = c("G1","G2","G3","G21","G22","G23"),
  logFC = c(1,1,1,-1,-1,-1)
)
 
index <- list(set1=index1, set2=index2, set3=index3, set4=index4)
camera(y, index, design, inter.gene.cor=NA, sort=FALSE)
camera(y, index, design, inter.gene.cor=0.01, sort=FALSE)

# Pre-ranked version
fit <- eBayes(lmFit(y, design))
logFC <- fit$coefficient[,2]
cameraPR(logFC, index, sort=FALSE)
}

\concept{gene set tests}
