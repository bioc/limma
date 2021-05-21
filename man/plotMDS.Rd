\title{Multidimensional scaling plot of distances between gene expression profiles}
\name{plotMDS}
\alias{plotMDS}
\alias{plotMDS.MDS}
\alias{plotMDS.default}
\alias{MDS-class}
\alias{show,MDS-method}

\description{
Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples.
}

\usage{
\method{plotMDS}{default}(x, top = 500, labels = NULL, pch = NULL, cex = 1,
     dim.plot = c(1,2), gene.selection = "pairwise",
     xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE, \dots)
\method{plotMDS}{MDS}(x, labels = NULL, pch = NULL, cex = 1, dim.plot = NULL,
     xlab = NULL, ylab = NULL, var.explained = TRUE, \dots)
}

\arguments{
  \item{x}{any data object which can be coerced to a matrix, for example an \code{ExpressionSet} or an \code{EList}. Rows represent genes or genomic features while columns represent samples.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. Defaults to \code{colnames(x)}.}
  \item{pch}{plotting symbol or symbols. See \code{\link{points}} for possible values. Ignored if \code{labels} is non-\code{NULL}.}
  \item{cex}{numeric vector of plot symbol expansions.}
  \item{dim.plot}{integer vector of length two specifying which principal components should be plotted.}
  \item{gene.selection}{character, \code{"pairwise"} to choose the top genes separately for each pairwise comparison between the samples or \code{"common"} to select the same genes for all comparisons.}
  \item{xlab}{title for the x-axis.}
  \item{ylab}{title for the y-axis.}
  \item{plot}{logical. If \code{TRUE} then a plot is created on the current graphics device.}
  \item{var.explained}{logical. If \code{TRUE} then the percentage variation explained is included in the axis labels.}
  \item{\dots}{any other arguments are passed to \code{plot}, and also to \code{text} (if \code{pch} is \code{NULL}).}
}

\details{
This function uses multidimensional scaling (MDS) to produce a principal coordinate (PCoA) or principal component (PCA) plot showing the relationships between the expression profiles represented by the columns of \code{x}.
If \code{gene.selection = "common"}, or if the \code{top} is equal to or greater than the number of rows of \code{x}, then a PCA plot is constructed from the \code{top} genes with largest standard deviations across the samples.

If \code{gene.section = "pairwise"} and \code{top} is less than \code{nrow(x)} then a PCoA plot is produced and distances on the plot represent the \emph{leading log2-fold-changes}.
The leading log-fold-change between a pair of samples is defined as the root-mean-square average of the \code{top} largest log2-fold-changes between those two samples.
The PCA and PCoA plots produced by \code{gene.selection="common"} and \code{gene.selection="pairwise"}, respectively, use similar distance measures but the PCA plot uses the same genes throughout whereas the PCoA plot potentially selects different genes to distinguish each pair of samples.
The pairwise choice is the default.
It potentially gives better resolution than a PCA plot if different molecular pathways are relevant for distinguishing different pairs of samples.

If \code{pch=NULL}, then each sample is represented by a text label, defaulting to the column names of \code{x}.
If \code{pch} is not \code{NULL}, then plotting symbols are used.

See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.
}

\value{
If \code{plot=TRUE} or if \code{x} is an object of class \code{"MDS"}, then a plot is created on the current graphics device.

An object of class \code{"MDS"} is also invisibly returned.
This is a list containing the following components:
\item{eigen.values}{eigen values}
\item{eigen.vectors}{eigen vectors}
\item{var.explained}{proportion of variance explained by each dimension}
\item{distance.matrix.squared}{numeric matrix of squared pairwise distances between columns of \code{x}}
\item{dim.plot}{dimensions plotted}
\item{x}{x-xordinates of plotted points}
\item{y}{y-cordinates of plotted points}
\item{gene.selection}{gene selection method}
}

\author{Di Wu and Gordon Smyth}

\references{
Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, and Smyth GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\seealso{

An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}

\examples{
# Simulate gene expression data for 1000 probes and 6 microarrays.
# Samples are in two groups
# First 50 probes are differentially expressed in second group
sd <- 0.3*sqrt(4/rchisq(1000,df=4))
x <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(x) <- paste("Gene",1:1000)
x[1:50,4:6] <- x[1:50,4:6] + 2
# without labels, indexes of samples are plotted.
mds <- plotMDS(x,  col=c(rep("black",3), rep("red",3)) )
# or labels can be provided, here group indicators:
plotMDS(mds,  col=c(rep("black",3), rep("red",3)), labels= c(rep("Grp1",3), rep("Grp2",3)))
}

\concept{plots}
