\name{cumOverlap}
\alias{cumOverlap}
\title{Cumulative Overlap Analysis of Ordered Lists}
\description{
Test whether the leading members of ordered lists significantly overlap.
}
\usage{
cumOverlap(ol1, ol2)
}

\arguments{
  \item{ol1}{vector containing first ordered list.}
  \item{ol2}{vector containing second ordered list.}
}

\value{
List containing the following components:
\item{n.min}{integer, top table length leading to smallest adjusted p-value.}
\item{p.min}{smallest adjusted p-value.}
\item{n.overlap}{integer, number of overlapping IDs in first \code{n.min}.}
\item{id.overlap}{vector giving the overlapping IDs in first \code{n.min}.}
\item{p.value}{numeric, vector of p-values for each possible top table length.}
\item{adj.p.value}{numeric, vector of Bonferroni adjusted p-values for each possible top table length.}
}

\details{
The function compares the top \code{n} members of each list, for every possible \code{n}, and conducts an hypergeometric test for overlap.
The function returns the value of \code{n} giving the smallest Bonferroni adjusted p-value.

This method was described in Chapter 4 of Wu (2011).
}

\author{Gordon Smyth and Di Wu}

\references{
Wu, D (2011).
Finding hidden relationships between gene expression profiles with application to breast cancer biology.
PhD thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/36278}
}

\examples{
ol1 <- letters[1:26]
ol2 <- letters[sample(26)]
coa <- cumOverlap(ol1, ol2)
coa$p.min
}

\concept{gene set enrichment analysis}
