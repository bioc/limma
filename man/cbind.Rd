\name{cbind}
\alias{cbind.RGList}
\alias{cbind.MAList}
\title{Combine RGList or MAList Objects}
\description{
Combine a series of \code{RGList} objects or combine a series of \code{MAList} objects.
}
\usage{
\method{cbind}{RGList}(\dots, deparse.level=1)
}
\arguments{
  \item{\dots}{\code{RGList} objects or \code{MAList} objects}
  \item{deparse.level}{not currently used, see \code{\link[base]{cbind}} in the base package}
}
\details{
The matrices of expression data from the individual objects are cbinded.
The data.frames of target information, if they exist, are rbinded.
The combined data object will preserve any additional components or attributes found in the first object to be combined.
}
\value{
An \code{\link[limma:rglist]{RGList}} or \code{\link[limma:malist]{MAList}} object holding data from all the arrays from the individual objects.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{cbind}} in the base package.
  
  \link{3.ReadingData} gives an overview of data input and manipulation functions in LIMMA.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("A1","A2")
MA1 <- new("MAList",list(M=M,A=A))

M <- A <- matrix(21:24,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("B1","B2")
MA2 <- new("MAList",list(M=M,A=A))

cbind(MA1,MA2)
}
\keyword{manip}