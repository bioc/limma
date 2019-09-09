\name{head}
\alias{head.RGList}
\alias{head.MAList}
\alias{head.EListRaw}
\alias{head.EList}
\alias{head.MArrayLM}
\alias{head.TestResults}
\alias{tail.RGList}
\alias{tail.MAList}
\alias{tail.EListRaw}
\alias{tail.EList}
\alias{tail.MArrayLM}
\alias{tail.TestResults}

\title{Return the First to Last Part of a Data Object}

\description{
Retrieve the first or last parts of an RGList, MAList, EListRaw, EList, MArrayLM or TestResults object.
}

\usage{
\method{head}{EList}(x, n = 6L, \dots)
\method{tail}{EList}(x, n = 6L, \dots)
}

\arguments{
  \item{x}{an object of class \code{RGList}, \code{MAList}, \code{EListRaw}, \code{EList}, \code{MArrayLM} or \code{TestResults}.}
  \item{n}{
    a single integer.
    If positive or zero, number rows of resulting object.
    If negative, all but the \code{n} last/first rows of \code{x}.
  }
  \item{\dots}{other arguments are not currently used.}
}

\details{
\code{head} (\code{tail}) returns the first (last) \code{n} rows when \code{n >= 0} or all but the last (first) \code{n} rows when \code{n < 0}.
}

\value{
An object like \code{x} but generally with fewer rows.
}

\author{Gordon Smyth}

\seealso{
  \code{\link{head}} in the utils package.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}

\examples{
E <- matrix(rnorm(40),20,2)
rownames(E) <- paste0("Gene",1:20)
colnames(E) <- c("A","B")
y <- new("EList",list(E=E))
head(y)
tail(y)
}

\keyword{array}
