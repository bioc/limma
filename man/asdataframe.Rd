\name{as.data.frame}
\alias{as.data.frame.MArrayLM}
\title{Turn a Microarray Linear Model Object into a Dataframe}
\description{
Turn a \code{MArrayLM} object into a \code{data.frame}.
}
\usage{
\method{as.data.frame}{MArrayLM}(x, row.names = NULL, optional = FALSE)
}
\arguments{
  \item{x}{an object of class \code{MArrayLM}}
  \item{row.names}{\code{NULL} or a character vector giving the row names for the
          data frame.  Missing values are not allowed.}
  \item{optional}{logical. If \code{TRUE}, setting row names and converting column
          names (to syntactic names) is optional.}
}
\details{
This method combines all the components of \code{x} which have a row for each probe on the array into a \code{data.frame}.
}
\value{
A data.frame.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{as.data.frame}} in the base package.
  
  \link{2.Classes} gives an overview of data classes used in LIMMA.
}
\keyword{array}