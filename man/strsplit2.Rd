\name{strsplit2}
\alias{strsplit2}

\title{Split Composite Names}
\description{
Split a vector of composite names into a matrix of simple names.}

\usage{
strsplit2(x, split, fixed = TRUE, \dots)
}

\arguments{
\item{x}{character vector}
\item{split}{character to split each element of vector on, see \code{strsplit}}
\item{fixed}{logical. If \code{FALSE}, then \code{split} will be treated as a regular expression, if \code{TRUE} then as a literal character string.}
\item{\dots}{other arguments are passed to \code{strsplit}}
}

\value{
A list containing components
\item{Name}{character vector of the same length as \code{x} containing first splits of each element}
\item{Annotation}{character vector of the same length as \code{x} containing second splits of each element}
}

\details{
This function is the same as \code{\link[base]{strsplit}} except that the output value is a matrix instead of a list.
The elements of \code{x} are treated as delimited text with fields separated by \code{split}.
The first column of the output matrix contains the first field from each element of \code{x}, the second column contains the second field etc.
The number of columns is equal to the maximum number of fields for any element of \code{x}.
If different rows have different numbers of fields, then trailing fields are filled out with \code{""}.

The motivation for this function in the limma package is to handle input columns which are composites of two or more annotation fields.
}

\seealso{
\code{\link[base]{strsplit}}.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\examples{
x <- c("AA196000;actinin, alpha 3",
"AA464163;acyl-Coenzyme A dehydrogenase, very long chain",
"3E7;W15277;No Annotation")
strsplit2(x,split=";")
}

\author{Gordon Smyth}

\keyword{character}
