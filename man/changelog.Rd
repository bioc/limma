\name{changeLog}
\alias{changeLog}

\title{Change Log}
\description{Show the most recent changes from a package change log or NEWS file.}

\usage{
changeLog(n = 30, package = "limma")
}

\arguments{
  \item{n}{integer, number of lines to write of change log.}
  \item{package}{character string giving name of package.}
}

\details{
The function will look for a file \code{changelog.txt} or \code{ChangeLog} in the top-level or \code{doc} directories of the installed package.
Failing that, it will look for \code{NEWS} or \code{NEWS.md} in the top-level directory.

Note that \code{changeLog} does not write the content of \code{NEWS.Rd}, which is a structured file.
Use \code{news(package="limma")} for that instead.
}

\value{No value is produced, but a number of lines of text are written to standard output.}

\author{Gordon Smyth}

\examples{
changeLog()
changeLog(package="statmod")
}

\seealso{
\link{01.Introduction}, \code{\link{news}}.
}

\keyword{documentation}
