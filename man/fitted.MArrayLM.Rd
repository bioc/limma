\name{fitted.MArrayLM}
\alias{fitted.MArrayLM}
\concept{regression}
\title{Fitted Values Method for MArrayLM Fits}
\usage{
\S3method{fitted}{MArrayLM}(object, \dots)
}
\arguments{
  \item{object}{a fitted object of class inheriting from \code{"MArrayLM"}.}
  \item{\dots}{other arguments are not currently used.}
}
\description{
  Obtain fitted values from a fitted linear model object.
}
\value{
  A numeric matrix of fitted values.
}
\seealso{
  \code{\link{residuals.MArrayLM}}.
  See also \code{\link{fitted}} in the stats package.
}

\author{Gordon Smyth}
\keyword{models}
\keyword{regression}
