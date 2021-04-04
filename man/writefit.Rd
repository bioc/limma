\name{write.fit}
\alias{write.fit}

\title{Write MArrayLM Object to a File}

\description{
Write a microarray linear model fit to a file.
}

\usage{
write.fit(fit, results = NULL, file, digits = NULL,
          adjust = "none", method = "separate", F.adjust = "none",
          quote = FALSE, sep = "\t", row.names = TRUE, \dots)
}

\arguments{
\item{fit}{object of class \code{MArrayLM} containing the results of a linear model fit.}
\item{results}{object of class \code{TestResults}.}
\item{file}{character string giving path name for the output file.}
\item{digits}{integer indicating rounding precision for output values. If \code{NULL}, then no rounding is done.}
\item{adjust}{character string specifying multiple-testing adjustment method for the t-statistic P-values, e.g., \code{"BH"}. See \code{\link[stats]{p.adjust}} for the available options. If \code{NULL} or \code{"none"} then the P-values are not adjusted.}
\item{method}{character string, should the P-value adjustment be \code{"global"} or \code{"separate"} for each contrast. Ignored if \code{adjust = "none"}.}
\item{F.adjust}{character string specifying adjustment method for the F-statistic P-values.}
\item{quote}{logical value. If \code{TRUE}, any character or factor columns will be surrounded by double quotes.}
\item{sep}{the field separator string. Values in the output file will be separated by this string.}
\item{row.names}{logical value, whether to include row names in the output file.}
\item{\dots}{other arguments are passed to \code{write.table}.}
}

\value{
No value is produced but a file is written to the current working directory.
}

\details{
This function writes a delimited text file containing for each gene (1) the average log2-intensity (\code{AveExpr}), (2) the coefficients or contrasts (log2-fold-changes, \code{Coef}), (3) moderated t-statistics, (4) t-statistic P-values, (5) F-statistic if available, (6) F-statistic P-values if available, (7) decideTests results if available and (8) gene names and annotation.

The \code{results} argument is optional. If supplied, it should be the output from \code{decideTests} for the same fit object, which indicates whether each contrast for each gene is considered statistically significant or not (coded 1 or -1 for positive or negative significant differences and 0 for non-significant values).

If \code{fit} contains row names and \code{row.names=TRUE}, then the row names will be the first column of the output file  with a blank column heading.
This behaviour is analogous to that of \code{write.csv} or to \code{write.table} with \code{col.names=NA}.
}

\author{Gordon Smyth}
\seealso{
\code{\link{write.table}} or \code{\link{write.csv}} in the base library.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\examples{
\dontrun{
# The following three alternatives are equivalent:
write.fit(fit, file = "temp.csv", sep = ",")
write.csv(fit, file = "temp.csv")
a <- as.data.frame(fit)
write.csv(fit, file = "temp.csv")
}
}

\keyword{IO}
