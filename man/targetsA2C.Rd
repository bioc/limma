\name{targetsA2C}
\alias{targetsA2C}
\alias{array2channel}
\title{Convert Two-Color Targets Dataframe from One-Row-Per-Array to One-Row-Per-Channel}
\description{
Convert a two-color targets dataframe with one row per array to one with one row per channel.
}
\usage{
targetsA2C(targets, channel.codes = c(1,2), channel.columns = list(Target=c("Cy3","Cy5")),
           grep = FALSE)
}
\arguments{
  \item{targets}{data.frame with one row per array giving information about target samples associated covariates.}
  \item{channel.codes}{numeric or character vector of length 2 giving codes for the channels}
  \item{channel.columns}{named list of character vectors of length 2.
  Each entry gives a pair of names of columns in \code{targets} which contain channel-specific information.
  This pair of columns should be assembled into one column in the output.}
  \item{grep}{logical, if \code{TRUE} then the channel column names are found by \code{grep}ing, i.e., the actual column names need only contain the names given by \code{channel.columns} as substrings}
}
\details{
The \code{targets} dataframe holds information about the RNA samples used as targets in the microarray experiment.
It is often read from a file using \code{\link{readTargets}}.
This function is used to convert the dataframe from an array-orientated format with one row for each array and two columns for the two channels into a channel-orientated format with one row for each individual channel observations.
In statistical terms, the first format treats the arrays as cases and treats the channels as repeated measurements.
The second format treats the individual channel observations as cases.
The second format may be more appropriate if the data is to be analyzed in terms of individual log-intensities.
}
\value{
data.frame with twice as many rows as \code{targets}.
Any pair of columns named by \code{channel.columns} will now be one column.
}
\references{
Smyth, GK, and Altman, NS (2013).
Separate-channel analysis of two-channel microarrays: recovering inter-spot information.
\emph{BMC Bioinformatics} 14, 165.
\doi{10.1186/1471-2105-14-165}
}
\seealso{
\code{targetsA2C} is used by the \code{coerce} method from \code{RGList} to \code{ExpressionSet} in the convert package.

An overview of methods for single channel analysis in limma is given by \link{07.SingleChannel}.
}
\examples{
targets <- data.frame(FileName=c("file1.gpr","file2.gpr"),Cy3=c("WT","KO"),Cy5=c("KO","WT"))
targetsA2C(targets)
}
\author{Gordon Smyth}

\concept{separate channel analysis}
