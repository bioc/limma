\name{sampleInfoFromGEO}
\alias{sampleInfoFromGEO}

\title{Read Sample Characteristics From GEO Series Matrix File}
\description{
Read sample annotation from a GEO Series Matrix File into data.frames.}

\usage{
sampleInfoFromGEO(file, remove.constant.columns = TRUE)
}

\arguments{
\item{file}{file name or path of GEO series matrix file.}
\item{remove.constant.columns}{logical, if \code{TRUE} then columns that have the same entry for every sample are removed.}
}

\value{
A list with three components:
\item{SampleInfo}{character matrix of sample annotation.}
\item{CharacteristicsCh1}{character matrix of sample characteristics associated with expression channel 1.}
\item{CharacteristicsCh2}{character matrix of sample characteristics associated with expression channel 2.}
Each sample corresponds to one row.
}

\details{
This function parses a GEO series matrix file.
Sample characteristics associated with expression channels 1 and 2 are separated into separate character matrices.
The function particularly allows for the fact that not every sample characteristic will have an entry for every sample.
}

\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\examples{
# This example downloads a series matrix file of about 33MB
\dontrun{
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40115/matrix/GSE40115-GPL15931_series_matrix.txt.gz"
download.file(url, "GSE40115.txt.gz")
a <- sampleInfoFromGEO("GSE40115.txt.gz")
colnames(a$SampleInfo)
colnames(a$CharacteristicsCh1)
colnames(a$CharacteristicsCh2)
}
}

\author{Gordon Smyth}

\keyword{reading data}
