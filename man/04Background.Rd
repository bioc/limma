\name{04.Background}
\alias{04.Background}
\title{Background Correction}

\description{
This page deals with background correction methods for two-color microarray data.

Usually one doesn't need to explicitly ask for background correction of the intensities because this is done by default by \code{\link{normalizeWithinArrays}},
which subtracts the background from the foreground intensities before applying the normalization method.
This default background correction method can be over-ridden by using \code{\link{backgroundCorrect}} which offers a number of alternative
background correct methods to simple subtraction.
The function \code{backgroundCorrect} is used to correct the \code{RGList} before applying \code{normalizeWithinArrays}.

\code{backgroundCorrect} uses utility functions \code{\link{ma3x3.matrix}}, \code{\link{ma3x3.spottedarray}}, \code{\link{normexp.fit}}, \code{\link{normexp.signal}} and \code{\link{normexp.m2loglik}}, \code{\link{normexp.grad}}.

\code{\link{kooperberg}} is a Bayesian background correction tool designed specifically for GenePix data.
\code{kooperberg} is not currently used as the default method for GenePix data because it is computationally intensive.
It requires several columns of the GenePix data files which are not read in by \code{read.maimages}, so you will need to use \code{read.series} instead of \code{read.maimages} if you wish to use \code{kooperberg}.
}

\author{Gordon Smyth}
\keyword{documentation}