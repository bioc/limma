\name{06.LinearModels}
\alias{06.LinearModels}
\title{Topic: Linear Models for Microarrays}

\description{
This page gives an overview of the LIMMA functions available to fit linear models and to interpret the results.
This page covers models for two color arrays in terms of log-ratios or for single-channel arrays in terms of log-intensities.
If you wish to fit models to the individual channel log-intensities from two colour arrays, see \link{07.SingleChannel}.

The core of this package is the fitting of gene-wise linear models to microarray data.
The basic idea is to estimate log-ratios between two or more target RNA samples simultaneously.
See the LIMMA User's Guide for several case studies.
}

\section{Fitting Models}{

The main function for model fitting is \code{\link{lmFit}}.
This is recommended interface for most users.
\code{lmFit} produces a fitted model object of class \code{\link[limma:marraylm]{MArrayLM}} containing coefficients, standard errors and residual standard errors for each gene.
\code{lmFit} calls one of the following three functions to do the actual computations:

\describe{
\item{ \code{\link{lm.series}} }{
	Straightforward least squares fitting of a linear model for each gene.}

\item{ \code{\link{mrlm}} }{
	An alternative to \code{lm.series} using robust regression as implemented by the \code{rlm} function in the MASS package.}

\item{ \code{\link{gls.series}} }{
	Generalized least squares taking into account correlations between duplicate spots (i.e., replicate spots on the same array) or related arrays.
	The function \code{\link{duplicateCorrelation}} is used to estimate the inter-duplicate or inter-block correlation before using \code{gls.series}.}
}

All the functions which fit linear models use \code{link{getEAW}} to extract data from microarray data objects, and \code{\link{unwrapdups}} which provides an unified method for handling duplicate spots.
}

\section{Forming the Design Matrix}{

\code{lmFit} has two main arguments, the expression data and the design matrix.
The design matrix is essentially an indicator matrix which specifies which target RNA samples were applied to each channel on each array.
There is considerable freedom in choosing the design matrix - there is always more than one choice which is correct provided it is interpreted correctly.

Design matrices for Affymetrix or single-color arrays can be created using the function \code{\link[stats]{model.matrix}} which is part of the R base package.
The function \code{\link{modelMatrix}} is provided to assist with creation of an appropriate design matrix for two-color microarray experiments.
For direct two-color designs, without a common reference, the design matrix often needs to be created by hand.
}

\section{Making Comparisons of Interest}{

Once a linear model has been fit using an appropriate design matrix, the command \code{\link{makeContrasts}} may be used to form a contrast matrix to make comparisons of interest.
The fit and the contrast matrix are used by \code{\link{contrasts.fit}} to compute fold changes and t-statistics for the contrasts of interest.
This is a way to compute all possible pairwise comparisons between treatments for example in an experiment which compares many treatments to a common reference.
}

\section{Assessing Differential Expression}{

After fitting a linear model, the standard errors are moderated using a simple empirical Bayes model using \code{\link{eBayes}} or \code{\link{treat}}.
A moderated t-statistic and a log-odds of differential expression is computed for each contrast for each gene.
\code{treat} tests whether log-fold-changes are greater than a threshold rather than merely different to zero.

\code{\link{eBayes}} and \code{\link{treat}} use internal functions \code{\link{squeezeVar}}, \code{\link{fitFDist}}, \code{\link{fitFDistRobustly}}, \code{\link{fitFDistUnequalDF1}}, \code{\link{tmixture.matrix}} and \code{\link{tmixture.vector}}.
}

\section{Summarizing Model Fits}{

After the above steps the results may be displayed or further processed using:
\describe{
\item{ \code{\link{topTable}} }{
	Presents a list of the genes most likely to be differentially expressed for a given contrast or set of contrasts.}

\item{ \code{\link{topTableF}} }{
	Presents a list of the genes most likely to be differentially expressed for a given set of contrasts.
	Equivalent to \code{topTable} with \code{coef} set to all the coefficients, \code{coef=1:ncol(fit)}.
}

\item{ \code{\link{volcanoplot}}}{
	Volcano plot of fold change versus the B-statistic for any fitted coefficient.}

\item{ \code{\link{plotlines}}}{
	Plots fitted coefficients or log-intensity values for time-course data.}

\item{ \code{\link{genas}}}{
	Estimates and plots biological correlation between two coefficients.}

\item{ \code{\link{write.fit}} }{
	Writes an \code{MarrayLM} object to a file.
	Note that if \code{fit} is an \code{MArrayLM} object, either \code{write.fit} or \code{write.table} can be used to write the results to a delimited text file.}
}
For multiple testing functions which operate on linear model fits, see \link{08.Tests}.
}

\section{Model Selection}{
\code{\link{selectModel}} provides a means to choose between alternative linear models using AIC or BIC information criteria.
}

\author{Gordon Smyth}
\references{
Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\url{http://projecteuclid.org/euclid.aoas/1469199900}

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology}, \bold{3}, No. 1, Article 3.
\url{https://gksmyth.github.io/pubs/ebayes.pdf}

Smyth, G. K., Michaud, J., and Scott, H. (2005). The use of within-array replicate spots for assessing differential expression in microarray experiments. \emph{Bioinformatics} 21(9), 2067-2075.
}

\seealso{
\link{01.Introduction},
\link{02.Classes},
\link{03.ReadingData},
\link{04.Background},
\link{05.Normalization},
\link{06.LinearModels},
\link{07.SingleChannel},
\link{08.Tests},
\link{09.Diagnostics},
\link{10.GeneSetTests},
\link{11.RNAseq}
}

\keyword{documentation}
