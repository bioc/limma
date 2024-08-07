\name{squeezeVar}
\alias{squeezeVar}

\title{Squeeze Sample Variances}
\description{
Squeeze a set of sample variances together by computing empirical Bayes posterior means.
}

\usage{
squeezeVar(var, df, covariate = NULL, robust = FALSE, winsor.tail.p = c(0.05,0.1),
           legacy = NULL)
}

\arguments{
\item{var}{numeric vector of independent sample variances.}
\item{df}{numeric vector of degrees of freedom for the sample variances. Can be a unit vector or of same length as \code{var}.}
\item{covariate}{numeric covariate of same length as \code{var} for estimating a trended prior variance. If \code{NULL}, then the prior variance \code{var.prior} is constant.}
\item{robust}{logical, should the estimation of \code{df.prior} and \code{var.prior} be robustified against outlier sample variances?}
\item{winsor.tail.p}{numeric vector of length 1 or 2, giving left and right tail proportions of \code{x} to Winsorize when \code{robust=TRUE}.}
\item{legacy}{logical. If \code{FALSE} then the new function \code{fitFDistUnequalDF1} will be called internally, otherwise the legacy functions \code{fitFDist} or \code{fitFDistRobustly} will be used. If \code{NULL}, then \code{fitFDistUnequalDF1} will be used whenever the degrees of freedom \code{df} are not all equal.}
}

\details{
This function implements empirical Bayes algorithms proposed by Smyth (2004) and Phipson et al (2016).

A conjugate Bayesian hierarchical model is assumed for a set of sample variances.
The hyperparameters are estimated by fitting a scaled F-distribution to the sample variances.
The function returns the posterior variances and the estimated hyperparameters.

Specifically, the sample variances \code{var} are assumed to follow scaled chi-squared distributions, conditional on the true variances,
and an scaled inverse chi-squared prior is assumed for the true variances.
The scale and degrees of freedom of this prior distribution are estimated from the values of \code{var}.

The effect of this function is to squeeze the variances towards a common value, or to a global trend if a \code{covariate} is provided.
The squeezed variances have a smaller expected mean square error to the true variances than do the sample variances themselves.

The amount of squeezing is controlled by the \code{prior.df}.
Both the global trend and the prior df are estimated internally but fitting an F-distribution to the sample variances, using either \code{fitFDist()} or \code{fitFDistRobustly()} or \code{fitFDistUnequalDF1()}.

If \code{covariate} is non-null, then the scale parameter of the prior distribution is assumed to depend on the covariate.
If the covariate is average log-expression, then the effect is an intensity-dependent trend similar to that in Sartor et al (2006).

\code{robust=TRUE} implements the robust empirical Bayes procedure of Phipson et al (2016), which allows some of the \code{var} values to be outliers.

The \code{legacy} argument was added in limma version 3.61.8 (August 2024).
If \code{legacy=FALSE}, then the new function \code{fitFDistUnequalDF1()} provides improved estimation of the global trend and prior df hyperparameters, especially when the \code{df} values are unequal.
\code{legacy=TRUE} provides legacy behavior for backward compatibility.
}

\note{
This function is called by \code{eBayes}, but beware a possible confusion with the output from that function.
The values \code{var.prior} and \code{var.post} output by \code{squeezeVar} correspond to the quantities \code{s2.prior} and \code{s2.post} output by \code{eBayes}, whereas \code{var.prior} output by \code{eBayes} relates to a different parameter.
}

\value{
A list with components
\item{var.post}{numeric vector of posterior variances. Of same length as \code{var}.}
\item{var.prior}{location or scale of prior distribution. A vector of same length as \code{var} if \code{covariate} is non-\code{NULL}, otherwise a single value.}
\item{df.prior}{degrees of freedom of prior distribution. A vector of same length as \code{var} if \code{robust=TRUE}, otherwise a single value.}
}

\author{Gordon Smyth}

\references{
Phipson B, Lee S, Majewski IJ, Alexander WS, and Smyth GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\doi{10.1214/16-AOAS920}

Sartor MA, Tomlinson CR, Wesselkamper SC, Sivaganesan S, Leikauf GD, Medvedovic M (2006).
Intensity-based hierarchical Bayes method improves testing for differentially expressed genes in microarray experiments.
\emph{BMC bioinformatics} 7, 538.

Smyth GK (2004).
Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology} Volume 3, Issue 1, Article 3.
\doi{10.2202/1544-6115.1027}.
See also the Preprint Version
\url{https://gksmyth.github.io/pubs/ebayes.pdf}
incorporating corrections to 30 June 2009.
}

\seealso{
This function is called by \code{\link{eBayes}}.

This function calls \code{\link{fitFDist}}, \code{\link{fitFDistRobustly}} or \code{\link{fitFDistUnequalDF1}} .

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\examples{
s2 <- rchisq(20,df=5)/5
squeezeVar(s2, df=5)
}

\concept{empirical Bayes}
