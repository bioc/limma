\name{eBayes}
\alias{eBayes}
\alias{treat}
\title{Empirical Bayes Statistics for Differential Expression}

\description{Given a linear model fit from \code{lmFit}, compute moderated t-statistics, moderated F-statistic, and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a global value.}

\usage{
eBayes(fit, proportion = 0.01, stdev.coef.lim = c(0.1,4),
       trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05,0.1))
treat(fit, fc = 1.2, lfc = NULL, trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05,0.1))
}

\arguments{
  \item{fit}{an \code{MArrayLM} fitted model object produced by \code{lmFit} or \code{contrasts.fit}.
  For \code{ebayes} only, \code{fit} can alternatively be an unclassed list produced by \code{lm.series}, \code{gls.series} or \code{mrlm} containing components \code{coefficients}, \code{stdev.unscaled}, \code{sigma} and \code{df.residual}.}
  \item{proportion}{numeric value between 0 and 1, assumed proportion of genes which are differentially expressed}
  \item{stdev.coef.lim}{numeric vector of length 2, assumed lower and upper limits for the standard deviation of log2-fold-changes for differentially expressed genes}
  \item{trend}{logical, should an intensity-trend be allowed for the prior variance? Default is that the prior variance is constant.}
  \item{robust}{logical, should the estimation of \code{df.prior} and \code{var.prior} be robustified against outlier sample variances?}
  \item{winsor.tail.p}{numeric vector of length 1 or 2, giving left and right tail proportions of \code{x} to Winsorize. Used only when \code{robust=TRUE}.}
  \item{fc}{a minimum fold-change below which changes are not considered scientifically meaningful.}
  \item{lfc}{a minimum log2-fold-change below which changes not considered scientifically meaningful. Defaults to \code{log2(fc)}.}
}

\value{
\code{eBayes} produces an object of class \code{MArrayLM} (see \code{\link{MArrayLM-class}}) containing everything found in \code{fit} plus the following added components:
  \item{t}{numeric matrix of moderated t-statistics.}
  \item{p.value}{numeric matrix of two-sided p-values corresponding to the t-statistics.}
  \item{lods}{numeric matrix giving the log-odds of differential expression (on the natural log scale).}
  \item{s2.prior}{estimated prior value for \code{sigma^2}. A row-wise vector if \code{covariate} is non-\code{NULL}, otherwise a single value.}
  \item{df.prior}{degrees of freedom associated with \code{s2.prior}. A row-wise vector if \code{robust=TRUE}, otherwise a single value.}
  \item{df.total}{row-wise numeric vector giving the total degrees of freedom associated with the t-statistics for each gene. Equal to \code{df.prior+df.residual} or \code{sum(df.residual)}, whichever is smaller.}
  \item{s2.post}{row-wise numeric vector giving the posterior values for \code{sigma^2}.}
  \item{var.prior}{column-wise numeric vector giving estimated prior values for the variance of the log2-fold-changes for differentially expressed gene for each constrast. Used for evaluating \code{lods}.}
  \item{F}{row-wise numeric vector of moderated F-statistics for testing all contrasts defined by the columns of \code{fit} simultaneously equal to zero.}
  \item{F.p.value}{row-wise numeric vector giving p-values corresponding to \code{F}.}
The matrices \code{t}, \code{p.value} and \code{lods} have the same dimensions as the input object \code{fit}, with rows corresponding to genes and columns to coefficients or contrasts.
The vectors \code{s2.prior}, \code{df.prior}, \code{df.total}, \code{F} and \code{F.p.value} correspond to rows, with length equal to the number of genes.
The vector \code{var.prior} corresponds to columns, with length equal to the number of contrasts.
If \code{s2.prior} or \code{df.prior} have length 1, then the same value applies to all genes.

\code{s2.prior}, \code{df.prior} and \code{var.prior} contain empirical Bayes hyperparameters used to obtain \code{df.total}, \code{s2.post} and \code{lods}.

\code{treat} a produces an \code{MArrayLM} object similar to that from \code{eBayes} but without \code{lods}, \code{var.prior}, \code{F} or \code{F.p.value}.
}

\details{
These functions are used to rank genes in order of evidence for differential expression.
They use an empirical Bayes method to squeeze the genewise-wise residual variances towards a common value (or towards a global trend) (Smyth, 2004; Phipson et al, 2016).
The degrees of freedom for the individual variances are increased to reflect the extra information gained from the empirical Bayes moderation, resulting in increased statistical power to detect differential expression.

Theese functions accept as input an \code{MArrayLM} fitted model object \code{fit} produced by \code{lmFit}.
The columns of \code{fit} define a set of contrasts which are to be tested equal to zero.
The fitted model object may have been processed by \code{contrasts.fit} before being passed to \code{eBayes} to convert the coefficients of the original design matrix into an arbitrary number of contrasts.

The empirical Bayes moderated t-statistics test each individual contrast equal to zero.
For each gene (row), the moderated F-statistic tests whether all the contrasts are zero.
The F-statistic is an overall test computed from the set of t-statistics for that probe.
This is exactly analogous the relationship between t-tests and F-statistics in conventional anova, except that the residual mean squares have been moderated between genes.

The estimates \code{s2.prior} and \code{df.prior} are computed by \code{fitFDist}.
\code{s2.post} is the weighted average of \code{s2.prior} and \code{sigma^2} with weights proportional to \code{df.prior} and \code{df.residual} respectively.
The log-odds of differential expression \code{lods} was called the \emph{B-statistic} by Loennstedt and Speed (2002).
The F-statistics \code{F} are computed by \code{classifyTestsF} with \code{fstat.only=TRUE}.

\code{eBayes} does not compute ordinary t-statistics because they always have worse performance than the moderated versions.
The ordinary (unmoderated) t-statistics can, however, can be easily extracted from the linear model output for comparison purposes---see the example code below.

\code{treat} computes empirical Bayes moderated-t p-values relative to a minimum fold-change threshold.
Instead of testing for genes that have true log-fold-changes different from zero, it tests whether the true log2-fold-change is greater than \code{lfc} in absolute value (McCarthy and Smyth, 2009).
In other words, it uses an interval null hypothesis, where the interval is [-lfc,lfc].
When the number of DE genes is large, \code{treat} is often useful for giving preference to larger fold-changes and for prioritizing genes that are biologically important.
\code{treat} is concerned with p-values rather than posterior odds, so it does not compute the B-statistic \code{lods}.
The idea of thresholding doesn't apply to F-statistics in a straightforward way, so moderated F-statistics are also not computed.
When \code{fc=1} and \code{lfc=0}, \code{treat} is identical to \code{eBayes}, except that F-statistics and B-statistics are not computed.
The \code{fc} threshold is usually chosen relatively small, because genes need to have fold changes substantially greater than the testing threshold in order to be considered statistically significant.
Typical values for \code{fc} are \code{1.1}, \code{1.2} or \code{1.5}.
The top genes chosen by \code{treat} can be examined using \code{\link{topTreat}}.

The \code{treat} threshold can be specified either as a fold-change via \code{fc} or as a log2-fold-change via \code{lfc}, with \code{lfc = log2(fc)}.
Note that the \code{treat} testing procedure is considerably more rigorous and conservative than simply applying same \code{fc} values as a fold-change cutoff to the list of differentially expressed genes.
Indeed, the observed log2-fold-change needs to substantially larger than \code{lfc} for a gene to be called as statistically significant by \code{treat}.
The threshold should be chosen as a small value below which results should be ignored rather than as a target fold-change.
In practice, modest values for \code{fc} such as \code{1.1}, \code{1.2} or \code{1.5} are usually the most useful.
Setting \code{fc=1.2} or \code{fc=1.5} will usually cause most differentially expressed genes to have estimated fold-changes of 2-fold or greater, depending on the sample size and precision of the experiment.
Larger thresholds are usually overly conservative and counter productive.
In general, the \code{fc} threshold should be chosen sufficiently small so that a worthwhile number of DE genes remain, otherwise the purpose of prioritizing genes with larger fold-changes will be defeated.

The use of \code{eBayes} or \code{treat} with \code{trend=TRUE} is known as the \emph{limma-trend} method (Law et al, 2014; Phipson et al, 2016).
With this option, an intensity-dependent trend is fitted to the prior variances \code{s2.prior}.
Specifically, \code{\link{squeezeVar}} is called with the \code{covariate} equal to \code{Amean}, the average log2-intensity for each gene.
The trend that is fitted can be examined by \code{\link{plotSA}}.
limma-trend is useful for processing expression values that show a mean-variance relationship.
This is often useful for microarray data, and it can also be applied to RNA-seq counts that have been converted to log2-counts per million (logCPM) values (Law et al, 2014).
When applied to RNA-seq logCPM values, limma-trend give similar results to the \code{\link{voom}} method.
The voom method incorporates the mean-variance trend into the precision weights, whereas limma-trend incorporates the trend into the empirical Bayes moderation.
limma-trend is somewhat simpler than \code{voom} because it assumes that the sequencing depths (library sizes) are not wildly different between the samples and it applies the mean-variance trend on a genewise basis instead to individual observations.
limma-trend is recommended for RNA-seq analysis when the library sizes are reasonably consistent (less than 3-fold difference from smallest to largest) because of its simplicity and speed.

If \code{robust=TRUE} then the robust empirical Bayes procedure of Phipson et al (2016) is used.
This is frequently useful to protect the empirical Bayes procedure against hyper-variable or hypo-variable genes, especially when analysing RNA-seq data.
See \code{\link{squeezeVar}} for more details.
}

\note{
The algorithm used by \code{eBayes} and \code{treat} with \code{robust=TRUE} was revised slightly in limma 3.27.6.
The minimum \code{df.prior} returned may be slightly smaller than previously.
}

\seealso{
\code{\link{squeezeVar}}, \code{\link{fitFDist}}, \code{\link{tmixture.matrix}}, \code{\link{plotSA}}.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\author{Gordon Smyth and Davis McCarthy}

\references{
Law, CW, Chen, Y, Shi, W, Smyth, GK (2014).
Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
\emph{Genome Biology} 15, R29.
\url{http://genomebiology.com/2014/15/2/R29}

Loennstedt, I., and Speed, T. P. (2002). Replicated microarray data. \emph{Statistica Sinica} \bold{12}, 31-46.

McCarthy, D. J., and Smyth, G. K. (2009).
Testing significance relative to a fold-change threshold is a TREAT.
\emph{Bioinformatics} 25, 765-771.
\url{http://bioinformatics.oxfordjournals.org/content/25/6/765}

Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016).
Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression.
\emph{Annals of Applied Statistics} 10, 946-963.
\url{http://projecteuclid.org/euclid.aoas/1469199900}

Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology} 3, Article 3.
\url{http://www.statsci.org/smyth/pubs/ebayes.pdf}
}

\examples{
#  See also lmFit examples

#  Simulate gene expression data,
#  6 microarrays and 100 genes with one gene differentially expressed
set.seed(2016)
sigma2 <- 0.05 / rchisq(100, df=10) * 10
y <- matrix(rnorm(100*6,sd=sqrt(sigma2)),100,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
y[1,4:6] <- y[1,4:6] + 1
fit <- lmFit(y,design)

#  Moderated t-statistic
fit <- eBayes(fit)
topTable(fit,coef=2)

#  Ordinary t-statistic
ordinary.t <- fit$coef[,2] / fit$stdev.unscaled[,2] / fit$sigma

# Treat relative to a 10% fold-change
tfit <- treat(fit, fc=1.1)
topTreat(tfit,coef=2)
}

\keyword{empirical Bayes}
