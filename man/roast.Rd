\name{roast}
\alias{roast}
\alias{roast.default}
\alias{mroast}
\alias{mroast.default}
\alias{Roast-class}
\alias{show,Roast-method}
\alias{fry}
\alias{fry.default}
\title{Rotation Gene Set Tests}
\description{
Rotation gene set testing for linear models.
}

\usage{
\method{roast}{default}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
      set.statistic = "mean", gene.weights = NULL, var.prior = NULL, df.prior = NULL,
      nrot = 1999, approx.zscore = TRUE, legacy = FALSE, \dots)
\method{mroast}{default}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
       set.statistic = "mean", gene.weights = NULL, var.prior = NULL, df.prior = NULL,
       nrot = 1999, approx.zscore = TRUE, legacy = FALSE, adjust.method = "BH",
       midp = TRUE, sort = "directional", \dots)
\method{fry}{default}(y, index = NULL, design = NULL, contrast = ncol(design), geneid = NULL,
      gene.weights = NULL, standardize = "posterior.sd", sort = "directional", \dots)
}

\arguments{
  \item{y}{numeric matrix giving log-expression or log-ratio values for a series of microarrays, or any object that can coerced to a matrix including \code{ExpressionSet}, \code{MAList}, \code{EList} or \code{PLMSet} objects.
        Rows correspond to probes and columns to samples.
        \code{NA} or infinite values are not allowed.
        If either \code{var.prior} or \code{df.prior} are \code{NULL}, then \code{y} should contain values for all genes on the arrays. If both prior parameters are given, then only \code{y} values for the test set are required.}
  \item{index}{index vector specifying which rows (probes) of \code{y} are in the test set.
        Can be a vector of integer indices, or a logical vector of length \code{nrow(y)}, or a vector of gene IDs corresponding to entries in \code{geneid}.
        Alternatively it can be a data.frame with the first column containing the index vector and the second column containing directional gene contribution weights.
        For \code{mroast} or \code{fry}, \code{index} is a list of index vectors or a list of data.frames. }
  \item{design}{design matrix}
  \item{contrast}{contrast for which the test is required.
        Can be an integer specifying a column of \code{design}, or the name of a column of \code{design}, or a numeric contrast vector of length equal to the number of columns of \code{design}.}
  \item{geneid}{gene identifiers corresponding to the rows of \code{y}.
        Can be either a vector of length \code{nrow(y)} or the name of the column of \code{y$genes} containing the gene identifiers.
        Defaults to \code{rownames(y)}.}
  \item{set.statistic}{summary set statistic. Possibilities are \code{"mean"},\code{"floormean"},\code{"mean50"} or \code{"msq"}.}
  \item{gene.weights}{numeric vector of directional (positive or negative) contribution weights specifying the size and direction of the contribution of each probe to the gene set statistics.
        For \code{mroast} or \code{fry}, this vector must have length equal to \code{nrow(y)}.
        For \code{roast}, can be of length \code{nrow(y)} or of length equal to the number of genes in the test set.} 
  \item{var.prior}{prior value for residual variances. If not provided, this is estimated from all the data using \code{squeezeVar}.}
  \item{df.prior}{prior degrees of freedom for residual variances. If not provided, this is estimated using \code{squeezeVar}.}
  \item{nrot}{number of rotations used to compute the p-values. Low values like 999 are suitable for testing but higher values such as 9999 or more are recommended for publication purposes.}
  \item{approx.zscore}{logical, if \code{TRUE} then a fast approximation is used to convert t-statistics into z-scores prior to computing set statistics. If \code{FALSE}, z-scores will be exact.}
  \item{legacy}{logical. See Note below for usage.}
  \item{adjust.method}{method used to adjust the p-values for multiple testing. See \code{\link{p.adjust}} for possible values.}
  \item{midp}{logical, should mid-p-values be used in instead of ordinary p-values when adjusting for multiple testing?}
  \item{sort}{character, whether to sort output table by directional p-value (\code{"directional"}), non-directional p-value (\code{"mixed"}), or not at all (\code{"none"}).}
  \item{standardize}{how to standardize for unequal probewise variances. Possibilities are \code{"residual.sd"}, \code{"posterior.sd"} or \code{"none"}.}
  \item{\dots}{any argument that would be suitable for \code{\link{lmFit}} or \code{\link{eBayes}} can be included.}
}

\value{
\code{roast} produces an object of class \code{"Roast"}.
This consists of a list with the following components:
  \item{p.value}{data.frame with columns \code{Active.Prop} and \code{P.Value}, giving the proportion of genes in the set contributing materially to significance and estimated p-values, respectively.
Rows correspond to the alternative hypotheses Down, Up, UpOrDown (two-sided) and Mixed.}
  \item{var.prior}{prior value for residual variances.}
  \item{df.prior}{prior degrees of freedom for residual variances.}

\code{mroast} produces a data.frame with a row for each set and the following columns:
  \item{NGenes}{number of genes in set}
	\item{PropDown}{proportion of genes in set with \code{z < -sqrt(2)}}
	\item{PropUp}{proportion of genes in set with \code{z > sqrt(2)}}
	\item{Direction}{direction of change, \code{"Up"} or \code{"Down"}}
	\item{PValue}{two-sided directional p-value}
	\item{FDR}{two-sided directional false discovery rate}
	\item{PValue.Mixed}{non-directional p-value}
	\item{FDR.Mixed}{non-directional false discovery rate}

\code{fry} produces the same output format as \code{mroast} but without the columns \code{PropDown} and \code{ProbUp}.
}

\details{
These functions implement rotation gene set tests proposed by Wu et al (2010).
They perform \emph{self-contained} gene set tests in the sense defined by Goeman and Buhlmann (2007).
For \emph{competitive} gene set tests, see \code{\link{camera}}.
For a gene set enrichment analysis (GSEA) style analysis using a database of gene sets, see \code{\link{romer}}.

\code{roast} and \code{mroast} test whether any of the genes in the set are differentially expressed.
They can be used for any microarray experiment that can be represented by a linear model.
The design matrix for the experiment is specified as for the \code{\link{lmFit}} function, and the contrast of interest is specified as for the \code{\link{contrasts.fit}} function.
This allows users to focus on differential expression for any coefficient or contrast in a linear model.
If \code{contrast} is not specified, then the last coefficient in the linear model will be tested.

The argument \code{index} is often made using \link{ids2indices} but does not have to be.
Each set to be tested is represented by a vector of row numbers or a vector of gene IDs.
Gene IDs should correspond to either the rownames of \code{y} or the entries of \code{geneid}.

All three functions support directional contribution gene weights, which can be entered either through the \code{gene.weights} argument or via \code{index}.
Directional gene weights allow each gene to be flagged as to its direction and magnitude of change based on prior experimentation.
A typical use is to make the \code{gene.weights} \code{1} or \code{-1} depending on whether the gene is up or down-regulated in the pathway under consideration.
Probes with directional weights of opposite signs are expected to have expression changes in opposite directions.
Gene with larger gene weights in absolute size will have more weight in the set statistic calculation.

Gene weights can be either genome-wide or set-specific.
Genome-wide weights can be entered via the \code{gene.weights} argument.
Set specific weights can be input by including the gene weights as part of the set's entry in \code{index}.
If any of the components of \code{index} are data.frames, then the second column will be assumed to be gene contribution weights for that set.
All three functions (\code{roast}, \code{mroast} and \code{fry}) support set-specific gene contribution weights as part of an \code{index} data.frame.

Set-specific directional gene weights are used to represent \emph{expression signatures} assembled from previous experiments, from gene annotation or from prior hypotheses.
In the output from \code{roast}, \code{mroast} or \code{fry}, a significant \code{"Up"} p-value means that the differential expression results found in \code{y} are positively correlated with the expression signature coded by the gene weights.
Conversely, a significant \code{"Down"} p-value means that the differential expression log-fold-changes are negatively correlated with the expression signature.

Note that the contribution weights set by \code{gene.weights} are different in nature and purpose to the precision weights set by the \code{weights} argument of \code{lmFit}.
\code{gene.weights} control the contribution of each gene to the formation of the gene set statistics and are directional, i.e., can be positive or negative.
\code{weights} indicate the precision of the expression measurements and should be positive.
The \code{weights} are used to construct genewise test statistics whereas \code{gene.weights} are used to combine the genewise test statistics.

The arguments \code{df.prior} and \code{var.prior} have the same meaning as in the output of the \code{\link{eBayes}} function.
If these arguments are not supplied, then they are estimated exactly as is done by \code{eBayes}.

The gene set statistics \code{"mean"}, \code{"floormean"}, \code{"mean50"} and \code{msq} are defined by Wu et al (2010).
The different gene set statistics have different sensitivities when only some of the genes in a set are differentially expressed.
If \code{set.statistic="mean"} then the set will be statistically significantly only when the majority of the genes are differentially expressed.
\code{"floormean"} and \code{"mean50"} will detect as few as 25\% differentially expressed in a set.
\code{"msq"} is sensitive to even smaller proportions of differentially expressed genes, if the effects are reasonably large.
Overall, the \code{"msq"} statistic gives the best power for rejecting the null hypothesis of no differentially expressed genes, but the significance can be driven by a small number of genes.
In many genomic applications it is appropriate to limit results to gene sets for which most of the genes response in a concordance direction, so the relatively conservative \code{"mean"} statistic is the default choice.

The output gives p-values three possible alternative hypotheses, 
\code{"Up"} to test whether the genes in the set tend to be up-regulated, with positive t-statistics,
\code{"Down"} to test whether the genes in the set tend to be down-regulated, with negative t-statistics,
and \code{"Mixed"} to test whether the genes in the set tend to be differentially expressed, without regard for direction.

\code{roast} estimates p-values by simulation, specifically by random rotations of the orthogonalized residuals (Langsrud, 2005), so p-values will vary slightly from run to run.
The p-value is computed as \code{(b+1)/(nrot+1)} where \code{b} is the number of rotations giving a more extreme statistic than that observed (Phipson and Smyth, 2010).
This means that the smallest possible mixed or two-sided p-values are \code{1/(nrot+1)}.
The function uses a symmetry argument to double the effective number of rotations for the one-sided tests, so the smallest possible \code{"Up"} or \code{"Down"} p-value is \code{1/(2*nrot+1)}.

The number of rotations \code{nrot} can (and should) be increased tTo get more precise p-values from \code{roast} or \code{mroast},
The default \code{nrot} is set fairly low to facilitate quick testing and experimentation but the smallest possible two-sided p-value is \code{1/(nrot+1)}.
To get definitive p-values for publication, at least \code{nrot=9999} or higher is recommended.

\code{mroast} does roast tests for multiple sets, including adjustment for multiple testing.
By default, \code{mroast} reports ordinary p-values but uses mid-p-values (Routledge, 1994) at the multiple testing stage.
Mid-p-values are probably a good choice when using false discovery rates (\code{adjust.method="BH"}) but not when controlling the family-wise type I error rate (\code{adjust.method="holm"}).

To improve the performance of the gene set statistics, \code{roast} and \code{mroast} transform the genewise moderated t-statistics to normality using \code{zscoreT}.
By default, an approximate closed-form transformation is used (\code{approx.zscore=TRUE}), which is very much faster than the exact transformation and performs just as well.
In Bioconductor 2.10, the transformation used has been changed from Hill's (1970) approximation to Bailey's (1980) formula because the latter is faster and gives more even accuracy; see \code{\link{zscoreT}} for more details.

\code{fry} is a fast alternative designed to approximate what \code{mroast} with \code{set.stat="mean"} would give for a very large (infinite) number of rotations.
In the special case that \code{df.prior} is very large and \code{set.statistic="mean"}, \code{fry} gives the same directional p-values that \code{mroast} would give if an infinite number of rotations could be performed.
In other circumstances, when genes have different variances, \code{fry} uses a standardization strategy to approximate the \code{mroast} results.
Using \code{fry} is recommended when performing tests for a large number of sets because it is fast and because it returns higher resolution p-values that are not limited by the number of rotations performed.
Note, the close approximation of \code{fry} to \code{mroast} is only for the directional p-values.
The \code{fry} mixed p-values are computed by a different method and will not necessarily be very close to those from \code{mroast}.
}

\note{
For Bioconductor 3.10, \code{roast} and \code{mroast} have been revised to use much less memory by conducting the rotations in chunks and to be about twice as fast by updating the normalizing transformation used when \code{approx.zscore=TRUE}.
For a limited time, users wishing to reproduce Bioconductor 3.9 results exactly can set \code{legacy=TRUE} to turn these revisions off.

\code{approx.score=TRUE} become the default in Bioconductor 3.0 (October 2014).

The default set statistic was changed from \code{"msq"} to \code{"mean"} in Bioconductor 2.7 (October 2010).
}

\seealso{
See \link{10.GeneSetTests} for a description of other functions used for gene set testing.
}
\author{Gordon Smyth and Di Wu}

\references{
Goeman, JJ, and Buhlmann, P (2007).
Analyzing gene expression data in terms of gene sets: methodological issues.
\emph{Bioinformatics} 23, 980-987. 

Langsrud, O (2005).
Rotation tests.
\emph{Statistics and Computing} 15, 53-60.

Phipson B, and Smyth GK (2010).
Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
\emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Article 39.
\url{http://www.statsci.org/smyth/pubs/PermPValuesPreprint.pdf}

Routledge, RD (1994).
Practicing safe statistics with the mid-p.
\emph{Canadian Journal of Statistics} 22, 103-110.

Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010). ROAST: rotation gene set tests for complex microarray experiments.
\emph{Bioinformatics} 26, 2176-2182.
\url{http://bioinformatics.oxfordjournals.org/content/26/17/2176}
}

\examples{
y <- matrix(rnorm(100*4,sd=0.3),100,4)
design <- cbind(Intercept=1,Group=c(0,0,1,1))

# First set of 5 genes are all up-regulated
index1 <- 1:5
y[index1,3:4] <- y[index1,3:4]+3
roast(y,index1,design,contrast=2)

# Second set of 5 genes contains none that are DE
index2 <- 6:10
mroast(y,list(set1=index1,set2=index2),design,contrast=2)
fry(y,list(set1=index1,set2=index2),design,contrast=2)

# Third set of 6 genes contains three down-regulated genes and three up-regulated genes
index3 <- 11:16
y[index3[1:3],3:4] <- y[index3[1:3],3:4]-3
y[index3[4:6],3:4] <- y[index3[4:6],3:4]+3

# Without gene weights
# Mixed p-value is significant for set3 but not the directional p-values
mroast(y,list(set1=index1,set2=index2,set3=index3),design,contrast=2)
fry(y,list(set1=index1,set2=index2,set3=index3),design,contrast=2)

# With gene weights
# Set3 is significantly up (i.e., positively correlated with the weights)
index3 <- data.frame(Gene=11:16,Weight=c(-1,-1,-1,1,1,1))
mroast(y,list(set1=index1,set2=index2,set3=index3),design,contrast=2)
fry(y,list(set1=index1,set2=index2,set3=index3),design,contrast=2)
}

\concept{gene set tests}
