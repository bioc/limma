\name{limmanews}
\title{limma News}
\encoding{UTF-8}


\section{version 4.0.0 (2026-10-30)}{\itemize{
  \item
  voomLmFit() transferred to limma from the edgeR package.
  Previously, voomLmFit() needed to call edgeR::glmFit() in order to
  identify exact zero fitted values, but voomLmFit() now calls a
  new limma C routine 'poisfit' for the same purpose, removing the
  dependence on edgeR.

  \item
  For the purpose of estimating sample weights, voomLmFit() calls
  arrayWeights() with 'method="reml"' if there is no loss of df due to
  exact zero fitted values and 'method="genebygene"' otherwise.

  \item
  The voomLmFit help page has been expanded, while the voom and
  voomWithQualityWeights help pages now have preambles that point
  users to voomLmFit.

  \item
  voomLmFit() now gives priority to 'offset' over 'offset.prior' or
  'lib.size' and no longer forces the rows of 'offset.prior' to add
  to zero.

  \item
  New argument 'contrasts' for voomLmFit(), lmFit(), lm.series() and
  gls.series(), which converts coefficient estimates and standard
  deviations into contrast estimates and standard deviations. This
  argument provides the  same functionality as calling 'contrasts.fit'
  on the output object, but with exact rather than approximate
  standard deviations.

  \item
  New argument 'nthreads' for voomLmFit(), lmFit(), lm.series(),
  gls.series(), duplicateCorrelation(), and
  arrayWeights() to enable OpenMP parallelization of genewise fits.

  \item
  C code backends for lm.series() and gls.series().

  \item
  C code backend for duplicateCorrelation(). duplicateCorrelation()
  previously called statmod::mixedModel2() to estimate genewise
  correlations, but this functionality is now absorbed into the C
  code. The new code is able to evaluate the genewise SVDs needed by
  mixedModel2() more efficiently than was possible for the univariate
  R function, resulting in substantial speed and memory improvements.

  \item
  C code backend for .arrayWeightsPrWtsREML(), which is called by
  arrayWeights() when 'method="reml"' and prior weights
  are set.

  \item
  Minor changes to weightedLowess() C code
  to fit in with the new C code elsewhere.

  \item
  New function normalizeBetweenSamples(). Similar to
  normalizeBetweenArrays() but works on ELists.

  \item
  camera() and cameraPR() now support gene sets as expression
  signatures, whereby a logFC or similar expected directional change
  is provided for each gene in the gene set.

  \item
  avearrays() methods for EList and MAList objects now sum the weights
  instead of averaging them. Other matrices are still averaged.

  \item
  Fix filtering error in Agilent rat case study in User's Guide.

  \item
  All the Rd pages have been checked for typos and inconsistent
  formating, especially regarding the reference lists. Reference URLs
  converted to DOIs where possible and some new DOIs added.
  Old statsci.org/smyth/pubs preprint URLs updated to 
  gksmyth.github.io/pubs URLs.
  Author lists converted to compact PubMed style without unnecessary
  punctuation. A few instances of British spellings converted to
  American spellings.

  \item
  Expand changelog.txt entries for limma 1.0 and earlier.

  \item
  Minor code update to strsplit2().

  \item
  Add test files dupcor.R, gls-series-c.R, lmfit-contrasts.R,
  lm-series-c.R, and voomlmfit-contrasts.R.
}}


\section{version 3.68.0 (2026-04-28)}{\itemize{
  \item
  New arguments 'offset' and 'offset.prior' for voom(). voom() now
  supports offsets in DGEList objects.

  \item
  New argument 'treat.lfc' for topSplice() to allow TREAT fold-change
  thresholding with diffSplice() output.

  \item
  In normalizeCyclicLoess(), default for 'adaptive.span' becomes 'TRUE',
  and the 'span' is chosen using the default settings of 
  chooseLoessSpan().
  
  \item
  The function topTableF(), which has been deprecated for many years,
  is now formally removed.
}}


\section{version 3.66.0 (2025-10-28)}{\itemize{
  \item
  New argument 'upshot' for treat().

  \item
  New argument 'span' for fitFDistUnequalDF1(), squeezeVar() and eBayes().

  \item
  The default for 'adaptive.span' in voom() and
  voomWithQualityWeights() is now 'TRUE'. All the voom functions now
  use chooseLowessSpan() by default to choose the 'span'.

  \item
  New argument 'xlab' for plotSplice() to allow labeling in terms of 
  transcripts or other splicing events instead of exons.

  \item
  The help pages for diffSplice(), topSplice() and plotSplice() have
  been updated to emphasize differential transcript usage instead of
  differential exon usage or differential splicing.
  "differential usage" concept tags have been added to the
  appropriate help pages, creating an entry in the reference manual
  index.

  \item
  makeContrasts() now uses the contrast names if 'contrasts' is a
  named vector.

  \item
  Fix typo in error message from makeContrasts() when the variable
  names are not syntactically valid.

  \item
  Bug fix to contrastAsCoef() when the 'contrast' matrix is not of
  full column rank.

  \item
  The columns of the 'covariates' matrix input to removeBatchEffect()
  are now mean-corrected.

  \item
  removeBatchEffect() now gives a message if 'design' is not
  specified.

  \item
  Bug fix for topTable() when 'confint=TRUE' and either 'p.value' is
  less than 1 or 'lfc' is greater than zero.

  \item
  Revise help page for lmFit(), removing out-of-date mention of
  default correlation.

  \item
  Minor edits to poolVar() and update of help page.

  \item
  Complete remaining dates for release versions in NEWS.Rd.
}}


\section{version 3.64.0 (2025-04-16)}{\itemize{
  \item
  In topSplice(), the default for 'test' is now "F" instead of
  "simes".

  \item
  diffSplice() is now an S3 generic function with a method 
  for MArrayLM objects.

  \item
  New argument 'legacy' for diffSplice().

  \item
  diffSplice() now adjusts transcript-level p-values using the full
  Simes method, instead of a modified Simes method leaving out the
  largest p-value for each gene that has been used until now. This 
  will cause the Simes p-values returned by diffSplice() to be
  slightly more conservative than before.

  \item
  diffSplice() now returns Bonferroni adjusted p-values as an
  alternative to Simes method.

  \item
  The 'gene.genes' data.frame returned by diffSplice() now has row
  names set to the gene IDs. Previously the last exon ID for each gene
  was used.

  \item
  diffSplice() now allows for exons with df.residual=0 and
  sigma=NA.

  \item
  Bug fix to diffSplice() when fitFDistUnequalDF1() was
  called with 'robust=TRUE', i.e., whenever 'robust=TRUE' and the
  number of transcripts differs between genes.

  \item
  Bug fix to genewise variance pooling in diffSplice() when the
  residual df are unequal. Previously the exon-wise variances were
  not weighted by df when pooling.

  \item
  New argument 'directional' for camera() and cameraPR().

  \item
  vooma() now allows missing values in both 'y' and
  'predictor', provided there are no entirely missing rows and
  provided all non-NA y values have a non-NA predictor.

  \item
  voomaLmFit() now allows missing values in 'y' and
  'predictor', provided there are no entirely missing rows and
  provided all non-NA y values have a non-NA predictor.

  \item
  voomaLmFit() now allows 'predictor' to be either a matrix or a
  column vector with genewise values.

  \item
  Edit voomaLmFit() help page. Add Li (2024) PhD Thesis as
  reference.
}}


\section{version 3.62.0 (2024-10-30)}{\itemize{
  \item
  New function fitFDistUnequalDF1() for improved empirical Bayes
  hyperparameter estimation. The new function is an alternative to
  the long-standing limma functions fitFDist() and fitFDistRobustly().
  The new function gives special
  attention to contexts where the df1 argument is unequal between cases
  or can include small fractional values, such as those resulting from
  edgeR's voomLmFit() or from the edgeR 4.0 quasi-likelihood pipeline.
  The new function uses profile
  maximum likelihood instead of moment estimation to estimate df.prior,
  and so may be more accurate than fitFDist() or fitFDistRobustly()
  even when the df1 values are all equal. The new hyperpameter estimation
  will be used by default by eBayes(), treat() and squeezeVar() whenever
  the residual degrees of freedom are unequal, but users can request
  legacy behavior for those functions by setting 'legacy=TRUE'.

  \item
  New argument 'adaptive.span' for voomWithQualityWeights() and
  normalizeCyclicLoess(), with the same purpose as the same argument for voom(). If 'TRUE' then
  an optimal value for the span is chosen based on the number of genes
  in the dataset. voom() and voomWithQualityWeights() now return the chosen 'span' as part of the
  output object if 'adaptive.span=TRUE'.

  \item
  contrasts.fit() now converts the coefficient name "(Intercept)"
  created automatically by model.matrix() to "Intercept", same as is
  done by makeContrasts(). This should avoid inconsistent coefficient
  names between the fitted model object and contrast matrices created
  by makeContrasts().

  \item
  New argument 'keep.EList' for voomaLmFit(), similar to the same
  argument for edgeR::voomLmFit(). The observation weights
  generated by voomaLmFit() are now stored only if 'keep.EList=TRUE'.

  \item
  voomaByGroup() now uses the overall mean-variance trend for groups
  with only one sample, for which a separate trend is not estimable.

  \item
  Add Ravindra et al (2023) as a reference and Mengbo Li as an author
  to the voomaByGroup() help page.
  Emphasize in the help details that only simple
  design matrices are supported by this function.

  \item
  Revise the eBayes/treat help page to further clarify the 'legacy'
  and 'fc' arguments.

  \item
  Bug fix to voomaByGroup() to save variance plot parameters
  correctly when 'save.plot=TRUE'.

  \item
  Remove argument 'zero.weights' from the MArrayLM method for
  plotMD() and plotMA().

  \item
  Fix bug in the MArrayLM methods for plotMD() and plotMA(), which
  were attempting to access a 'weights' component of the fitted model
  object, analogously to what plotMD() does with the EList or MAList
  objects, even though MArrayLM objects do not contain a 'weights' component.

  \item
  Use of the stats package digamma() function 
  replaced with statmod's logmdigamma() in fitFDist(), fitFDistUnequalDF1() and intraspotCorrelation().
  Slight improvements to speed, accuracy and code simplicity,
  but results should remain the same to at least 13 significant figures.
}}


\section{version 3.60.0 (2024-05-01)}{\itemize{
  \item
  The default settings for the 'small.n' and 'min.span' arguments of
  chooseLowessSpan() have been increased, increasing the chosen span
  value.

  \item
  New argument 'adaptive.span' for voom(). If 'TRUE', then
  chooseLowessSpan() is used to select an optimal 'span' value
  depending on the number of genes, same as is done for vooma(),
  voomaByGroup() and voomaLmFit().

  \item
  The plot information saved by voom() when 'save.plot=TRUE' now
  includes 'pch' and 'cex' plotting character settings.

  \item
  The default 'span' set by vooma() is increased slightly to the value given by
  'chooseLowessSpan(ngenes, small.n=50, min.span=0.3, power=1/3)'. A new
  argument 'legacy.span' had been added to optionally restore the old
  default settings for 'span' for users who want backward
  compatibility.

  \item
  vooma() argument 'covariate' renamed to 'predictor'.

  \item
  vooma() has been revised in several other ways to match the
  behavior of voom() more closely. The order of arguments has been
  adjusted to match voom() and lmFit(). A new argument 'save.plot'
  has been added similar to the same argument for voom(). Changes to
  have also been made to the title and xlab for the variance trend
  plot when 'predictor' is non-NULL.

  \item
  New function voomaLmFit() with the same functionality as vooma() but
  which automates the estimation of sample weights and intrablock
  correlation from vooma(). It produces an MArrayLM fit object instead
  of an EList object. It is the analogous to edgeR::voomLmFit() but
  for continuous microarray-like data.

  \item
  Edits to the help pages for decideTests(), voom(),
  voomWithQualityWeights(), vooma() and "11RNAseq". Example added to
  the vooma() help page.

  \item
  Add RNA-seq to package description.

  \item
  Add Charity Law, Goknur Giner and Mengbo Li to author list in
  package description.
}}


\section{version 3.58.0 (2023-10-25)}{\itemize{
  \item
  New argument 'covariate' for vooma() that allows the variance
  trend to depend on an extra covariate in addition to average
  log-expression. This allows vooma() to take account of any
  predictor that is correlated with the precision of each expression
  value.

  \item
  New argument 'group' to simplify the typical calling sequence for removeBatchEffect().
  The 'group' argument provides an alternative and possibly simpler way to specify the experimental conditions to be preserved.

  \item
  treat() now works when there are NA coefficients. Previously it gave
  an error if any estimated coefficient was NA.

  \item
  Fix behaviour of topTable() with 'adjust.method=NULL' so that it
  gives "BH" adjustment as documented. Previously the argument was
  passed to p.adjust(), which gave "holm" adjustment.

  \item
  Convert Introduction vignette from Sweave/pdf from to Rmarkdown/html.
  The Introduction vignette is also somewhat expanded.

  \item
  cameraPR() has been addeed to the '10GeneSetTests.Rd' help topic.

  \item
  The pdf limma User's Guide is now formally listed as a vignette.

  \item
  Non-user-visible changes:
  limma now imports statmod rather than suggesting it.
  BiocStyle has been added to the Suggests field of DESCRIPTION.
}}


\section{version 3.56.0 (2023-04-25)}{\itemize{
  \item
  Rename readSampleInfoFromGEO() to sampleInfoFromGEO().

  \item
  Add new argument 'p.value' for topGO() and topKEGG().

  \item
  Fix a problem with inconsistent pathway names from KEGG. Due to a
  change on the KEGG website, pathway names from getGeneKEGGLinks()
  have "path:" prefixes but pathway names from getKEGGPathwayNames()
  do not. To make the pathway names consistent, getGeneKEGGLinks() now
  removes the "path:" prefix.

  \item
  eBayes() now checks whether 'fit' is a list object before
  undertaking other tests.

  \item
  Revise voom.Rd to specify that y-values in voom plot are sqrt
  residual standard deviations.

  \item
  Update out-of-date URLs in the help pages and the User's Guide.
  Convert reference URLs to DOIs.

  \item
  Convert CITATION to use the newer bibentry() instead of citEntry().
}}


\section{version 3.54.0 (2022-11-02)}{\itemize{
  \item
  New function goanaTrend() to estimate a covariate-dependent trend in
  the probability of differential expression for use in goana() and
  kegga() with 'trend=TRUE'. The trend estimated by goanaTrend() is
  squeezed slightly towards constancy to provide stability when the
  number of DE genes is small. The amount of squeezing decreases with
  the number of DE genes.

  goana() and kegga() now call goanaTrend() when the 'trend' argument
  is used. When 'plot=TRUE' the plot is now created by goanaTrend()
  instead of by barcodeplot().

  \item
  Rename argument 'prior.prob' to 'null.prob' in goana() and kegga().

  \item
  Update goana() and kegga() code to catch NA gene IDs or covariate
  values.

  \item
  goana() and kegga() no longer give an error when a trend is
  estimated but there are no DE genes.

  \item
  kegga() now uses https instead of http links when reading from rest.kegg.jp.

  \item
  Add new argument 'fc' to topTable() and topTableF().
  Change default for 'lfc' to NULL instead of 0.

  \item
  lmFit() now checks explicitly for NAs in the design matrix.

  \item
  Update references listed on help pages to use DOIs instead of URLs.
}}


\section{version 3.52.0 (2022-04-27)}{\itemize{
  \item
  New function readSampleInfoFromGEO().

  \item
  Allow the 'trend' argument of eBayes() to specify a general variance
  covariate.

  \item

  More detailed checking and error messages when the 'object' input to
  lmFit() is a data.frame.

  \item
  Improve checking for incorrectly specified 'contrasts' argument in
  contrasts.fit().
}}


\section{version 3.50.0 (2021-04-27)}{\itemize{
  \item
  Improve help pages for fry() and barcodeplot().

  \item
  Revise Section 18.1.10 of User's Guide so that coloring of
  X and Y genes is consistent between plots.

  \item
  Fix bug in the MDS class method of plotMDS() (introduced in the previous release) when new 'dim.plot'
  values are set.

  \item
  Fix bug in read.idat() when 'annotation' contains NA values.
}}


\section{version 3.48.0 (2021-05-20)}{

  \subsection{New functionality}{\itemize{

    \item
    Explicitly setting 'weights=NULL' in a call to lmFit() no longer
    over-rides the 'weights' value found in 'object'.
    Default settings in lmFit() changed from 'ndups=' and 'spacing=1'
    to 'ndups=NULL' and 'spacing=NULL', although this doesn't change function behavior
    from a user point of view.

    \item
    A number of improvements to duplicateCorrelation() to make the results more robust
    and to make the interface consistent with lmFit().
    duplicateCorrelation() sets 'weights' same as lmFit(). Setting
    'weights=NULL' in the function call no longer overwrites weights
    found in 'object'.
    duplicateCorrelation() now checks whether the block factor is
    spanned by the design matrix. If so, it returns intrablock
    correlations of zero with a warning. Previously this usage error
    was not specifically trapped and could lead to correlations that
    were or NA or close to 1 depending on floating point errors.
    duplicateCorrelation() now issues a simplified message when design
    is not of full rank and uses message() to do so instead of cat().
    In terms of output, duplicateCorrelation() now bounds the genewise correlations away
    from the upper and lower bounds by 0.01 so that the correlation
    matrix will always be positive-definite.
    There is also a fix to the value returned by duplicateCorrelation() when no blocks or
    duplicates are present.

    \item
    New argument 'fc' for treat() so that the fold-change threshold can
	optionally be specified on the fold-change scale rather than as a
	log2-fold-change.

    \item
    plotMDS() no longer calls cmdscale() but instead performs the necessary
    eigenvector computations directly. Proportion of variance explained
    by each dimension is now computed and is optionally added to the dimension labels.
    The 'ndim' argument is now
    removed. All eigenvectors are now stored so that plotMDS.MDS does
    not need to recompute them when different dimensions are plotted.

    \item
    New arguments 'path' and 'bgxpath' for read.idat().
    read.idat() now checks for gzipped IDAT files and, if detected,
    gives an informative error message.
    read.idat() now checks for existence of input files before calling
    illuminaio read functions.

  }}

  \subsection{Other code improvements}{\itemize{

    \item
    The average log-expression column written by write.fit() now has
	column heading "AveExpr" to match the output from topTable().
	The column was previously called "A".
 
  }}

  \subsection{Documentation}{\itemize{

    \item
    Additional documentation for the 'design' argument of 'lmFit'
    using the term "samples" instead of "arrays" and mentioning that
    the design matrix defaults to 'object$design' when that component
    is not NULL.
 
    \item
    duplicateCorrelation() help page revised including new code example.
 
    \item
    Help page for voom() now explains that the design matrix will be
    set from the 'group' factor of the DGEList object if available.
 
    \item
    coolmap() help page now clarifies which heatmap.2() arguments are
    reserved and which can be included in the coolmap call.
 
  }}

  \subsection{Bug fixes}{\itemize{
 
    \item
    Fix typo in voom() warning when negative counts are detected.
 
  }}
}


\section{version 3.46.0 (2020-10-28)}{

  \subsection{New functionality}{\itemize{

    \item
    New function chooseLowessSpan().

    \item
    fitFDist() with a non-NULL 'covariate' now fits a smoother trend
    (with fewer spline knots) than before unless there are at least
    30 observations. This also affects squeezeVar() and eBayes() with
    'trend=TRUE'.

    \item
    New argument 'output.style' to weightedLowess().

    \item
    write.fit() now outputs a blank column heading for the row names
    so that the number of column names in the delimited output file
    will agree with the number of columns (similar to write.csv in the
    base package).

  }}

  \subsection{Code improvements}{\itemize{

    \item
    Subsetting and contrasts.fit() now work on MArrayLM objects even
    if the cov.coefficients component is absent.

    \item
    volcanoplot() now prompts user to run eBayes() if the object
    doesn't contain p-values.

    \item
    Add checks and more informative error messages to getEAWP() and
    lmFit() when the data object is missing, NULL or has zero rows.

    \item
    Subsetting for RGList, MAList, EListRaw, EList, MArrayLM or
    TestResults objects now accepts arguments other than 'i' or 'j'
    but ignores them without an error. This is relevant if a user
    adds a 'drop' argument by analogy with matrix subsetting. Previously
    that produced an error.

    \item
    Subsetting of TestResults now requires two arguments, same as all
    of the other data object classes listed above. Previously single
    index subsetting (for example object[i]) was allowed.

    \item
    Improved treatment of zero weights by weighted.median().

    \item
    More consistent use of rep.int(), rep_len() and rep() throughout
    the package.

    \item
    Replace NA with NA_real_ where appropriate.

    \item
    Add a message to topTableF() to warn that the function is
    obsolete and will be removed in a future version of limma.
    Remove usages of topTableF() from the User's Guide.

    \item
    Remove obsolete function toptable(), which has been officially 
    deprecated in favor of topTable() since 1 Feb 2018.

  }}

  \subsection{Documentation}{\itemize{

    \item
    Revise and expand the weightedLowess help page.

    \item
    Add more explanation about "heirarchical" method in the details
    section of the decideTests() help page.

    \item
    Add Note to voom help page to emphasise that voom is not an appropriate
    tool to analyse expresssion measures that have been normalized for library size.

    \item
    Revise the help page for plotMDS() to clarify that the function can
    produce PCA plots as well as PCoA plots.

    \item
    Add example to EList-class help page.

    \item
    Edits to the eBayes help page.

    \item
    Edit concept entries for help files (for reference manual index).

  }}

  \subsection{Bug fixes}{\itemize{

    \item
    Revise write.fit() so that column names of the output file are
    correct when 'digits=NULL' and the 'fit' object contains only one
    coefficient column. Previously the same contrast name was repeated
    as the column name for the coefficients, t-statistics and p.values,
    making them hard to distinguish.
 
    \item
    Bug fix to arrayWeights() when 'y' contains NAs, the design matrix
    has several columns and some but not all genes have no residual df.

  }}
}


\section{version 3.44.0 (2020-04-28)}{

  \subsection{New functionality}{\itemize{

    \item
    Add new argument 'package' to changeLog(). changeLog() now works
    for any package rather than only for limma.

  }}

  \subsection{Code improvements}{\itemize{

    \item
    Improved treatment of NA coefficients by contrasts.fit().
    contrasts.fit() has always removed coefficients that were NA because
    of singularities in the design matrix, but any extra NA coefficients
    caused by NA expression values would cause all the stdev.unscaled
    values returned by contrasts.fit() to be NA for those genes.
    contrasts.fit() now returns non-NA coefficients and
    non-NA stdev.unscaled values if all the NA coefficients are
    multiplied by zero contrast multipliers.

    \item
    Speed improvement for contrasts.fit() if some input coefficients are
    not involved (have zero multipliers) in any of the contrasts.

    \item
    contrasts.fit() now preserves the 'pivot' component of the fit
    object instead of removing it.

    \item
    voom now checks explicitly for NA or negative counts and gives an
    informative error message.

    \item
    fitFDist() now allows a covariate trend even with fewer than 4
    useable observations. (Useable means positive residual df and
    finite variance and covariate values.)

    \item
    Internal code cleaned up to avoid partial matching of function
    arguments, attributes or list component names. The automatic package
    tests are now run with the warnPartialMatchArgs,
    warnPartialMatchAttr and warnPartialMatchDollar options all set to
    TRUE.

  }}

  \subsection{Documentation}{\itemize{

    \item
    More explanation of 'covariates' argument to removeBatchEffects().

  }}

  \subsection{Bug fixes}{\itemize{
 
    \item
    Fix bug in contrasts.fit() when some coefficients are NA due to NA
    expression values and the corresponding contrast multipliers are
    zero. The new code fixes the improvement originally introduced in
    limma 3.35.9. The new code will return non-NA contrasts if the NA
    coefficients all get 0 weight in the contrast. Previous code was
    correct if all genes had the same NA coefficients but could give an
    avoidable NA result if some genes had NA coefficients but the
    first gene did not.

    \item
    Improvements to subsetting MArrayLM objects by column. Previously
    the presence of non-estimable coefficients sometimes caused the
    cov.coefficient matrix to be subsetted incorrectly; that is now
    fixed. Subsetting zero columns is now allowed even when F-statistics
    are present.

    \item
    Fix to classifyTestsF() when one or more of the coefficients have
    zero variance. Previously the zero variance will introduce NAs into
    the correlation matrix and hence cause an error. This is now
    avoided. This also fixes eBayes() when one or more of the
    coefficients are identically zero (usually caused by an all-zero
    contrast).

  }}
}


\section{version 3.42.0 (2019-10-30)}{

  \subsection{New functionality}{\itemize{

    \item
     New head() and tail() methods for all limma data classes.

    \item
    New unique() and subsetting methods for TestResults objects.
    Previously a
    subsetted TestResults object became an ordinary numeric matrix, but
    now the TestResults class is preserved. Single index subsetting is
    also allowed but produces a numeric vector.

    \item
    New argument 'gene.weights' for fry() so as to match the arguments
    and behavior of mroast(). Previously gene weights could be input to
    fry() only through the 'index' argument.
    fry's mixed p-values now take account of gene weights.
    Previously gene weights were not used in fry's mixed
    p-value calculation.

    \item
    Add arguments 'block', 'correlation' and 'weights' to voom() and
    remove the '...' argument.

    \item
    Update the arguments of voomWithQualityWeights() to match the
    revisions to arrayWeights() in limma 3.40.0.
    Add new argument 'var.group'.

    \item
    New arguments 'xlab', 'ylab' and 'main' for plotFB(). The default
    plot title now includes the column name from the data object.
    Default for 'pch' increased to 0.3.

    \item
    plotWithHighlights and plotMD.MArrayLM now gives special treatment
    to the case that 'status' is a TestResults object.

  }}

  \subsection{Code improvements}{\itemize{

    \item
    roast() and mroast() now use less memory. Memory usage now
    remains bounded regardless of the number of rotations used.
    Both functions are faster than before when \code{approx.zscore=TRUE}.
    A new argument 'legacy' is added to allow users to turn these improvements
    off if they wish to reproduce numerical results from limma 3.40.0 exactly.
    The default number of rotations 'nrot' is increased  
    from 999 to 1999.
    A bug has been fixed for roast() and mroast() when 'index' is a character vector of geneids or
    rownames. Previously this usage failed.

    \item
    Complete rewrite of zscoreT(), which is now somewhat faster and
    offer two new options for the normalizing transformation. The new
    argument 'method' indicates which transformation is used when
    'approx=TRUE'. The default approximation method is changed from
    "hill" to "bailey".
    zscoreT() now allows NA values when 'approx=TRUE' and works
    correctly for all valid df values. Previously 'approx=TRUE'
    returned NA results for 'df' infinite or <= 0.05. A simple robust
    approximation from Wallace (1959) has been introduced to handle
    very large or very small 'df' values. Previously 'approx=FALSE'
    treated any 'df' greater than 10000 as infinite; this threshold is
    now raised to 1e300.

    \item
    tZScore() is slightly faster.

    \item
    All uses of the approx() function now set
    'ties=list("ordered",mean)', except for fitFDistRobustly()
    which uses 'ties=mean'.
    The change was prompted by a change in
    R 3.6.0 whereby using the default for 'ties' generates a warning
    message whenever ties are present in 'x'. The new code gives the
    same results but is very slightly faster and avoids the warning
    message.
    limma now depends on R >= 3.6.0 because of this change.

  }}

  \subsection{Documentation}{\itemize{

    \item
    Update Yoruba case study in User's Guide to call voom() and
    duplicateCorrelation() twice each instead of once.

    \item
    Update Users Guide: update Rsubread reference, correct remarks
    about prior distribution for 1/sigma^2 in Section 13.2 and rerun
    the Yoruba and Pasilla case studies.

    \item
    Use "RNA-seq" and "ChIP-seq" consistently in the documentation instead of
    "RNA-Seq" and "ChIP-Seq".

  }}

 \subsection{Bug fixes}{\itemize{

   \item
   Fix bug in barcodeplot() when 'index' or 'index2' are character
   vectors.

  }}
}


\section{version 3.40.0 (2019-05-03)}{\itemize{

\item
Major rewrite of arrayWeights() to improve speed and stability.
  The arrayWeightsSimple() function has been removed and its
  functionality incorporated into arrayWeights().
A new 'var.group' argument has been added to simplify specification of the 
  variance design matrix.
The weights are now squeezed slightly towards unit and a new argument
  'prior.n' has been added to control the prior weight with which the array weights are squeezed.
arrayWeights() now chooses between the REML and gene-by-gene algorithms automatically by default. REML is chosen
  when there are no prior weights or missing values and gene-by-gene
  is used otherwise.
arrayWeights() now checks for and skips over any genes with zero
  residual variances.

\item
voom() now checks for experiments that have no replication. In this
  case it now returns weights that are all 1, with an informative
  warning but no error.
Two references have been added to the voom help page.

\item
classifyTests.Rd now focuses on nested F-tests.
The function FStat() is now removed. It was a convenience wrapper
  for classifyTestsF() with fstat.only=TRUE but is not used often
  enough to justify itself.
classifyTestsP() is no longer exported as it is not intended to be
  called by users.
The function classifyTestsT() is now removed, having been
  superceded by more sophisticated multiple testing methods.

\item
contrasts.fit() now removes any test statistics or p-values found
  in the fit object. This avoids any potential mis-match between coefficients and
  test statistics if eBayes() has been run on the object previously.

More details have been added to the warning about approximations in the
  contrasts.fit() help page.

\item
The old function ebayes(), which was deprecated a year ago in favor
  of eBayes(), is now removed.

\item
The fitted() and residuals() methods for MArrayLM objects now give
an informative error message if the object contains contrasts
  instead of the original coefficients.

\item
Rewrite Filtering chapter of User's Guide.

\item
Add Research Fields to the biocViews entry in the package
  DESCRIPTION file.
}}


\section{version 3.38.0 (2018-10-31)}{\itemize{

\item
New function plotExonJunc() to plot results from diffSplice().

\item
New function logsumexp().

\item
New argument hl.col for volcanoplot(), allowing users to specify
  the color for the gene names when highlight > 0.

\item
barcodeplot() no longer assumes that 'statistic'
  has unique names. Previously it returned an error if
  names(statistic) contained any duplicated values.

\item
The colors "blue", "red" and "yellow" used by coolmap() changed to
  "blue2", "red2" and "yellow2" when used in a color panel with white.

\item
goana.Rd now explains more explicitly that p-values are unadjusted
  for multiple testing.

\item
arrayWeights.Rd now mentions minimum dimensions for expression
  object.

\item
More advice on how to choose 'lfc' added to the treat() help page.

\item
Minor bug fix to the mixed p-value from roast() and mroast() when
  set.statistic="floormean".

\item
Bug fix for cumOverlap(), which was under-counting overlaps in some
  cases.

}}


\section{version 3.36.0 (2018-05-01)}{\itemize{

\item
New arguments 'quote' and 'row.names' for write.fit(). Write.fit()
  now outputs row.names by default, which previously were suppressed.

\item
coolmap() will now accept an arbitrary vector of colors. A new
  preset "whitered" panel is now also supported.

\item
In the past, contrasts.fit() always returned NA contrasts for genes
  with any NA fitted coefficients. contrasts.fit() will now ignore an
  NA coefficient if the contrast multiplier is zero for that
  coefficient.

\item
The treat() default for 'lfc' changed from lfc=0 to lfc=log2(1.2).

\item
kegga() and goana() now check whether a data.frame has been input
  by mistake and generates an error. Previously a data.frame value
  for 'de' was interpreted as a list of gene sets without any error.

\item
voomWithQualityWeights() now returns 'sample.weights' as a column
  of the 'targets' data.frame instead of as a separate vector.

\item
"TestResults" objects now include a 'labels' attribute, defaulting
  to c("Down","NotSig","Up").

\item
Functions ebayes() and toptable(), long ago replaced by eBayes() and
topTable(), are now formally deprecated.

\item
Update the User's Guide case study that analyses Lrp- E. Coli
  samples and profiled by Affymetrix arrays.

\item
Update the Agilent single-channel case study in Section 17.4 of the
  Users' Guide to use the Agilent gIsWellAboveBG detection column.

\item
Comments added to voom.Rd and eBayes.Rd to clarify the relationship
  between limma-trend and voom.

\item
In the kegga help page, clarify default gene ID system used by
kegga() when species="dme".

\item
Update wsva() help page to describe correct number of columns for output.

\item
bug fix for goana() when 'covariate' is specified and some of the
  'universe' genes don't have GO annotation.

\item
Fix to how roast(), mroast() and fry() handle weights. The 'weight' argument
can now be a matrix, or a vector of length nrow(y), or a vector of length ncol(y).
This allows the functions to accept precision weights, or gene weights, or sample weights.

}}


\section{version 3.34.0 (2017-10-31)}{\itemize{

\item 
read.idat() can now handle idat files in SNP format.

\item
New argument 'sort.by' for topSplice().

\item
diffSplice() now treats any NA geneid as a separate gene with one
  exon. Previously all NA geneids were treated as the same gene.

\item 
Bug fixes for plotSplice() which, in some circumstances, was not
  highlight significant exons.

\item
Default value for 'style' in volcanoplot() changed to "p-value"
  instead of "p". This doesn't change the behavior of the function.

\item
volcanoplot() didn't work with MArrayLM objects produced by
  treat(). Now fixed.

\item 
The gene.info.file argument of alias2SymbolUsingNCBI() will now
  optionally accept a data.frame instead of a file name.

\item
alias2SymbolUsingNCBI() now always produces a data.frame.
  Previously a single-column data.frame was returned as a vector.

\item
Bug fixes to camera() and cameraPR() when 'index' contains character row names.

\item
New argument 'restrict.universe' for the default kegga() method.
  The default is now not to the restrict the universe to genes with
  KEGG annotation.

\item
goana.default() code is rewritten to more closely match
  kegga.default. Slight speed improvement. Prior probabilities are now
  computed using the restricted universe with GO annotation instead of
  the whole universe.

\item
Bug fix for kegga() when trend=TRUE or a prior.prob or covariate is
  specified. Previously the Down p-values were substantially too small
  and the Up p-values were too large.

\item
plotWithHighlights() checks whether 'status' is a factor and
  converts to character.

\item
Bug fixes for beadCountWeights(). Default is now set correct for
  'design' and the function now works correctly when 'y' is an EList
  object contain bead standard errors but not standard deviations.

}}


\section{version 3.32.0 (2017-04-25)}{\itemize{

\item
New function cameraPR(), which implemented a pre-ranked version
of camera().

\item
New function alias2SymbolUsingNCBI(), which converts gene aliases
  or synonyms into official gene symbols using an NCBI gene-info
  file.

\item
New function wsva() for weighted surrogate variable analysis.

\item
New function coolmap(). This is essentially a wrapper for the
  heatmap.2() function in the ggplots package, but with sensible
  default settings for genomic log-expression data.

\item
decideTests() is now an S3 generic function with a default method
  and a method for MArrayLM objects. decideTests() now selects all
  null hypotheses as rejected if p.value=1.

\item
length() methods removed all limma data objects (objects of class
  EList, EListRaw, RGList, MAList or MArrayLM). length(x) will now
  return the number of list components in the object rather than the
  number of elements in the expression matrix.

\item
New argument 'style' for volcanoplot(). The default is now to use
  -log10(p-value) for the y-axis instead of the B-statistic.

\item
New argument 'xlab' for barcodeplot().

\item
New argument 'col' for plotSA(). plotSA() now longer plots a
  lowess curve trend, but if appropriate both high and low outlier
  variances are highlighted in a different color.

\item
Argument 'replace.weights' removed from voomWithQualityWeights().
  The function now always produces an EList, similar to voom(). The
  default behavior of the function is unchanged.

\item
barcodeplot() now ranks statistics from low to high, instead of
  from high to low, following the usual style of axes in R plots.
  This means that left and right are now interchanged.

\item
plotSA() now plots quarter-root variances instead of
  log2(variances).

\item
Default for 'legend' argument of plotWithHighlights() changed from
  "topleft" to "topright".

\item
fitFDist() now estimates the scale by mean(x) when df2 is estimated
  to be Inf. This will make the results from eBayes() less
  conservative than before when df.prior=Inf.

\item
plotSA() now indicates, by way of an open plotting symbol, any
  points that have low robust df.prior values.

\item
Clearer error message from fitFDistRobustly() when some variances
  are zero.

\item
C functions are now registered using R_registerRoutines.

\item
Bug fix for contrastAsCoef() when there is more than one contrast.
  Previously the coefficients for the transformed design matrix were
  correct only for the first contrast.

\item
Bug fix for kegga() when the universe is explicitly specified.

\item
Bug fix for fitFDistRobustly() when there is an extreme outlier.
  Previously floating point underflow for the outlier p-value could
  cause an error.

\item
Bug fix to mroast(), which was ignoring 'geneid' argument.

\item
Bug fix to printHead() for arrays with 1 column.
}}


\section{version 3.30.0 (2016-10-18)}{\itemize{

\item
  New function cumOverlap() to analyse the overlap between two
  ordered lists.

\item
  New function detectionPValues() to compute detection p-values from
  negative control probes.

\item
  New function fitmixture() to estimate genewise fold changes and
  expression values from mixture experiments. Previously this
  funtion was only available as part of the illumina package available
  from http://bioinf.wehi.edu.au/illumina

\item
  New function logcosh() to compute log(cosh(x)) accurately without
  floating underflow or overflow.

\item
  The default settings for the 'inter.gene.cor' and 'use.neg.cor'
  arguments of camera() have been changed. camera() now uses by
  default a preset value for the inter-gene correlation. This has the
  effect that it tends to rank highly co-regulated,
  biologically-interpretable more highly than before.

\item
  New flexibility for the roast() and mroast() functions, similar to
  that previously implemented for fry(). The index vector for each
  gene set can now be a data.frame, allowing each gene set to have
  its own set of gene weights. The indices can now optionally be a
  vector of gene names instead of a vector of indices. roast() and
  mroast() now support the robust empirical Bayes option of
  squeezeVar(). roast() and mroast() can now accept, via '...', any
  argument that would be normally passed to lmFit() or eBayes().

\item
  Slight change to the standardize="posterior.sd" method for fry().

\item
  goana(), alias2Symbol() and alias2SymbolTable() now work for any species for which an Entrez Gene based
  organism package exists.

\item
  The 'species' argument of kegga() can now accept any Bioconductor
  species abbreviation.

\item
  topGO() now breaks ties by the number of genes in the GO Term and
  by the name of the Term if the p-values are equal. (This is the
  same behavior as topKEGG.)

\item
  New argument 'plot' for plotMDS(), to optionally allow an MDS object
  to be returned without making a plot.

\item
  New arguments 'annotation' and 'verbose' for read.idat().
  The first of these to allows users to read
  any required columns from the manifest file.

\item
  New arguments 'pch.y' and 'pch.z' for plotRLDF.

\item
  Unnecessary argument 'design' removed from fitted.MArrayLM.

\item
  normalizeBetweenArrays() now checks whether the input 'object' is a
  data.frame, and converts to a matrix if possible.

\item
  duplicateCorrelation() now expands weights using expandAsWeights(),
  making it consistent with lmFit().

\item
  The lowess line drawn by plotSA() is now more robust with respect
  to NA variances.

\item
  More informative error message from voom() when there is only one
  row of data.

\item
  More informative error message from getEAWP() and lmFit() when
  'object' is a non-normalized data object.

\item
  Update Phipson et al (2016) reference for robust empirical Bayes.

\item
  Bug fix to fitFDistRobustly(), which in some circumstances was
  trying to save extra (undocumented) results that had not been
  computed.

\item
  Bug fixes for fry() with robust=TRUE or when 'index' has NULL names.

\item
  Bug fix to propTrueNull() when method="hist" and all the p-values
  are less than 1/nbins.

\item
  Bug fix to alias2Symbol() with expand.symbol=TRUE.

}}


\section{version 3.28.0 (2016-05-04)}{\itemize{

\item
  Improved capabilities and performance for fry(). fry() has two
  new arguments, 'geneid' and 'standardize'. The index argument of
  fry() can now be a list of data.frames containing identifiers and
  weights for each set. The options introduced by the standardize
  argument allows fry() to be more statistically powerful when genes have unequal
  variances. fry() now accepts any arguments that would be suitable
  for lmFit() or eBayes().
  The 'sort' argument of fry() is now the same as for mroast().

\item
  roast(), mroast(), fry() and camera() now require the 'design'
  matrix to be set. Previously the design matrix defaulted to a
  single intercept column.

\item
  Two changes to barcodeplot(): new argument 'alpha' to set
  semitransparency of positional bars in some circumstances;
  the default value for 'quantiles' is increased slightly.

\item
  kegga() has several new arguments and now supports any species
  supported by KEGG. kegga() can now accept annotation as a
  data.frame, meaning that it can work from user-supplied annotation.

\item
  New functions getGeneKEGGLinks() and getKEGGPathwayNames() to get
  pathway annotation from the rest.kegg website.

\item
  topKEGG() now breaks tied p-values by number of genes in pathway and
  by name of pathway.

\item
  goana() now supports species="Pt" (chimpanzee).

\item
  plotRLDF() now produces more complete output and the output is
  more completely documented. It can now be used as a complete LDF
  analysis for classification based on training data. The argument
  'main' has been removed as an explicit argument because it and
  other plotting arguments are better passed using the ... facility.
  New arguments 'ndim' and 'plot' have been added. The first allows
  all possible discriminant functions to be completed even if only
  two are plotted. The second allows dicriminant functions to be
  computed without a plot.
  
  plotRLDF() also now uses squeezeVar() to estimate how by how much the
  within-group covariance matrix is moderated. It has new arguments
  'trend' and 'robust' as options to the empirical Bayes estimation
  step. It also has new argument 'var.prior' to allow the prior
  variances to be optionally supplied. Argument 'df.prior' can now be
  a vector of genewise values.

\item
  New argument 'save.plot' for voom().

\item
  diffSplice() now returns genewise residual variances.

\item
  removeExt() has new 'sep' argument.

\item
  Slightly improved algorithm for producing the df2.shrunk values
  returned by fitFDistRobustly(). fitFDistRobustly() now returns
  tail.p.value, prob.outlier and df2.outlier as well as previously
  returned values. The minimum df.prior value returned by
  eBayes(fit, robust=TRUE) may be slightly smaller than previously.

\item
  tmixture.vector() now handles unequal df in a better way. This
  will be seen in slightly changed B-statistics from topTable when
  the residual or prior df are not all equal.

\item
  If a targets data.frame is provided to read.maimages() through the
  'files' argument, then read.maimages() now overwrites its rownames
  on output to agree with the column names of the RGList object.

\item
  More explicit handling of namespaces.
  Functions needed from the grDevices, graphics, stats and utils
  packages are now imported explicitly into the NAMESPACE, avoiding any possible masking by other packages.
  goana(), alias2Symbol(), alias2SymbolTable() now load the name spaces of GO.db
  the relevant organism package org.XX.eg.db instead of loading the
  packages.
  This keeps the user search path cleaner.

\item
  Various minor bug fixes and improvements to error messages.

}}


\section{version 3.26.0 (2015-10-14)}{\itemize{

\item
  New functions kegga() and topKEGG() to conduct KEGG pathway over-representation analyses.
  kegga() is an S3 generic function with a method for MArrayLM fitted model objects as well a default method.
  It uses the KEGGREST package to access up-to-date KEGG pathways using an internet connection.

\item
  goana() with trend now chooses the span of the tricube smoother
  based on the number of differentially expressed genes. Larger
  spans give smoother trends and are used for smaller numbers of
  genes. Arguments 'species' and 'plot' are no longer part of
  goana.MArrayLM() and are instead passed if present in a call to goana.default().

\item
  topGO() has a new argument 'truncate.term' to limit number of
  characters appearing in a column of the output data.frame.

\item
  Improved power for the romer() function. The empirical Bayes
  shrinkage of the residuals, removed on 4 April, has been
  re-instated in a corrected form. It can optionally be turned off
  using the 'shrink.resid' argument.

\item
  camera() has a new argument inter.gene.cor. This allows a preset
  inter gene correlation to be set for each test test, resulting in
  potentially less conservative tests.
  camera() no longer gives NA results when a gene has a residual
  variance exactly zero.

\item
  mroast() now sorts significant gene sets by proportion of genes
  changing as well as by p-value and set size. When p-values are 
  equal, sets with a higher proportion of changing genes will now
  rank higher in ordered results.
  The contrast argument of roast() can now optionally be a character
  string giving the name of a column of the design matrix.

\item
  fry() now produces "mixed" (non-directional) p-values and FDRs as
  well as directional p-values and FDRs.
  Bug fix for fry() when the gene set has a small number of genes.

\item
  camera(), fry(), mroast() and romer() now give more informative
  error messages when the index list is empty.

\item
  Speedup for ids2indices().

\item
  tricubeMovingAverage() has new argument 'power' and old argument
  'full.length' is removed.
  tricubeMovingAverage() now checks whether span is correctly
  specified. It now gives identical results for any span>=1 or any
  span<=0.

\item
  lmFit() and friends now always produce an output component 'rank'
  that returns the column rank of the design matrix. Previously it
  depended on the options used whether this component was given.

\item
  lmFit() now assumes that if the 'object' of expression values is a simple vector, that it represents
  a single gene rather than a single sample.

\item
  lmFit() now always treats infinite expression values as NA. In the past,
  this was done somewhat inconsistently.

\item
  Checking for NA coefficients in lmFit() is now done more
  efficiently, meaning that lmFit() is now faster on large datasets.

\item
  More informative error message from lmscFit() when correlation is
  NA.

\item
  getEAWP() now works on all eSet objects, not just ExpressionSet,
  provided there is an data element called "exprs".
  This allows all the linear modelling functions to handle general eSet objects.

\item
  plotMD() can now accept arguments 'array' or 'coef' as appropriate
  as synonyms for 'column'. This may help users transitioning from
  plotMA ot plotMD.

\item
  Arguments hi.pch, hi.col and hi.cex of plotWithHighlights() renamed
  to hl.pch, hl.col and hl.cex. ('hl' is short for 'highlight'.)

\item
  New 'tolerance' argument for read.idat() to allow manifest and idat
  to differ by a certain number of probes

\item
  The 'genes' argument of controlStatus() now accepts EListRaw or
  Elist objects.

\item
  All data classes (EListRaw, EList, RGList, MAList, MArrayLM) now 
  require two arguments for subsetting. This restores the behavior
  in limma version 3.19.10 and earlier. In limma 3.19.11 through
  3.25.13, subsetting by one argument y[i] was equivalent to y[i,].
  This will now give an error message.

\item
  The dimnames<- method for EListRaw objects now sets rownames for
  the background matrix Eb as well as for the foreground matrix E.

\item
  avearrays() now gives a more informative error message if input
  argument 'x' is not a matrix.

\item
  auROC() now allows for tied stat values.

\item
  Bug fix to write.fit() when method="global" and fit contains more
  than one column.
}}


\section{version 3.24.0 (2015-04-17)}{\itemize{
\item
Limma review article (Ritchie et al, Nucleic Acids Research, 2015) is now published, and becomes the primary citation for limma.
References in User's Guide, help files and CITATION have been updated to refer to this article.

\item
New function plotMD() for mean-difference plots. This will take
  over the functionality of the plotMA() function.

\item
mdplot() has new arguments 'columns' and 'main'.

\item
Arguments 'pch', 'cex' and 'col' of plotWithHighlights() renamed to
  'hi.pch', 'hi.cex' and 'hi.col'. Improved help page for
  plotWithHighlights().

\item
plot() methods removed for EList, MAList and MArrayLM objects. Use plotMD() instead.

\item
New function fry() provides a fast version of mroast()
  when there is little heteroscedasticity between genes.

\item 
Minor change to mroast() output so that FDR values are never
  smaller than p-values, even when mid-p-values are used.
 
\item
camera(), roast(), mroast() and romer() now give a warning if
  any arguments found in ... are not used.

\item
romer() is now an S3 generic function.
Speed improvements for romer() when the "mean" or "floormean" set
  statistics are used.
The empirical Bayes shrinkage of the residuals, introduced 5 October 2009, has been removed --
this means that romer() results will revert to being somewhat more conservative.

\item
topGO() can now sort gene ontology results by several p-value
  columns at once, using the minimum of the columns. The default is
  now to sort results using all the p-value columns found in the
  results instead of just by the first column.

\item
goana() has new arguments 'plot' and 'covariate'.
The plot argument allows users to view any estimated trend in the genewise significance values.
The help files goana.Rd and goana.MArrayLM.Rd are consolidated into one file.

\item
plotDensities() has a new argument 'legend', allowing user to
  reposition the legend.

\item
voomWithQualityWeights() has a new argument 'col', allowing the barplot of 
  sample weights to be colour-coded.

\item
Improvements and bug fixes to plotMDS().
plotMDS now stores 'axislabel', depending on the method used to
  compute distances, and uses this to make more informative axis labels.
plotMDS now traps problems when the distances are all zero or if
  one or more of the dimensions are degenerate.
Also, a bug when gene.selection="common" and top=nrow(x) has been fixed.

\item
The arguments for predFCm() have been renamed. Old argument
  'VarRel' replaced by logical argument 'var.indep.of.fc'. Old
  argument 'prop' replaced by 'all.de'. New argument
  'prop.true.null.method' passes through method to propTrueNull().

\item
vennDiagram() now fills circles with specified colors.

\item
plotMA() now has a method for EListRaw objects.

\item
voomWithQualityWeights() now stores the sample weights in the output EList object.
  The values stored in $weights are unchanged, they are still the combined 
  observation and sample weights.

\item
nec() and neqc() now remove the background Eb values if they exist in the
  EListRaw object.

\item
Acknowledgements edited in User's Guide.

\item
fix URLs for a number of references in help pages.

\item
Improvements to lowess C code.

\item
When loading suggested packages within function calls, usages of
  require() have been converted to requireNamespace() wherever
  possible. Functions from suggested packages are called using '::'.

\item
topTable() now gives an information error message if eBayes() or
  treat() have not already been run on the fitted model object.

\item
genas() now gives an error message when applied to fitted model for
  which the coefficient correlations are gene-specific.

\item
Fix a bug in the MAList method for plotDensities that was introduced with
  the use of NextMethod() in version 3.21.16.
}}


\section{version 3.22.0 (2014-10-14)}{\itemize{
\item
New functions goana() and topGO() provide gene ontology analyses of differentially genes from a linear model fit.
The tests include the ability to adjust for gene length or abundance biases in differential expression detection, similar to the goseq package.

\item
Improvements to diffSplice.
diffSplice() now calculates Simes adjusted p-values for gene level inferences, in addition to the exon level t-tests and gene level F-tests.
topSplice() now has three ranking methods ("simes", "F" or "t"), with "simes" now becoming the default.
diffSplice() also has a new argument 'robust' giving access to robust empirical Bayes variance moderation.

\item
New function plotExons() to plot log-fold-changes by exon for a given gene.

\item
New function voomWithQualityWeights() allows users to estimate sample quality weights or allow for heteroscedasticity between treatment groups when doing an RNA-seq analysis.

\item
Improvement to arrayQualightyWeights().
It now has a new argument 'var.design' which allows users to model variability by treatment group or other covariates.

\item
Improved plotting for voomaByGroup().

\item
barcodeplot() can now plot different weights for different genes in the set.

\item
Improvements to roast() and mroast().
The directional (up and down) tests done by roast() now use both the original rotations and their opposite signs, effectively doubling the number of effective rotations for no additional computational cost.
The two-sided tests are now done explicitly by rotation instead of doubling the smallest one-sided p-value.
The two-sided p-value is now called "UpOrDown" in the roast() output.
Both functions now use a fast approximation to convert t-statistics into z-scores,
making the functions much faster when the number of rotations or the number of genes is large.
The contrast argument can now optionally be a character string giving a column name of the design matrix.

\item
zscoreT() can optionally use a fast approximation instead of the slower exact calculation.

\item
symbols2indices() renamed to ids2indices().

\item
Improvements to removeBatchEffect().
It can now take into account weights and other arguments that will affect the linear model fit.
It can now accept any arguments that would be acceptable for lmFit().
The behavior of removeBatchEffect() with design supplied has also changed so that it is now consistent with that of lmFit() when modelling batches as additive effects.
Previously batch adjustments were made only within the treatment levels defined by the design matrix.

\item
New function plotWithHighlights(), which is now used as the low-level function for plotMA() and plot() methods for limma data objects.

\item
The definition of the M and A axes for an MA-plot of single channel
  data is changed slightly. Previously the A-axis was the average of
  all arrays in the dataset -- this has been definition since MA-plots
  were introduced for single channel data in April 2003. Now an
  artificial array is formed by averaging all arrays other than the
  one to be plotted. Then a mean-difference plot is formed from the
  specified array and the artificial array. This change ensures the 
  specified and artificial arrays are computed from independent data,
  and ensures the MA-plot will reduce to a correct mean-difference
  plot when there are just two arrays in the dataset.

\item
plotMDS() can now optionally plot samples using symbols instead of
  text labels.
It no longer has a 'col' argument, which instead is handled by ....

\item
vennDiagram() now supports circles of different colors for any
  number of circles. Previously this was supported only up to three
  sets. 

\item
getEAWP() will now find a weights matrix in an ExpressionSet object
  if it exists.

\item
update to helpMethods().

\item
Substantial updates to the two RNA-seq case studies in the User's Guide.
In both cases, the short read data has been realigned and resummarized.

\item
Improvements to many Rd files.
Many keyword entries have been revised.
Many usage and example lines been reformated to avoid over long lines.

\item
biocViews keywords updated.

\item
Subsetting columns of a MArrayLM object no longer subsets the design matrix.

\item
Bug fix for read.maimages: default value for 'quote' was not being
  set correctly for source="agilent.mean" or source="agilent.median".

\item
Bug fix to topTableF() and topTable().
The ordering of Amean values was sometimes incorrect when sorting by F-statistic and a lfc or p.value filter had been set.

\item
Bug fix to read.ilmn() when sep=",".
}}


\section{version 3.20.0 (2014-04-12)}{\itemize{

\item
New functions diffSplice(), topSplice() and plotSplice() provide functionality to analyse differential splicing using exon-level expression data from either microarrays or RNA-seq.

\item
New Pasilla case study added to User's Guide, demonstrating differential splicing analysis of RNA-seq data.

\item
new function weightedLowess() which fits a lowess curve with prior weights.
Unlike previous implementations of lowess or loess, the weights are used in calculating which neighbouring points to include in each local regression as well as in the local regression itself.

\item
weightedLoess() now becomes the default method used by loessFit() to fit the loess curve when there are weights.
The previous locfit and loess() methods are offered as options.

\item
linear model fit functions lm.series(), mrlm.series() and gls.series() no longer drop the dimensions of the components of the fitted object when there is just coefficient or just one gene.
Previously this was done inconsistently in some cases but not others.
Now the matrix components always keep dimensions.

\item
The functions lmFit(), eBayes() and tmixture.vector() now work even when there is just one gene (one row of data).

\item
New function subsetListOfArrays(), which is used to simplify the subsetting code for RGList, MAList, EList, EListRaw and MArrayLM objects.

\item
new function tricubeMovingAverage() for smoothing a time series.

\item
barcodeplot() has a new option to add enrichment worms to the plot, making use of tricubeMovingAverage().

\item
New plot() methods for RGList, MAList, EList and MArrayLM class objects. In each case, this produces a similar result to plotMA().
When using plot() or plotMA() on an MArrayLM object, the column is now specified by the 'coef' argument instead of by 'array'.

\item
plotMA3by2() now works on single channel data objects as well as on MAList objects.

\item
New function read.idat() to read files from Illumina expression beadarrays in IDAT format.

\item
The ctrlpath argument of read.ilmn() now defaults to the same as path for regular probes.
This means that only one path setting is required if the regular and control probe profiles are in the same directory.

\item
read.ilmn() now sets the same probe IDs as rownames for both the expression matrix E and the annotation data.frame genes, providing that the probe IDs are unique.

\item
beadCountWeights() can now work with either probe-wise standard errors or probe-wise standard deviations.

\item
treat() has new arguments robust and winsor.tail.p which are passed through to robust empirical Bayes estimation.

\item
topTreat() now includes ... argument which is passed to topTable().

\item
topTable() with confint=TRUE now produces confidence intervals based on the t-distribution instead of on the normal distribution.
It also now accepts a numeric value for the confint argument to specify a confidence level other the default of 0.95.

\item
topTable() will now work on an MArrayLM fit object that is missing the lods component, for example as produced by treat().

\item
roast() and mroast() now permit array weights and observation weights to both be specified.

\item
camera(), roast() and mroast() now use getEAWP() to interpret the data object. This means that they now work on any class of data object that lmFit() will.

\item
romer() now uses propTrueNull(method="lfdr") instead of convest().
This makes it substantially faster when the number of genes is large.

\item
genas() now uses fit$df.total from the MArrayLM object.
This prevents df.total from exceeding the total pooled residual df for the dataset.
The genas() results will change slightly for datasets for which df.prior was very lage.

\item
plotDensities() is now an S3 generic function with methods for RGList, MAList, EListRaw and EList objects.

\item
plotFB is now an S3 generic function with methods for RGList and EList data objects.

\item
New topic help pages 10GeneSetTests.Rd and 11RNAseq.Rd.
The page 10Other.Rd is deleted.
All topic help pages are now listed under 'See also' in the package introduction page accessible by ?limma.

\item
avereps() was never intended to be applied to RGList or EListRaw objects.
It now gives an error when applied to these objects instead of returning a matrix of questionable value.

\item
Bug fix: fitFDistRobustly() was failing when there were missing values or zero df values and covariate was NULL.

\item
Bug fix: vennDiagram() wasn't passing extra arguments (...) to plot() when the number of sets was greater than 3.

\item
Bug fix to topTreat(). Rownames were incorrectly ordered when p<1.

\item
bug fix to genas(), which was not handling vector df.prior correctly when the fit object was generated using robust=TRUE.

\item
bug fix to squeezeVar().
Previously there was an error when robust=TRUE and trend=FALSE and some of the estimated df.prior were infinite.

\item
bug fix to topTable() and topTableF() when sorting by F-statistic combined with p-value or lfc cutoffs.
}}


\section{version 3.18.0 (2013-10-15)}{\itemize{

\item
new function beadCountWeights() to estimate quantitative weights
from the bead counts for each probe for Illumina BeadArrays.

\item
New function contrastAsCoef(), which reforms a design matrix so
that one or more specified contrasts become coefficients. This
function is called by roast().

\item
plotMA() is now an S3 generic function.

\item
The legend argument to plotMA() can now take a character value
giving the position to place the legend. 

\item
toptable(), topTable() and topTreat() now preserve the rownames of
the fit object, unless the fit object has duplicated names, in
which case the rownames are copied to the ID column. Empty
rownames are replaced with 1:nrow(fit).

\item 
read.ilmn() no longer adds the probe IDs to the gene annotation
data.frame, leaving them instead as rownames of the expression
matrix. It longer creates a targets file since the sample names
are already preserved as column names of the expression matrix.

\item
loessFit() now uses the locfit.raw in the locfit package when
weights are provided instead of loess in the stats package. The
function now runs very efficiently even on very long data vectors.
The output results will change slightly when weights are provided.

\item
voom() now outputs lib.size as a column of targets instead of as a
separate component.

\item
cbind for EList and EListRaw objects now recognizes a design matrix
if it is present.

\item
plotMDS() now checks explicitly that there are at least 3 samples
to plot.

\item
normexp.fit.detection.p() now tolerates some non-monotonicity in
the detection p-pvalues as a function of expression.

\item
fitFDistRobustly() now uses a smoother for the smallest df.prior
values. This may result in smaller tail values than before when
a group of input x values appear to be outliers but the largest
value is not individually a stand-out value.

\item
New merge methods for EList and EListRaw objects.

\item
topTable() and treat() now give more informative error messages
when the argument fit is not a valid MArrayLM object.

\item
roast() now calls mroast() if the index vector is a list.
Bug fix to mroast(), which had been ignoring the weights.

\item
Updates to genas() function.
argument chooseMethod renamed to subset and option "n" renamed to "all". 
Function now returns NA results and a message when no genes satisfy
the criterion for inclusion in the analysis.
Some editing of help page and streamlining of code.

\item
Roles of contributors now specified in author field of DESCRIPTION
file using standard codes.

\item
Additions and updates to references in the help pages.
Removed defunct Berkeley Press links to published Smyth (2004)
article in several Rd files. Replacing with link to Preprint.
Added link to Phipson (2013) thesis in two Rd files.
Add Majewski et al reference to genas.Rd.
Add Phipson et al and Sartor et al references to squeezeVar.Rd.
Add Phipson et al reference to eBayes.Rd.
Update lmscFit and voom references.

\item
Update mammmary stem cell case study in User's Guide. As well
as reflecting changes to read.ilmn() and topTable(), this now
demonstrates how to find signature genes for particular cell type.

\item
documentation about rownames and column names and the use of
rownames(fit) and colnames(fit) added to lmFit.Rd.

\item
improvements to help pages for data classes.

\item
Edits to normalizeBetweenArrays help page (i) to further clarify
which normalization methods are available for single-channel data
and which are available for two-color data and (ii) to give a cross
citation to the neqc() function for Illumina BeadChips.

\item
Edits to voomaByGroup help page.

\item
duplicateCorrelation() now uses the weights matrix when block is
set. Previously the weights were ignored when block was used.

\item
Bug fix to subsetting for MArrayLM objects: the df.total component
was not being subsetted.

\item
bug fix to eBayes(robust=TRUE) when some of the df.prior values are
infinite.

\item
Bug fix to ebayes(), which was not passing the 'robust' argument
correctly on to squeezeVar().

\item
Bug fix to fitFDistRobustly(), which affected the estimated scale
when df2 is estimated to be Inf.
}}


\section{version 3.16.0 (2013-04-04)}{\itemize{

\item
New section in User's Guide on time course experiments with many time points.
The RNA-seq case study in User's Guide has also been revised.

\item
Improvements to various help pages including read.maimages.Rd, squeezeVar.Rd, fitFDist.Rd, trigammaInverse.Rd, normalizeRobustSpline.Rd, genas.Rd and roast.Rd.
Previously the meaning of source="agilent" was mis-stated in read.maimages.Rd.

\item
New robust method for estimating the empirical Bayes prior, called
by specifying robust=TRUE in the call to eBayes(). When this is
TRUE the output df.prior is now a vector instead of a scalar.

\item  
New function fitFDistRobustly() estimates the parameters of a
scaled F-distribution robustly using Winsorized values. Outlier
observations receive smaller values for df.prior than non-outliers.
This permits robust methods for squeezeVar(), ebayes() and eBayes(),
all of which now have a new argument wins.tail.p to specify the
tail proportions for Winsorizing.

\item
fitFDist() now permits infinite values for the covariate.
It also gracefully handles cases where the covariate takes only a small number of distinct values.
Similarly for eBayes() and squeezeVar() that call fitFDist(). 

\item
All the functions that perform gene set tests have been revised to make the input and output formats more consistent.

roast(), mroast() and camera() are now S3 generic functions, with methods for EList and MAList objects.

The order of arguments has been changed for roast(), mroast() and camera() so that the first argument is now y.

All functions that perform gene sets now use the argument 'index' to specify which genes are included in the test set.
Previously this argument was called 'iset' for roast() and romer() and 'indices' for camera().

camera() and mroast() now produce a data.frames.
Instead of separate up and down p-value columns, there is now a two-sided p-value and a column indicating direction of change. 
There are new columns giving FDR values and the number of genes in each set.
There is a new argument 'sort' to indicate whether output results should be sorted by p-value.

mroast() has a new argument 'weights' for observational weights, to bring it into line with roast(),

\item
vennDiagram() can now plot up to five sets (previously limited to three).

\item
genas() now optionally draws a plot in which ellipses are used to represent the technical and biological components of correlation.
It also now has the ability to automatically select which probes are used for the correlation analysis, and a new argument controls the method used for this selection.

\item
New options for the method argument of propTrueNull().

\item
New functions vooma() and voomaByGroup() for computing precision weights based on a mean-variance trend.
vooma() is similar to voom() but for microarray data instead of RNA-Seq.
voomaByGroup() allows different groups to have systematically different variances.

\item
New function predFCm() to compute predictive (shrunk) log fold changes.

\item
New function fitGammaIntercept() for estimating the intercept of a gamma glm with an offset.
Used by genas().

\item
New function zscoreHyper() for computing z-score equivalents of deviates from a hypergeometric distribution.

\item
New function qqf() for qq-plots relative to an F-distribution.

\item
normalizeWithinArrays() with method="robustspline" now longer requires the layout argument to be set.
The layout argument for normalizeRobustSpline() now defaults to a single print-tip group.

\item
fitFDist() now coerces degrees of freedom df1 below 1e-15 to zero.

\item
Due to changes in R, loessFit() no longer makes direct calls to foreign language code in the stats package, and instead calls R functions.
Unfortunately, this makes loessFit() about 25-30\% slower than previously when weights are used.

\item
Bug fix to read.maimages(), which was not accepting source="agilent.mean".

\item
Bug fix for contrasts.fit() when the covariance matrix of the coefficients (cov.coefficients) is not found in the fitted model object.
This situation doesn't arise using any of the standard limma analysis pipelines.

\item
Bug fix to lmscFit() when the residual df = 1.

\item
Bug fix to readTargets() to avoid warning message when targets$Label is used to set row names but targets$Label contains duplicated entries.

}}


\section{version 3.14.0 (2012-10-02)}{\itemize{

\item
limma license upgraded to GPL (>=2) instead of LGPL to match R itself.

\item
Many updates to the User's Guide. 
Sections have been added on reading single channel Agilent and Illumina data.
The chapter on experimental designs has been split into three chapters on single-channel, common reference and
two-color designs respectively. The material on the fixed effect approach
  to technical replication has been deleted.
There are new sections on nested interactions for factorial designs
and on multi-level designs.

\item
The links to the Apoa1, Weaver and Bob1 datasets in the User's Guide have been updated
to help users download the data themselves if they wish to repeat the case study analyses.

\item
The help page for camera() now cites published paper Wu and Smyth (NAR, 2012).
In view of the results of this paper, the claim is no longer made on help page for geneSetTest()
  that genes might be treated as independent when the experimental
  units are genetically identical mice.

\item
Minor edits to CITATION file.

\item
New function propTrueNull() for fast estimation of the proportion of true null hypotheses from a vector of p-values.

\item
New function zscore() to compute z-score equivalents for deviates from any continuous distribution.
Includes the functionality of the older functions zscoreGamma() and zscoreT() as special cases.

\item
roast() now accepts observation level weights, through a new
  argument 'weights'.

\item
loessFit() now applies minimum and maximum bounds by default to avoid zero or
  infinite weights. Equal weights are now treated as if the weights were NULL, even all
  zero weights, so that the lowess code is called instead of the loess code.

\item
When there are no weights, loessFit() now extracts residuals
directly from the C code output instead of computing in R.

\item
  fitFDist() now permits missing values for x or zero values for df1
  even when there is a covariate. This means that squeezeVar() and
  eBayes() now work with trends even when not all the data values
  are informative.

\item
New argument 'file' for convest(), implementing edits contributed
  by Marcus Davy. Arguments doplot and dereport renamed to 'plot'
  and 'report'.

\item
Two improvements for plotMDS().
It now coerces labels to be character, and now makes extra room on the plot when the text labels are wide.

\item
plotMDS() no longer gives an error when the requested number of top
  genes is greater than the total number of rows of data.

\item
Code speed-up for alias2SymbolTable()

\item
any(duplicated()) replaced by anyDuplicated() in several functions.

\item
Fix to voom() so that it computes weights correctly even when the
  design matrix is not of full rank.

\item
Bug fix for roast() when the fitted model has only one coefficient.

}}


\section{version 3.12.0 (2012-03-30)}{\itemize{

\item
read.maimages() with source="agilent" now reads median foreground
estimates instead of mean foreground. New option source=
"agilent.mean" preserves earlier meaning of source="agilent".

\item
Agilent single-channel case study added to User's Guide.

\item
removeBatchEffect() now corrects for continuous covariates
as well as qualitative factors.

\item
new function camera() performs competitive
gene set tests while adjusting for inter-gene correlation.

\item
new function interGeneCorrelation() estimates the average
intergene correlation for a set of genes.

\item
columns in output from roast() have been re-ordered.

\item
arguments 'selected' and 'selected2' renamed to 'index'
and 'index2' in functions barcodeplot(), geneSetTest()
and wilcoxGST().

\item
default labels for barcodeplot() are now somewhat more
explicit.

\item
new function rankSumTestWithCorrelation extends the
Wilcoxon-Mann-Whitney test to allow for correlation between
cases in one of the groups. geneSetTest() now calls this
function instead of wilcox.test, with a consequence improvement
in speed.

\item
The lfc (log-fold-change) cutoff argument of topTable() is now
applied to the minimum absolute logFC when ranking by F-statistic.
Previously lfc was only used when ranking by t-statistic.

\item
new methods "fast" and "affy" for normalizeCyclicLoess(),
with "fast" becoming the default method.
New argument 'cyclic.method' for normalizeBetweenArrays() gives
access to the different cyclic loess methods.

\item
There were problems with using the argument gene.weights in
mroast(). This argument is now permitted to be of the same
length as the number of probes in the data set. It is then
automatically subsetted for each gene set.

\item
mroast() now uses mid-p-values by default when adjusting for multiple testing.

\item
neqc(), nec() and normexp.fit.control() now give user-friendly error messages
when no negative control probes or no regular probes are found.

}}

\section{version 3.10.0 (2011-11-01)}{\itemize{

\item
New function voom() allows RNA-Seq experiments to be analysed using the standard limma pipeline.
An RNA-Seq case study is added to User's Guide.

\item
treat(), roast() and mroast() can now estimate and work with a trend on the prior variance, bringing them into line with eBayes().

\item
barcodeplot() and barcodeplot2() merged into one function.

\item
removeBatchEffect() can now correct for two batch factors.

\item
plotMDS is now an S3 generic function. This allows MDS plots to be redrawn with new labels without needing to repeat the distance or scaling calculations.
New S4 class "MDS" to hold the multidimensional scaling information output from plotMDS.

\item
getEAWP() now gets probe annotation from the expression rownames of an EList object, if no other probe annotation is available.

\item
topRomer() now ranks gene sets by secondary columns as well the primary criterion specified, to give a more meaningful ranking when the p-values are tied.

\item
wilcoxGST() now accepts signed or unsigned test statistics.
Change to p-value calculation in geneSetTest() when rank.only=FALSE to avoid zero p-values and follow recommendation of Phipson and Smyth (SAGMB, 2010).

\item
plotMA() now recognizes ElistRaw and EList objects appropriately.

\item
Default span for normalizeCyclicLoess increased from 0.4 to 0.7.
Speed improved when weights=NULL.

\item
Weaver case study (Section 11.5) in User's Guide is updated and rewritten.
Data classes ElistRaw and Elist now described in the quick start section of the User's Guide. Other minor updates to User's Guide.

\item
Bug fix for normalizeBetweenArrays() when object is an EListRaw and method="cyclicloess". Previously this function was applying cyclicloess to the raw intensities, then logging. Now it logs first, then applies cyclicloess.

\item
Bug fix to avereps() for EList objects x when x$other is not empty.

}}


\section{version 3.8.0 (2011-04-14)}{\itemize{

\item
New method "cyclicloess" for normalizeBetweenArrays().

\item
New function normalizeMedianValues().
This is now used by normalizeBetweenArrays() for method="scale" when the
data object is a matrix or EListRaw object.

\item
ebayes() and eBayes() now include an option to allow the prior variance
to be intensity-dependent. Associated with this is a new optional
argument 'covariate' for fitFDist() and squeezeVar().

\item
The total degrees of freedom used by ebayes() and eBayes() to compute
p-values, and therefore also by topTable() and so on, is now not allowed
to be larger than the total residual degrees of freedom pooled over all
probes. This limit will normally become active only when the prior
degrees of freedom is estimated to be infinite. Previously the total df
was also allowed to be infinite in this case, so that p-values were
computed as if the t-statistics followed a standard normal distribution.
The change will not affect microarray analyses with thousands of genes,
but can come into play with simulated data or with data sets containing
only a few probes. In these cases, the ebayes p-values now agree with
those from a pooled analysis in which the genewise variances are assumed
equal.

\item
df.total is now included as a component of the MArrayLM object output by
ebayes() and eBayes().

\item
topTable() and toptable() now optionally output confidence intervals for
the log-fold-changes.

\item
neqc() now supports normexp by control background correction using
negative control values inferred from detection p-values. This is mainly
applicable to Illumina BeadChips. The default value of the
'other.columns' argument of read.ilmn() is changed to "Detection".

\item
plotSA() now adds the prior variance estimated by eBayes(), if available,
to the sigma vs A plot.

\item
symbols2indices() now returns the indices of all occurrences of each
symbol on the array, not just the first occurrence. The argument names
have been changed to be more easily understood, and a new argument
remove.empty has been added to optionally remove gene sets of length
zero. Previously this was always done.

\item
normalizeWithinArrays() no longer checks for printer layout information
before doing "control" normalization.

\item
read.columns() and read.ilmn() now allow column names with leading or
trailing blanks, although the white space will be trimmed.

\item
Gene symbol annotation has been added to the Estrogen case study in the
User's Guide. The Introduction to the User's Guide now reflects limma's
capabilities for pre-processing data from single-channel and Illumina
platforms. A reference to Shi et al (2010) has been added for neqc(),
the installation instructions have been updated, and use of the sma
package has been removed from the code examples.

\item
The help pages for nec() and neqc() have been combined into one file.

\item
Bug fix to squeezeVar(): the output var.posterior was too long when there
was a non-null covariate and the prior variance was estimated to be
infinite.

\item
Bug fixes to avereps() and avearrays() for matrices, which had been
treating NA expression values as if they were zero. Averages are now
computed correctly over the non-missing values.

\item
Bug fix to mroast(), which had been assigning prior.df incorrectly for
the second and subsequent gene sets.

\item
Fix to gls.series() so that it works correctly even when there is only
one column of microarray data.

\item
Bug fix to plotMDS() when gene.selection="common".

\item
Bug fixes to rbind() for EList and EListRaw objects, and to subsetting,
cbind() and rbind() for EListRaw objects, which were not operating
correctly on the Eb component.

\item
read.ilmn() is now careful to return the gene annotation as a data.frame,
not a character vector, even when it has only one column.
}}


\section{version 3.6.0 (2010-10-18)}{\itemize{

\item
New function avearrays() to average over technical replicate arrays.

\item
New function normalizeVSN() to apply vsn normalization to limma data
objects. normalizeBetweenArrays() no longer supports method="vsn".

\item
New function mroast(), which organizes roast() tests for multiple gene
sets.

\item
New function topRomer() to summarize the output from romer().

\item
romer() and romer2() have been consolidated into one function romer().
There is a new argument 'set.statistic' to choose the set summary
statistic, and the old argument 'floor' has been removed as no longer
needed. romer() no longer outputs p-values for the "either" alternative
hypothesis; p-values are output for "Up", "Down" and "Mixed".

\item
New argument 'set.statistic' for roast() to specify which gene set
summary statistic should be used. The new default is "mean" and the
original setting is now called "msq". The output from roast() no longer
includes a column for the "either" alternative hypothesis.

\item
New S4 class "Roast" for the output from roast().

\item
The backgroundCorrect() code has been simplified and a new function
backgroundCorrect.matrix() introduced for matrix arguments. Background
intensities are now used correctly when the argument RG is an EListRaw
object.

\item
neqc() now subtracts background intensities Rb if they are found in the
data object.

\item
read.maimages() now supports ImaGene Version 9.0.

\item
read.ilmn() now uses readGenericHeader() to skip over header information
in BeadStudio files, and no longer issues a warning if control probe
files are not provided.

\item
marray is no longer a suggested package.

\item
Updates to the 04.Background help page, which gives an overview of the
background correction functions, and more complete documentation for
roast() and romer(). Citation details have been updated for roast() and
propexpr(), and the Pomelo II reference has been added to the User's
Guide.

\item
avereps() for MAList objects was not averaging the "other" data matrices
correctly. Now fixed.

\item
Bug fix to normexp.fit.control() when the argument 'x' is a matrix.

\item
Bug fixes to roast(), affecting the case when set.statistic="msq" and the
gene weights change sign, and the estimation of the proportion of active
genes when there are gene weights.
}}


\section{version 3.4.0 (2010-04-23)}{\itemize{

\item
New function normalizeCyclicLoess().

\item
read.maimages() has a new argument green.only, which replaces the older
argument channels. read.maimages() now outputs an EListRaw object
whenever single channel data is read, and an RGList object when two color
data is read.

\item
read.maimages() now provides the option of reading Agilent Feature
Extraction data using median foreground summaries.

\item
New argument 'robust' for normexp.fit.control() and neqc(), giving the
option of robust estimation of the background mean and standard
deviation. Robustness uses the huber() function from the MASS package
and is judged on the log-scale. Computation of the background mean and
standard deviation is also faster.

\item
The default value for the offset in neqc() has been reduced from 50 to
16.

\item
New argument gene.selection for plotMDS().

\item
avedups() and avereps() now have methods for EList objects.

\item
New as.data.frame() methods for MAList, EList and EListRaw objects.

\item
makeContrasts() requires coefficient names (levels) to be syntactically
valid variable names in R. The intercept column name "(Intercept)",
which is often produced by the model.matrix() function, is now
automatically converted, and a more detailed diagnostic message is output
in the case of non-valid names.

\item
Speed improvement to avereps(), achieved by calling the base package's
rowsum() function, following a suggestion of Michael Lawrence.

\item
MArrayLM objects now re-calculate F-statistics if they are subsetted by
columns. This ensures that the F-statistics always match the contrast
columns in the object.

\item
The output from treat() now includes a component df.prior.

\item
topTreat() now allows genelist to be a vector as well as a data.frame.
Previously it failed when a vector was input.

\item
The truth argument of auROC() can now be logical.

\item
Output from vennDiagram() improved to draw complete unbroken circles.

\item
New "quote" argument for read.ilmn().

\item
The raw data for the IlluminaWG-6 case study in the User's Guide is now
available online, with the URL given in the User's Guide. The probe
filtering in the Illumina case study has been corrected and the results
updated.

\item
Bug fix to backgroundCorrect() to ensure that normexp.method is always
correctly passed to normexp.fit(). Previously normexp.method was being
ignored, in favour of the default, when RG was an RGList object.

\item
Bug fix to the p-value computation in treat(). This affected only
p-values greater than 0.5.

\item
Bug fix to the row numbers returned by toptable() and topTable() when the
arguments p or lfc are specified.

\item
Bug fix to the sigma calculation in normexp.fit.control().
}}


\section{version 3.2.0 (2009-10-28)}{\itemize{

\item
New functions read.ilmn() and read.ilmn.targets() for reading Illumina
Whole Genome Expression BeadChip data.

\item
New functions normexp.fit.control() and neqc(), which perform normexp
background correction using negative control information. This is
especially useful for Illumina BeadChip data. neqc() undertakes
background correction and quantile normalization in one step for EListRaw
data objects with negative control probes.

\item
New function propexpr() for estimating the proportion of expressed probes
in a microarray, using negative control probe information. Especially
useful for Illumina BeadChip data.

\item
New function genas(), which tests for associations between contrasts in a
linear model.

\item
New function topTreat() for listing the top results from treat().
treat() now also outputs thresholded t-statistics.

\item
New function alias2SymbolTable() for converting a set of alias names to a
vector of official gene symbols of the same length.

\item
New function wilcoxGST() introduced as a synonym for geneSetTest() with
ranks.only=TRUE.

\item
New function barcodeplot2() for displaying directional gene sets with up
and down genes shown separately.

\item
New options and improved power for romer(). Empirical Bayes shrinkage of
the residuals is now used to improve statistical power, and there are two
new options for the set summary statistic. The original mean50 statistic
is provided in a new function romer2(). The new romer() uses an ordinary
mean set statistic by default, with a floor-mean set statistic as an
option. Both functions give better statistical power than the previous
version, and romer() is also computationally quicker. romer() now
outputs the number of genes in each gene set.

\item
Rotation p-values in roast() and romer() are now computed as
(x+1)/(nrot+1) instead of x/nrot, where x is the number of random
rotations giving a more extreme statistic than that observed.

\item
backgroundCorrect() and normalizeWithinArrays() now accept EListRaw
objects.

\item
propexpr() now accepts EListRaw objects, or any other data object with an
as.matrix method.

\item
getEAWP() and lmFit() now return Amean values for EList objects.

\item
arrayWeights() now displays convergence information when TRACE=TRUE for
method="genebygene". The convergence criterion for arrayWeightsSimple()
is now more lenient for large data sets, while still achieving about 5
decimal place accuracy, and failed convergence from floating point errors
is now detected and reported with a diagnostic message.

\item
strsplit2() now uses ... to pass arguments to strsplit().

\item
Speed improvements to avereps(), following suggestions from Axel Klenk.

\item
read.ilmn() now accepts more variations in the name of the ProbeID
column.

\item
A case study with Illumina BeadChip microarrays has been added to the
User's Guide. The Oshlack et al (2007) reference on normalization of
boutique arrays has been added to the normalizeWithinArrays help page and
to the User's Guide.

\item
The sma package has been removed from the list of suggested packages.

\item
Bug fix to removeBatchEffect(), which gave incorrect results when
design=NULL.

\item
Two fixes to contrasts.fit() when there are non-estimable coefficients in
the linear model. The function previously gave a warning about coercing
double to logical, and also failed in this situation when the contrast
was a single vector instead of a matrix.
}}


\section{version 2.18.0 (2009-04-21)}{\itemize{

\item
New class EList for one-channel microarray data and new class EListRaw
for non-normalized one-channel data. Basic matrix-like methods,
including cbind, rbind and subsetting, are provided for EListRaw objects.

\item
New argument 'channels' for read.maimages(), which allows the function to
read one-channel data.

\item
New function plotMDS(), which makes a multidimensional scaling plot using
a particular distance measure.

\item
New function romer() for rotation gene set enrichment analysis. This can
be viewed as a version of GSEA for linear models.

\item
New function symbols2indices(), which matches a list of gene sets as
symbols against a vector of gene symbols.

\item
New function plotSA(), which produces a sigma vs A plot for a fitted
microarray linear model. The color of the trend line drawn by plotSA()
has been changed to red.

\item
New argument expand.symbols for alias2Symbol(). This can be used to
prevent aliases from being expanded if they are already official symbols.

\item
ebayes() now returns lods as a matrix. This ensures that MArrayLM
objects can always be subsetted by column names, even when there is only
one column.

\item
lmFit() now looks for a design matrix in an MAList data object, if one
exists. The lmFit() warning message on NA coefficients is now more
refined, appearing only when the coefficient vector is partially NA for
some probes.

\item
printHead() now prints compactly for recursive lists of arbitrary depth.

\item
The C code dynamic library is now loaded when the package is loaded
rather than when it is attached, to enable other packages to import
limma.

\item
The biocViews entry has been changed from Statistics to Bioinformatics.

\item
Documentation updates for backgroundCorrect(), treat(), roast(),
geneSetTest(), topTable(), topTableF() and alias2Symbol().

\item
Bug fixes to the calculation of the up, down, mixed and either p-values
in romer().

\item
Bug fix to toptable(), which caused t-statistics to be out of order when
p<1 or lfc>0, and another which caused an error when p<1 and lfc>0.

\item
Bug fixes to the way that topTable() and topTableF() handle sort="none".

\item
Bug fix to rowname assignment for MAList and RGList objects containing
matrices in the 'other' component.

\item
Fix to gls.series() so that it does not fail when the weight matrix
contains NA values.

\item
Bug fix to printHead() to handle object components which are neither
atomic nor recursive.
}}


\section{version 2.16.0 (2008-10-22)}{\itemize{

\item
New function roast(), which does rotation gene set testing for linear
models.

\item
New function treat(), which computes empirical Bayes moderated-t p-values
relative to a minimum required fold-change threshold.

\item
New function removeBatchEffect(), which removes batch effects from
microarray data.

\item
New function barcodeplot(), which gives a graphical representation of
geneSetTest() using ranks.

\item
New function selectModel(), which implements Akaike Information Criterion
(AIC) or Bayesian Information Criterion (BIC) model selection between
alternative linear models. A new method "mallowscp" has also been added.

\item
New function alias2Symbol() to convert gene aliases to official symbols,
using the org.Hs.eg.db and similar annotation packages.

\item
New S3 generic function avereps() for averaging over irregularly spaced
replicate spots, with a method for MAList objects and a default method
intended for matrices.

\item
New function getDupSpacing(), which determines duplicate spacing from a
vector of gene IDs. getLayout() with guessdups=TRUE can now handle
values for ndups greater than 2, and getSpacing() recognizes a new value
"subarrays", meaning that all the probes on an array are repeated in two
or more subarrays.

\item
Substantial tidying up of the normexp.fit() code. The optimization step
is now done in C, speeding up the computation several fold, and the old
pure-R normexp.fit() has been removed. The old method "neldermead" is
now replaced by "saddle", and the old methods "bfgs", "nlminb" and
"nlminblog" are now replaced by "mle". All the user-visible functions
previously called by normexp.fit() have been removed, including
dnormexp(), dnormexp.saddle(), normexp.m2loglik() and
normexp.m2loglik.saddle().

\item
backgroundCorrect() has a new argument normexp.method to specify the
parameter estimation method for normexp background correction. The
options are "mle", "saddle", "rma" and "rma75".

\item
Better treatment of array weights in linear modelling. lm.series() and
gls.series() now treat array weights as a special case when there are no
missing values in the expression data and no probe-specific weights.
This means that the array weights are now incorporated into
fit$cov.coefficients, so that contrasts.fit() correctly takes array
weights into account when computing the contrast standard errors, which
it did not do previously. It also means that lmFit() is much faster when
array weighting is used. There is also a speed improvement to
gls.series() when there are no missing values and no probe-specific
weights.

\item
New capabilities and arguments for write.fit(). The argument 'method'
allows p-values to be adjusted globally instead of separately, and the
argument F.adjust allows the F-statistic p-values to be adjusted.

\item
New arguments sort.by and p.value for topTableF(), making it more
consistent with topTable(). topTableF() now also allows the input
argument 'number' to be larger than the number of rows of the fitted
model object, in which case the requested number is set to the number of
available rows.

\item
New option "none" for the sort argument of topTable().

\item
In read.maimages(), support has been added for ArrayVision data files
with 'Artifact-removed' (ARM) density values.

\item
gls.series() now preserves the rownames (probe IDs) of the expression
object. This affects lmFit() when ndups>1 or block is non-NULL, and
makes gls.series() consistent with lm.series() and mrlm().

\item
The following obsolete functions, deprecated in December 2006 in favour
of read.maimages(), have been removed: m.spot, a.spot, read.matrix,
read.series, rg.series.spot, rg.genepix, rg.quantarray and rg.spot.

\item
The criterion for the "Active" column output by roast() has been changed
to be more discriminating.

\item
Bug fixes to plotMA(), which was ignoring the array argument for RGList
objects and was not working for MArrayLM objects.

\item
Bug fixes to getEAWP() for PLMset objects, for which the probe names,
weights and Amean components were not being set correctly.

\item
Bug fix to removeBatchEffect(), which previously returned an error.
}}


\section{version 2.14.0 (2008-04-30)}{\itemize{

\item
New function getEAWP(), introduced to unify the treatment of the
different classes of microarray data objects that can be input to the
linear model functions lmFit(), arrayWeights() and printtipWeights().

\item
Streamlining of the plotMA() code. plotMA() should now work correctly
for any object class that can be coerced to a matrix, for example
LumiBatch objects.

\item
lmFit() now tries to guess whether a data matrix contains log-ratios or
log-intensities. When the input data object is a matrix, it is not
possible to know for certain which is the case, and so lmFit() did not
previously set the Amean component of the output MArrayLM object.
lmFit() now sets Amean if the input data values are all positive, because
the data are almost certainly log-intensities in that case.

\item
lmFit() now warns if coefficients become inestimable due to missing
values.

\item
All usage of Biobase's exprSet class has been removed or replaced with
ExpressionSet. The class exprSet2 has been removed.

\item
normalizeBetweenArrays(method="vsn") has been updated to reflect changes
in the vsn package, now uses vsnMatrix() when the data object is a
matrix, and gives an error if the object is an MAList. The results are
again on the log2 scale, as they were before December 2007.

\item
The method option "neldermean" for normexp.fit() has been corrected to
"neldermead".

\item
The function rlm.series(), deprecated in 2004, and the deprecated
function splitName() have been removed.

\item
arrayWeights() now gives a nice error message if the number of arrays is
too few to give statistically useful results.

\item
All uses of I() in data.frame() arguments have been removed and replaced
with stringsAsFactors=FALSE. This ensures that the appropriate columns
in data.frames have class "character" instead of "AsIs".

\item
Comments have been added to the geneSetTest() help page to warn that the
test assumes genewise independence.

\item
A section on Multiple Testing Across Contrasts has been added to the
User's Guide, and the help page for backgroundCorrect() has been edited
for improved readability.

\item
The BibTeX entry returned by citation("limma") now uses incollection
style instead of inbook, thanks to Joern Toedling.

\item
Bug fix to getEAWP(), which was not getting probe names correctly from
the rownames of the data object when the object is a matrix.

\item
Bug fix when a matrix data object containing missing values is passed to
lmFit(). The check for log-intensities vs log-ratios did not remove
missing values and hence produced an error.
}}


\section{version 2.12.0 (2007-10-08)}{\itemize{

\item
topTable() and toptable() have new arguments p.value and lfc, allowing
users to cut the gene list on adjusted p-value and fold change, rather
than just taking a set number of genes as before.

\item
backgroundCorrect() now accepts a matrix argument.

\item
New argument n.pts for normexp.fit(), which becomes the default for
method="normexp" in backgroundCorrect(). This decreases computation time
and introduces some robustness. The rule used to choose the quantiles
has also been improved to give more nearly unbiased parameter estimators.

\item
Improvements to the numeric computations of dnormexp.saddle(), which is
used by backgroundCorrect(method="normexp"). New function
normexp.m2loglik.saddle(), which is the same as normexp.m2loglik() but
uses the saddle-point approximation. normexp.fit() has a new argument
'methods'.

\item
New as.matrix() methods for vsn, ExpressionSet and LumiBatch objects, so
that these objects can be used directly in lmFit().

\item
decideTests(method="hierarchical") now passes adjust.method down to the
second (t-test) level of testing. Previously the second level testing
was always "holm". The spelling has also been corrected from
"heirarchical" to "hierarchical".

\item
contrasts.fit() now warns if the row names of the contrast matrix do not
match the column names of the coefficients.

\item
plotMA3by2() has a new argument 'device' to specify the graphic format.

\item
The plot character and size can now be passed to volcanoplot() as
arguments.

\item
The Ritchie et al (2007) reference has been added to the User's Guide and
to the backgroundCorrect and normexp help pages, and loessfit.Rd has been
edited to better acknowledge the lowess and loess authors.

\item
Correction to loessFit() because the lowess() function is now in stats
instead of base.

\item
Bug fix to read.columns() to stop a spurious warning message when
text.to.search has length greater than one.

\item
Bug fix for as.matrix.vsn(), which allows vsn objects from vsn2() to be
used directly in lmFit().
}}


\section{version 2.10.0 (2007-04-26)}{\itemize{

\item
A NAMESPACE has been introduced to the package. The functions matvec()
and vecmat() have been renamed to .matvec() and .vecmat() and are no
longer exported. All objects with names not starting with "." are
exported.

\item
Changes to the column headings from topTable(). 'M' is now 'logFC' and 'A' is
now 'AveExpr'.

\item
lmFit() now supports ExpressionSet objects, and there is a new function
as.matrix.ExpressionSet().

\item
New function strsplit2(), which replaces splitName(). splitName() is
deprecated.

\item
New function read.columns(), which is similar to read.delim() but reads
specified columns only. read.columns() replaces the functions
getColClasses() and namesInFun().

\item
New function mergeScansRG().

\item
New argument 'other.columns' for read.imagene().

\item
The old functions designed to read image data files into data.frames are
now deprecated. This affects read.matrix, rg.series.spot, read.series,
m.spot, a.spot, rg.spot, rg.quantarray and rg.genepix.

\item
wtVariables() has been renamed to namesInFun().

\item
normalizeRobustSpline() now works even when there is only one print-tip group,
although it is designed to handle multiple print-tip groups.

\item
The function .onAttach() replaces the old .First.lib().

\item
normalizeForPrintorder(plot=TRUE) now produces no output, to agree with
the documentation.

\item
The default number of digits in write.fit() has been changed to 3.

\item
The User's Guide sections on spot quality weights and quality assessment
have been expanded, a section on array quality weights has been added,
and the case studies have been updated for the current version of limma.

\item
Bug fix to read.maimages() for source="arrayvision". read.images() now reads
ArrayVision files again, which it has not done since version limma 2.3.1.

\item
Bug fixes to read.maimages() for source="quantarray", which was reading
too few lines when searching for "Begin Data" in the second and following
files and was not setting the annotation columns, and for the Agilent and
generic sources.

\item
Bug fix to read.columns(), to read SPOT files with a leading column of
row numbers.

\item
Bug fix to topTable(). The values sort="M" and resort.by="M" needed to
be updated to sort="logFC" and resort.by="logFC", and the values
sort="AveExpr" and resort.by="AveExpr" are now allowed.

\item
Bug fix to lmFit(), which was giving an error with marrayNorm objects
because of incorrect setting of Amean.

\item
Bug fix to write.fit(), which was writing ambiguous column names for the
coefficients, t-statistics and p-values when the MArrayLM object had only
one column.

\item
Bug fix to lm.series() and gls.series(), which were returning errors when
the design matrix was not of full rank and the columns were named.

\item
volcanoplot() now allows a NULL value for the 'names' argument.
}}


\section{version 2.8.0 (2006-10-04)}{\itemize{

\item
New function topTableF() to rank genes on the basis of F-statistics.
topTable() now accepts a vector argument coef and calls topTableF() to
produce a ranking by F-statistic if length(coef)>1. topTable() with
coef=NULL ranks genes by the F-statistic for all the contrasts.

\item
New argument 'coefficients' for contrasts.fit(). This provides a simple
way to specify the contrast matrix when the contrasts are just a subset
of the original coefficients.

\item
New argument 'contrasts' for makeContrasts(), which accepts a character
vector of contrasts. The matrix produced by makeContrasts() now has
named dimensions.

\item
New argument 'lfc' for decideTests() to require a minimum fold-change.

\item
normalizeWithinArrays() now uses asMatrixWeights() to process the weights
argument, so that the weights can be a vector of probe weights, a vector
of array weights, or a matrix of spot weights. The weights argument of
modifyWeights() now has a default value.

\item
New function exprs.MA() to extract the log-expression matrix from an
MAList object.

\item
New function plotlines().

\item
New function plotFG() for foreground-background plots.

\item
New function as.matrix.PLMset(), which allows lmFit() to work correctly
on PLMset objects.

\item
New arguments for vennDiagram() to allow colors and the printing of
multiple sets of counts on the same diagram.

\item
wtflags() has a new argument cutoff, and the argument w has been renamed
to weight.

\item
In geneSetTest(), the values "two.sided", "less" and "greater" are now
permitted as synonyms for "either", "down" and "up". geneSetTest() no
longer issues a warning with type="t" and all positive values for the
statistic.

\item
read.imagene() now uses mean background instead of median if the auto
segmentation method has been set, and readImaGeneHeader() now allows for
extraneous trailing white space in the lines of the data files.

\item
In read.maimages(), flush=TRUE has been added to all instances of
read.table(). This makes read.maimages() tolerant of spurious columns
added to the end of data lines.

\item
Improvements to lm.series(), gls.series() and rlm.series() to better
preserve column names.

\item
Improvement to normalizeQuantiles() to reduce memory requirements.

\item
decideTests() now preserves probe row names.

\item
Subsetting for MArrayLM objects has been revised so that the contrast
matrix is subsetted rather than the design matrix.

\item
makeContrasts() now gives an informative error message if the parameter
names are not syntactically valid.

\item
The previously deprecated function array2channel() has been removed.

\item
A new section on Fitted Model Objects has been added to the User's Guide,
together with a table of image analysis programs, an updated Bob case
study, and the CARMAweb reference. The warning that lmscFit() is
"experimental" has been removed.

\item
normalizeWithinArrays() was incorrect with method="composite" and
weights=NULL. Now fixed.

\item
Bug fixes to kooperberg(), plotlines(), plotMA3by2() and to
read.imagene(), which was trying to read too many data rows when there
were multiple ImaGene fields.

\item
Fix to read.maimages(): using readGenericHeader() to get column names was
not working when sep or quote were different from the defaults.
}}


\section{version 2.6.0 (2006-04-27)}{\itemize{

\item
New function arrayWeightsSimple(), with a new argument 'maxratio' which
stops the iteration when the ratio of largest to smallest weights becomes
large. arrayWeights() now returns a vector instead of a matrix, and is
faster and more accurate.

\item
New function asMatrixWeights(), which is used by all the linear
model functions so that they
can accept vectors of probe or array weights in place of a full matrix
of weights.

\item
New generic function avedups(), with default and MAList methods.

\item
New function plotFB().

\item
New option "rma" for the method argument of backgroundCorrect().

\item
All the normexp functions now use simpler arguments: 'x' instead of
'foreground' and 'background', and a vector argument 'par' instead of
separate arguments for each parameter. normexp.fit() now uses a
saddlepoint approximation to the normexp density and reverts to the
Nelder-Mead optimisation method, giving numerically reliable
computations. New functions dnormexp() and dnormexp.saddle().

\item
fitFDist() now treats very small and zero variances in a unified way, and
the bounding of standard deviations away from zero now uses a smaller
bound. This slightly changes the hyperparameter estimates and
significance results produced by eBayes().

\item
topTable() now includes separate columns for raw and adjusted p-values.

\item
New S3 methods fitted() and residuals() for the MArrayLM class, and new
S3 methods for 'dimnames<-' for RGList and MAList objects.

\item
Subsetting of RGList, MAList and MArrayLM objects by rownames now works
correctly even if not all components have rownames. When an MArrayLM
object is subsetted, the targets frame is no longer subsetted and the
design matrix subsets only on columns.

\item
duplicateCorrelation() and intraspotCorrelation() now return
atanh-correlations, and the corresponding output components have been
renamed to 'atanh.correlations'.

\item
The correlation argument of lmFit() no longer has a default value.

\item
Support for ScanArray Express has been added to read.maimages(), and
read.maimages() now reads BlueFuse v3.2 files correctly. The files
argument of read.maimages() can now be a data.frame, assumed to be the
targets frame. read.maimages() now uses check.names=FALSE to preserve
annotation column names as input, uses the X position column to determine
printer layout information for source="genepix", and no longer tries to
guess ndups for GenePix data, which was found to be unreliable for small
boutique arrays.

\item
readGenericHeader() replaces readBlueFuseHeader() and is now used by
read.maimages() to get the skip value and column names for generic source
files. readGenericHeader() now uses protectMetachar() so that special
characters can be used in column headings.

\item
Changes to the argument list and defaults for geneSetTest() to make it
more intuitive.

\item
decideTests() now supports adjust.method="BY".

\item
normalizeMedianDeviations() and normalizeMedians() have been replaced by
normalizeMedianAbsValues(), which is resistant to floating underflow.
The default method for normalizeBetweenArrays() has been changed to
"Aquantile".

\item
The deprecated functions dupcor.series() and
kooperberg_old() have been removed.

\item
readTargets() now checks that the labels are unique before making
row.names.

\item
MAList subsetting now supports the 'other' component, and cbind() and rbind()
now support the 'other' component in RGList and MAList objects.

\item
limma now formally depends on R version 2.2.0 or higher.

\item
A section on Paired Samples has been added to the User's Guide and the
help page for normalizeWithinArrays() has been edited and expanded. The
help pages for lm.series(), gls.series() and mrlm() now make clear that
they are utility functions called by lmFit(), and the help page for
topTable() better explains the differences between toptable() and
topTable(). The definition of interaction in the Weaver case study has
been corrected.

\item
Bug fix to decideTests(), which failed with method="nestedF" and
adjust.method="BH".

\item
Bug fix to ensure that makeContrasts() works correctly when the level
names happen to be the same as local variable names.

\item
Bug fixes to merge.RGList() and merge.MAList(), which used to break when
there were any components of the data objects other than intensity
matrices.

\item
Bug fixes to topTable() when some p-values are missing and to
as.data.frame.MArrayLM() when the row subset is of length 1.

\item
Bug fixes to read.maimages() for source="genepix" when there was only one
block of spots, for source="imagene" when a data.frame of file names was
confused with a targets file, and for source="genepix" when the meta-grid
layout could not be determined.

\item
kooperberg() now gives an informative error message if the input is not
an RGList.
}}


\section{version 2.2.0 (2005-10-13)}{\itemize{

\item
kooperberg() now accepts an RGList object instead of requiring each gpr
file to be read into a separate data.frame. The old function has been
renamed kooperberg_old().

\item
Extensive clean-up of the background correction functions based on the
normexp distribution. fit.normexp0() and m2loglik.normexp() have been
removed and the functions are now called normexp.fit(),
normexp.m2loglik(), normexp.grad() and normexp.signal().

\item
The deprecated functions designMatrix() and classifyTests() have been
removed.

\item
makeContrasts() now accepts expression- or character-valued variable
arguments.

\item
"fdr" has been changed to "BH" throughout as the p.adjust() method, and
the default adjust.method for topTable() has been changed to "BH".

\item
New options source="genepix.custom" and source="genepix.median" for
read.maimages(). The default value for the source argument is now NULL.
The columns argument can now be arbitrary, provided that the two
foreground channels are included, and the inclusion of background columns
is now optional.

\item
Support for BlueFuse files has been added to read.maimages(), with a new
function readBlueFuseHeader().

\item
readImaGeneHeader() has been entirely re-written and now supports ImaGene
6.0 and multiple sub-array fields.

\item
New method "control" for normalizeWithinArrays().

\item
normalizeWithinArrays(method="median") now uses a weighted median where
appropriate, and there is a new function weighted.median().
normalizeWithinArrays() now checks whether M or A are vectors, and gives
an informative error for method="printtiploess" if the layout is
incompatible with the data dimensions.

\item
arrayWeights() now includes a full "reml" option, and works correctly on
the chip level coefficients for PLMset objects.

\item
eBayes() now returns the F.p.value component as a vector rather than a
matrix.

\item
Point size 20 and a new path argument have been added to imageplot3by2().

\item
New introductory chapters 04.Background and 08.Tests for the online help
system, together with a general cleaning up of the chapter system.

\item
convest() now checks that the input p-values are between 0 and 1.

\item
Changes to write.fit() to make it tolerant of not finding expected
components such as Amean in the fit object.

\item
The User's Guide has been converted to LaTeX report style, with new
sections on getting help and on citing limma, an updated Weaver case
study, and the dendritic cell direct comparison example now using
modelMatrix(). A note has been added to the duplicateCorrelation() help
page about the minimum number of arrays needed for the function to return
statistically useful results.

\item
Bug fix to plotDensities(), where the background was being subtracted
twice.

\item
Work-around in eBayes() for the pf(q,df1,df2=Inf) bug in R 2.1.0.

\item
Bug fix to normexp.signal() to allow NAs in the signal.

\item
Bug fix to RGList subsetting, where the targets component was subsetted
incorrectly with j but missing i.

\item
Bug fix to as.data.frame.MArrayLM() when the fit contains only one
contrast and there are matrices stored as vectors.

\item
modelMatrix() now gives a helpful error message if 'targets' is missing.
}}


\section{version 1.9.6 (2005-05-18)}{\itemize{

\item
New function arrayWeights(), contributed by Matt Ritchie, which estimates
relative quality weights for each array in a multi-array experiment. The
previous version of the function has been renamed arrayWeightsQuick().

\item
New function convest() to estimate the proportion of true null hypotheses
from a vector of p-values.

\item
New function volcanoplot().

\item
New as.matrix() methods for all the limma data classes, and a new
as.data.frame() method for MArrayLM objects.

\item
New BFGS version of fit.normexp(), contributed by Jeremy Silver, with
several new utility functions to compute the likelihood and gradient.
The old Nelder-Mead version has been renamed fit.normexp0().
signal.normexp() is now more resistant to floating overflow and
underflow, and fit.normexp() now catches and reports optimization errors.

\item
lm.series() and lmFit() now treat a linear model with no weights and no
missing values as a special case, giving a great increase in speed for
typical Affymetrix analyses.

\item
decideTests() now allows for NA p-values for method="separate" or
"global", where it was previously giving wrong answers, and
summary.TestResults() also allows for NA p-values. decideTests() now
runs eBayes() if the MArrayLM object is found not to contain
t-statistics.

\item
The log and unlog arguments have been removed from MA.RG() and RG.MA().
RGList objects are now always assumed to contain unlogged intensities.
There are corresponding changes to plotDensities(), for which
log.transform=TRUE is now the default, and RG.MA() has a new argument
unlog.

\item
limma now officially requires R version 2.0 or higher.

\item
limmaUsersGuide() replaces limmaguideURL() and now outputs the file name
rather than a URL. The html version of the User's Guide has been
removed.

\item
The default column name for gene names in topTable() output is now
'ProbeID', and lmFit() now captures rownames if the data object is a
matrix.

\item
readTargets() has a new '...' argument and now looks for columns "Label"
or "FileName" to make row.names if the row.names argument is NULL.
removeExt() now coerces its argument to character.

\item
duplicateCorrelation() now traps the error condition when all the
non-missing observations for a gene belong to the same block.

\item
readImaGeneHeader() now allows for headers up to 100 lines, previously
60.

\item
The 'columns' argument no longer overwrites the 'source' argument of
read.maimages().

\item
New sections in the User's Guide on Time Course Experiments, Separate
Channel Analyses, Several Groups, Adjusting for Multiple Testing, and
software packages using limma. The User's Guide now uses BibTeX for its
reference list. The CITATION file has been updated to the Bioconductor
book chapter.

\item
Bug fixes to lmFit() when there are no weights or missing values and the
residual degrees of freedom are 0 or 1, and when there are non-estimable
coefficients.

\item
Bug fix to uniquegenelist(), which now preserves the data.frame or matrix
class. This fixes a bug when fitting a model to an exprSet object with
duplicate probes.

\item
Bug fix to normalizeRobustSpline() when Sigma0 is zero, that is, when the
method reduces to a global spline.

\item
Bug fix to normalizeBetweenArrays() with method="Gquantile" or
"Rquantile", which now changes only the A-values.

\item
Bug fix to plotDensities() when the object is an MAList and
log.transform=TRUE, which should give a log rather than a log-log scale.

\item
Bug fix so that lmFit() calls mrlm() rather than the deprecated
rlm.series().

\item
Bug fixes to changeLog(), to getLayout2() when the GAL file has no
header, to readGAL() for the stricter regular expression syntax of R
2.1.0, and to subsetting of MArrayLM objects when the components are
vectors rather than matrices.

\item
intraspotCorrelation() now uses try() to tolerate some fit errors from
remlscore().
}}


\section{version 1.8.5 (2004-10-25)}{\itemize{

\item
All the linear modelling functions now return the covariance matrix of
the estimated coefficients. This improves support for models with
covariance structures, such as random effects models for technical
replicates, mixed models for log-intensity analysis of two color data,
and duplicate spots. contrasts.fit() and classifyTestsF() now use the
covariance matrix stored in the data object.

\item
New function decideTests() for global multiple testing across probes and
contrasts.

\item
New function geneSetTest().

\item
New function nonEstimable(). lmFit() now reports non-estimable
coefficients, and contrasts.fit() stops with an error message if an
attempt is made to take a contrast of a non-estimable parameter.

\item
New background correction method "normexp" for backgroundCorrect(),
together with associated functions. backgroundCorrect() also has a new
argument 'offset', and MA.RG() and normalizeWithinArrays() now accept
arguments 'bc.method' and 'offset'.

\item
New functions imageplot3by2() and plotMA3by2() to write plots to files.
imageplot() now shows a z-range legend, has a new argument 'mar', and
traps infinite z values with an informative error message.

\item
plotMA() now does single-channel MA-plots and accepts exprSet, matrix and
MArrayLM objects. Points with zero weight are no longer plotted by
default.

\item
eBayes() now returns F-statistics and the corresponding p-values, and
checks that the design matrix is of full rank before computing the
F-statistic.

\item
New summary methods for RGList, MAList and MArrayLM objects, and a
slicker summary method for TestResults objects.

\item
New argument other.columns for read.maimages(). Subsetting of RGList
objects and MA.RG() now accommodate an 'other' component.

\item
read.maimages() now handles the new format SMD files, tries to set the
printer layout for Agilent data, and calls readGPRHeader() to get the
skip records for each gpr file. New functions readSMDHeader() and
getLayout2(), the latter reading GAL file headers.

\item
New utility functions protectMetachar(), trimWhiteSpace() and tZscore().

\item
New function changeLog() to display the limma change log, which has been
moved to the inst/doc directory of the package.

\item
lmFit() now accepts PLMset objects and better recognises when a
marrayNorm object has a non-empty list of gene names.

\item
rlm.series() is deprecated in favour of mrlm().

\item
kooperberg() now returns an RGList.

\item
heatDiagram() now chooses the color range more carefully, so that red is
associated only with positive fold changes and green with negative,
separates the primary up and down genes with white space, and has a new
argument 'mar'.

\item
The install.R file has been removed and replaced with LazyLoad in the
DESCRIPTION file. limma now passes R CMD check under R 2.0 and was for
the first time posted to CRAN as well as to Bioconductor.

\item
anova.MAList() now checks whether the design matrix has more than one
column, and write.fit() now works on MArrayLM objects when eBayes() has
not been used.

\item
New sections in the User's Guide on Installation and on Special Designs,
the Estrogen case study redone, and many edits throughout. A citation
entry has been added for limma. The Swirl and ApoAI case studies have
also been updated.

\item
Bug fixes to contrasts.fit() when there are non-estimable coefficients,
to decideTests() for method="hierarchical" and method="nestedF", to
kooperberg() when the background is zero, and to read.imagene().

\item
Bug fixes to read.maimages() for Agilent data, to readTargets(), and to
dim() for RGList, MAList and MArrayLM objects, which was returning NULL
when there was just one column.

\item
Bug fix to vennCounts() when there are no significant results for one or
more contrasts.

\item
normalizeWithinArrays(method="composite") now allows for NAs, and calls
to rlm() in normalizeRobustSpline() now explicitly set prior weights.
}}


\section{version 1.6.6 (2004-05-17)}{\itemize{

\item
New functions lmscFit() and intraspotCorrelation() for separate channel
analysis of two color data, together with a new documentation entry on
single channel analysis. lmscFit() was originally introduced as
lmsc.series() and now produces an MArrayLM object.

\item
New function write.fit() for outputting the results of a linear model
fit.

\item
New function modelMatrix() to construct design matrices for microarray
experiments with common references. modelMatrix() replaces
designMatrix(), which is deprecated. The use of 'ref' is now more
general and does not have to occur on every array.

\item
New function array2channel() to convert targets files from a by-array to
a by-channel representation.

\item
New function auROC() to compute the exact area under the ROC curve.

\item
New function poolVar() for Welch-style pooling of sample variances with
unequal variances.

\item
New function squeezeVar(), originally introduced as smoothVar().

\item
New function FStat(). classifyTests() has been renamed to
classifyTestsF(), and classifyTestsP() now extracts the degrees of
freedom from the data object.

\item
The S3 class "classification" has been replaced by the S4 class
"TestResults", including show and summary methods. classifyTestsF(),
classifyTestsT() and classifyTestsP() now produce TestResults objects.

\item
duplicateCorrelation() now uses randomizedBlockFit() from the statmod
package instead of gls() from the nlme package. This results in a
several-fold increase in speed for duplicateCorrelation() and
dupcor.series(). dupcor.series() is now deprecated and limma depends on
statmod.

\item
New argument 'block' for duplicateCorrelation(), gls.series() and lmFit()
to treat technical replicates as randomized blocks.

\item
New methods Rquantile, Gquantile and Tquantile for
normalizeBetweenArrays(), and a new method "vsn".

\item
New function getSpacing(), and getLayout() now has an option to guess the
spacing between duplicate spots. New function blockDiag() to form block
diagonal matrices. printorder() has been upgraded to accommodate spacing
by rows or topbottom.

\item
New function modifyWeights().

\item
New argument 'resort.by' for toptable() and topTable().

\item
New arguments 'mar' and 'cex' for vennDiagram(), which reduce white-space
and control text size.

\item
New argument 'annotation' for read.maimages(), and source="agilent" added
to read.maimages().

\item
readGAL() now sorts GAL files into standard order if they are not already
sorted, and read.maimages() now sets the gene list directly from the gpr
files, so readGAL() will usually no longer be required.

\item
New merge method for MAList objects, and new rbind methods for RGList and
MAList objects.

\item
lmFit() now extracts probe IDs from exprSet objects.

\item
vennCounts() now uses table() for the major computation and is therefore
much faster.

\item
backgroundCorrect(method="edwards") now handles NAs.

\item
A Sweave-style vignette has been added, as well as the User's Guide, and
new functions usersguideURL() and .First.lib() implement a limma entry in
the Windows drop-down Vignettes menu.

\item
Estrogen and Ecoli-lrp case studies have been added to the User's Guide,
and internal links have been added to the contents section.

\item
Substantial upgrade to the functions which identify and track control
spots from the gene list, including readTargets(), readSpotTypes(),
controlStatus() and plotMA(). New 'path' and '...' arguments for
readSpotTypes() and readTargets().

\item
Bug fixes to designMatrix(), which was giving the wrong matrix when the
common reference was Cy5 for the first array, and to printorder() with
spacing="topbottom", which was returning vectors of the wrong length.

\item
Bug fixes to imageplot() and heatdiagram(), which caused over-printing
when using par(mfrow).

\item
Bug fix to normalizeRobustSpline() to prevent an error when the spline
basis becomes singular.

\item
Changes to loessFit() to avoid errors with small numbers of observations.

\item
Bug fix to duplicateCorrelation(), which was not finding the M matrix
when the object was not an MArrayLM, and a check has been added for
correlations less than -1.

\item
Bug fixes to write.fit(), which was failing when adjusting p-values with
multiple contrasts, and further fixes to accommodate the TestResults
class.

\item
Refinements to the estimation of var.prior in ebayes() and related
functions.

\item
Bug fix to gls.series() to ensure that the spacing and weights arguments
are passed to duplicateCorrelation().

\item
Bug fixes to backgroundCorrect(method="edwards") and to plotPrintTipLoess().
}}


\section{version 1.2.8 (2003-10-30)}{\itemize{

\item
New function lmFit(), an object-orientated interface to linear modelling
which accepts RGList, MAList, marrayNorm and exprSet objects. New
MArrayLM class to hold the fitted model, and a new virtual class
LargeDataObject with an associated show method.

\item
Object-orientated versions of several functions are being gradually
introduced: eBayes() for ebayes(), topTable() for toptable() and
heatDiagram() for heatdiagram().

\item
New function makeContrasts() to construct contrast matrices.

\item
New function duplicateCorrelation() as an object-orientated wrapper for
dupcor.series().

\item
New function designMatrix() to construct design matrices for microarray
experiments with common references.

\item
New functions vennCounts() and vennDiagram() for Venn diagram plots, with
a new argument 'include'.

\item
New function loessFit(), giving substantial speed improvements for loess
normalization. loessFit() has a new argument 'bin' and the interpolation
settings have been fine tuned.

\item
New functions for Kooperberg Bayesian background correction, with the
internal functions now hidden from the user.

\item
New function plotDensities(), contributed by Natalie Thorne.

\item
New backgroundCorrect() methods "half", "minimum", "movingmin", "edwards"
and "none".

\item
read.maimages() now supports ImaGene, ArrayVision and SMD data files, and
reading of QuantArray files has been fixed. New functions
readGPRHeaders() and readImageneHeaders(), and read.imagene() now allows
for changes of header information between arrays.

\item
classifyTests() now supports MArrayLM objects, accepts design and
contrast matrices, and returns an F-statistic.

\item
New class exprSet2, and new function topTable().

\item
Support for correlated design matrices has been added to contrasts.fit().

\item
New cbind, rbind and subsetting methods for RGList and MAList objects,
and support for much extra structure in these objects.

\item
New argument 'ties' for normalizeBetweenArrays() and new argument 'mar'
for imageplot().

\item
normalizeWithinArrays(), normalizeForPrintorder() and plotPrintorder()
are no longer generic, which simplifies the documentation. The
documentation for normalizeBetweenArrays() has also been consolidated.

\item
New section in the User's Guide on single channel normalization,
contributed by Natalie Thorne, expanded sections on reading data, and a
new affy example.

\item
Bug fixes to classifyTests(), to loessFit(), to toptable(), which was
displaying incorrect A-values, to wtIgnore.Filter() with GenePix data,
and to backgroundCorrect() with method="minimum".
}}


\section{version 1.0.8 (2003-05-29)}{\itemize{

\item
limma was included in the Bioconductor project for the first time, as
part of Bioconductor release 1.2.

\item
The package was originally called smawehi and was made public for the
first time on 11 November 2002. It was renamed limma, with the title
"Linear Models for Microarray Data", in April 2003. The early versions
of the package were intended to be used in conjunction with the sma
package from the University of California at Berkeley.

\item
New S4 classes RGList and MAList for two color microarray data, including
setValidity() methods. Previously the data objects had the same structure
but were unclassed.

\item
New function MA.RG() to make MAList objects from RGList objects or lists,
and a new merge() method for RGList objects.

\item
New function backgroundCorrect().

\item
New S4 generic functions normalizeWithinArrays() and
normalizeBetweenArrays(), together with normalizeMedians(),
normalizeQuantiles(), normalizeMedianDeviations() and
normalizeForPrintorder(). Several of these functions were renamed from
earlier versions: normalize.quantile() to normalizeQuantiles(),
norm.scale.ma() to normalizeMedianDeviations() and normalize.printorder()
to normalizeForPrintorder().

\item
New function normalizeRobustSpline(), with a corresponding option for
normalizeWithinArrays(), which implements print-tip normalization by
fitting spline curves robustly to each print-tip group and then applying
empirical Bayes moderation across the array.

\item
New functions wtarea() and wtflags() to make quality weights for
microarray data from image analysis statistics, and a new function
wtIgnore.Filter() for QuantArray data.

\item
New functions classifyTests() and classifyTests.vector() for multiple
testing across genes and more than one contrast.

\item
New function plotMA() to make MA plots for both two color and
single-channel data, and a new function plotPrintTipLoess() to display
lowess normalization curves for printed two color microarrays.
plotPrintorder() was renamed from plot.printorder().

\item
Other functions renamed at this time include anova.series() to
anova.MAList(), fdist.fit() to fitFDist() and trigamma.inverse() to
trigammaInverse().

\item
New utility functions isNumeric() and helpMethods().

\item
A printable pdf version of the User's Guide is now provided.

\item
Bug fixes to ebayes() when there are more than two coefficients, to
read.maimages() for QuantArray data, and to classifyTests().
}}


\section{version 0.3 (2002-11-25)}{\itemize{

\item
The limma package was originally called the smawehi package, because
it was intended to be used in conjunction with the sma package from
the University of California at Berkeley, and was
made public via the repository http://bioinf.wehi.edu.au/smawehi.
The package originally focused on the analysis of two-color microarrays.

Version 0.3 of the package formed the basis of the One Day Course
"Statistics for Microarray Data Analysis" presented at the Australian Genstat
Conference in Perth, Australia, on 2 December 2002 by
Terry Speed, Gordon Smyth and Natalie Thorne.

\item
The package contained the following functions: a.spot, anova.series, bwss,
bwss.matrix, contrasts.fit, dl2norm, dnorml2norm, dupcor.series,
ebayes, fdist.fit, gls.series, gridc, gridr, heatdiagram, image.array,
lm.series, m.spot, matvec, norm.scale.ma, normalize.printorder,
normalize.quantile, plot.printorder, qqt, read.series, rg.genepix,
rg.quantarray, rg.series.spot, rg.spot, rlm.series, spotc, spotr,
txmixture.matrix, txmixture.vector, toptable, trigamma.inverse,
uniquegenelist, unwrapdups, vecmat.
}}
