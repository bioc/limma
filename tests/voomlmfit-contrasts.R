library(limma)

#	Tests for the contrasts argument of voomLmFit(). Because voom precision
#	weights are probe-specific, the contrast is applied in the final weighted fit
#	and the contrast stdev.unscaled is exact per gene -- equal to a direct weighted
#	refit and different from the post-hoc contrasts.fit() approximation.

set.seed(2024)
maxabs <- function(a,b) max(abs(as.matrix(a)-as.matrix(b)))

ngenes <- 200L
group <- factor(rep(1:3, each=4))
#	Treatment coding gives correlated coefficients, so a multi-coefficient contrast
#	exercises the per-gene-vs-shared correlation difference under voom weights.
design <- model.matrix(~group)
colnames(design) <- c("Int","grp2","grp3")
#	Counts with gene-specific baselines
counts <- matrix(rpois(ngenes*nrow(design), lambda=rep(sample(20:200,ngenes,replace=TRUE),nrow(design))), ngenes, nrow(design))
rownames(counts) <- paste0("G",1:ngenes)
contrast.matrix <- makeContrasts(grp2, grp3-grp2, levels=design)   # 2nd contrast mixes coefficients

fp <- suppressMessages(voomLmFit(counts, design))                       # coefficient space
fc <- suppressMessages(voomLmFit(counts, design, contrasts=contrast.matrix))  # contrast space

#	1. Contrast-space output of the correct shape.
stopifnot(!is.null(fc$contrasts))
stopifnot(ncol(fc$coefficients)==ncol(contrast.matrix))
stopifnot(identical(colnames(fc$coefficients), colnames(contrast.matrix)))
stopifnot(identical(dim(fc$stdev.unscaled), dim(fc$coefficients)))
#	df.residual and sigma are gene-level and unchanged by the contrast.
stopifnot(isTRUE(all.equal(fc$sigma, fp$sigma, check.attributes=FALSE)))
stopifnot(identical(fc$df.residual, fp$df.residual))

#	2. Coefficients equal the post-hoc route exactly (contrasts are linear).
fpost <- contrasts.fit(fp, contrast.matrix)
stopifnot(maxabs(fc$coefficients, fpost$coefficients) < 1e-9)

#	3. stdev.unscaled is exact: equals the final weighted fit with contrasts, and a
#	direct per-gene gold standard sqrt(diag(C' (X'W_gX)^-1 C)).
W <- fp$EList$weights
ydat <- fp$EList$E
fref <- lmFit(ydat, design, weights=W, contrasts=contrast.matrix)
stopifnot(maxabs(fc$stdev.unscaled, fref$stdev.unscaled) < 1e-10)
gold <- matrix(NA_real_, ngenes, ncol(contrast.matrix))
for(i in 1:ngenes) {
	M <- chol2inv(chol(crossprod(design*sqrt(W[i,]))))
	gold[i,] <- sqrt(diag(t(contrast.matrix) %*% M %*% contrast.matrix))
}
stopifnot(maxabs(fc$stdev.unscaled, gold) < 1e-10)

#	4. The exact result genuinely differs from the post-hoc approximation.
stopifnot(maxabs(fc$stdev.unscaled, fpost$stdev.unscaled) > 1e-6)

#	5. eBayes/topTable run on the contrast fit.
eb <- eBayes(fc)
tt <- topTable(eb, coef=1, n=5)
stopifnot(nrow(tt)==5L)

#	6. Blocked path returns a contrast-space fit without error.
fb <- suppressMessages(voomLmFit(counts, design, block=rep(1:6,each=2), contrasts=contrast.matrix))
stopifnot(!is.null(fb$contrasts), ncol(fb$coefficients)==ncol(contrast.matrix))
