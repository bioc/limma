library(limma)

#	Tests for the contrasts argument of lmFit() and the lm.series/gls.series
#	behaviour: the sub-functions return a contrast-space fit in all cases (the
#	no-weights fast path applies contrasts.fit() directly, the genewise weighted
#	path applies the contrast exactly per gene), and method="robust" rejects contrasts.

set.seed(2024)
maxabs <- function(a,b) max(abs(as.matrix(a)-as.matrix(b)))

n <- 100L
narrays <- 8L
design <- cbind(Int=1, A=rep(c(0,1),each=4), B=rep(c(0,1),4))
y <- matrix(rnorm(n*narrays), n, narrays, dimnames=list(paste0("G",1:n),NULL))
contrast.matrix <- makeContrasts(A, B, AmB=A-B, levels=design)

#	1. Fast path (no weights): the sub-functions are self-contained and return a
#	contrast-space fit, equal to contrasts.fit() applied to the coefficient fit.
fastlm <- lm.series(y, design, contrasts=contrast.matrix)
stopifnot(!is.null(fastlm$contrasts))
stopifnot(ncol(fastlm$coefficients)==ncol(contrast.matrix))
reflm <- contrasts.fit(lm.series(y, design), contrast.matrix)
stopifnot(isTRUE(all.equal(fastlm$coefficients, reflm$coefficients, check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(fastlm$stdev.unscaled, reflm$stdev.unscaled, check.attributes=FALSE)))
fastgls <- gls.series(y, design, block=rep(1:4,each=2), correlation=0.3, contrasts=contrast.matrix)
stopifnot(!is.null(fastgls$contrasts))
stopifnot(ncol(fastgls$coefficients)==ncol(contrast.matrix))
refgls <- contrasts.fit(gls.series(y, design, block=rep(1:4,each=2), correlation=0.3), contrast.matrix)
stopifnot(isTRUE(all.equal(fastgls$coefficients, refgls$coefficients, check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(fastgls$stdev.unscaled, refgls$stdev.unscaled, check.attributes=FALSE)))

#	2. lmFit(contrasts=) on the no-weights path equals contrasts.fit(lmFit()) exactly.
f1 <- eBayes(lmFit(y, design, contrasts=contrast.matrix))
f2 <- eBayes(contrasts.fit(lmFit(y, design), contrast.matrix))
stopifnot(!is.null(f1$contrasts))
stopifnot(isTRUE(all.equal(f1$coefficients, f2$coefficients, check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(f1$stdev.unscaled, f2$stdev.unscaled, check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(f1$cov.coefficients, f2$cov.coefficients, check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(f1$t, f2$t, check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(f1$F, f2$F, check.attributes=FALSE)))

#	3. Weighted slow path applies the contrast in the C kernel (sets contrasts) and
#	matches the exact per-gene gold standard sqrt(diag(C' (X'WX)^-1 C)).
W <- matrix(runif(n*narrays, 0.3, 3), n, narrays)
fw <- lmFit(y, design, weights=W, contrasts=contrast.matrix)
stopifnot(!is.null(fw$contrasts))
stopifnot(ncol(fw$coefficients)==ncol(contrast.matrix))
gold.sd <- gold.coef <- matrix(NA_real_, n, ncol(contrast.matrix))
for(i in 1:n) {
	wi <- W[i,]
	M <- chol2inv(chol(crossprod(design*sqrt(wi))))
	b <- M %*% crossprod(design*wi, y[i,])
	gold.coef[i,] <- t(contrast.matrix) %*% b
	gold.sd[i,] <- sqrt(diag(t(contrast.matrix) %*% M %*% contrast.matrix))
}
stopifnot(maxabs(fw$coefficients, gold.coef) < 1e-10)
stopifnot(maxabs(fw$stdev.unscaled, gold.sd) < 1e-10)

#	4. Thread invariance on the weighted contrast path.
fw1 <- lmFit(y, design, weights=W, contrasts=contrast.matrix, nthreads=1L)
fw2 <- lmFit(y, design, weights=W, contrasts=contrast.matrix, nthreads=2L)
stopifnot(identical(fw1$coefficients, fw2$coefficients))
stopifnot(identical(fw1$stdev.unscaled, fw2$stdev.unscaled))

#	5. method="robust" rejects contrasts with an informative error.
err <- try(lmFit(y, design, method="robust", contrasts=contrast.matrix), silent=TRUE)
stopifnot(inherits(err, "try-error"))
stopifnot(grepl("robust", conditionMessage(attr(err,"condition"))))
#	robust without contrasts still works.
stopifnot(is(lmFit(y, design, method="robust"), "MArrayLM"))
