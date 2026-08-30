#ifndef LIMMA_R_EXPORT_H
#define LIMMA_R_EXPORT_H

#include <Rinternals.h>

SEXP dupcorfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP awremlfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP glsfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP poisfit(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP lmfit(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP weighted_lowess(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif
