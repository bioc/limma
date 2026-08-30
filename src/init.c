#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "limma.h"
#include "R_exports.h"

/* .C calls */
extern void fit_saddle_nelder_mead(void *, void *, void *, void *, void *, void *);
extern void normexp_gm2loglik(void *, void *, void *, void *, void *, void *);
extern void normexp_hm2loglik(void *, void *, void *, void *, void *, void *);
extern void normexp_m2loglik(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"fit_saddle_nelder_mead", (DL_FUNC) &fit_saddle_nelder_mead, 6},
    {"normexp_gm2loglik",      (DL_FUNC) &normexp_gm2loglik,      6},
    {"normexp_hm2loglik",      (DL_FUNC) &normexp_hm2loglik,      6},
    {"normexp_m2loglik",       (DL_FUNC) &normexp_m2loglik,       6},
    {NULL, NULL, 0}
};

/* .Call calls */
static const R_CallMethodDef CallEntries[] = {
    {"awremlfit", (DL_FUNC) &awremlfit, 9},
    {"dupcorfit", (DL_FUNC) &dupcorfit, 6},
    {"glsfit", (DL_FUNC) &glsfit, 6},
    {"lmfit",  (DL_FUNC) &lmfit,  5},
    {"poisfit", (DL_FUNC) &poisfit, 5},
    {"weighted_lowess", (DL_FUNC) &weighted_lowess, 6},
    {NULL, NULL, 0}
};

void R_init_limma(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
