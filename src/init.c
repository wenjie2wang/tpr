#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP bootsSample_rap(SEXP, SEXP, SEXP);
extern SEXP myinterp(SEXP, SEXP);
extern SEXP pfEst_rap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"bootsSample_rap", (DL_FUNC) &bootsSample_rap,  3},
    {"myinterp",        (DL_FUNC) &myinterp,         2},
    {"pfEst_rap",       (DL_FUNC) &pfEst_rap,       14},
    {NULL, NULL, 0}
};

void R_init_tpr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
