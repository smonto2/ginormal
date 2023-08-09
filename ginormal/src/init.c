#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(pbdv)(double *v, double *x, double *dv, double *dp, double *pdf, double *pdd);

static const R_FortranMethodDef FortranEntries[] = {
  {"pbdv", (DL_FUNC) &F77_NAME(pbdv), 6},
  {NULL, NULL, 0}
};

void R_init_ginormal(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
