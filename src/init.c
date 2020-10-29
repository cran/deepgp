#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>

/* .C function calls */
extern void alcGP_R(void *, void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *, void *, void *);
extern void inv_det_R(void *, void *, void *, void *);
extern void distance_R(void *, void *, void *, void *, void *, void *);
extern void distance_symm_R(void *, void *, void *, void *);
extern void Wij_R(void *, void *, void *, void *, void *, void *, void *,
                  void *, void *);

static const R_CMethodDef CEntries[] = {
  {"alcGP_R",         (DL_FUNC) &alcGP_R,        13},
  {"inv_det_R",       (DL_FUNC) &inv_det_R,       4},
  {"distance_R",      (DL_FUNC) &distance_R,      6},
  {"distance_symm_R", (DL_FUNC) &distance_symm_R, 4},
  {"Wij_R",           (DL_FUNC) &Wij_R,           9},
  {NULL, NULL, 0}
};


void R_init_deepgp(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

