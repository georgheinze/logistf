#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void linpack_choleski(void *, void *);
extern void linpack_inv_det(void *, void *, void *);
extern void logistffit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void logistffit_revised(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void logistffit_IRLS(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void logistplfit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"linpack_choleski", (DL_FUNC) &linpack_choleski,  2},
    {"linpack_inv_det",  (DL_FUNC) &linpack_inv_det,   3},
    {"logistffit",       (DL_FUNC) &logistffit,       25},
    {"logistffit_revised",       (DL_FUNC) &logistffit,       25},
    {"logistffit_IRLS",  (DL_FUNC) &logistffit_IRLS,  24},
    {"logistplfit",      (DL_FUNC) &logistplfit,      21},
    {NULL, NULL, 0}
};

void R_init_logistf(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}