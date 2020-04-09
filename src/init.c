#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void linpack_choleski(void *, void *);
extern void linpack_inv_det(void *, void *, void *);
extern void logistffit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void logistpl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"linpack_choleski", (DL_FUNC) &linpack_choleski,  2},
    {"linpack_inv_det",  (DL_FUNC) &linpack_inv_det,   3},
    {"logistffit",       (DL_FUNC) &logistffit,       24},
    {"logistpl",         (DL_FUNC) &logistpl,         20},
    {NULL, NULL, 0}
};

void R_init_logistf(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}