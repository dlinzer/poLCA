#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void d2lldbeta2(void *, void *, void *, void *, void *, void *, void *, void *);
extern void postclass(void *, void *, void *, void *, void *, void *, void *, void *);
extern void probhat(void *, void *, void *, void *, void *, void *, void *);
extern void ylik(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"d2lldbeta2", (DL_FUNC) &d2lldbeta2, 8},
    {"postclass",  (DL_FUNC) &postclass,  8},
    {"probhat",    (DL_FUNC) &probhat,    7},
    {"ylik",       (DL_FUNC) &ylik,       7},
    {NULL, NULL, 0}
};

void R_init_poLCA(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
