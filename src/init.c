#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "multLPM.h"

static R_CallMethodDef callMethods[] = {
  {"loglikJointMult", (DL_FUNC) &loglikJointMult, 26},
  {"predcondY", (DL_FUNC) &predcondY, 18},
  {NULL, NULL, 0}
};



void R_init_multLPM(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
