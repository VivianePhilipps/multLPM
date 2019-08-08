#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "multLPM.h"

static R_CallMethodDef callMethods[] = {
  {"loglikUACV", (DL_FUNC) &loglikUACV, 26},
  {"loglikJointMult", (DL_FUNC) &loglikJointMult, 26},
  {"predcondY", (DL_FUNC) &predcondY, 18},
  {"probaEC", (DL_FUNC) &probaEC, 20},
  {"probaED", (DL_FUNC) &probaED, 20},
  {NULL, NULL, 0}
};



void R_init_multLPM(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
