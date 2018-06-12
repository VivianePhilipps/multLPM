#include <R_ext/RS.h>


SEXP loglikJointMult(SEXP b0,
		     SEXP bfix0,
		     SEXP fix0,
		     SEXP Y0,
		     SEXP X0,
		     SEXP Xd0,
		     SEXP idD0,
		     SEXP D0,
		     SEXP Xseuil0,
		     SEXP nmes0,
		     SEXP nv0,
		     SEXP idx0,
		     SEXP idiag0,
		     SEXP ncor0,
		     SEXP ny0,
		     SEXP nalea0,
		     SEXP ntr0,
		     SEXP link0,
		     SEXP nodes0,
		     SEXP epsY0,
		     SEXP nRE0,
		     SEXP nBM0,
		     SEXP nbevt0,
		     SEXP idcause0,
		     SEXP entreRetard0,
		     SEXP nbTronc0);


SEXP precondY(SEXP b0,
	      SEXP bfix0,
	      SEXP fix0,
	      SEXP Y0,
	      SEXP X0,
	      SEXP nmes0,
	      SEXP nv0,
	      SEXP idx0,
	      SEXP idiag0,
	      SEXP ncor0,
	      SEXP ny0,
	      SEXP nalea0,
	      SEXP ntr0,
	      SEXP link0,
	      SEXP nodes,
	      SEXP nRE0,
	      SEXP nBM0,
	      SEXP chol0);
