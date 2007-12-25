using namespace std;

#include "tntsupp.h"
#include "geese.h"
#include "utils.h"
#include "inter.h"
#include "lgtdl.h"
#include "famstr.h"


extern "C" {
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
}

/***************************************************************\
 Boostrapping significance and goodness-of-fit test for constant
 based on influence functions for covariate effect
\***************************************************************/
Vector<DVector>  bootsSample (DMatrix &infl, DVector &sig, int nsim) {
  int N = infl.num_rows(), T = infl.num_cols();
  int i, j;
  Vector<DVector> samp(nsim);
  DVector S0(T); S0 = 0.0; samp = S0;
  DVector G(N);
  for (i = 1; i <= nsim; i++) {
    for (j = 1; j <= N; j++) {
      G(j) = rnorm(0.0, 1.0) ;
    }
    DMatrix foo = SMult(G, infl);
    for (j = 1; j <= T; j++) {
      samp(i)(j) = sum(asVec(MatCol(foo, j))) / sig(j) / (double) N;
    }
  }
  return samp;
}


extern "C" {
  SEXP bootsSample_rap (SEXP infl_R, SEXP sig_R, SEXP nsim_R) {
    DMatrix infl = asDMatrix(infl_R);
    DVector sig = asDVector(sig_R);
    int nsim = INTEGER(nsim_R)[0];
    //    printf("nsim = %d\n", nsim);
    GetRNGstate();
    Vector<DVector> samp = bootsSample(infl, sig, nsim);
    PutRNGstate();
    //wrapp up
    SEXP ans = asSEXP(samp);
    return ans;
  }
}
