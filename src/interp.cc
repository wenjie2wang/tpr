// extern "C" {
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
// }

#include "tntsupp.h"
#include "geese.h"
#include "utils.h"
#include "inter.h"
#include "lgtdl.h"
#include "famstr.h"

using namespace std;

/******************************************************************\
  faster interpprev
\******************************************************************/
extern "C" {
  SEXP myinterp (SEXP Y_R, SEXP Tis_R ) {
    Vector<Lgtdl> Y = asVLgtdl(Y_R);
    DVector Tis = asDVector(Tis_R);
    int N = Y.size();
    DVector Yt(N);
    for (int i = 1; i <= N; i++) {
      Yt(i) = Y(i).interpprev(Tis(i));
    }
    SEXP ans;
    ans = asSEXP(Yt);
    return ans;
  }
}
