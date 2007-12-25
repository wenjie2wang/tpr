#include "tntsupp.h"

extern "C"{
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
}

#include "famstr.h"
#include "param.h"
#include "inter.h"
#include "lgtdl.h"

DVector approx(const DVector &x, const DVector &y, 
	       const DVector &xout, int method) {
  int nx = x.size(), lxout = xout.size();
  DVector ans(xout);
  double yleft = y(1), yright = y(nx), f = 0.0;
  R_approx(x.begin(), y.begin(), &nx, ans.begin(), &lxout,
	   &method, &yleft, &yright, &f);
  return ans;
}

double approx(const DVector &x, const DVector &y, 
	      double t, int method) {
  int nx = x.size(), lxout = 1;
  double ans = t;
  double yleft = y(1), yright = y(nx), f = 0.0;
  R_approx(x.begin(), y.begin(), &nx, &ans, &lxout,
	   &method, &yleft, &yright, &f);
  return ans;
}

Vector<Lgtdl> asVLgtdl(SEXP a) {
  int len = GET_LENGTH(a);
  Vector<Lgtdl> ans(len);
  SEXP dati;
  for (int i = 0; i < len; i++) {
    dati = VECTOR_ELT(a, i);
    DVector time = asDVector(VECTOR_ELT(dati, 0));
    DVector cov = asDVector(VECTOR_ELT(dati, 1));
    Lgtdl Lg(time, cov);
    ans(i + 1) = Lg;
  }
  return ans;
}

DVector interpprev(double t, Vector<Lgtdl> &Yt) {
  int n = Yt.size();
  DVector ans(n);
  for (int i = 1; i <= n; i++) 
    ans(i) = Yt(i).interpprev(t);
  return ans;
}


Vector<VLgtdl> asVVLgtdl (SEXP a) {
  int len = GET_LENGTH(a);
  Vector<VLgtdl> ans(len);
  SEXP ai;
  for (int i = 0; i < len; i++) {
    ai = VECTOR_ELT(a, i);
    ans(i + 1) = asVLgtdl(ai);
  }
  return ans;
}
