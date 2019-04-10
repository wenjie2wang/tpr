#ifndef LGTDL_H
#define LGTDL_H

#include "tntsupp.h"
#include "geese.h"
#include <R.h>
#include <Rdefines.h>

extern "C" {
extern void R_approx(double *, double *, int *, double *, int *,
		     int *, double *, double *, double *);

}



DVector approx(const DVector &x, const DVector &y, 
	       const DVector &xout, int method);

double approx(const DVector &x, const DVector &y, 
	      double t, int method);

class Lgtdl{
protected:
  DVector time_;
  DVector cov_;
public:
  Lgtdl(DVector &time, DVector &cov) : time_(time), cov_(cov) {}
  Lgtdl() : time_(0), cov_(0) {}
  ~Lgtdl() {}
  //  Lgtdl(DMatrix &dat);
  DVector time() {return time_;}
  DVector cov() {return cov_;}
  double interpprev(double time) {
    return approx(time_, cov_, time, 2);
  }
  double interplinear(double time) {
    return approx(time_, cov_, time, 1);
  }
  DVector interpprev(DVector &time) {
    return approx(time_, cov_, time, 2);
  }
  DVector interplinear(DVector &time) {
    return approx(time_, cov_, time, 1);
  }
};

typedef Vector<Lgtdl> VLgtdl;

Vector<Lgtdl> asVLgtdl(SEXP a);

DVector interpprev(double t, Vector<Lgtdl> &Yt);

Vector<VLgtdl> asVVLgtdl (SEXP a);

#endif //LGTDL_H
