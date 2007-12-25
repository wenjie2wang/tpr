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

// #define MAXIT 50
// #define EPSILON 0.0001

// global
int MAXIT, SMOOTH, INTERCEPTSMOOTH;
double EPSILON;
DVector DelAlpha; 
IVector UsedIt(1); // number of iterations used

// time-varying covariates
Vector<VLgtdl> Xtv, Ztv;
int pv, qv; // number of time-varying covariates in X and Z

DVector getTvCov(DMatrix &V, Vector<VLgtdl> &TV, int i, double t)
{
  DVector v = asVec ( MatRow (V, i) );
  int vd = v.dim(), tvd = TV.dim();
  if (tvd == 0) return v;
  for (int k = 1; k <= tvd; k++) {
    v(vd - tvd + k) = TV(k)(i).interpprev(t);
  }
  return v;
}


void getVtMat (double t, DMatrix &V, Vector<VLgtdl> &TV) 
{
  int n = V.num_rows(), vd = V.num_cols(), tvd = TV.dim();
  for (int k = 1; k <= tvd; k++) {
    int col = vd - tvd + k;
    for (int j = 1; j <= n; j++) {
      V(j, col) = TV(k)(j).interpprev(t);
    }
  }
}

/********************************************************************\
 Define class EVStr, expectation / variance structure
\********************************************************************/
class EVStr 
{
protected:
  Link link_;
  Variance var_;
public:
  // constructor
  EVStr(int link, int v) {
    link_ = Link(link);
    var_ = Variance(v);
  }
  // destructor
  ~EVStr() {}
  double linkfun ( double mu ) { return link_.linkfun ( mu ); }
  double linkinv ( double eta) { return link_.linkinv ( eta ); }
  double mu_eta ( double eta) { return link_.mu_eta (eta ); }
  double var (double mu) { return var_.v( mu ); }
}; // this semi-colon is important

EVStr asEVStr ( SEXP s )
{ 
  //s is a list of link and v
  int link = INTEGER(VECTOR_ELT(s, 0))[0];
  int v = INTEGER(VECTOR_ELT(s, 1))[0];
  EVStr str(link, v);
  return str;
}

/******************************************************************\
 kernel functions and kernel structure
\******************************************************************/
double EpanechnikovKernel ( double u )
{
  return ( fabs(u) <= 1) ? 0.75 * (1 - u * u) : 0.0;
}

double triangularKernel ( double u ) 
{
  return ( fabs(u) <= 1 ) ? 1.0 - fabs(u) : 0.0;
}
      

double uniformKernel ( double u ) 
{
  return ( fabs(u) <= 1 ) ? 1.0 : 0.0;
}

class KernStr {
protected:
  double band_;
  int kern_;
  int poly_;
public:
  // constructor
  KernStr(int kern, int poly, double band) {
    kern_ = kern;
    band_ = band;
    poly_ = poly;
  }
  // destructor
  ~KernStr() {}
  int poly () { return poly_; }
  double band () { return band_; }
  double kernfun ( double u ) { 
    double ku;
    switch ( kern_ ) {
    case 1:  
      return EpanechnikovKernel ( u );
    case 2:
      return triangularKernel ( u );
    default:
      return uniformKernel ( u );
    }
  }
}; // this semi-colon is important

KernStr asKernStr ( SEXP s ) 
{ 
  int kern = INTEGER(VECTOR_ELT(s, 0))[0];
  int poly = INTEGER(VECTOR_ELT(s, 1))[0];
  double band = REAL(VECTOR_ELT(s, 2))[0];
  KernStr str(kern, poly, band);
  return str;
}


/********************************************************************\ 
  General newton-raphson iteration for estimating equaiton G = 0
  with derivative H 
\********************************************************************/

double newton_raphson_1step ( DMatrix &H, DVector &G, DVector &theta ) 
{
  DVector Del = solve ( H, G );
  theta = theta + Del;
  double del = fmax ( fabs (Del) );
  return del;
}

/*********************************************************************\
  Model:
  $$ E[Y(t) | X(t), Z(t)] = g[ X(t)\alpha(t) + Z(t) \theta ] $$

\*********************************************************************/

/********************************************************************\
  Given $\alpha(t)$, EE for \theta is
  $$\sum_{i=1}^n \sum_{l=1}^L 
    \delta_i(t) \dot g [X(t)\alpha(t) + Z(t) \theta]
     { Y_i(t_{l + 1} - g[ X(t)\alpha(t) + Z(t) \theta ] } * Z_i = 0
\********************************************************************/

void prepTheta ( Vector<Lgtdl> &Y, 
		 Vector<Lgtdl> &Delta,
		 DMatrix &X, DMatrix &Z,
		 DVector &W,
		 DVector &Tis, 
		 Vector<DVector> &Alpha, DVector &theta,
		 EVStr &str,
		 //output
		 DMatrix &H, DVector &G ) 
{
  int T = Tis.size(), N = Y.size(), p = X.num_cols(), q = Z.num_cols();
  double t, delta, eta, mu, mu_eta, v, w, y;

  H = 0.0; G = 0.0;

  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= T - 1; j++) {
      t = Tis(j);
      delta = Delta(i).interpprev ( t ); //changed j to i
      if ( delta == 0 ) continue;
      DVector alpha = Alpha(j);

      DVector x = getTvCov ( X, Xtv, i, t);
      DVector z = getTvCov ( Z, Ztv, i, t);

      
      eta = dot_prod (x, alpha) + dot_prod (z, theta);
      mu = str.linkinv ( eta );
      v = str.var ( mu );
      mu_eta = str.mu_eta ( eta );
      w = W(j) / v * mu_eta;
      y = Y(i).interpprev ( t );

      G = G + ( w * ( y - mu ) ) * z;
      H = H + ( w * mu_eta ) * outerprod ( z );
    }
  }
}

/*********************************************************************\
  Given $\theta$, EE for $\alpha(t)$ is
  \sum_{i=1}^n \delta_i(t) \dot g [X(t)\alpha(t) + Z(t) \theta]
    { Y_i(t) - g [ X(t)\alpha(t) + Z(t) \theta ] } = 0
\*********************************************************************/

void prepAlpha_j ( Vector<Lgtdl> &Y, 
		   Vector<Lgtdl> &Delta,
		   DMatrix &X, 
		   DMatrix &Z, 
		   DVector &Tis, 
		   int j,
		   Vector<DVector> &Alpha, DVector &theta,
		   EVStr &str,
		   //output
		   DMatrix &H, DVector &G ) 
{
  int N = Y.size();
  double t, delta, eta, mu, mu_eta, v, w, y;
  
  H = 0.0; G = 0.0;
  t = Tis(j);
  DVector alpha = Alpha(j);

  for (int i = 1; i <= N; i++) {
    delta = Delta(i).interpprev(t);
    if (delta == 0) continue;

    DVector x = getTvCov ( X, Xtv, i, t);
    eta = dot_prod ( x, alpha );
    if (Z.num_cols() > 0) {
      DVector z = getTvCov ( Z, Ztv, i, t);
      eta = eta + dot_prod ( z, theta );
    }
    mu = str.linkinv ( eta );
    v = str.var ( mu );
    mu_eta = str.mu_eta ( eta );
    w = mu_eta / v;
    y = Y(i).interpprev ( t );
    
    G = G + ( w * ( y - mu ) ) * x;
    H = H + ( w * mu_eta ) * outerprod ( x );
  }
}


/******************************************************************\
 Newton-Raphson algorithm 1-step update
\******************************************************************/
double updateTheta ( Vector<Lgtdl> &Y, 
		     Vector<Lgtdl> &Delta,
		     DMatrix &X, DMatrix &Z,
		     DVector &W, //weight over time
		     DVector &Tis, 
		     Vector<DVector> &Alpha, DVector &theta,
		     EVStr &str ) 
{
  double del = 0.0;
  int q = Z.num_cols(); //assuming q > 0
  DMatrix H(q, q); DVector G(q);
  prepTheta ( Y, Delta, X, Z, W, Tis, Alpha, theta, str, H, G);
  //update theta
  del = newton_raphson_1step ( H, G, theta );
  return del;
}

DVector updateAlpha ( Vector<Lgtdl> &Y, 
		      Vector<Lgtdl> &Delta,
		      DMatrix &X, DMatrix &Z,
		      //DVector &W,
		      DVector &Tis, 
		      Vector<DVector> &Alpha, DVector &theta,
		      EVStr &str, KernStr &kern ) 
{
  int T = Tis.size(), N = Y.size();
  int p = X.num_cols(), q = Z.num_cols();
  DMatrix Ht(p, p); DVector Gt(p);
  DVector Del(T);

  //DVector Eta(N); Eta = 0.0;
  //if (q > 0) Eta = Z * theta;
  //assuming p > 0
  for ( int j = 1; j <= Tis.size(); j++ ) {
    Ht = 0.0; Gt = 0.0;
    prepAlpha_j (Y, Delta, X, Z, Tis, j, Alpha, theta, str, Ht, Gt);
  
    //update alpha_j
    Del(j) = newton_raphson_1step ( Ht, Gt, Alpha(j) );
  }
  return Del;
}

void pfEst ( Vector<Lgtdl> &Y, 
	     Vector<Lgtdl> &Delta,
	     DMatrix &X, DMatrix &Z,
	     DVector &W,
	     DVector &Tis, 
	     Vector<DVector> &Alpha, DVector &theta,
	     EVStr &str, KernStr &kern) 
{
  int p = X.num_cols(), q = Z.num_cols();
  double del_alpha = 0.0, del_theta = 0.0;
  for (int iter = 0; iter < MAXIT; iter++) {
    UsedIt(1) = iter;
    //update alpha(t), p >= 1
    DelAlpha = updateAlpha(Y, Delta, X, Z, Tis, Alpha, theta, str, kern);
    del_alpha = fmax(fabs(DelAlpha));
    //cerr << "del_alpha = " << del_alpha << endl;

    //update theta
    if (q > 0) {
      del_theta = updateTheta(Y, Delta, X, Z, W, Tis, Alpha, theta, str);
      //cerr << "del_theta = " << del_theta << endl;
    }
    
    if (del_theta < EPSILON && del_alpha < EPSILON) break;
  }
}

/********************************************************************\
 Local polynomial for time-varying coefficients
\********************************************************************/

DVector prepXAug_it(DMatrix &X, int i, double t, double khu, KernStr &kern) {
  int p = X.num_cols();
  DVector xi = getTvCov ( X, Xtv, i, t);
  /* x needs to be changed to accomodate local polynomial */
  int xdim = (INTERCEPTSMOOTH == 1) ? p * kern.poly() + p : p * kern.poly() + 1;
  DVector x ( xdim );
  if (INTERCEPTSMOOTH == 1) {
    for (int l = 0; l <= kern.poly(); l++) {
      Index1D I( l * p + 1, (l + 1) * p);
      VecSubs(x, I) = pow(khu, l) * xi;
    }
  } else {
    Index1D I(1, 1);
    VecSubs(x, I) = asVec ( VecSubs(xi, I) );
    if (p > 1) {
      Index1D J(2, p);
      DVector xinoint = asVec ( VecSubs ( xi, J ) );
      for (int l = 0; l <= kern.poly(); l++) {
	Index1D I( l * (p - 1) + 2, (l + 1) * (p - 1) + 1 );
	VecSubs(x, I) = pow(khu, l) * xinoint;	    
      }	  
    }	
  }
  return x;
}

void prepAlphaLP_j ( Vector<Lgtdl> &Y,
		     Vector<Lgtdl> &Delta,
		     DMatrix &X,
		     DMatrix &Z,
		     DVector &Tis,
		     int j, 
		     DVector &AlphaAug_j, DVector &theta,
		     EVStr &str,
		     KernStr &kern,
		     //output
		     DMatrix &H, DVector &G ) 
{
  int N = Y.size(), T = Tis.size();
  int p = X.num_cols();
  double t, delta, eta, mu, mu_eta, v, w, y;

  H = 0.0; G = 0.0;
  t = Tis(j);

  double band = kern.band(), kernwei;

  for (int i = 1; i <= N; i++) {
    for (int k = 1; k <= T; k++) {
      double u = Tis(k);
      delta = Delta(i).interpprev(u);
      if (delta == 0) continue;

      double khu = (u - t) / band;
      kernwei = kern.kernfun( khu ) / band;
      if (kernwei == 0.0) continue;

      /* construct augumented x */
      DVector x = prepXAug_it ( X, i, u, khu, kern );

      eta = dot_prod ( x, AlphaAug_j );
      if (Z.num_cols() > 0) {
	DVector z = getTvCov ( Z, Ztv, i, u);
	eta = eta + dot_prod ( z, theta );
      }
      mu = str.linkinv ( eta );
      v = str.var ( mu );
      mu_eta = str.mu_eta ( eta );
      w = mu_eta / v;
      y = Y(i).interpprev ( u );
      
      G = G + ( kernwei * w * ( y - mu ) ) * x;
      H = H + ( kernwei * w * mu_eta ) * outerprod ( x );
    }
  }
}


DVector updateAlphaLP ( Vector<Lgtdl> &Y, 
			Vector<Lgtdl> &Delta,
			DMatrix &X, DMatrix &Z,
			DVector &Tis, 
			Vector<DVector> &AlphaAug, 
			DVector &theta,
			EVStr &str, KernStr &kern ) 
{
  int T = Tis.size(), N = Y.size();
  int p = X.num_cols(), q = Z.num_cols();
  int pAug = AlphaAug(1).size();
  DMatrix Ht(pAug, pAug); DVector Gt(pAug);
  DVector Del(T); 
  int k;

  UsedIt(1) = 0;
  for ( int j = 1; j <= Tis.size(); j++ ) {
    Ht = 0.0; Gt = 0.0;

    DVector AlphaAug_j = AlphaAug(j);
    //if (j == 1) cerr << "j = " << j << endl << Alpha(j) << Alpha_j;
    
    prepAlphaLP_j (Y, Delta, X, Z, Tis, j, AlphaAug_j, theta, str, kern, Ht, Gt);
      
    //update Alpha_j
    Del(j) = newton_raphson_1step ( Ht, Gt, AlphaAug_j );
  }
  return Del;
}


Vector<DVector> prepAlphaAugInit ( Vector<DVector> &Alpha, 
				   DVector &Tis,
				   KernStr &kern )
{
  int p = Alpha(1).size();
  int xdim = (INTERCEPTSMOOTH == 1) ? p * kern.poly() + p : p * kern.poly() + 1;
  DVector a0(xdim); a0 = 0.0;
  int T = Tis.size();
  Vector<DVector> AlphaAug(T); AlphaAug = a0;

  for (int j = 1; j <= T; j++) {
    DVector Alpha_j(xdim);    Alpha_j = 0.0;
    Index1D I(1, p);
    //get initial value
    DVector alpha(p); alpha = 0.0;
    double wei = 0.0;
    double t = Tis(j);
    for ( int l = 1; l <= T; l++) {
      double u = Tis(l);
      double band = kern.band();
      double khu = (u - t) / band;
      double kernwei = kern.kernfun( khu ) / band;
      if (kernwei == 0.0) continue;
      alpha = alpha + kernwei * Alpha(l);
      //if (j == 1) cerr << "l = " << l << "\n alpha" << alpha;
      wei = wei + kernwei;
    }
    VecSubs(Alpha_j, I) = (1.0 / wei) * alpha;
    AlphaAug(j) = Alpha_j;
  }
  return AlphaAug;
}

void pfEstLP ( Vector<Lgtdl> &Y, 
	       Vector<Lgtdl> &Delta,
	       DMatrix &X, DMatrix &Z,
	       DVector &W,
	       DVector &Tis, 
	       Vector<DVector> &Alpha, DVector &theta,
	       EVStr &str, KernStr &kern ) 
{
  int p = X.num_cols(), q = Z.num_cols();
  double del_alpha = 0.0, del_theta = 0.0;
  
  Vector<DVector> AlphaAug = prepAlphaAugInit ( Alpha, Tis, kern );
  for (int iter = 0; iter < MAXIT; iter++) {
    UsedIt(1) = iter;

    DelAlpha = updateAlphaLP(Y, Delta, X, Z, Tis, AlphaAug, theta, str, kern);
    del_alpha = fmax(fabs(DelAlpha));

    // update Alpha from AlphaAug
    Index1D I(1, p);
    for (int j = 1; j <= Tis.size(); j++) {
      Alpha(j) = asVec ( VecSubs(AlphaAug(j), I) );      
    }

    if (q > 0) {
      del_theta = updateTheta(Y, Delta, X, Z, W, Tis, Alpha, theta, str);
    }

    if (del_theta < EPSILON && del_alpha < EPSILON) break;
  }
}


/**************************************************************\
 quantities needed for variance estimator 
\**************************************************************/
DVector getG2v_t (double t,
		  DVector &delta, 
		  DMatrix &X, DVector &alpha, 
		  DMatrix &Z, DVector &theta,
		  EVStr &str) {
  int N = X.num_rows();
  double eta, mu, v, mu_eta;

  DVector g2v(N); g2v = 0.0;
  for (int i = 1; i <= N; i++) {
    if (delta(i) == 0) continue;
    DVector x = getTvCov ( X, Xtv, i, t);
    eta = dot_prod ( x, alpha );// + Eta(i);
    if (Z.num_cols() > 0) {
      DVector z = getTvCov ( Z, Ztv, i, t);
      eta = eta + dot_prod(z, theta);
    }
    mu = str.linkinv ( eta );
    v = str.var ( mu );
    mu_eta = str.mu_eta ( eta );
    g2v(i) = mu_eta * mu_eta / v;
  }
  return g2v;
}

DMatrix getG2vxx_t (double t, 
		    DVector &delta, DMatrix &X, DVector &g2v) {
  int N = delta.size(), p = X.num_cols();

  DMatrix g2vxx(p, p); g2vxx = 0.0;
  for (int i = 1; i <= N; i++) {
    if (delta(i) == 0) continue;
    DVector x = getTvCov ( X, Xtv, i, t);
    g2vxx = g2vxx + g2v(i) * outerprod(x);
  }
  g2vxx = (1 / (double) N) * g2vxx;
  return g2vxx;
}

DMatrix getG2vxz_t (double t,
		    DVector &delta, DMatrix &X, DMatrix &Z, DVector &g2v) {
  int N = delta.size(), p = X.num_cols(), q = Z.num_cols();

  DMatrix g2vxz(p, q);
  for (int i = 1; i <= N; i++) {
    if (delta(i) == 0) continue;
    DVector x = getTvCov ( X, Xtv, i, t);
    //assuming q > 0
    DVector z = getTvCov ( Z, Ztv, i, t);
    g2vxz = g2vxz + g2v(i) * outerprod(x, z);
  }
  g2vxz = (1 / (double) N) * g2vxz;
  return g2vxz;
}

DVector getM_t ( double t,
		 DVector &delta,
		 DVector &y, 
		 DMatrix &X, DVector &alpha, 
		 DMatrix &Z, DVector &theta,
		 EVStr &str ) {
  int N = delta.size();
  double eta, mu, v, mu_eta, gv;
  DVector M(N); M = 0.0;
  for (int i = 1; i <= N; i++) {
    if (delta(i) == 0) continue;
    DVector x = getTvCov ( X, Xtv, i, t);
    eta = dot_prod ( x, alpha );// + Eta(i);
    if (Z.num_cols() > 0) {
      DVector z = getTvCov ( Z, Ztv, i, t);
      eta = eta + dot_prod(z, theta);
    }
    mu = str.linkinv ( eta );
    v = str.var ( mu );
    mu_eta = str.mu_eta ( eta );
    gv = mu_eta / v;
    M(i) = gv * (y(i) - mu);
  }
  return M;
}

void prepInfl (Vector<Lgtdl> &Y,
	       Vector<Lgtdl> &Delta,
	       DMatrix &X, DMatrix &Z,
	       DVector &Tis,
	       EVStr &str,
	       //parameters
	       Vector<DVector> &Alpha, DVector &theta,
	       //output
	       Vector<DVector> &g2v,
	       Vector<DMatrix> &g2vxx, // p by p
	       Vector<DMatrix> &g2vxz, // p by q
	       Vector<DVector> &M,    //N by 1
	       Vector<DMatrix> &ZbarT //p by q
	       ) {
  int T = Tis.size();
  int p = X.num_cols(), q = Z.num_cols();

  for (int i = 1; i <= T; i++) {
    DVector alpha_t = Alpha(i);
    DVector delta_t = interpprev(Tis(i), Delta);
    DVector y_t = interpprev(Tis(i), Y);
    g2v(i) = getG2v_t ( Tis(i), delta_t, X, alpha_t, Z, theta, str);
    g2vxx(i) = getG2vxx_t ( Tis(i), delta_t, X, g2v(i) );
    M(i) = getM_t( Tis(i), delta_t, y_t, X, alpha_t, Z, theta, str);
    if (q > 0) {
      g2vxz(i) = getG2vxz_t ( Tis(i), delta_t, X, Z, g2v(i) );
      ZbarT(i) = solve(g2vxx(i), g2vxz(i));
    }
    //cerr << "M(" << t << ")=" << M(t);
  }
}

/**************************************************************\
 influence functions
\**************************************************************/

// returns a q x N matrix
DMatrix getInflTheta (DVector &Tis, DVector &W,
		      DMatrix &X, DMatrix &Z,
		      //intermediate quantities
		      Vector<DVector> &g2v,
		      Vector<DVector> &M,
		      Vector<DMatrix> &ZbarT ) {
  int T = ZbarT.size(), N = X.num_rows(), q = Z.num_cols(), p = X.num_cols();

  DMatrix A(q, q), inflB(N, q); 
  A = 0.0; inflB = 0.0;

  for (int t = 1; t <= T; t++) {
    //cerr << "t = " << t << "@@ZbatT(t) is " << ZbarT(t);
    // wrong: A = A + Transpose_View<DMatrix>(Z) * ( g2v(t) * W(t) * X * ZbarT(t) );
    if (pv > 0) getVtMat(Tis(t), X, Xtv); //reference passing
    if (qv > 0) getVtMat(Tis(t), Z, Ztv);

    DMatrix ZtT = transpose(Z);
    A = A + ZtT * SMult(W(t) * g2v(t), Z - X * ZbarT(t) );
    //cerr << " M(t) = " <<  M(t);
    inflB = inflB + SMult(W(t) * M(t), Z - X * ZbarT(t) ); // N by q
  }
  A = (1 / (double) N) * A;
  DMatrix inflB_trans = transpose ( inflB ); // q by N
  DMatrix inflTheta = solve(A, inflB_trans); // q by N
  return inflTheta;
}

//returns a vector of p by N matrix
Vector<DMatrix> getInflAlpha (DVector &Tis, DMatrix &X, DMatrix &Z,
			      Vector<DMatrix> &g2vxx,
			      Vector<DVector> &M,
			      Vector<DMatrix> &ZbarT,
			      DMatrix &inflTheta // N by q
			      ) {
  int T = ZbarT.size(), p = X.num_cols(), 
    q = inflTheta.num_cols(), N = X.num_rows();
  Vector<DMatrix> inflAlpha(T);
  DMatrix i0(p, N); i0 = 0.0; inflAlpha = i0;

  for (int t = 1; t <= T; t++) {
    if (pv > 0) getVtMat(Tis(t), X, Xtv); //reference passing
    if (qv > 0) getVtMat(Tis(t), Z, Ztv);

    DMatrix XM_t = SMult ( M(t), X ); // N by p
    DMatrix XM_t_trans = transpose(XM_t);
    inflAlpha(t) =  solve(g2vxx(t), XM_t_trans);
    if (q > 0) inflAlpha(t) = inflAlpha(t) - ZbarT(t) * inflTheta;
  }
  return inflAlpha; // p by N
}

/**************************************************************\
 variance estimator; for unsmoothed version only now.
\**************************************************************/

void pfVar(Vector<Lgtdl> &Y,
	   Vector<Lgtdl> &Delta,
	   DMatrix &X, DMatrix &Z, DVector &W,
	   DVector &Tis,
	   EVStr &str,
	   //parameters
	   Vector<DVector> &Alpha,
	   DVector &Theta,
	   //output
	   Vector<DMatrix> &VarAlpha,
	   DMatrix &VarTheta,
	   Vector<DMatrix> &inflAlpha,
	   DMatrix &inflTheta) {
  int T = Tis.size(), N = Y.size(), p = X.num_cols(), q = Z.num_cols();
  //DVector Eta = Z * Theta;
  //preparation
  Vector<DVector> g2v(T), M(T);
  Vector<DMatrix> g2vxx(T), g2vxz(T), ZbarT(T);
  DVector M0(N); M0 = 0.0; g2v = M0; M = M0;
  DMatrix XZ0(p, q); XZ0 = 0.0; ZbarT = XZ0; g2vxz = XZ0;
  DMatrix XX0(p, p); XX0 = 0.0; g2vxx = XX0;
  
  prepInfl (Y, Delta, X, Z, Tis, str, Alpha, Theta,
	    g2v, g2vxx, g2vxz, M, ZbarT);
  //influence functions and variances
  if (q > 0) {
    inflTheta = getInflTheta(Tis, W, X, Z, g2v, M, ZbarT);
    DMatrix inflThetaT = transpose(inflTheta);
    VarTheta = (1.0 / N / N ) * (inflTheta * inflThetaT);
  }
  
  inflAlpha = getInflAlpha(Tis, X, Z, g2vxx, M, ZbarT, inflTheta);
  for (int t = 1; t <= T; t++) {
    DMatrix inflAlpha_tT = transpose(inflAlpha(t));
    VarAlpha(t) = (1.0 / N / N ) * (inflAlpha(t) * inflAlpha_tT );
  }
}


/**************************************************************\
 R Wrappers
\**************************************************************/
extern "C" {
  SEXP pfEst_rap (SEXP Y_R, SEXP Delta_R, 
		  SEXP X_R, SEXP Xtv_R,
		  SEXP Z_R, SEXP Ztv_R,
		  SEXP W_R, SEXP Tis_R, 
		  SEXP Alpha_R, 
		  SEXP Theta_R,
		  SEXP Err_R,
		  SEXP str_R,
		  SEXP kern_R,
		  SEXP contr_R) {
    Vector<Lgtdl> Y = asVLgtdl(Y_R);
    Vector<Lgtdl> Delta = asVLgtdl(Delta_R);
    DMatrix X = asDMatrix(X_R), Z = asDMatrix(Z_R);

    pv = GET_LENGTH(Xtv_R); qv = GET_LENGTH(Ztv_R);
    if (pv > 0) Xtv = asVVLgtdl(Xtv_R); 
    if (qv > 0) Ztv = asVVLgtdl(Ztv_R);
    //cerr << "pv = " << pv << " qv = " << qv << endl;
   

    DVector W = asDVector(W_R), Tis = asDVector(Tis_R);
    Vector<DVector> Alpha = asVDVector(Alpha_R);
    DVector Theta = asDVector(Theta_R);
    DVector Err = asDVector(Err_R);
    EVStr str = asEVStr(str_R);
    KernStr kern = asKernStr(kern_R);

    DelAlpha = Err;

    MAXIT = INTEGER(VECTOR_ELT(contr_R, 0))[0];
    EPSILON =  REAL(VECTOR_ELT(contr_R, 1))[0];
    SMOOTH = INTEGER(VECTOR_ELT(contr_R, 2))[0];
    INTERCEPTSMOOTH = INTEGER(VECTOR_ELT(contr_R, 3))[0];

    //estimate
    if (SMOOTH == 0) pfEst(Y, Delta, X, Z, W, Tis, Alpha, Theta, str, kern);
    else pfEstLP(Y, Delta, X, Z, W, Tis, Alpha, Theta, str, kern);
    //variance
    int T = Tis.size(), N = Y.size(), p = X.num_cols(), q = Z.num_cols();
    Vector<DMatrix> VarAlpha(T); 
    DMatrix VA0(p, p); VA0 = 0.0; VarAlpha = VA0;
    DMatrix VarTheta(q, q); VarTheta = 0.0;
    //influence functions
    Vector<DMatrix> InflAlpha(T);
    DMatrix i0(p, N); i0 = 0.0; InflAlpha = i0;
    DMatrix InflTheta(q, N); InflTheta = 0.0;

    //cerr << "parfunVar start here\n";
    pfVar(Y, Delta, X, Z, W, Tis, str, Alpha, Theta, 
	  VarAlpha, VarTheta,
	  InflAlpha, InflTheta);

    //cout << InflAlpha(2) << VarAlpha(2);
    //wrapp up
    SEXP ans;
    PROTECT ( ans = NEW_LIST(8) );
    SET_VECTOR_ELT ( ans, 0, asSEXP(Alpha) );
    SET_VECTOR_ELT ( ans, 1, asSEXP(Theta) );
    SET_VECTOR_ELT ( ans, 2, asSEXP(VarAlpha) );
    SET_VECTOR_ELT ( ans, 3, asSEXP(VarTheta) );
    SET_VECTOR_ELT ( ans, 4, asSEXP(UsedIt) );
    SET_VECTOR_ELT ( ans, 5, asSEXP(InflAlpha) );
    SET_VECTOR_ELT ( ans, 6, asSEXP(InflTheta) );
    SET_VECTOR_ELT ( ans, 7, asSEXP(DelAlpha) );
    UNPROTECT(1);
    return ans;
  }
}
