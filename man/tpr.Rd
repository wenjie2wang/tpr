\name{tpr}
\alias{tpr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Temporal Process Regression }
\description{
  Regression for temporal process responses and time-independent
  covariate. Some covariates have time-varying coefficients while others
  have time-independent coefficients.  
}
\usage{
tpr(y, delta, x, xtv=list(), z, ztv=list(), w, tis,
    family = poisson(),
    evstr = list(link = 5, v = 3),
    alpha = NULL, theta = NULL,
    tidx = 1:length(tis),
    kernstr = list(kern=1, poly=1, band=range(tis)/50),
    control = list(maxit=25, tol=0.0001, smooth=0, intsmooth=0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ Response, a list of "lgtdl" objects. }
  \item{delta}{ Data availability indicator, a list of "lgtdl" objects.}
  \item{x}{ Covariate matrix for time-varying coefficients.}
  \item{xtv}{ A list of list of "lgtdl" for time-varying covariates with
    time-varying coefficients.}
  \item{z}{ NOT READY YET; Covariate matrix for time-independent coefficients.}
  \item{ztv}{ NOT READY YET; A list of list of "lgtdl" for time-varying
    covariates with time-independent coefficients.}
  \item{w}{ Weight vector with the same length of \code{tis}. }
  \item{tis}{ A vector of time points at which the model is to be fitted.}
  \item{family}{ Specification of the response distribution; see
    \code{family} for \code{glm}; this argument is used in getting
    initial estimates.}
  \item{evstr}{ A list of two named components, link function and
    variance function.
    link: 1 = identity, 2 = logit, 3 = probit, 4 = cloglog, 5 = log;
    v: 1 = gaussian, 2 = binomial, 3 = poisson
  }
  \item{alpha}{ A matrix supplying initial values of alpha. }
  \item{theta}{ A numeric vector supplying initial values of theta. }
  \item{tidx}{ indices for time points used to get initial values. }
  \item{kernstr}{ A list of two names components:
    kern: 1 = Epanechnikov, 2 = triangular, 0 = uniform;
    band: bandwidth
  }
  \item{control}{ A list of named components:
    maxit: maximum number of iterations;
    tol: tolerance level of iterations.
    smooth: 1 = smoothing; 0 = no smoothing.
  }
}
\details{
  This rapper function can be made more user-friendly in the future. For
  example, \code{evstr} can be determined from the \code{family} argument.
}
\value{
  An object of class "tpr":
  \item{tis }{same as the input argument}
  \item{alpha }{estimate of time-varying coefficients}
  \item{beta }{estimate of time-independent coefficients}
  \item{valpha}{a matrix of variance of alpha at tis}
  \item{vbeta}{a matrix of variance of beta at tis}
  \item{niter}{the number of iterations used}
  \item{infAlpha}{a list of influence functions for alpha}
  \item{infBeta}{a matrix of influence functions for beta}
}
\references{
  Fine, Yan, and Kosorok (2004). Temporal Process
  Regression. Biometrika.

  Yan and Huang (2009). Partly Functional Temporal Process Regression
  with Semiparametric Profile Estimating Functions. Biometrics.
}
\author{ Jun Yan <jun.yan@uconn.edu> }
\keyword{ regression }
\keyword{ robust }
