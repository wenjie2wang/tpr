\name{tpr.test}
\alias{tpr.test}
\alias{sig.test.int.ff}
\alias{sig.test.boots.ff}
\alias{gof.test.int.ff}
\alias{gof.test.boots.ff}
\alias{gof.test.boots.pf}
\title{ Significance and Goodness-of-fit Test of TPR }
\description{
  Two kinds of tests are provided for inference on the coefficients in a
  fully functional TRP model: integral test and bootstrap test.
}
\usage{
sig.test.int.ff(fit, chypo = 0, idx, weight = TRUE, ncut = 2)
sig.test.boots.ff(fit, chypo = 0, idx, nsim = 1000, plot = FALSE)
gof.test.int.ff(fit, cfitList = NULL, idx, weight = TRUE, ncut = 2)
gof.test.boots.ff(fit, cfitList = NULL, idx, nsim = 1000, plot = FALSE)
gof.test.boots.pf(fit1, fit2, nsim, p = NULL, q = 1)
}
\arguments{
  \item{fit}{ a fitted object from \code{tpr} }
  \item{chypo}{ hypothesized value of coefficients }
  \item{idx}{ the index of the coefficients to be tested }
  \item{weight}{ whether or not use inverse variation weight }
  \item{ncut}{ the number of cuts of the interval of interest in
    integral test }
  \item{cfitList}{ a list of fitted object from \code{cst.fit.ff} }
  \item{nsim}{ the number of bootstrap samples in bootstrap test }
  \item{plot}{ whether or not plot }
  \item{fit1}{ fit of H0 model (reduced) }
  \item{fit2}{ fit of H1 model (full) }
  \item{p}{ the index of the time-varying estimation in fit2 }
  \item{q}{ the index of the time-independent estimation in fit1 }
}
\value{
  Test statistics and their p-values.
}
\references{
  Fine, Yan, and Kosorok (2004). Temporal Process Regression. Biometrika.
}
\author{ Jun Yan \email{jun.yan@uconn.edu} }
\seealso{
  \code{\link{tpr}}
}
\examples{
## see ?tpr
}
\keyword{htest}
\keyword{multivariate}
