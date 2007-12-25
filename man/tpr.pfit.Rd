\name{tpr.pfit}
\alias{cst.fit.ff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Constant fit of coefficients in a TPR model }
\description{
  Weighted least square estimate of a constant model for time-varying
  coefficients in a TPR model.
}
\usage{
cst.fit.ff(fit, idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{ a fitted object from \code{tpr} }
  \item{idx}{ the index of the }
}
\value{
  The estimated constant fit, standard error, z-value and p-value.
}
\references{
  Fine, Yan, and Kosorok (2004). Temporal Process Regression. Biometrika.
}
\author{ Jun Yan \email{jyan@stat.uiowa.edu} }
\seealso{
  \code{\link{tpr.test} }
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ robust }
