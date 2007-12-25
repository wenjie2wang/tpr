\name{ci.plot}
\alias{ci.plot}
\title{ Confidence Interval Plot }
\description{
  Plotting time-varying coefficient with pointwise confidence.
}
\usage{
ci.plot(x, y, se, level = 0.95, ylim = NULL, newplot = TRUE, fun = gaussian()$linkinv, dfun = gaussian()$mu.eta, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the x coordinate }
  \item{y}{ the y coordinate }
  \item{se}{ the standard error of y }
  \item{level}{ confidence level }
  \item{ylim}{ the range of y axis }
  \item{newplot}{ if TRUE, draw a new plot }
  \item{fun}{ a transform function }
  \item{dfun}{ the derivative of the tranform function }
  \item{\dots}{ arguments to be passed to plot }
}
\author{Jun Yan \email{jyan@stat.uiowa.edu} }
\keyword{ hplot }
\keyword{ aplot }
