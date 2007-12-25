\name{dnase}
\alias{dnase}
\docType{data}
\title{ rhDNase Data }
\description{
  Randomized trial of rhDNase for treatment of cystic fibrosis
}
\usage{data(dnase)}
\format{
  A data frame with 767 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{rx}}{treatment arm: 0 = placebo, 1 = rhDNase}
    \item{\code{fev}}{forced expirotary volume, a measure of lung capacity}
    \item{\code{futime}}{follow time}
    \item{\code{iv1}}{IV start time}
    \item{\code{iv2}}{IV stop time}
  }
}
\details{
  During an exacerbation, patients received intravenous (IV) 
  antibiotics and were considered unsusceptible until seven
  exacerbation-free days beyond the end of IV therapy.

  A few subjects were infected at the time of enrollment,
  for instance a subject has a first infection interval of
  -21 to 7. We do not count this first infection as an "event",
  and the subject first enters the risk set at day 7.
}
\source{
  Therneau and Grambsch (2000). Modeling Survival Data: Extending the
  Cox model. Springer. 
  http://www.mayo.edu/hsr/people/therneau/book/data/dnase.html
}
\references{
  Yan and Fine (2007). Analysis of Episodic Data with Application to
  Recurrent Pulmonary Exacerbations in Cystic Fibrosis Patients. JASA.
}
\examples{
## This example steps through how to set up for the tpr function.
## Three objects are needed:
##     1) response process (an object of "lgtdl")
##     2) data availability process (an object of "lgtdl")
##     3) a time-independent covariate matrix

data(dnase)

## extracting the unique id and subject level information
dat <- unique(dnase[,c("id", "futime", "fev", "rx")])

## construct temporal process response for recurrent enent
rec <- lapply(split(dnase[,c("id", "iv1", "futime")], dnase$id),
              function(x) {
                v <- x$iv1
                maxfu <- max(x$futime)
                ## iv1 may be negative!!!
                if (is.na(v[1])) c(0, maxfu + 1.0)
                else if (v[1] < 0) c(v[1] - 1, v[!is.na(v)], maxfu + 1.0)
                else c(0, v[!is.na(v)], maxfu + 1.0)
              })

yrec <- lapply(rec,
               function(x) {
                 dat <- data.frame(time=x, cov=1:length(x)-1)
                 len <- length(x)
                 dat$cov[len] <- dat$cov[len - 1]
                 as.lgtdl(dat)
               })

## construct temporal process response for accumulative days exacerbation
do1.acc <- function(x) {
  gap <- x$iv2 - x$iv1 + 1
  if (all(is.na(gap))) yy <- tt <- NULL
  else {
    gap <- na.omit(gap)
    yy <- cumsum(rep(1, sum(gap)))
    tt <- unlist(sapply(1:length(gap), function(i)
                        seq(x$iv1[i], x$iv2[i], by=1.0)))
  }
  yy <- c(0, yy, rev(yy)[1])
  if (!is.null(tt[1]) && tt[1] < 0)
    tt <- c(tt[1] - 1, tt, max(x$futime) + 1.0)
  else tt <- c(0, tt, max(x$futime) + 1.0)
  as.lgtdl(data.frame(time=tt, cov=yy))
}

yacc <- lapply(split(dnase[,c("id", "iv1", "iv2", "futime")], dnase$id),
              do1.acc)

## construct data availability (or at risk) indicator process
tu <- max(dat$futime) + 0.001
rt <- lapply(1:nrow(dat),
             function(i) {
               x <- dat[i, "futime"]
               time <- c(0, x, tu)
               cov <- c(1, 0, 0)
               as.lgtdl(data.frame(time=time, cov=cov))
             })

## time-independent covariate matrix
xmat <- model.matrix(~ rx + fev, data=dat)
## time-window in days
tlim <- c(10, 168)
good <- unlist(lapply(yrec, function(x) x$time[1] == 0))

## fully functional temporal process regression
## for recurrent event
m.rec <- tpr(yrec, rt, xmat[,1:3], list(), xmat[,-(1:3),drop=FALSE], list(),
             tis=10:160, w = rep(1, 151), family = poisson(),
             evstr = list(link = 5, v = 3))
par(mfrow=c(1,3), mgp=c(2,1,0), mar=c(4,2,1,0), oma=c(0,2,0,0))
for(i in 1:3) ci.plot(m.rec$tis, m.rec$alpha[,i], sqrt(m.rec$valpha[,i]))
## hypothesis test of significance
## integral test, covariate index 2 and 3
sig.test.int.ff(m.rec, idx=2:3, ncut=2)
sig.test.boots.ff(m.rec, idx=2:3, nsim=1000)
## constant fit
cfit <- cst.fit.ff(m.rec, idx=2:3)

## goodness-of-fit test for constant fit
gof.test.int.ff(m.rec, idx=2:3, ncut=2)
gof.test.boots.ff(m.rec, idx=2:3, nsim=1000)





## for cumulative days in exacerbation
m.acc <- tpr(yacc, rt, xmat[,1:3], list(), xmat[,-(1:3),drop=FALSE], list(),
             tis=10:160, w = rep(1, 151), family = gaussian(),
             evstr = list(link = 1, v = 1))
par(mfrow=c(1,3), mgp=c(2,1,0), mar=c(4,2,1,0), oma=c(0,2,0,0))
for(i in 1:3) ci.plot(m.acc$tis, m.acc$alpha[,i], sqrt(m.acc$valpha[,i]))
}
\keyword{datasets}
