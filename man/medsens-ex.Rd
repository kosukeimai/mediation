\name{medsens.cont} 
\alias{medsens.cont} 
\title{Causal Mediation Analysis - Sensitivity Analysis} 
\description{ 
Function to perform sensitivity analysis with continuous mediator and outcome.
} 
\usage{

medsens.cont(mmodel, ymodel, INT=FALSE, T="T", M="M",DETAIL=TRUE)

} 

\arguments{ 
\item{mmodel}{R model object for mediator.  Can be of class lm, polr, glm, or gam.} 
\item{ymodel}{R model object for outcome.  Can be of class lm, glm, gam, or rq.} 

\item{INT}{If TRUE this indicates that treatment is interacted with mediator in ymodel object. Default is FALSE} 
\item{T}{Name of binary treatment indicator.}
\item{M}{Name of mediator variable.}
\item{DETAIL}{If FALSE estimation is on .1 grid for rho vs. a grid of .01). Default is TRUE}
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu} }
 
\seealso{See also \code{\link{medsens.cont}}, \code{\link{medsens.binary}} }

\examples{ 

#Example with JOBS II Field experiment
#Required packages needed for below examples: MASS, AER, mgcv, quantreg
rm(list=ls())
library(MASS)
library(foreign)
library(AER)
library(Zelig)
library(mgcv)
library(quantreg)
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/data")
load("job.RData")
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/R")
source("mediate.R")
source("medsens.R")
source("medsens.cont.R")

#data(job)
attach(job)

#Estimate outcome and mediator models for continuous treatment and outcome
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=job)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
sens.cont <- medsens.cont(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE)
sens.cont <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE)

summary(sens.cont)
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))

#Estimate outcome and mediator models for continuous treatment and outcome
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
model.y <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)

sens.dich <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE)
summary(sens.dich)
plot(sens.dich, main="JOBS", ylim=c(-.2,.2))


} 
