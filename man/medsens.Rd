\name{medsens} 
\alias{medsens} 
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
load("sim.RData")
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/R")
source("mediate.R")
source("medsens.R")
source("medsens.cont.R")

#data(job)
attach(job)
attach(sim)

#continuous mediator and continuous outcome
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=job)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
#sens.cont <- medsens.cont(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE)
sens.cont <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE, sims=1000)
summary(sens.cont)
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))





#binary outcome and continuous mediator
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
model.y <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)
sens.dichO <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE)
summary(sens.dichO)
plot(sens.dichO, main="JOBS", ylim=c(-.2,.2))



#continuous outcome and binary mediator w no interaction
model.m <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, family=binomial(link="probit"))
model.y <- lm(depress2 ~ treat + job_dich  + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
sens.dichM <- medsens(model.m, model.y, T="treat", M="job_dich", INT=FALSE,  DETAIL=FALSE)
summary(sens.dichM)
plot(sens.dichM, main="JOBS", ylim=c(-.2,.2))
#continuous outcome and binary mediator w interaction
model.m <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, family=binomial(link="probit"))
model.y <- lm(depress2 ~ treat + job_dich + treat:job_dich + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
sens.dichM.int <- medsens(model.m, model.y, T="treat", M="job_dich", INT=TRUE,  DETAIL=FALSE)
summary(sens.dichM.int)
plot(sens.dichM.int, main="JOBS", ylim=c(-.2,.2))


#Simulation Data
model.m <- lm(M.cont ~ T+X.1, data=sim)
model.y <- glm(Y.dich ~ T+M.cont+X.1, family=binomial(link="probit"), data=sim)
binary.outcome.sens <- medsens(b, c, sims=1000, T="T", M="M.cont")
binary.outcome <- mediate(b, c, sims=1000, T="T", M="M.cont")
summary(binary.outcome.sens)
summary(binary.outcome)

model.m  <- glm(M.dich ~ T+X.1, , family=binomial(link="probit"), data=sim)
model.y  <- lm(Y.cont.dich.int ~ T+M.dich+T:M.dich+X.1)
sens.dichM.int <- medsens(model.m, model.y, T="T", M="M.dich", INT=TRUE,  DETAIL=FALSE)
binary.mediator <-mediate(b, c, sims=1000, T="T", M="M.dich", INT=TRUE)
summary(binary.mediator.int)


b <- glm(M.dich ~ T+X.1, , family=binomial(link="probit"), data=sim)
c <- lm(Y.cont.dich.int ~ T+M.dich+T:M.dich+X.1)
sens.dichM.int <- medsens(b, c, T="T", M="M.dich", INT=TRUE,  DETAIL=FALSE)


} 
