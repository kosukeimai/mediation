\name{medsens} 
\alias{medsens} 
\title{Causal Mediation Analysis - Sensitivity Analysis} 
\description{ 
Function to perform sensitivity analysis with continuous mediator and outcome.
} 
\usage{

medsens(mmodel, ymodel, INT=FALSE, T="T", M="M",DETAIL=TRUE)

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
 
\seealso{See also \code{\link{mediate}} }

\examples{ 

#Example with JOBS II Field experiment
#Required packages needed for below examples: MASS, AER, mgcv, 
#For illustration purposes simulations set to low number.

#Example with JOBS II Field experiment
data(jobs)
attach(jobs)
#Example with simulation data
data(sim)
attach(sim)



#continuous mediator and continuous outcome
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=job)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
#sens.cont <- medsens.cont(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE)
sens.cont <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE, sims=10)
summary(sens.cont)
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))





#binary outcome and continuous mediator
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
model.y <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)
sens.dichO <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE, sims=10)
summary(sens.dichO)
plot(sens.dichO, main="JOBS", ylim=c(-.2,.2))




#continuous outcome and binary mediator w no interaction
model.m <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, family=binomial(link="probit"))
model.y <- lm(depress2 ~ treat + job_dich  + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
sens.dichM <- medsens(model.m, model.y, T="treat", M="job_dich", INT=FALSE,  DETAIL=FALSE, sims=10)
summary(sens.dichM)
plot(sens.dichM, main="JOBS", ylim=c(-.2,.2))
#continuous outcome and binary mediator w interaction
model.m <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, family=binomial(link="probit"))
model.y <- lm(depress2 ~ treat + job_dich + treat:job_dich + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
sens.dichM.int <- medsens(model.m, model.y, T="treat", M="job_dich", INT=TRUE,  DETAIL=FALSE, sims=10)
summary(sens.dichM.int)
plot(sens.dichM.int, main="JOBS", ylim=c(-.2,.2))

############################3
#Simulation Data
############################

#No interaction-only binmed sens does not work and match analytical
model.m <- lm(M.cont ~ T+X.1, data=sim)
model.y <- lm(Y.cont ~ T+M.cont+X.1, data=sim)
cont.cont.sens <- medsens(model.m, model.y,  T="T", M="M.cont", DETAIL=FALSE, sims=10)
cont.cont <- mediate(model.m, model.y, sims=1000, T="T", M="M.cont")
summary(cont.cont.sens)
summary(cont.cont)


model.m <- lm(M.cont ~ T+X.1, data=sim)
model.y <- glm(Y.dich ~ T+M.cont+X.1, family=binomial(link="probit"), data=sim)
binary.outcome.sens <- medsens(model.m, model.y,  T="T", M="M.cont", INT=FALSE,  DETAIL=FALSE, sims=10)
summary(binary.outcome.sens)


model.m  <- glm(M.dich ~ T+X.1, , family=binomial(link="probit"), data=sim)
model.y  <- lm(Y.cont.dich ~ T+M.dich+X.1, data=sim)
binary.mediator.sens <- medsens(model.m , model.y,  T="T", M="M.dich", INT=FALSE,  DETAIL=FALSE, sims=10)
summary(binary.mediator.sens)



#interaction
model.m  <- glm(M.dich ~ T+X.1, , family=binomial(link="probit"), data=sim)
model.y  <- lm(Y.cont.dich.int ~ T+M.dich+T:M.dich+X.1)
binary.mediator.sens.int <- medsens(model.m, model.y,  T="T", M="M.dich", INT=TRUE,  DETAIL=FALSE, sims=10)
summary(binary.mediator.sens.int)




} 
