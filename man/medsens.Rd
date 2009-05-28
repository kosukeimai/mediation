\name{medsens} 
\alias{medsens} 
\title{Causal Mediation Analysis - Sensitivity Analysis} 
\description{ 
Function to perform sensitivity analysis with continuous mediator and outcome.
} 
\usage{

medsens(model.m, model.y, T="treat.name", M="med.name", INT=FALSE, DETAIL=FALSE, sims=1000, pr.plot=FALSE)

} 

\arguments{ 
\item{model.m}{R model object for mediator.  Can be of class lm, polr, glm, or gam.} 
\item{model.y}{R model object for outcome.  Can be of class lm, glm, gam, or rq.} 
\item{INT}{If TRUE this indicates that treatment is interacted with mediator in ymodel object. Default is FALSE} 
\item{T}{Name of binary treatment indicator.}
\item{M}{Name of mediator variable.}
\item{DETAIL}{If FALSE estimation is on .1 grid for rho vs. a grid of .01).}
\item{sims}{Number of draws for parametric bootstrap.}
\item{pr.plot}{Plot proportion mediated for binary mediator.}
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

#Example with simulation data
data(sim)




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
#model.m <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, family=binomial(link="probit"))
#model.y <- lm(depress2 ~ treat + job_dich  + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
#sens.dichM <- medsens(model.m, model.y, T="treat", M="job_dich", INT=FALSE,  DETAIL=FALSE, sims=10)
#summary(sens.dichM)
#plot(sens.dichM, main="JOBS", ylim=c(-.2,.2))
#continuous outcome and binary mediator w interaction
#model.m <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, family=binomial(link="probit"))
#model.y <- lm(depress2 ~ treat + job_dich + treat:job_dich + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
#sens.dichM.int <- medsens(model.m, model.y, T="treat", M="job_dich", INT=TRUE,  DETAIL=FALSE, sims=10)
#summary(sens.dichM.int)
#plot(sens.dichM.int, main="JOBS", ylim=c(-.2,.2))


} 
