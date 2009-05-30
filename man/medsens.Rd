\name{medsens} 
\alias{medsens} 
\title{Causal Mediation Analysis - Sensitivity Analysis} 
\description{ 
Function to perform sensitivity analysis on mediation effect for violations of sequential ignorability assumption. This allows for a correlation between the error terms of the outcome model and the mediator model. Sensitivity analysis is possible with continuous mediator and continuous outcome and binary outcome and continuous mediator.
Future versions of the function will also permit sensitivity analysis for a continuous outcome and binary mediator. Output from the function can be passed through summary or plot functions which display estimated mediation effects for given values of rho.
} 
\usage{

medsens(model.m, model.y, T="treat.name", M="med.name", INT=FALSE, DETAIL=FALSE, sims=1000, pr.plot=FALSE)

} 

\arguments{ 
\item{model.m}{R model object for mediator.  Can be of class lm or glm using a probit link function.} 
\item{model.y}{R model object for outcome.  Can be of class lm or glm using a probit link function.} 
\item{INT}{If TRUE this indicates that treatment is interacted with mediator in ymodel object. Default is FALSE. Interaction implemented by including a T:M interaction term in the outcome model. This option is currently only available for continuous outcome continuous mediator case.} 
\item{T}{Name of binary treatment indicator, 1 for treatment and 0 if control. Variable name listed in quotations.}
\item{M}{Name of mediator. Variable name listed in quotations.}
\item{DETAIL}{If FALSE estimation is on .1 grid for rho (sensitivity parameter) vs. a grid of .01). Due to computational time FALSE is recommended.}
\item{sims}{Number of draws for parametric bootstrap used to calculate 95 percent confidence intervals.}
\item{pr.plot}{Plot proportion mediated for binary outcome variable.}
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{imai@princeton.edu} }
 
\seealso{See also \code{\link{mediate}} }

\examples{ 

#Example with JOBS II Field experiment
#For illustration purposes simulations set to low number.

#Example with JOBS II Field experiment
data(job)

#Example with simulation data
data(sim)



#########################################
#continuous mediator and continuous outcome
#########################################

#fit parametric model
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=job)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
#pass model objects through medsens function
sens.cont <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE, sims=10)
#use summary function to display values of rho where 95% confidence interval overlaps with 0.
summary(sens.cont)
#plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))




#########################################
#binary outcome and continuous mediator
#########################################

model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
model.y <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)
sens.dichO <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE, sims=100)
summary(sens.dichO)
plot(sens.dichO, main="JOBS", ylim=c(-.2,.2))



#########################################
#continuous outcome and binary mediator w no interaction
#NOT YET IMPLEMENTED...DO NOT DRAW INFERENCES FROM RESULTS
#########################################

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
