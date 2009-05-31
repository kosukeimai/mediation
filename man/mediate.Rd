\name{mediate} 
\alias{mediate} 
\title{Causal Mediation Analysis} 
\description{ 
Function to estimate causal mediation effects.
} 
\usage{

#Default Method 
mediate(model.m, model.y, sims=1000, boot=FALSE, INT=FALSE, T="treat.name", M="med.name", C=NULL)

} 

\arguments{ 
\item{model.m}{R model object for mediator.  Can be of class lm, polr, glm, or gam.} 
\item{model.y}{R model object for outcome.  Can be of class lm, glm, gam, or rq.} 
\item{sims}{Number of draws for bootstrap.} 
\item{boot}{If FALSE parametric bootstrap is used for confidence interval, if TRUE nonparametric bootstrap will be used.}
\item{INT}{If true this indicates that treatment is interacted with mediator in ymodel object.} 
\item{T}{Name of binary treatment indicator.}
\item{M}{Name of mediator variable.}
\item{C}{Name of binary treatment indicator for control. Only necessary for gam() object with interaction.}
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{imai@princeton.edu} }
 
\seealso{See also \code{\link{medsens}} }

\examples{ 

#Required packages needed for below examples: MASS, AER, mgcv, quantreg
#For illustration purposes simulations set to low number.

#Example with JOBS II Field experiment
data(job)


####################################################
#Continuous Outcome - Continuous Mediator--Equivalent to Baron-Kenny Method
####################################################

b <- lm(job_seek ~ treat + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=job)
c <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)

#Calculates quantities using parametric bootstrap
#The general format of the function is to record two model objects, the first for the mediator the second for the outcome variable 
continuous <- mediate(b, c, sims=10, T="treat", M="job_seek")
summary(continuous)
#Calculates quantities using the non-parametric bootstrap - This takes several minutes
continuous_boot <- mediate(b, c, sims=10, boot=TRUE, T="treat", M="job_seek")
summary(continuous_boot)

#INTERACTION TERM BTWN TREATEMENT AND MEDIATOR IS INCLUDED IN THE MODEL
b <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=job)
d <- lm(depress2 ~ treat + job_seek + treat:job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
cont.int <- mediate(b, d, sims=10, INT=TRUE, T="treat", M="job_seek")
summary(cont.int)


######################################################
#Binary Outcome
######################################################

b <- lm(job_seek ~ treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp+income, data=job)
c <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)
binary <- mediate(b, c, sims=10, T="treat", M="job_seek")
summary(binary)



} 
