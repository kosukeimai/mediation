\name{mediate} 
\alias{mediate} 
\title{Causal Mediation Analysis} 
\description{ 
Function to estimate average causal mediation effects.
} 
\usage{

#Default Method 
mediate(mmodel, ymodel, sims=1000, boot=FALSE, INT=FALSE, T="T", M="M")

} 

\arguments{ 
\item{mmodel}{R model object for mediator.  Can be of class lm, polr, glm, or gam.} 
\item{ymodel}{R model object for outcome.  Can be of class lm, glm, gam, or rq.} 
\item{sims}{Number of draws for bootstrap.} 
\item{boot}{If FALSE parametric bootstrap is used for confidence interval, if TRUE nonparametric bootstrap will be used.}
\item{INT}{If true this indicates that treatment is interacted with mediator in ymodel object.} 
\item{T}{Name of binary treatment indicator.}
\item{M}{Name of mediator variable.}
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} }
 
\seealso{See also \code{\link{medsens.cont}}, \code{\link{medsens.binary}} }

\examples{ 

#Example with JOBS II Field experiment
data(jobs)
attach(jobs)

#Continuous Mediator and Outcome
m.model <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=job)
y.model <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)

#Nonparametric Boostrap - NB this takes several minutes
med.1 <- mediate(m.model, y.model, sims=1000, boot=TRUE, T="treat", M="job_seek")
med.2 <- mediate(y.model, y.model, sims=1000, T="treat", M="job_seek")

summary(med.1)
summary(med.2)

#Relax No-Interaction Assumption

y.model <- lm(depress2 ~ treat + job_seek + treat:job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)

med.3 <- mediate(m.model, y.model, sims=1000, boot=TRUE, INT=TRUE, T="treat", M="job_seek")

summary(med.3)

#Generalized additive model (GAM) for outcome allowing smoothing on the mediator
a<-lm(depress2 ~ treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp+income, data=job)
b <- lm(job_seek ~ treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp+income, data=job)
c <- gam(depress2 ~ treat +s(job_seek, bs="cr")+ depress1+age+educ+sex+nonwhite+econ_hard+marital+occp+income , data=job)
gam <- mediate.binary(b, c, sims=sims, boot=TRUE, T="treat", M="job_seek")
summary(gam)



#Quantile regression for outcome variable
b <- lm(job_seek ~ treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp, data=job)
c <- rq(depress2 ~ treat +job_seek+ depress1+age+educ+sex+nonwhite+econ_hard+marital+occp +income, data=job, tau=.5)



#Binary Outcome
y.model <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)

quantile <- mediate(m.model, y.model, sims=1000, boot=TRUE, T="treat", M="job_seek")

summary(quantile)


#Ordinal Mediator

b <- polr(job_disc ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, method="probit", Hess=TRUE)
c <- lm(depress2 ~ treat + job_disc + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)

discrete<- mediate.binary(b, c, sims=1000, T="treat", M="job_disc")
summary(discrete)

} 
