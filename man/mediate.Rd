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

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu}, Dustin Tingley, Princeton University, \email{dtingley@princeton.edu} }
 
\seealso{See also \code{\link{medsens.cont}}, \code{\link{medsens.binary}} }

\examples{ 

#Example with JOBS II Field experiment
data(jobs)
attach(jobs)


#Required packages needed for below examples: MASS, AER, mgcv, quantreg
#rm(list=ls())
#library(MASS)
#library(foreign)
#library(AER)
#library(Zelig)
#library(mgcv)
#library(quantreg)
#setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/data")
#load("job.RData")
#setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/R")
#source("mediate.R")
#source("medsens.cont.R")



####################################################
#Continuous Outcome - Continuous Mediator--Equivalent to Baron-Kenny
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
cont.int.boot  <- mediate(b, d, sims=10, boot=TRUE, INT=TRUE, T="treat", M="job_seek")
summary(cont.int.boot)




######################################################
#GAM: Estimation of quantities of interest using generalized additive models, GAM, for the outcome variable while smoothing over the mediator
######################################################
#Effect With No Interaction Assumption
b <- lm(job_seek ~ treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp+income, data=job)
c <- gam(depress2 ~ treat +s(job_seek, bs="cr")+ depress1+age+educ+sex+nonwhite+econ_hard+marital+occp+income , data=job)
gam <- mediate(b, c, sims=10, boot=TRUE, T="treat", M="job_seek")
summary(gam)

#Effect Without No Interaction Assumption
#Note, to make GAM play nice you have to create separate dummies for the treatment and control variables.
b <- lm(job_seek ~ treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp+income, data=job)
d3 <- gam(depress2 ~treat+ s(job_seek, by=control)+ s(job_seek, by=treat) + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp+income, data=job)
gam.int <- mediate(b, d3, sims=10, boot=TRUE, T="treat", M="job_seek", C="control", INT=TRUE)
summary(gam.int)

#This creates nice plots of the effects as calculated in the GAM model
#pdf("GAM_PLOT_Interaction.pdf", width=8.5, height=3, onefile=FALSE, paper="special")
#par(mfrow=c(1,3),mai=c(.55,.55,.4,.1), cex.lab=1.2, cex.axis=1.2 , cex.main=1.2, font.main=1)
#plot(c,   main="All Observations", xlab="Job-search Self-efficacy (M)", ylab="Expected Level of Depression s(M)", ylim=c(-1,3))
#mtext("No Interaction", side=3, line=2, font=2)
#plot(d3,  select=1, main="Control Group", xlab="Job-search Self-efficacy (M)", ylab=expression(s[0]), ylim=c(-1,3))
#mtext("With Interaction", side=3, line=2, adj=1.9, font=2)
#plot(d3,  select=2, main="Treatment Group", xlab="Job-search Self-efficacy (M)", ylab=expression(s[1]), ylim=c(-1,3))
#dev.off()


######################################################
#Quantile regression for outcome variable
#Calculates mediation effects at user selected quantiles of the outcome variable, here the median
######################################################

b <- lm(job_seek ~ treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp+income, data=job)
c <- rq(depress2 ~ treat +job_seek+ depress1+age+educ+sex+nonwhite+econ_hard+marital+occp +income, data=job, tau=.5)
quantile<- mediate(b, c, sims=10, boot=TRUE, T="treat", M="job_seek")
summary(quantile)



######################################################
#Ordinal Mediator
######################################################
#Use polr to estimate ordered probit

b <- polr(job_disc ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, method="probit", Hess=TRUE)
c <- lm(depress2 ~ treat + job_disc + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)

discrete<- mediate(b, c, sims=100, T="treat", M="job_disc")
summary(discrete)


######################################################
#Binary Outcome
######################################################

b <- lm(job_seek ~ treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp+income, data=job)
c <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)
binary <- mediate(b, c, sims=1000, T="treat", M="job_seek")
summary(binary)


#EXAMPLES USING SIMULATED ANSWERS
H:\imai_methods\Mediation\CheckoutCVS\mediation\mediation\data
load("sim.RData")
load("sim.answers.RData")
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/R")
source("mediate.R")

b <- lm(M.cont ~ T+X.1, data=sim)
c <- lm(Y.cont ~T+M.Cont+X.1, data=sim)
cont-cont<-mediate(b,c)
summary(cont-cont)

b <- lm(M.cont ~ T+X.1, data=sim)
c <- glm(Y.dich ~ T+M+X.1, family=binomial(link="probit"), data=sim)
binary <- mediate(b, c, sims=1000, T="T", M="M.cont")
summary(binary)

sim.answers

} 
