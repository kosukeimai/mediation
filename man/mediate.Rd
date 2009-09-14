\name{mediate} 
\alias{mediate} 
\title{Causal Mediation Analysis} 
\description{ 
Function to estimate causal mediation effects.
} 
\usage{

mediate(model.m, model.y, sims=1000, boot=FALSE, INT=FALSE, treat="treat.name",
mediator="med.name", control=NULL)

} 

\arguments{ 
\item{model.m}{ R model object for mediator.  Can be of class lm, polr, glm, or gam.} 
\item{model.y}{ R model object for outcome.  Can be of class lm, glm, gam, or rq.} 
\item{sims}{ Number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.} 
\item{boot}{ If FALSE a quasi-Bayesian approximation is used for confidence interval, if TRUE nonparametric bootstrap will be used.}
\item{INT}{ If TRUE this indicates that treatment is interacted with mediator in \code{model.y} object. Default is FALSE.} 
\item{treat}{ Name of binary treatment indicator.}
\item{mediator}{ Name of mediator variable.}
\item{control}{ Name of binary treatment indicator for control. Only necessary for gam() object with interaction.}
}

\details{This is the workhorse function for estimating mediation effects for a variety of data types. For a continuous mediator variable and a continuous outcome variable, the results will be identical to the usual Baron and Kenny method.  The function can, however, accomodate other data types including binary outcomes and mediators and discrete mediators.  Continuous variables can also be modeled nonparametrically or semiparametrically.

In fitting the model.y object with an interaction between treatment and mediator, the user must input the interation in the format treatment:mediator.  Furthermore, the user must make sure to set the INT option to TRUE since the current function does not sense the existence of the interaction term automatically. 

Users should note that use of the nonparametric bootstrap requires several minutes of computing time.}

\section{Warning}{These functions assume that all missing values have been removed from the data set.  This can be done using the \code{na.omit()} command before the outcome and mediation models are fitted.}

\value{\code{mediate} returns an object of class "\code{mediate}".  The function \code{summary} is used to obtain a table of the results.  The object \code{mediate} is a list that contains the following components.  Some of these elements are not available depending on whether INT is specified as TRUE or FALSE by the user.
\item{d0}{point estimate for mediation effect under control.}
\item{d1}{point estimate for mediation effect under treatment.}
\item{d0.ci}{confidence interval for mediation effect under control. For all confidence intervals, the confidence level is set to 95 percent.}
\item{d1.ci}{confidence interval for mediation effect under control.}
\item{tau.coef}{point estimate for total effect.}
\item{tau.ci}{confidence interval for total effect.}
\item{z0}{point estimate for direct effect under control.}
\item{z1}{point estimate for direct effect under treatment.}
\item{z0.ci}{confidence interval for direct effect under control.}
\item{z1.ci}{confidence interval for direct effect under control.}
\item{pct.coef}{percentage of total effect due to mediation.}
\item{pct.ci}{confidence interval for percentage of total effect due to mediation.}
\item{boot}{the boot argument used.}
\item{INT}{the INT argument used.}
\item{treat}{name of the 'treat' variable used.}
\item{mediator}{name of the 'mediator' variable used.}
\item{call.m}{call of the mediator model used.}
\item{call.y}{call of the outcome model used.}
\item{env.m}{environment in which the mediator model was fitted.}
\item{env.y}{environment in which the outcome model was fitted.}
}

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.
Imai, Kosuke, Luke Keele, Dustin Tingley, and Teppei Yamamoto. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer. Forthcoming.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{See also \code{\link{medsens}} }

\examples{ 
#Example with JOBS II Field experiment

#**For illustration purposes simulations set to low number**

data(jobs)

####################################################
#Continuous Outcome and Continuous Mediator
####################################################

b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
c <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age, data=jobs)

#Calculates quantities using quasi-Bayesian approximation
#The general format of the function is to record two model objects
#The first for the mediator the second for the outcome variable 
continuous <- mediate(b, c, sims=10, treat="treat", mediator="job_seek")
summary(continuous)

#Calculates quantities using the non-parametric bootstrap - This takes several minutes
continuous_boot <- mediate(b, c, boot=TRUE, sims=10, treat="treat", mediator="job_seek")
summary(continuous_boot)

#Interaction term between treatment and mediator is included in the model
b <- lm(job_seek ~ treat + depress1, data=jobs)
d <- lm(depress2 ~ treat + job_seek + treat:job_seek + depress1, data=jobs)
cont.int <- mediate(b, d, sims=10, INT=TRUE, treat="treat", mediator="job_seek")
summary(cont.int)

######################################################
#Binary Outcome
######################################################

b <- lm(job_seek ~ treat + depress1, data=jobs)
c <- glm(work1 ~ treat + job_seek + depress1, family=binomial(link="probit"), data=jobs)
binary <- mediate(b, c, sims=10, treat="treat", mediator="job_seek")
summary(binary)

} 
