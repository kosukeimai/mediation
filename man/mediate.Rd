\name{mediate} 
\alias{mediate} 
\title{Causal Mediation Analysis} 
\description{ 
Function to estimate causal mediation effects.
} 
\usage{

mediate(model.m, model.y, sims=1000, boot=FALSE, 
        treat="treat.name", mediator="med.name", 
        control=NULL, conf.level=.95, 
        control.value=0, treat.value=1, 
        long=TRUE, dropobs=TRUE, robustSE=FALSE, 
        ...)

} 

\arguments{ 
\item{model.m}{ a fitted model object for mediator.  Can be of class 'lm', 'polr', 'glm', 'gam', or 'rq'.} 
\item{model.y}{ a fitted model object for outcome.  Can be of class 'lm', 'polr', 'glm', 'gam', 'vglm', or 'rq'.} 
\item{sims}{ number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.} 
\item{boot}{ if 'FALSE' a quasi-Bayesian approximation is used for confidence interval; if TRUE nonparametric bootstrap will be used.}
\item{conf.level}{ level of the returned two-sided confidence intervals.  Default is to return
the 2.5 and 97.5 percentiles of the simulated quantities.} 
\item{treat}{ a character string indicating the name of treatment variable used in the models.  The treatment can be either binary (integer or a two-valued factor) or continuous (numeric).}
\item{mediator}{ a character string indicating the name of mediator variable used in the models.}
\item{control}{ a character string indicating the name of control group indicator. Only relevant if 'model.y' is of class 'gam'. If provided, 'd0', 'z0' and 'n0' are allowed to differ from 'd1', 'z1' and 'n1', respectively.}
\item{control.value}{ value of the treatment variable in the "control" condition. Only relevant when 'treat' is continuous.  Default is 0.}
\item{treat.value}{ value of the treatment variable in the "treatment" condition. Only relevant when 'treat' is continuous. Default is 1.}
\item{long}{ Logical. If 'TRUE', the output contains the entire sets of simulation draws of the the ACME, direct, proportion mediated, and total effects. Default is 'TRUE'.}
\item{dropobs}{ Logical. If 'TRUE', }
\item{robustSE}{ A logical value. If 'TRUE', Heteroskedasticisy-consistent standard errors will be used in quasi-Bayesian simulations. Ignored if 'boot' is 'TRUE' or neither 'model.m' nor 'model.y' has a method for \code\link{{vcovHC}} in the \code\link{{sandwich}} package.}
\item{...}{ Other arguments passed to \code\link{{vcovHC}} in the \code{\link{sandwich}} package: typically the 'type' argument. Ignored if 'robustSE' is 'FALSE'.}

}

\details{This is the workhorse function for estimating mediation effects for a variety of data types. For a continuous mediator variable and a continuous outcome variable, the results will be identical to the usual Baron and Kenny method.  The function can, however, accomodate other data types including binary outcomes and mediators and discrete mediators.  Continuous variables can also be modeled nonparametrically or semiparametrically.

In fitting the model.y object with an interaction between treatment and mediator, the user must input the interation in the format treatment:mediator.  Furthermore, the user must make sure to set the INT option to TRUE since the current function does not sense the existence of the interaction term automatically. 

Users should note that use of the nonparametric bootstrap requires several minutes of computing time. The current version does not require specifying the INT option as the function senses whether there is a treatment by mediator interaction specified.}

\section{Warning}{These functions assume that all missing values have been removed from the data set.  This can be done using the \code{na.omit()} command before the outcome and mediation models are fitted.}

\value{\code{mediate} returns an object of class "\code{mediate}".  The function \code{summary} is used to obtain a table of the results.  The object \code{mediate} is a list that contains the following components.  Some of these elements are not available depending on whether INT is specified as TRUE or FALSE by the user.
\item{d0}{point estimate for mediation effect under control.}
\item{d1}{point estimate for mediation effect under treatment.}
\item{d0.ci}{confidence interval for mediation effect under control. For all confidence intervals, the confidence level is set to 95 percent.}
\item{d1.ci}{confidence interval for mediation effect under control.}
\item{d0.sims}{Vector for mediation effect under control.}
\item{d1.sims}{Vector for mediation effect under treatment.}
\item{z0}{point estimate for direct effect under control.}
\item{z1}{point estimate for direct effect under treatment.}
\item{z0.ci}{confidence interval for direct effect under control.}
\item{z1.ci}{confidence interval for direct effect under treatment.}
\item{z0.sims}{Vector for direct effect under control.}
\item{z1.sims}{Vector for direct effect under treatment.}
\item{n0}{percentage of total effect due to mediation under control.}
\item{n1}{percentage of total effect due to mediation under treatment.}
\item{n0.ci}{confidence interval for percentage of total effect due to mediation under control.}
\item{n1.ci}{confidence interval for percentage of total effect due to mediation under treatment.}
\item{n0.sims}{Vector for percentage of total effect due to mediation under control.}
\item{n1.sims}{Vector for percentage of total effect due to mediation under treatment.}
\item{tau.coef}{point estimate for total effect.}
\item{tau.ci}{confidence interval for total effect.}
\item{tau.sims}{Vector for total effect.}
\item{boot}{the boot argument used.}
\item{INT}{the INT argument used for presence of interaction term.}
\item{treat}{name of the 'treat' variable used.}
\item{mediator}{name of the 'mediator' variable used.}
\item{conf.level}{Confidence interval used.}
\item{model.m}{call of the mediator model used.}
\item{model.y}{call of the outcome model used.}
\item{treat.value}{Value of treatment variable, default is 1.}
\item{control.value}{Value of control variable, default is 0.}
\item{env.m}{environment in which the mediator model was fitted.}
\item{env.y}{environment in which the outcome model was fitted.}
}

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.
Imai, Kosuke, Luke Keele, Dustin Tingley, and Teppei Yamamoto. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer. Forthcoming.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{\code{\link{medsens}}, \code{\link{plot.mediate}}, \code{\link{summary.mediate}}, \code{\link{mediations}}, \code{\link{vcovHC}}}

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
cont.int <- mediate(b, d, sims=10, treat="treat", mediator="job_seek")
summary(cont.int)

######################################################
#Binary Outcome
######################################################

b <- lm(job_seek ~ treat + depress1, data=jobs)
c <- glm(work1 ~ treat + job_seek + depress1, family=binomial(link="probit"), data=jobs)
binary <- mediate(b, c, sims=10, treat="treat", mediator="job_seek")
summary(binary)

} 
