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
\item{boot}{ a logical value. if 'FALSE' a quasi-Bayesian approximation is used for confidence interval; if TRUE nonparametric bootstrap will be used.}
\item{conf.level}{ level of the returned two-sided confidence intervals.  Default is to return
the 2.5 and 97.5 percentiles of the simulated quantities.} 
\item{treat}{ a character string indicating the name of treatment variable used in the models.  The treatment can be either binary (integer or a two-valued factor) or continuous (numeric).}
\item{mediator}{ a character string indicating the name of mediator variable used in the models.}
\item{control}{ a character string indicating the name of control group indicator. Only relevant if 'model.y' is of class 'gam'. If provided, 'd0', 'z0' and 'n0' are allowed to differ from 'd1', 'z1' and 'n1', respectively.}
\item{control.value}{ value of the treatment variable used as the control condition. Only relevant when 'treat' is continuous.  Default is 0.}
\item{treat.value}{ value of the treatment variable used as the treatment condition. Only relevant when 'treat' is continuous. Default is 1.}
\item{long}{ a logical value. If 'TRUE', the output will contain the entire sets of simulation draws of the the average causal mediation effects, direct effects, proportions mediated, and total effect. Default is 'TRUE'.}
\item{dropobs}{ a logical value indicating the behavior when the model frames of 'model.m' and 'model.y' are composed of different observations. If 'TRUE', models will be re-fitted using the intersection of the two data frames. If 'FALSE', error is returned. Default is 'FALSE'.}
\item{robustSE}{ a logical value. If 'TRUE', Heteroskedasticisy-consistent standard errors will be used in quasi-Bayesian simulations. Ignored if 'boot' is 'TRUE' or neither 'model.m' nor 'model.y' has a method for \code{\link{vcovHC}} in the \code{\link{sandwich}} package. Default is 'FALSE'.}
\item{...}{ other arguments passed to \code{\link{vcovHC}} in the \code{\link{sandwich}} package: typically the 'type' argument. Ignored if 'robustSE' is 'FALSE'.}
}

\details{This is the workhorse function for estimating causal mediation effects for a variety of data types. When both the mediator and outcome models are linear, the results will be identical to the usual Baron and Kenny method.  The function can, however, accommodate other data types including binary, ordered and count outcomes and mediators as well as censored outcomes.  Continuous variables can also be modeled nonparametrically, semiparametrically, or using quantile regressions.

Users should note that use of the nonparametric bootstrap often requires significant computing time. }


\value{\code{mediate} returns an object of class "\code{mediate}", a list that contains the components listed below.  Some of these elements are not available if 'long' is set to 'FALSE' by the user.

The function \code{summary} (i.e., \code{summary.mediate}) can be used to obtain a table of the results.  The function \code{plot} (i.e., \code{plot.mediate}) can be used to produce a plot of the estimated average causal mediation, average direct, and total effects along with their confidence intervals.

\item{d0, d1}{point estimates for average causal mediation effects under the control and treatment conditions.}
\item{d0.ci, d1.ci}{confidence intervals for average causal mediation effects. The confidence level is set at the value specified in 'conf.level'.}
\item{d0.sims, d1.sims}{vectors of length 'sims' containing simulation draws of average causal mediation effects.}
\item{z0, z1}{point estimates for average direct effect under the control and treatment conditions.}
\item{z0.ci, z1.ci}{confidence intervals for average direct effects.}
\item{z0.sims, z1.sims}{vectors of length 'sims' containing simulation draws of average direct effects.}
\item{n0, n1}{the "proportions mediated", or the percentages of total effect transmitted by the mediator under the control and treatment conditions.}
\item{n0.ci, n1.ci}{confidence intervals for the proportions mediated.}
\item{n0.sims, n1.sims}{vectors of length 'sims' containing simulation draws of the proportions mediated.}
\item{tau.coef}{point estimate for total effect.}
\item{tau.ci}{confidence interval for total effect.}
\item{tau.sims}{a vector of length 'sims' containing simulation draws of the total effect.}
\item{boot}{logical, the 'boot' argument used.}
\item{INT}{a logical value indicating whether the model specification allows the effects to differ between the treatment and control conditions.}
\item{treat}{a character string indicating the name of the 'treat' variable used.}
\item{mediator}{a character string indicating the name of the 'mediator' variable used.}
\item{conf.level}{the confidence level used. }
\item{model.m}{the mediator model used.}
\item{model.y}{the outcome model used.}
\item{treat.value}{value of the treatment variable used as the treatment condition.}
\item{control.value}{value of the treatment variable used as the control condition.}
\item{nobs}{number of observations in the model frame for 'model.m' and 'model.y'.}
}

\references{Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.} 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, 
Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, 
Luke Keele, Ohio State University, \email{keele.4@osu.edu}, 
Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
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
