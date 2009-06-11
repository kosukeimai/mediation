\name{mediate} 
\alias{mediate} 
\title{Causal Mediation Analysis} 
\description{ 
Function to estimate causal mediation effects.
} 
\usage{

#Default Method 
mediate(model.m, model.y, sims=1000, boot=FALSE, INT=FALSE, T="treat.name",
M="med.name", C=NULL)

} 

\arguments{ 
\item{model.m}{ R model object for mediator.  Can be of class lm, polr, glm, or gam.} 
\item{model.y}{ R model object for outcome.  Can be of class lm, glm, gam, or rq.} 
\item{sims}{ Number of draws for bootstrap.} 
\item{boot}{ If FALSE parametric bootstrap is used for confidence interval, if TRUE nonparametric bootstrap will be used.}
\item{INT}{ If TRUE this indicates that treatment is interacted with mediator in \code{model.y} object. Default is FALSE.} 
\item{T}{ Name of binary treatment indicator.}
\item{M}{ Name of mediator variable.}
\item{C}{ Name of binary treatment indicator for control. Only necessary for gam() object with interaction.}
}

\details{This is the workhorse function for estimating mediation effects for a variety of data types. For a continuous mediator variable and a continuous outcome variable, the results will be identical to the usual Baron and Kenny method.  The function can, however, accomodate other data types including binary outcomes and mediators and discrete mediators.  Continuous variables can also be modeled nonparametrically or semiparametrically.

Users should note that use of the nonparametric bootstrap requires several minutes of computing time.}

\section{Warning}{These functions assume that all missing values have been removed from the data set.  This can be done using the \code{na.omit()} command before the outcome and mediation models are fitted.}

\value{\code{mediate} returns an object of class "\code{mediate}".  The function \code{summary} is used to obtain a table of the results.  The object \code{mediate} is a list that contains the following components.  Some of these elements are not available depending on whether INT is specified as TRUE or FALSE by the user.
\item{d0}{Point estimate for mediation effect under control.}
\item{d1}{Point estimate for mediation effect under treatment.}
\item{d0.ci}{Confidence interval for mediation effect under control.}
\item{d1.ci}{Confidence interval for mediation effect under control.}
\item{tau.coef}{Point estimate for total effect.}
\item{tau.ci}{Confidence interval for total effect.}
\item{z0}{Point estimate for direct effect under control.}
\item{z1}{Point estimate for direct effect under treatment.}
\item{z0.ci}{Confidence interval for direct effect under control.}
\item{z1.ci}{Confidence interval for direct effect under control.}
\item{pct.coef}{Percentage of total effect due to mediation.}
\item{pct.ci}{Confidence interval for percentage of total effect due to mediation.}
}

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.
Imai, Kosuke, Luke Keele, Dustin Tingley and Teppei Yamamoto (2009) Causal Mediation Analysis in R.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{See also \code{\link{medsens}} }

\examples{ 
#For illustration purposes simulations set to low number.

#Example with JOBS II Field experiment
data(jobs)

####################################################
#Continuous Outcome - Continuous Mediator--Equivalent to Baron-Kenny Method
####################################################

b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
c <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age, data=jobs)

#Calculates quantities using parametric bootstrap
#The general format of the function is to record two model objects
#The first for the mediator the second for the outcome variable 
continuous <- mediate(b, c, sims=1000, T="treat", M="job_seek")
summary(continuous)

\dontrun{
#Calculates quantities using the non-parametric bootstrap - This takes several minutes
continuous_boot <- mediate(b, c, boot=TRUE, T="treat", M="job_seek")
summary(continuous_boot)

#INTERACTION TERM BTWN TREATEMENT AND MEDIATOR IS INCLUDED IN THE MODEL
b <- lm(job_seek ~ treat + depress1, data=jobs)
d <- lm(depress2 ~ treat + job_seek + treat:job_seek + depress1, data=jobs)
cont.int <- mediate(b, d, sims=10, INT=TRUE, T="treat", M="job_seek")
summary(cont.int)


######################################################
#Binary Outcome
######################################################

b <- lm(job_seek ~ treat + depress1, data=jobs)
c <- glm(work1 ~ treat + job_seek + depress1, family=binomial(link="probit"), data=jobs)
binary <- mediate(b, c, T="treat", M="job_seek")
summary(binary)
}


} 
