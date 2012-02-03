\name{mediate.npdiscrete} 
\alias{mediate.npdiscrete} 
\title{Causal Mediation Analysis} 
\description{ 
'mediate.npdiscrete' is used to estimate various average causal mediation effects (indirect effect) nonparametrically when the mediator is discrete
} 
\usage{

mediate.npdiscrete(M, Y, T, boot=FALSE, sims=1000, conf.level = 0.95)

} 

\arguments{ 
\item{M}{ The mediator variable, must be discrete.} 
\item{Y}{ The outcome variable.} 
\item{T}{ The treatment indicator, must be binary.} 
\item{boot}{ a logical value. if 'FALSE' a large sample Delta method approximation is used for confidence intervals; if 'TRUE' nonparametric bootstrap will be used. Default is 'FALSE'.}
\item{sims}{ number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.} 
\item{conf.level}{ level of the returned two-sided confidence intervals.  Default is to return
the 2.5 and 97.5 percentiles of the simulated quantities.} 
}

\details{This function estimates mediation effects nonparametrically when the mediator is discrete.  This function basically estimates the mediation effects from a set of cross-tabulations as defined by the mediator and treatment variables.  This function is not meant to be the main function for estimating mediation effects.  It simply allows users one method for nonparametric estimation.  This function does not accommodate conditioning on pretreatment covariates or allow for a sensitivity analysis. In general, users can achieve similar or better functionality using the \code{mediate} function. This function will produce errors if for some combinations of treatment and mediator, there are no outcomes observed.  If this happens use the mediate function, but be aware that the estimates are a function of a model based interpolation.
}


\value{\code{mediate.npdiscrete} returns an object of class "\code{mediate.npdiscrete}", a list that contains the components listed below.  

The function \code{summary} (i.e., \code{summary.mediate.npdiscrete}) can be used to obtain a table of the results. 
\item{d0, d1}{point estimates for average causal mediation effects under the control and treatment conditions.}
\item{d0.ci, d1.ci}{confidence intervals for average causal mediation effects. The confidence level is set at the value specified in 'conf.level'.}
\item{boot}{logical, the 'boot' argument used.}
\item{conf.level}{the confidence level used. }
\item{nobs}{number of observations in the data.}
}

\references{Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.} 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, 
Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, 
Luke Keele, Penn State University, \email{ljk20@psu.edu}, 
Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{\code{\link{mediate}}}

\examples{ 
# Examples with JOBS II Field Experiment

data(jobs)

#Use the discrete version of the mediator.

foo.1 <- mediate.npdiscrete(job_disc, depress2, treat)
summary(foo.1)

foo.2 <- mediate.npdiscrete(job_disc, depress2, treat, boot=TRUE)
summary(foo.2)

} 
