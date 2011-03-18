\name{plot.mediations} 
\alias{plot.mediations} 
\title{Plotting Indirect, Direct, and Total Effects from Multiple Mediation Analyses} 
\description{ 
Plots results from mediate function. Y axis lists indirect, direct, and total effects and x-axis plots the respective magnitudes. Most standard options for plot function available.
} 
\usage{

\method{plot}{mediations}(model,ask=FALSE,...) {

} 

\arguments{ 
\item{model}{Output from function mediations.} 

}

\details{This function plots the indirect, direct, and total effects from the mediate function using horizontal bars to represent confidence intervals. Because mediations conducts multiple mediation analyses, this plot commands allows users to cycle through the results.
}

\arguments{ 
\item{x}{object of class \code{mediate} as produced by \code{mediate} or \code{mediations} function.}
\item{ask}{Whether to cycle through the plots one at a time. If TRUE then user controls advancing through each step.}
}

\section{Warning}{Users wishing additional functionalities can easily update the function.}

\value{\code{mediations} returns an object of class "\code{mediations}".  The function \code{summary} is used to obtain a table of the results. The plot function plots these quantities.  

}

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.
Imai, Kosuke, Luke Keele, Dustin Tingley, and Teppei Yamamoto. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer. Forthcoming.} 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu} }
 
\seealso{See also \code{\link{mediate}} }

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
plot.mediate(continuous)


} 
