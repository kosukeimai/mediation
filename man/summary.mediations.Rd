\name{summary.mediations} 
\alias{summary.mediations} 
\alias{print.summary.mediations} 

\title{Summarizing Output from Multiple Mediation Analyses} 
\description{ 
Function to report results from mediation analysis. Reported categories are mediation effect, direct effect, total effect, and proportion of total effect mediated. 
All quantities reported with confidence intervals specified in mediations. Prints summaries for all models estimated in mediations.
} 
\usage{

  \method{summary}{mediations}(object, ...) 

  \method{print}{summary.mediations}(x, ...)

} 

\arguments{ 
\item{object}{Output from mediations function.} 
\item{x}{Output from mediations function.} 
\item{...}{Additional arguments to be passed.}
}

\author{Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}}
 
\seealso{\code{\link{mediations}} }

\examples{ 

#Example with JOBS II Field experiment

#**For illustration purposes simulations set to low number**

data(jobs)

#########################################
#continuous mediator and continuous outcome
#########################################

model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=jobs)
continuous <- mediate(model.m, model.y, sims=10, treat="treat", mediator="job_seek")
summary(continuous)


} 
