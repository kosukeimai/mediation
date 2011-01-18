\name{summary.mediations} 
\alias{summary.mediations} 
\alias{print.summary.mediations} 

\title{Summarizing output from mediation analysis coming out of the wrapper function mediations.} 
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

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.
Imai, Kosuke, Luke Keele, Dustin Tingley, and Teppei Yamamoto. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer. Forthcoming.} 

\author{Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}}
 
\seealso{See also \code{\link{mediate}} and \code{\link{mediations}} }

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
