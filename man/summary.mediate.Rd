\name{summary.mediate} 
\alias{summary.mediate} 
\alias{print.summary.mediate} 

\title{Summarizing output from mediation analysis.} 
\description{ 
Function to report results from mediation analysis. Reported categories are mediation effect, direct effect, total effect, and proportion of total effect mediated. 
All quantities reported with 95 percent confidence intervals. If the \code{INT = TRUE} option is selected in the mediate function, which allows for interaction between treatment and mediator, then effects are reported separately for the treatment and control groups.
} 
\usage{

  \method{summary}{mediate}(object, ...) 

  \method{print}{summary.mediate}(x, ...)

} 

\arguments{ 
\item{object}{Output from mediate function.} 
\item{x}{Output from mediate function.} 
 \item{...}{Additional arguments to be passed.}
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.
Imai, Kosuke, Luke Keele, Dustin Tingley, and Teppei Yamamoto. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer. Forthcoming.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{See also \code{\link{mediate}} }

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
