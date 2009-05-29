\name{summary.mediation} 
\alias{summary.mediation} 
\alias{print.summary.mediation} 

\title{Summarizing output from mediation analysis.} 
\description{ 
Function to report summary statistics from mediation analysis. Reported categories are mediation effect, direct effect, total effect, and proportion of total effect mediated. 
All quantities reported with 95 percent confidence intervals. If the INT=TRUE option is selected in the mediate function, which allows for a treatment X mediator interaction, then effects are reported separetly for the treatment and control groups.
} 
\usage{

  \method{summary}{mediation}(object, ...) 

  \method{print}{summary.mediation}(x, ...)

} 

\arguments{ 
\item{object}{Output from mediate function.} 
\item{x}{Output from mediate function.} 
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{imai@princeton.edu} }
 
\seealso{See also \code{\link{mediate}} }

\examples{ 

#Example with JOBS II Field experiment
#For illustration purposes simulations set to low number.

#Example with JOBS II Field experiment
data(jobs)


#########################################
#continuous mediator and continuous outcome
#########################################

#fit parametric model
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=job)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
#pass model objects through medsens function
sens.cont <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE, sims=10)
#use summary function to display values of rho where 95% confidence interval overlaps with 0.
summary(sens.cont)
#plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))





} 
