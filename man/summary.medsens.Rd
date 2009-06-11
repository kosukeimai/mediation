\name{summary.medsens} 
\alias{summary.medsens}
\alias{print.summary.medsens} 
\title{Summarizing results from sensitivity analysis of mediation effects.} 
\description{ 
Function to perform sensitivity analysis on mediation effect for violations of sequential ignorability assumption. The procedure estimates the mediation effect allowing for a correlation between the error terms of the outcome model and the mediator model. The extent of this correlation is expressed in terms of the parameter rho. Sensitivity analysis is possible with continuous mediator and continuous outcome and binary outcome and continuous mediator.
Future versions of the function will also permit sensitivity analysis for a continuous outcome and binary mediator. Output from the function can be passed through summary or plot functions which display estimated mediation effects for given values of rho.
} 
\usage{


  \method{summary}{medsens}(object, ...) 

  \method{print}{summary.medsens}(x, ...)

} 

\arguments{ 
\item{x}{Output from medsens function.}
\item{object}{Output from medsens function.}
\item{...}{Additional arguments to be passed.}
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{See also \code{\link{medsens}} }

\examples{ 

#Example with JOBS II Field experiment
#For illustration purposes simulations set to low number.

#Example with JOBS II Field experiment
data(jobs)

\dontrun{
#########################################
#continuous mediator and continuous outcome
#########################################

#fit parametric model
model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age  + occp + marital + nonwhite + educ + income, data=jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=jobs)
#pass model objects through medsens function
sens.cont <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE, sims=1000)
#use summary function to display values of rho where 95% confidence interval overlaps with 0.
summary(sens.cont)
#plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))




#########################################
#binary outcome and continuous mediator
#########################################

model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=jobs)
model.y <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=jobs)
sens.dichO <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE)
summary(sens.dichO)
plot(sens.dichO, main="JOBS", ylim=c(-.2,.2))
}


} 
