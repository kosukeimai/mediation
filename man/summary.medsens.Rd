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
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.
Imai, Kosuke, Luke Keele, Dustin Tingley, and Teppei Yamamoto. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer. Forthcoming.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{See also \code{\link{medsens}} }

\examples{ 

\dontrun{
#########################################
#Continuous mediator and continuous outcome
#########################################

#Fit parametric models
model.m <- lm(job_seek ~ treat + depress1, data=jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1, data=jobs)

#Pass model objects through mediate function
med.cont <- mediate(model.m, model.y, treat="treat", mediator="job_seek", sims=1000)

#Pass mediate output through medsens function
sens.cont <- medsens(med.cont, sims=1000, rho.by=.1)

#Use summary function to display values of rho where 95% confidence interval overlaps with 0
summary(sens.cont)

#########################################
#binary outcome and continuous mediator
#########################################

model.m <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=jobs)
model.y <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=jobs)
med.bin <- mediate(model.m, model.y, treat="treat", mediator="job_seek", sims=1000)
sens.bin <- medsens(med.bin, sims=1000, rho.by=.1)
summary(sens.bin)

} 
}
