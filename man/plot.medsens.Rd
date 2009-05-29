\name{plot.medsens} 
\alias{plot.medsens} 


\title{Plotting results from sensitivity analysis of mediation effects.} 
\description{ 
Plots results from medsens function. Y axis plots mediation effect and x-axis plots the error correlation rho. Standard options for plot function available.
} 
\usage{

  \method{plot}{medsens}(x, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL,...) 

 # \method{print}{plot.medsens}(z)

} 

\arguments{ 
\item{x}{Output from medsens function.} 
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{imai@princeton.edu} }
 
\seealso{See also \code{\link{medsens}} }

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
#plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))



} 
