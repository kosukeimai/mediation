\name{medsens.cont} 
\alias{medsens.cont} 
\title{Causal Mediation Analysis - Sensitivity Analysis} 
\description{ 
Function to perform sensitivity analysis with continuous mediator and outcome.
} 
\usage{

medsens.cont(mmodel, ymodel, INT=FALSE, T="T", M="M")

} 

\arguments{ 
\item{mmodel}{R model object for mediator.  Can be of class lm, polr, glm, or gam.} 
\item{ymodel}{R model object for outcome.  Can be of class lm, glm, gam, or rq.} 

\item{INT}{If true this indicates that treatment is interacted with mediator in ymodel object.} 
\item{T}{Name of binary treatment indicator.}
\item{M}{Name of mediator variable.}
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} }
 
\seealso{See also \code{\link{medsens.cont}}, \code{\link{medsens.binary}} }

\examples{ 

#Example with JOBS II Field experiment
data(jobs)
attach(jobs)

model.m <- lm(job_seek ~ treat + depress1, data=job)
model.y <- lm(depress2 ~ treat + job_seek + depress1, data=job)

sens.1 <- medsens.cont(model.m, model.y, T="treat", M="job_seek")

summary(sens.1)
plot(sens.1, main="Sensitivity Analysis", ylim=c(-.2,.2))


} 
