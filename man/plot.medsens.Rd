\name{plot.medsens} 
\alias{plot.medsens} 


\title{Plotting results from sensitivity analysis of mediation effects.} 
\description{ 
Plots results from medsens function. Y axis plots mediation effect and x-axis plots the error correlation rho. Standard options for plot function available.
} 
\usage{

  \method{plot}{medsens}(x, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL, sens.par="rho", r.type=1, sign.prod=1, pr.plot=FALSE, ...) 

 # \method{print}{plot.medsens}(z)

} 

\arguments{ 
\item{x}{Output from medsens function.} 
\item{xlab}{x-axis label.} 
\item{ylab}{y-axis label.} 
\item{xlim}{range for x-axis.} 
\item{ylim}{range for y-axis.} 
\item{main}{main title for graph.}
\item{sens.par}{Type of sensitivity interpretation to be plotted. Default plots ACME in terms of sensitivity parameter rho. If "R2" is specified then in terms of unexplained variance.}
\item{r.type}{Type of R^2 interpretation to be used. If 1 is specified then percentage of previously unexplained variance is used. If 2 is specified then percentage of total unexplained variance is used. Only relevant if sens.par="R2" is specified.}
\item{sign.prod}{Whether the omitted variable has the same impact on both the mediator and outcome variable (1) or a different impact (2). Only relevant if sens.par="R2" is specified.}
\item{pr.plot}{If pr.plot=TRUE then proportion mediated will be plotted.}
\item{...}{Additional arguments to be passed.} 
} 

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{See also \code{\link{medsens}} }

\examples{ 

#Example with JOBS II Field experiment
#For illustration purposes simulations set to low number.

data(jobs)


\dontrun{
#########################################
#continuous mediator and continuous outcome
#########################################

#fit parametric model
model.m <- lm(job_seek ~ treat + depress1, data=jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1, data=jobs)
#pass model objects through mediate function
med.cont <- mediate(model.m, model.y, T="treat", M="job_seek", INT=FALSE, sims=1000)
#pass mediate output through medsens function
sens.cont<-medsens(med.cont, sims=1000, rho.by=.2)
#use summary function to display values of rho where 95% confidence interval overlaps with 0.
summary(sens.cont)
#plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))
#plot sensitivity analysis using R^2 method. See plot.medsens for additional detail
plot(sens.cont, sens.par="R2", r.type=2, sign.prod=1)
}


} 
