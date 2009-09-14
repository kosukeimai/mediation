\name{plot.medsens} 
\alias{plot.medsens} 


\title{Plotting results from sensitivity analysis of mediation effects.} 
\description{ 
Plots results from medsens function. Y axis plots mediation effect and x-axis plots the error correlation rho. Standard options for plot function available.
} 
\usage{

  \method{plot}{medsens}(x, sens.par="rho", r.type=1, sign.prod=1, pr.plot=FALSE, smooth.effect=FALSE, smooth.ci=FALSE, levels=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL, ...)

 # \method{print}{plot.medsens}(z)

} 

\arguments{ 
\item{x}{output from \code{medsens} function.} 
\item{sens.par}{type of sensitivity interpretation to be plotted. Default plots mediation effects in terms of sensitivity parameter rho. If "R2" is specified then in terms of variances explained by an unobserved pretreatment confounder.}
\item{r.type}{type of R square interpretation to be used. If 1 is specified then proportion of previously unexplained variance is used. If 2 is specified then proportion of total unexplained variance is used. Only relevant if sens.par is set to "R2".}
\item{sign.prod}{whether the omitted variable affects the mediator and outcome variable in the same direction (1) or different directions (-1). Only relevant if sens.par is set to "R2".}
\item{pr.plot}{if pr.plot=TRUE then proportion mediated will be plotted.}
\item{smooth.effect}{whether the estimated mediation effects are smoothed via a lowess smoother before being plotted.}
\item{smooth.ci}{whether the confidence bands are smoothed via a lowess smoother before being plotted.}
\item{levels}{vector of levels at which to draw contour lines. Only relevant if sens.par is set to "R2".}
\item{xlab}{x-axis label.} 
\item{ylab}{y-axis label.} 
\item{xlim}{range for x-axis.} 
\item{ylim}{range for y-axis.} 
\item{main}{main title for graph.}
\item{...}{additional arguments to be passed.} 
} 

\section{Warning}{The smooth.effect and smooth.ci options should be used with caution since the smoothing could affect substantive implications of the graphical analysis in a significant way.}

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.
Imai, Kosuke, Luke Keele, Dustin Tingley, and Teppei Yamamoto. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer. Forthcoming.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{See also \code{\link{medsens}} }

\examples{ 

#Example with JOBS II Field experiment

data(jobs)

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

#Plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))

#Plot sensitivity analysis using R^2 method. See plot.medsens for additional detail
plot(sens.cont, sens.par="R2", r.type=2, sign.prod=1)
}


} 
