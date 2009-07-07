\name{medsens} 
\alias{medsens} 
\title{Causal Mediation Analysis - Sensitivity Analysis} 
\description{ 
Function to perform sensitivity analysis on mediation effect for violations of sequential ignorability assumption. This allows for a correlation between the error terms of the outcome model and the mediator model. Sensitivity analysis is possible with 1) continuous mediator and continuous outcome, 2) binary outcome and continuous mediator, and 3) continuous outcome and binary mediator. Output from the function can be passed through summary or plot functions which display estimated mediation effects for given values of rho.
} 
\usage{

medsens(x,  rho.by=.1, sims=1000, eps=.Machine$double.eps)

} 

\arguments{ 
\item{x}{Output from \code{mediate} function} 
\item{rho.by}{Specifies the grid for rho (sensitivity parameter) on which estimation is calculated for. Due to computational time we do not recommend decreasing this.}
\item{sims}{Number of draws for parametric bootstrap used to calculate 95 percent confidence intervals.}
\item{eps}{Convergence parameter for FGLS estimation.}
} 

\details{This is the workhorse function for estimating sensitivity analyses for mediation effects. In fitting models for a binary variable a probit function must be used. In fitting the model.y object with a treatmentXmediator interaction the user must input the interation in the format treatment:mediator. Furthermore, the user must make sure to also include the main effect for the mediator in this case.  

Users should note that computation time is several minutes for these functions. Setting rho.by to smaller numbers significantly increases computational time, as does increasing the number of simulations.}

\section{Warning}{These functions assume that all missing values have been removed from the data set.  This can be done using the \code{na.omit()} command before the outcome and mediation models are fitted.}

\value{\code{medsens} returns an object of class "\code{medsens}". Some of the below elements are not available depending on whether INT is
     specified as TRUE or FALSE by the user or depending on the type of model fit   The function \code{summary} is used to obtain a table of the results. The function \code{plot} is used to plot the results.  

      \item{d0}{Point estimate for mediation effect under control.}
      \item{d1}{Point estimate for mediation effect under treatment.}
\item{upper.d0}{Upper confidence interval for mediation effect under control.}
\item{lower.d0}{Lower confidence interval for mediation effect under control.}
\item{upper.d1}{Upper confidence interval for mediation effect under
          treatment.}
\item{lower.d1}{Lower confidence interval for mediation effect under
          treatment.}
     \item{tau}{Point estimate for total effect.}
\item{upper.tau}{Upper confidence interval for total effect.}
\item{lower.tau}{Lower confidence interval for total effect.}
      \item{nu}{Proportion of total effect mediated.}
\item{upper.nu}{Upper confidence interval for proportion mediated.}
\item{lower.nu}{Upper confidence interval for proportion mediated.}
\item{R2star.prod}{Value of the product of R^2 terms for the percentage of additional variance left unexplained, corresponding to value of rho.}
\item{R2tilda.prod}{Value of the product of R^2 terms for the percentage of total variance left unexplained, corresponding to value of rho.}
\item{R2star.thresh}{Value of R^2 for percentage of additional variance left unexplained which corresponds to a ACME of 0.}
\item{R2tilda.thresh}{Value of R^2 for percentage of total variance left unexplained which corresponds to a ACME of 0.}
\item{r.square.y}{R^2 for the outcome model.}
\item{r.square.m}{R^2 for the mediator model.}
}

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.
Imai, Kosuke, Luke Keele, Dustin Tingley and Teppei Yamamoto (2009) Causal Mediation Analysis in R.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
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
model.m <- lm(job_seek ~ treat + depress1, data=jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1, data=jobs)
#pass model objects through mediate function
med.cont <- mediate(model.m, model.y, T="treat", M="job_seek", INT=FALSE, sims=1000)
#pass mediate output through medsens function
sens.cont<-medsens(med.cont, sims=100, rho.by=.2)
#use summary function to display values of rho where 95% confidence interval overlaps with 0.
summary(sens.cont)
#plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))
#plot sensitivity analysis using R^2 method. See plot.medsens for additional detail
plot(sens.cont, sens.par="R2", r.type=2, sign.prod=1)

} 
