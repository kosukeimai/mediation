\name{medsens} 
\alias{medsens} 
\title{Causal Mediation Analysis - Sensitivity Analysis} 
\description{ 
Function to perform sensitivity analysis on mediation effect for violations of sequential ignorability assumption. This allows for a correlation between the error terms of the outcome model and the mediator model. Sensitivity analysis is possible with 1) continuous mediator and continuous outcome, 2) binary outcome and continuous mediator, and 3) continuous outcome and binary mediator. Output from the function can be passed through summary or plot functions which display estimated mediation effects for given values of rho.
} 
\usage{

medsens(x, rho.by=.1, sims=1000, eps=.Machine$double.eps)

} 

\arguments{ 
\item{x}{Output from \code{mediate} function} 
\item{rho.by}{Specifies the grid for rho (sensitivity parameter) on which estimation is calculated for. Decreasing
this will result in more precise sensitivity estimates but significantly more computational time.}
\item{sims}{Number of Monte Carlo draws for calculating 95 percent confidence intervals. This argument
is used only in cases where either binary mediator or binary outcome is involved.}
\item{eps}{Convergence parameter for FGLS estimation. This argument is used only in the continuous-continuous case.}
} 

\details{This is the workhorse function for estimating sensitivity analyses for mediation effects. In fitting models for a binary variable a probit function must be used.

Users should note that computation time is several minutes for these functions. Setting rho.by to smaller numbers significantly increases computational time, as does increasing eps or the number of simulations.}

\section{Warning}{These functions assume that all missing values have been removed from the data set.  This can be done using the \code{na.omit()} command before the outcome and mediation models are fitted.}

\value{\code{medsens} returns an object of class "\code{medsens}". Some of the below elements are not available depending on whether INT is
     specified as TRUE or FALSE by the user or depending on the type of model fit. The function \code{summary} is used to obtain a table of the results. The function \code{plot} is used to plot the results.

      \item{d0}{vector of point estimates for mediation effect under control.}
      \item{d1}{vector of point estimates for mediation effect under treatment.}
\item{upper.d0}{vector of upper confidence intervals for mediation effect under control.}
\item{lower.d0}{vector of lower confidence intervals for mediation effect under control.}
\item{upper.d1}{vector of upper confidence intervals for mediation effect under
          treatment.}
\item{lower.d1}{vector of lower confidence intervals for mediation effect under
          treatment.}
     \item{tau}{vector of point estimates for total effect.}
\item{upper.tau}{vector of upper confidence intervals for total effect.}
\item{lower.tau}{vector of lower confidence intervals for total effect.}
      \item{nu}{vector of estimated proportions of total effect mediated.}
\item{upper.nu}{vector of upper confidence intervals for proportion mediated.}
\item{lower.nu}{vector of lower confidence intervals for proportion mediated.}
\item{rho}{vector of sensitivity parameter rho used.}
\item{rho.by}{increment of rho used.}
\item{err.cr}{value of rho corresponding to the mediation effect of zero.}
\item{R2star.prod}{value of the product of R square terms for the proportion of additional variance left unexplained, corresponding to value of rho.}
\item{R2tilda.prod}{value of the product of R square terms for the proportion of total variance left unexplained, corresponding to value of rho.}
\item{R2star.thresh}{value of R square for the proportion of additional variance left unexplained which corresponds to the mediation effect of zero.}
\item{R2tilda.thresh}{value of R square for the proportion of total variance left unexplained which corresponds to the mediation effect of zero.}
\item{r.square.y}{R square for the outcome model.}
\item{r.square.m}{R square for the mediator model.}
\item{INT}{the INT argument used in the original mediate object.}
\item{type}{type of the mediator and outcome models used. Currently either 'ct', 'bm' or 'bo'.}
}

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2009) A General Approach to Causal Mediation Analysis.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2009) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects.
Imai, Kosuke, Luke Keele, Dustin Tingley and Teppei Yamamoto (2009) Causal Mediation Analysis in R.} 

\author{Luke Keele, Ohio State University, \email{keele.4@osu.edu} , Dustin Tingley, Princeton University, \email{dtingley@princeton.edu}, Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{See also \code{\link{mediate}} }

\examples{ 

#Example with JOBS II Field experiment

#**For illustration purposes simulations set to low number**

data(jobs)

#########################################
#continuous mediator and continuous outcome
#########################################

#Fit parametric models
model.m <- lm(job_seek ~ treat + depress1, data=jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1, data=jobs)

#Pass model objects through mediate function
med.cont <- mediate(model.m, model.y, treat="treat", mediator="job_seek", sims=10)

#Pass mediate output through medsens function
sens.cont <- medsens(med.cont, rho.by=.2, eps=.01)

#Use summary function to display values of rho and R^2 where 95% confidence interval overlaps with 0.
summary(sens.cont)

#Plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))

#Plot sensitivity analysis using R^2 method. See plot.medsens for additional detail
plot(sens.cont, sens.par="R2", r.type=2, sign.prod=1)

} 
