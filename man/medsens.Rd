\name{medsens} 
\alias{medsens} 
\title{Causal Mediation Analysis - Sensitivity Analysis} 
\description{ 
Function to perform sensitivity analysis on mediation effect for violations of sequential ignorability assumption. This allows for a correlation between the error terms of the outcome model and the mediator model. Sensitivity analysis is possible with 1) continuous mediator and continuous outcome, 2) binary outcome and continuous mediator, and 3) continuous outcome and binary mediator. Output from the function can be passed through summary or plot functions which display estimated mediation effects for given values of rho.
} 
\usage{

medsens(model.m, model.y, T="treat.name", M="med.name", INT=FALSE, DETAIL=FALSE, sims=1000, eps=.Machine$double.eps)

} 

\arguments{ 
\item{model.m}{R model object for mediator.  Can be of class lm or glm using a probit link function.} 
\item{model.y}{R model object for outcome.  Can be of class lm or glm using a probit link function.} 
\item{INT}{If TRUE this indicates that treatment is interacted with mediator in ymodel object. Default is FALSE. Interaction implemented by including a T:M interaction term in the outcome model. This option is currently only available for continuous outcome continuous mediator case.} 
\item{T}{Name of binary treatment indicator, 1 for treatment and 0 if control. Variable name listed in quotations.}
\item{M}{Name of mediator. Variable name listed in quotations.}
\item{DETAIL}{If FALSE estimation is on .1 grid for rho (sensitivity parameter) vs. a grid of .01). Due to computational time FALSE is recommended.}
\item{sims}{Number of draws for parametric bootstrap used to calculate 95 percent confidence intervals.}
\item{eps}{Convergence parameter for FGLS estimation.}
} 

\details{This is the workhorse function for estimating sensitivity analyses for mediation effects. In fitting models for a binary variable a probit function must be used. In fitting the model.y object with a treatmentXmediator interaction the user must input the interation in the format treatment:mediator. Furthermore, the user must make sure to also include the main effect for the mediator in this case.  

Users should note that computation time is several minutes for these functions. Setting DETAIL=TRUE significantly increases computational time, as does increasing the number of simulations.}

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
#pass model objects through medsens function
sens.cont <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE, sims=1000)
#use summary function to display values of rho where 95% confidence interval overlaps with 0.
summary(sens.cont)
#plot mediation effect and 95% CI's for each value of rho
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))

\dontrun{
#########################################
#binary outcome and continuous mediator
#########################################

model.m <- lm(job_seek ~ treat + depress1, data=jobs)
model.y <- glm(work1 ~ treat + job_seek + depress1, family=binomial(link="probit"), data=jobs)
sens.dichO <- medsens(model.m, model.y, T="treat", M="job_seek", INT=FALSE,  DETAIL=FALSE)
summary(sens.dichO)
plot(sens.dichO, main="JOBS", ylim=c(-.2,.2))

#########################################
#continuous outcome and binary mediator w no interaction
#########################################

model.m <- glm(job_dich ~ treat + depress1, data=jobs, family=binomial(link="probit"))
model.y <- lm(depress2 ~ treat + job_dich + depress1, data=jobs)
sens.dichM <- medsens(model.m, model.y, T="treat", M="job_dich", INT=FALSE,  DETAIL=FALSE)
summary(sens.dichM)
plot(sens.dichM, main="JOBS", ylim=c(-.2,.2))

#########################################
#continuous outcome and binary mediator w interaction
#########################################

model.m <- glm(job_dich ~ treat + depress1, data=jobs, family=binomial(link="probit"))
model.y <- lm(depress2 ~ treat + job_dich + treat:job_dich + depress1, data=jobs)
sens.dichM.int <- medsens(model.m, model.y, T="treat", M="job_dich", INT=TRUE,  DETAIL=FALSE)
summary(sens.dichM.int)
plot(sens.dichM.int, main="JOBS", ylim=c(-.2,.2))
}

} 
