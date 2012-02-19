\name{mediate.pd} 
\alias{mediate.pd} 
\title{Parallel Design Bounds} 
\description{ 'mediate.pd' estimates the bounds on causal mediation effects for the parallel design. If a treatment/mediator interaction is allowed then the bounds are calculated. If a treatment/mediator interaction is not allowed then the mediation effect is point identified in which case the point estimate and bootstrapped confidence intervals are estimated. See Imai, Tingley, and Yamamoto JRSS-A for additional details.
} 

\usage{

mediate.pd(outcome, mediator, treat, manipulated, data,
           NINT = TRUE, sims = 1000, conf.level = 0.95)

}

\arguments{
 \item{outcome}{vector of experiment outcome variable}
 \item{mediator}{vector of experiment mediating variable}
 \item{treatment}{vector of the treatment variable}
 \item{manipulated}{vector of the indicator of whether observation received mediator manipulation.}
 \item{NINT}{logical for whether user is assuming no interaction between treatment and mediator. If TRUE then no interaction is assumed}
  \item{sims}{ number of Monte Carlo draws for nonparametric bootstrap.}
 \item{conf.level}{ level of the returned two-sided confidence intervals.  Default is to return
    the 2.5 and 97.5 percentiles of the simulated quantities.}
}

\details{This function calculates causal mediation effects (CME) for the parallel design described in Imai, Tingley, Yamamoto JRSS-A.
In the design two experimental arms are used. In one the treatment is randomized and the mediator and outcome variables are measured. In the second design the treatment is randomized, and the other the mediator is perfectly manipulated and the mediator and outcome variables measured.

The function requires the outcome, mediator, treatment, and manipulation variables to be binary (0/1). Users must specify whether or not to allow a treatment/mediator interaction.
}


\value{\code{mediate.pd} returns an object of class "\code{mediate.design}", a list that contains the components listed below.

The function \code{summary} (i.e., \code{summary.mediate.design}) can be used to obtain a table of the results.

\item{d0, d1}{point estimates or lower/upper bounds for causal mediation effects (no treratment/mediator interaction) or lower/upper bounds for causal mediation effects under the control and treatment conditions (allows treatment/mediator interaction)}
\item{d0.ci, d1.ci}{confidence intervals for average causal mediation effects (no treratment/mediator interaction). The confidence level is set at the value specified in 'conf.level'.}
\item{manip}{indicator for whether unit received experimental manipulation of mediator.}
\item{nobs}{number of observations used.}
\item{conf.level}{confidence level used for two-sided confidence intervals.}
\item{sims}{number of bootstrap simulations used for confidence interval calculation (no treratment/mediator interaction).}
\item{design}{"PD.NINT" if no interaction assumed, "PD" if interaction allowed.}
}

\references{
Imai, K., Tingley, D. and Yamamoto, T. (2012) Experimental Designs for Identifying Causal Mechanisms. Journal of the Royal Statistical Society, Series A (Statistics in Society)

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the Black Box of Causality: Learning about Causal Mechanisms from Experimental and Observational Studies, American Political Science Review, Vol. 105, No. 4 (November), pp. 765-789.

Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.}

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, Teppei Yamamoto, MIT, \email{teppei@mit.edu}, }

\seealso{\code{\link{mediate}}, \code{\link{medsens}}, \code{\link{plot.mediate}}, \code{\link{summary.mediate}}, \code{\link{mediations}} }

\examples{ # Examples with simulated data

#################################################### # Example 1: Parallel Design Experiment####################################################

data(boundsdata)
bound2<-mediate.pd(boundsdata$Y,boundsdata$M,boundsdata$T,boundsdata$D,NINT=TRUE,sims=5,conf.level=.95)
summary(bound2)


bound2.1<-mediate.pd(boundsdata$Y,boundsdata$M,boundsdata$T,boundsdata$D,NINT=FALSE)
summary(bound2.1)

}
