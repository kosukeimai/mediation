\name{mediate.pd} \alias{mediate.pd} \title{Parallel Design Bounds} \description{ 'mediate.pd' estimates the bounds on causal mediation effects for the parallel design. If a treatment/mediator interaction is allowed then the bounds are calculated. If a treatment/mediator interaction is not allowed then the mediation effect is point identified in which case the point estimate and bootstrapped confidence intervals are estimated. See Imai, Tingley, and Yamamoto JRSS-A for additional details.
} \usage{

mediate.ced(outcome,mediator,treatment,manipulated,NINT=TRUE,conf.level=.95) {

}

\arguments{
 \item{outcome}{vector of experiment outcome variable}
 \item{mediator}{vector of experiment mediating variable}
 \item{treatment}{vector of the treatment variable}
 \item{manipulated}{vector of the indicator of whether observation received manipulation.}
 \item{NINT}{logical for whether user is assuming no interaction between treatment and mediator. If TRUE then no interaction is assumed}
 \item{conf.level}{ level of the returned two-sided confidence intervals.  Default is to return
    the 2.5 and 97.5 percentiles of the simulated quantities.}
}

\details{This function calculates causal mediation effects (CME) for the parallel design described in Imai, Tingley, Yamamoto JRSS-A. In the design two experimental arms are used. In one the mediator is not manipulated and in the other the mediator is perfectly manipulated.

%Teppei, should we just refer people to the papers for the expressions. Doesn't seem too common for man files to have all the expressions.
%The CME represents the expected difference in the potential outcome when the mediator took the value that would realize under the treatment condition as opposed to the control condition, while the treatment status itself is held constant. \deqn{\delta(t) \ = \ E\{Y(t, M(t_1)) - Y(t, M(t_0))\},}{%
%    \delta(t) = E[Y(t, M(t1)) - Y(t, M(t0))],}
%where \eqn{t, t_1, t_0}{t, t1, t0} are particular values of the treatment \eqn{T} such that \eqn{t_1 \neq t_0}{t1 != t0}, \eqn{M(t)} is the potential mediator, and \eqn{Y(t,m)} is the potential outcome variable.

The function requires the outcome, mediator, treatment, and manipulation variables to be binary (0/1). Users must specify whether or not to allow a treatment/mediator interaction.
}


\value{\code{mediate.ced} returns an object of class "\code{mediate.design}", a list that contains the components listed below.

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
Imai, K., Tingley, D. and Yamamoto, T. (2012) Experimental Designs for Identifying Causal Mechanisms. Journal of the Royal Statistical Society, Series A (Statistics in Society)"

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the Black Box of Causality: Learning about Causal Mechanisms from Experimental and Observational Studies, American Political Science Review, Vol. 105, No. 4 (November), pp. 765-789.

Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.}

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, Teppei Yamamoto, MIT, \email{teppei@mit.edu}, }

\seealso{\code{\link{mediate}}, \code{\link{medsens}}, \code{\link{plot.mediate}}, \code{\link{summary.mediate}}, \code{\link{mediations}} }

\examples{ # Examples with simulated data

#################################################### # Example 1: Parallel Design Experiment####################################################

data(boundsdata)
bound2<-mediate.pd(data$Y,data$M,data$T,data$D)
summary(bound2)


}
