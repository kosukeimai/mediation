\name{mediate.sed} \alias{mediate.sed} \title{mediate.sed} \description{ 'mediate.sed' estimates causal mediation effects for the single experiment design. The two options are to use either the sequential ignorability (SI) assumption in which nonparametric estimates of the average causal mediation effect are produced, or, to relax the SI assumption and to calculate the bounds on the causal mediation effect. Assumptions for the bounds are discussed in Imai, Tingley, Yamamoto JRSS-A.
} \usage{

mediate.sed(outcome,mediator,treatment,SI=FALSE)

}

\arguments{
 \item{outcome}{vector of experiment outcome variable. Must be binary (0/1) if SI is not made.}
 \item{mediator}{vector of experiment mediating variable. Must be binary if SI assumption is not made. Must be discrete if SI assumption is made.}
 \item{treatment}{vector of the treatment variable. Must be binary (0/1).}
 \item{SI}{Logical indicator for whether sequential ignorability assumption is made.}
 \item{boot}{ a logical value. if 'FALSE' a large sample Delta method approximation is used for confidence intervals; if 'TRUE' nonparametric bootstrap will be used. Default is 'FALSE'.}
\item{sims}{ number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.}
\item{conf.level}{ level of the returned two-sided confidence intervals.  Default is to return
the 2.5 and 97.5 percentiles of the simulated quantities.}
}

\details{This function calculates causal mediation effects (CME) for the single experiment design when the user wants either the non-parametric estimates based on the SI assumption or the bounds when the SI assumption is not made. The CME represents the expected difference in the potential outcome when the mediator took the value that would realize under the treatment condition as opposed to the control condition, while the treatment status itself is held constant. \deqn{\delta(t) \ = \ E\{Y(t, M(t_1)) - Y(t, M(t_0))\},}{%
    \delta(t) = E[Y(t, M(t1)) - Y(t, M(t0))],}
where \eqn{t, t_1, t_0}{t, t1, t0} are particular values of the treatment \eqn{T} such that \eqn{t_1 \neq t_0}{t1 != t0}, \eqn{M(t)} is the potential mediator, and \eqn{Y(t,m)} is the potential outcome variable.

The experimental design covered by the function the single experiment design where only the treatment is randomized.
The function requires the treatment, mediator, and outcome variables to be binary (0/1).

}


\value{\code{mediate.sed} returns an object of class "\code{mediate.design}", a list that contains the components listed below.

The function \code{summary} (i.e., \code{summary.mediate.design} can be used to obtain a table of the results.

\item{d0, d1}{point estimates or lower/upper bounds for causal mediation effects under the control and treatment conditions for the single experiment design.}
\item{d0.ci, d1.ci}{confidence intervals for average causal mediation effects for the non-parametric estimates. The confidence level is set at the value specified in 'conf.level'.}
\item{boot}{logical, the 'boot' argument used.}
\item{conf.level}{the confidence level used. }
\item{sims}{number of bootstrap simulations used for confidence interval calculation.}
\item{nobs}{number of observations used.} }

\references{
Imai, K., Tingley, D. and Yamamoto, T. (2012) Experimental Designs for Identifying Causal Mechanisms. Journal of the Royal Statistical Society, Series A (Statistics in Society)"

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the Black Box of Causality: Learning about Causal Mechanisms from Experimental and Observational Studies, American Political Science Review, Vol. 105, No. 4 (November), pp. 765-789.

Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.}

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, Luke Keele, PSU, \email{ljk20@psu.edu}, Teppei Yamamoto, MIT, \email{teppei@mit.edu}, }

\seealso{\code{\link{mediate}}, \code{\link{medsens}}, \code{\link{plot.mediate}}, \code{\link{summary.mediate}}, \code{\link{mediations}} }

\examples{ #


#################################################### # Example 1: Bounds without SI assumption####################################################
data(boundsdata)

data.temp<-subset(boundsdata,boundsdata$D==0)#Only use data from no encouragement set.
bound1<-mediate.sed(data.temp$Y,data.temp$M,data.temp$T,SI=FALSE)
summary(bound1)


#################################################### # Example 2: Nonparametric estimate of ACME under SI assumption####################################################
# Example with JOBS II Field Experiment
data(jobs)
#Use the discrete version of the mediator.

foo.1 <- mediate.sed(job_disc, depress2, treat, SI=TRUE)
summary(foo.1)

foo.2 <- mediate.sed(job_disc, depress2, treat, SI=TRUE, boot=TRUE)
summary(foo.2)


}
