\name{mediate.ped} \alias{mediate.ped} \title{Parallel Encouragement Design Bounds} \description{ 'mediate.ped' estimates the bounds on causal mediation effects for the parallel design. If a treatment/mediator interaction is allowed then the bounds are calculated. If a treatment/mediator interaction is not allowed then the mediation effect is point identified in which case the point estimate and bootstrapped confidence intervals are estimated. See Imai, Tingley, and Yamamoto JRSS-A for additional details.
} \usage{

mediate.ped(outcome,mediator,treatment,encouragement) {

}

\arguments{
 \item{outcome}{vector of experiment outcome variable}
 \item{mediator}{vector of experiment mediating variable}
 \item{treatment}{vector of the treatment variable}
 \item{encouragement}{vector of the indicator of whether observation received manipulation.}
}

\details{This function calculates causal mediation effects (CME) for the parallel encouragment design described in Imai, Tingley, Yamamoto JRSS-A. In the design two experimental arms are used. In one the mediator is not manipulated and in the other the mediator is encouraged in either a positive or negative direction. 

%Teppei, should we just refer people to the papers for the expressions. Doesn't seem too common for man files to have all the expressions.
%The CME represents the expected difference in the potential outcome when the mediator took the value that would realize under the treatment condition as opposed to the control condition, while the treatment status itself is held constant. \deqn{\delta(t) \ = \ E\{Y(t, M(t_1)) - Y(t, M(t_0))\},}{%
%    \delta(t) = E[Y(t, M(t1)) - Y(t, M(t0))],}
%where \eqn{t, t_1, t_0}{t, t1, t0} are particular values of the treatment \eqn{T} such that \eqn{t_1 \neq t_0}{t1 != t0}, \eqn{M(t)} is the potential mediator, and \eqn{Y(t,m)} is the potential outcome variable.

The function requires the outcome, mediator, treatment, and manipulation variables to be binary (0/1). The value of the manipulation indicator must be trichotomous (-1 (discouragement)/0 (no encouragement) /1 (encouragement)) in the parallel encouragement design.
}


\value{\code{mediate.ped} returns an object of class "\code{mediate.design}", a list that contains the components listed below.

The function \code{summary} (i.e., \code{summary.mediate.design}) can be used to obtain a table of the results.

\item{BE.c, BE.t}{lower/upper bounds for population causal mediation effects under the control and treatment conditions for the parallel encouragement experiment design.}
\item{BET.c, BET.t}{lower/upper bounds for complier causal mediation effects under the control and treatment conditions for the parallel encouragement experiment design.}
\item{manip}{indicator for encouragement status of mediator.}
\item{nobs}{number of observations used.}
\item{design}{"PED"}
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

# **For illustration purposes a small number of simulations are used and a randomly generate data set.**

#Create fake data
D<-rep(seq(-1,1, by=1), 10)#encouragment condition under PED
n<-length(D)
Y<-rbinom(n,1,.5)
M1<-rbinom(n,1,.5)
M2<-rbinom(n,1,.5)
T<-rbinom(n,1,.5)
V<-rbinom(n,1,.5)

#################################################### # Example 1: Parallel Design Experiment#################################################### 
pd.exp<-mediate.ped(Y,M1,T,D)
summary(pd.exp)

} 
