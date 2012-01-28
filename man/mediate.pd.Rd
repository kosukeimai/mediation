\name{mediate.pd} \alias{mediate.pd} \title{Parallel Design Bounds} \description{ 'mediate.pd' estimates the bounds on causal mediation effects for the parallel design allowing for a treatment/mediator interaction. Additional assumptions for the bounds under each design are discussed in Imai, Tingley, and Yamamoto JRSS-A.
} \usage{

mediate.pd(outcome=outcome,mediator,treatment,manipulated) {

}

\arguments{
 \item{outcome}{vector of experiment outcome variable.}
 \item{mediator}{vector of experiment mediating variable. }
 \item{treatment}{vector of the treatment variable}
 \item{manipulated}{vector of the indicator of whether observation received manipulation. Null for single experiment design.}
}

\details{This function calculates causal mediation effects (CME) for three different types of experimental designs. The CME represents the expected difference in the potential outcome when the mediator took the value that would realize under the treatment condition as opposed to the control condition, while the treatment status itself is held constant. \deqn{\delta(t) \ = \ E\{Y(t, M(t_1)) - Y(t, M(t_0))\},}{%
    \delta(t) = E[Y(t, M(t1)) - Y(t, M(t0))],}
where \eqn{t, t_1, t_0}{t, t1, t0} are particular values of the treatment \eqn{T} such that \eqn{t_1 \neq t_0}{t1 != t0}, \eqn{M(t)} is the potential mediator, and \eqn{Y(t,m)} is the potential outcome variable.

The function requires the treatment, mediator, and outcome variables to be binary (0/1). The value of the manipulation indicator must be binary (0 (no encouragement)/ 1 (encouragement)).
}


\value{\code{mediate.pd} returns an object of class "\code{mediate.design}", a list that contains the components listed below.

The function \code{summary} (i.e., \code{summary.mediate.design} can be used to obtain a table of the results.

\item{d0, d1}{lower/upper bounds for causal mediation effects under the control and treatment conditions for the single experiment design.}
%\item{BE.c, BE.t}{lower/upper bounds for population causal mediation effects under the control and treatment conditions for the parallel encouragement experiment design.}
%\item{BET.c, BET.t}{lower/upper bounds for complier causal mediation effects under the control and treatment conditions for the parallel encouragement experiment design.}
\item{manip}{indicator for whether experimental manipulation of mediator.}
\item{nobs}{number of observations used for computation of bounds.} }

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
pd.exp<-mediate.pd(Y,M1,T,D)
summary(pd.exp)

} 