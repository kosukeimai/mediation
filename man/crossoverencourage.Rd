\name{crossoverencourage} \alias{crossoverencourage} \title{Crossover Encouragement Design} \description{ 'crossoverencourage' estimates the causal mediation effects for the crossover encouragement design. Assumptions for the point identification of the ACME are discussed in Imai, Tingley, Yamamoto JRSS-A 2012.
} \usage{

crossoverencourage(outcome=outcome,mediator1=mediator1,mediator2=mediator2,treatment=treatment,encouragement=encouragement,sims=1000, conf.level=.95) 

}

\arguments{
 \item{outcome}{vector of experiment outcome variable.}
 \item{mediator1}{vector of experiment mediating variable in the first stage. }
 \item{mediator2}{vector of experiment mediating variable in the second stage. }
 \item{treatment}{vector of the treatment variable}
 \item{encouragement}{vector of the indicator of whether observation was encouraged to take on the same value of the mediator as in the first period (1) or no encouragment (0).}
 \item{sims}{ number of Monte Carlo draws for nonparametric bootstrap.} 

}

\details{This function calculates causal mediation effects (CME) for the parallel encouragement design. The CME represents the expected difference 

\deqn{\delta(t) \ = \ E\{Y(t, M(t_1)) - Y(t, M(t_0))\},}
#in the potential outcome when the mediator took the value that would realize under the treatment condition as opposed to the control condition, while the treatment status itself is held constant. \deqn{\delta(t) \ = \ E\{Y(t, M(t_1)) - Y(t, M(t_0))\},}{%
#    \delta(t) = E[Y(t, M(t1)) - Y(t, M(t0))],}
#where \eqn{t, t_1, t_0}{t, t1, t0} are particular values of the treatment \eqn{T} such that \eqn{t_1 \neq t_0}{t1 != t0}, \eqn{M(t)} is the potential mediator, and \eqn{Y(t,m)} is the potential outcome variable.

The function requires the treatment, mediator, and outcome variables to be binary (0/1). The value of the encouragement indicator must be binary (0 (no encouragement)/ 1 (encouragement)).
}


\value{\code{crossoverencourage} returns an object of class "\code{crossoverencourage}", a list that contains the components listed below.

The function \code{summary} (i.e., \code{summary.crossoverencourage} can be used to obtain a table of the results.

\item{d.p.c.mu , d.p.t.mu }{point estimates for average causal mediation effects under the control and treatment conditions.}
\item{d.p.c.ci, d.p.t.ci}{confidence intervals for average causal mediation effects. The confidence level is set at the value specified in 'conf.level'.}
\item{d.p.c,d.p.t}{a vector of length 'sims' containing simulation draws of average causal mediation effects under the control and treatment conditions.}
\item{conf.level}{the confidence level used. }
\item{sims}{number of bootstrap simulations used.}
\item{nobs}{number of observations used for computation.}
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
#Crossover encouragement test
n<-1000
Y<-rbinom(n,1,.5)
M1<-rbinom(n,1,.5)
M2<-rbinom(n,1,.5)
T<-rbinom(n,1,.5)
V<-rbinom(n,1,.5)

####Calculate Quantities of Interest####
ce<-crossoverencourage(Y,M1,M2,T,V,sims=100)
summary(ce)



} 
