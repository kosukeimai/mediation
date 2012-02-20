\name{mediate.ced} 
\alias{mediate.ced} 
\title{Estimating Average Causal Mediation Effects under the Crossover Encouragement Design} 

\description{ 'mediate.ced' estimates the average causal mediation effects for the crossover encouragement design.
} 

\usage{

mediate.ced(outcome, med.1, med.2, treat, encourage, data,
            sims = 1000, conf.level = .95)

}

\arguments{

\item{outcome}{ name of the outcome variable in 'data'. The variable must be binary (factor or numeric 0/1).}
\item{med.1}{ variable name in 'data' containing the mediator values observed in the first experiment. The variable must be binary (factor or numeric 0/1).}
\item{med.2}{ variable name in 'data' containing the mediator values observed in the second experiment.}
\item{treat}{ name of the treatment variable in 'data'. Must be binary (factor or numeric 0/1).}
\item{encourage}{ name of the encouragement indicator in 'data'. Must be binary (factor or numeric 0/1).}
\item{data}{ a data frame containing all the above variables.}
\item{sims}{ number of bootstrap simulations.}
\item{conf.level}{ level of the returned two-sided confidence intervals.}

}

\details{This function calculates causal mediation effects (CME) for the crossover encouragement design described in Imai, Tingley, Yamamoto JRSS-A.
The design has two stages. In the first experiment the treatment is randomized and the mediator/outcome variables measured.
In the second the treatment is set to value opposite of first period and a randomly selected group of units receives encouragement to take on the mediator observed value in period 1.
The function requires the outcome, both mediator, treatment, and manipulation variables to be binary (0/1).
}


\value{\code{mediate.ced} returns an object of class "\code{mediate.design}", a list that contains the components listed below.

The function \code{summary} (i.e., \code{summary.mediate.design}) can be used to obtain a table of the results.

\item{d0, d1}{point estimates or lower/upper bounds for causal mediation effects under the control and treatment conditions.}
\item{d0.ci, d1.ci}{confidence intervals for average causal mediation effects. The confidence level is set at the value specified in 'conf.level'.}
\item{manip}{indicator for whether unit received experimental manipulation of mediator.}
\item{nobs}{number of observations used.}
\item{conf.level}{confidence level used for two-sided confidence intervals.}
\item{sims}{number of bootstrap simulations used for confidence interval calculation.}
\item{design}{"CED"}

}

\references{
Imai, K., Tingley, D. and Yamamoto, T. (2012) Experimental Designs for Identifying Causal Mechanisms. Journal of the Royal Statistical Society, Series A (Statistics in Society)"

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the Black Box of Causality: Learning about Causal Mechanisms from Experimental and Observational Studies, American Political Science Review, Vol. 105, No. 4 (November), pp. 765-789.

Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.}

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, Teppei Yamamoto, MIT, \email{teppei@mit.edu}, }

\seealso{\code{\link{mediate}}, \code{\link{medsens}}, \code{\link{plot.mediate}}, \code{\link{summary.mediate}}, \code{\link{mediations}} }

\examples{ 

# Generate fake data
n <- 1000
Y <- rbinom(n,1,.5)
M1 <- rbinom(n,1,.5)
M2 <- rbinom(n,1,.5)
T <- rbinom(n,1,.5)
V <- rbinom(n,1,.5)
DD <- as.data.frame(cbind(Y, M1, M2, T, V))

# Apply the function
ced <- mediate.ced("Y", "M1", "M2", "T", "V", DD, sims = 100)
summary(ced)

}
