\name{mediate.ped} 
\alias{mediate.ped} 
\title{Computing Bounds on Average Causal Mediation Effects under the Parallel Encouragement Design} 

\description{
'mediate.ped' computes the nonparametric bounds on the average causal mediation effects for the parallel encouragement design.
} 

\usage{

mediate.ped(outcome, mediator, treat, encourage, data)

}

\arguments{

\item{outcome}{ name of the outcome variable in 'data'.}
\item{mediator}{ name of the mediator in 'data'. The variable must be binary (factor or numeric 0/1).}
\item{treat}{ name of the treatment variable in 'data'. Must be binary (factor or numeric 0/1).}
\item{encourage}{ name of the encouragement variable in 'data'. The variable must be a numeric vector taking on either -1, 0, or 1.}
\item{data}{ a data frame containing all the above variables.}

}

\details{
This function calculates average causal mediation effects (ACME) for the parallel encouragement design.
In the design two experimental arms are used. In one the treatment is randomized and the mediator and outcome variables are measured. In the second arm the treatment is randomized, the mediator is randomly encouraged either up or down, and the outcome variable is measured.

Two type of causal quantities are estimated: the population ACME and the complier ACME. The latter refers to the subpopulation of the units for whom the encouragement has its intended effect, and the width of its bounds are tighter than that of the population ACME. See Imai, Tingley and Yamamoto (2012) for details.
}


\value{\code{mediate.pd} returns an object of class "\code{mediate.design}", a list that contains the components listed below.

The function \code{summary} (i.e., \code{summary.mediate.design}) can be used to obtain a table of the results.

\item{d0, d1}{ estimated nonparametric sharp bounds for the population ACME under the control and treatment conditions.}
\item{d0.p, d1.p}{ estimated nonparametric sharp bounds for the complier ACME under the control and treatment conditions.}
\item{nobs}{ number of observations used.}
\item{design}{ indicates the design. Always equals "PED".}
}

\references{
Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. (2014). "mediation: R package for Causal Mediation Analysis", Journal of Statistical Software, Vol. 59, No. 5, pp. 1-38. 

Imai, K., Tingley, D. and Yamamoto, T. (2012) Experimental Designs for Identifying Causal Mechanisms. Journal of the Royal Statistical Society, Series A (Statistics in Society)"

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the Black Box of Causality: Learning about Causal Mechanisms from Experimental and Observational Studies, American Political Science Review, Vol. 105, No. 4 (November), pp. 765-789.

Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.}

\author{
Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu};
Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}.}

\seealso{\code{\link{mediate}}, \code{\link{medsens}}, \code{\link{plot.mediate}}, \code{\link{summary.mediate}}, \code{\link{mediations}} }

\examples{

data(boundsdata)

bound3 <- mediate.ped("out.enc", "med.enc", "ttt", "enc", boundsdata)
summary(bound3)

}

