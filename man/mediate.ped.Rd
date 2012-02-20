\name{mediate.ped} 
\alias{mediate.ped} 
\title{Computing Bounds on Average Causal Mediation Effects under the Parallel Encouragement Design} 

\description{
'mediate.ped' computes the nonparametric bounds on average causal mediation effects for the parallel encouragement design. See Imai, Tingley, and Yamamoto JRSS-A for additional details.
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

\details{This function calculates causal mediation effects (CME) for the parallel encouragment design described in Imai, Tingley, Yamamoto JRSS-A.
In the design two experimental arms are used. In one the treatment is randomized and the mediator and outcome variables are measured. In the second design the treatment is randomized, and the other the mediator is encouraged either up, down, or no encouragement. The mediator and outcome variables are measured.

The function requires the outcome, mediator, and treatment to be binary (0/1). The value of the manipulation indicator must be trichotomous (-1 (discouragement)/0 (no encouragement) /1 (encouragement)) in the parallel encouragement design.
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

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, Teppei Yamamoto, MIT, \email{teppei@mit.edu} }

\seealso{\code{\link{mediate}}, \code{\link{medsens}}, \code{\link{plot.mediate}}, \code{\link{summary.mediate}}, \code{\link{mediations}} }

\examples{

data(boundsdata)

bound3 <- mediate.ped("out.enc", "med.enc", "ttt", "enc", boundsdata)
summary(bound3)

}

