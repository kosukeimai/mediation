\name{summary.ivmediate}
\alias{summary.ivmediate}
\alias{print.summary.ivmediate}

\title{Summarizing Output from Mediation Analysis with Treatment Noncompliance}
\description{
Function to report results from mediation analysis with treatment noncompliance. Reported categories are local average causal mediation effects, local average natural direct effects and local average treatment (total) effect.
}
\usage{

  \method{summary}{ivmediate}(object, conf.level = object$conf.level[1], ...)

  \method{print}{summary.ivmediate}(x, ...)

}

\arguments{
\item{object}{ output from mediate function.}
\item{conf.level}{ confidence level for the intervals reported in the summary table.}
\item{x}{ output from summary.mediate function.}
\item{...}{ additional arguments affecting the summary produced.}
}

\references{
Yamamoto, T. (2013). Identification and Estimation of Causal Mediation Effects with Treatment Noncompliance. Unpublished manuscript.
} 

\author{Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}.}
 
\seealso{\code{\link{ivmediate}}, \code{\link{plot.ivmediate}}}

\examples{ 
# Examples with JOBS II Field Experiment

# ** For illustration purposes a small number of simulations are used **

require(parallel)
require(MASS)

data(jobs)

a <- lm(comply ~ treat + sex + age + marital + nonwhite + educ + income, 
        data = jobs)
b <- glm(job_dich ~ comply + treat + sex + age + marital + nonwhite + educ + income, 
        data = jobs, family = binomial)
c <- lm(depress2 ~ job_dich * (comply + treat) + sex + age + marital + nonwhite + educ + income, 
        data = jobs)

out <- ivmediate(a, b, c, sims = 50, boot = FALSE,
                 enc = "treat", treat = "comply", mediator = "job_dich",
                 mc.cores = min(2, getOption("mc.cores", 2L)))
                 
summary(out)
plot(out)

} 

