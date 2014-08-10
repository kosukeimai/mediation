\name{test.TMint}
\alias{test.TMint}
\alias{test.TMint.default}
\alias{test.TMint.mediate}
\alias{test.TMint.mediate.order}
\alias{print.htest.order}
\title{Significance Test for Treatment-Mediator Interaction in Causal Mediation Analysis}
\description{
Function to test whether the average causal mediation effects and direct effects are significantly different between the treatment and control contitions.
}
\usage{
\method{test.TMint}{mediate}(x, conf.level = x$conf.level, ...)
}

\arguments{
\item{x}{ output from \code{mediate} function.}
\item{conf.level}{ level of the returned two-sided confidence intervals for the effect differences. By default it is set to the value used in the original mediate call.}
\item{...}{ additional arguments.}
}

\value{\code{test.TMint} returns an object of class "\code{htest}" when applied to a \code{mediate} object. See \code{\link{t.test}} for more explanations of the contents. The function returns an object of class "\code{htest.order}" which has its own \code{print} method included in this package.
}

\references{
Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. (2014). "mediation: R package for Causal Mediation Analysis", Journal of Statistical Software, Vol. 59, No. 5, pp. 1-38. 

Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.}

\author{Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}.}

\seealso{\code{\link{mediate}}}

\examples{ 

# Examples with JOBS II Field Experiment

# **For illustration purposes a small number of simulations are used**

data(jobs)

# Fit mediator and outcome models allowing for treatment-mediator interaction
b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
d <- lm(depress2 ~ treat*job_seek + econ_hard + sex + age, data=jobs)

# Test for significance of interaction
fit <- mediate(b, d, sims=50, treat="treat", mediator="job_seek")
test.TMint(fit)

}

