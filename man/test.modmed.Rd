\name{test.modmed}
\alias{test.modmed}
\alias{test.modmed.default}
\alias{test.modmed.mediate}
\alias{test.modmed.mediate.order}
\alias{print.test.modmed.mediate}
\alias{print.test.modmed.mediate.order}
\title{Significance Test for Moderated Mediation}
\description{
Function to test whether the average causal mediation effects and direct effects are significantly different between two moderator strata.
}
\usage{
\method{test.modmed}{mediate}(object, covariates.1, covariates.2,
                              sims = object$sims, conf.level = object$conf.level, ...)
\method{print}{test.modmed.mediate}(x, ...)
}

\arguments{
\item{object}{ output from \code{mediate} function.}
\item{covariates.1}{ first set of value(s) of covariates (moderators) included in the mediator and outcome models. See documentation for the \code{covariates} argument for the \code{\link{mediate}} function.}
\item{covariates.2}{ second set of value(s) of covariates (moderators) included in the mediator and outcome models.}
\item{sims}{ number of simulation draws the test will be based on. Defaults to the number used in the original mediate fit.}
\item{conf.level}{ level of the returned two-sided confidence intervals for the effect differences. By default it is set to the value used in the original mediate call.}
\item{x}{ output from \code{test.modmed} function.}
\item{...}{ additional arguments.}
}

\value{When applied to a \code{mediate} object, \code{test.modmed} returns an object of class "\code{test.modmed.mediate}", a list composed of "\code{htest}" objects. See \code{\link{t.test}} for more explanations of \code{htest} objects. When applied to a \code{mediate.order} object, the function returns an object of class "\code{test.modmed.mediate.order}" which is a list composed of "\code{htest.order}" objects.
}

\details{The function takes the original call to \code{mediate} and reruns the algorithm twice with the two sets of \code{covariates} values. It assumes that the objects in the environment in which the original mediate call was made also exist in the current environment under the same variable names, i.e., it evaluates the updated call in the current environment.
}

\references{Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.}

\author{Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}.}

\seealso{\code{\link{mediate}}, \code{\link{test.TMint}}}

\examples{ 

# Examples with JOBS II Field Experiment

# **For illustration purposes a small number of simulations are used**

data(jobs)

# Fit mediator and outcome models allowing for interaction with moderator
b.int <- lm(job_seek ~ treat*age + econ_hard + sex, data=jobs)
d.int <- lm(depress2 ~ treat*job_seek*age + econ_hard + sex, data=jobs)

# Initial mediate fit
fit <- mediate(b.int, d.int, sims=50, treat="treat", mediator="job_seek")

# Test for significance of moderated mediation
test.modmed(fit, list(age = 20), list(age = 70), sims = 100)

}

