\name{mediate} 
\alias{mediate} 
\title{Causal Mediation Analysis} 
\description{ 
'mediate' is used to estimate various quantities for causal mediation analysis, including average causal mediation effects (indirect effect), average direct effects, proportions mediated, and total effect.
} 

\usage{

mediate(model.m, model.y, sims = 1000, 
        boot = FALSE, boot.ci.type = "perc",
        treat = "treat.name", mediator = "med.name", 
        covariates = NULL, outcome = NULL,
        control = NULL, conf.level = .95, 
        control.value = 0, treat.value = 1, 
        long = TRUE, dropobs = FALSE,
        robustSE = FALSE, cluster = NULL, group.out = NULL, ...)

} 

\arguments{ 
\item{model.m}{ a fitted model object for mediator.  Can be of class 'lm', 'polr', 'bayespolr', 'glm', 'bayesglm', 'gam', 'rq', 'survreg', or 'merMod'.} 
\item{model.y}{ a fitted model object for outcome.  Can be of class 'lm', 'polr', 'bayespolr', 'glm', 'bayesglm', 'gam', 'vglm', 'rq', 'survreg', or 'merMod'.} 
\item{sims}{ number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.} 
\item{boot}{ a logical value. if 'FALSE' a quasi-Bayesian approximation is used for confidence intervals; if 'TRUE' nonparametric bootstrap will be used. Default is 'FALSE'.}
\item{boot.ci.type}{ a character string indicating the type of bootstrap confidence intervals. If "bca" and boot = TRUE, bias-corrected and accelerated (BCa) confidence intervals will be estimated. If "perc" and boot = TRUE, percentile confidence intervals will be estimated. Default is "perc".}
\item{conf.level}{ level of the returned two-sided confidence intervals.  Default is to return
the 2.5 and 97.5 percentiles of the simulated quantities.} 
\item{treat}{ a character string indicating the name of the treatment variable used in the models.  The treatment can be either binary (integer or a two-valued factor) or continuous (numeric).}
\item{mediator}{ a character string indicating the name of the mediator variable used in the models.}
\item{covariates}{ a list or data frame containing values for a subset of the pre-treatment covariates in 'model.m' and 'model.y'. If provided, the function will return the estimates conditional on those covariate values.}
\item{outcome}{ a character string indicating the name of the outcome variable in `model.y'. Only necessary if 'model.y' is of class 'survreg'; otherwise ignored.}
\item{control}{ a character string indicating the name of the control group indicator. Only relevant if 'model.y' is of class 'gam'. If provided, 'd0', 'z0' and 'n0' are allowed to differ from 'd1', 'z1' and 'n1', respectively.}
\item{control.value}{ value of the treatment variable used as the control condition. Default is 0.}
\item{treat.value}{ value of the treatment variable used as the treatment condition. Default is 1.}
\item{long}{ a logical value. If 'TRUE', the output will contain the entire sets of simulation draws of the the average causal mediation effects, direct effects, proportions mediated, and total effect. Default is 'TRUE'.}
\item{dropobs}{ a logical value indicating the behavior when the model frames of 'model.m' and 'model.y' (and the 'cluster' variable if included) are composed of different observations. If 'TRUE', models will be re-fitted using common data rows. If 'FALSE', error is returned. Default is 'FALSE'.}
\item{robustSE}{ a logical value. If 'TRUE', heteroskedasticity-consistent standard errors will be used in quasi-Bayesian simulations. Ignored if 'boot' is 'TRUE' or neither 'model.m' nor 'model.y' has a method for \code{vcovHC} in the \code{sandwich} package. Default is 'FALSE'.}
\item{cluster}{ a variable indicating clusters for standard errors. Note that this should be a vector of cluster indicators itself, not a character string for the name of the variable.}
\item{group.out}{ a character string indicating the name of the lmer/glmer group on which the mediate output is based. Can be used even when a merMod function is applied to only one of the mediator or the outcome. If merMod functions are applied to both the mediator and the outcome, default is the group name used in the outcome model; if the mediator group and the outcome group are different and the user is interested in the mediate output based on the mediator group, then set group.out to the group name used in the mediator merMod model. If a merMod function is applied to only one of the mediator or the outcome, group.out is automatically set to the group name used in the merMod model.}  
\item{...}{ other arguments passed to \code{vcovHC} in the \code{sandwich} package: typically the 'type' argument. Ignored if 'robustSE' is 'FALSE'.}
}

\details{This is the workhorse function for estimating causal mediation effects for a variety of data types. The average causal mediation effect (ACME) represents the expected difference in the potential outcome when the mediator took the value that would realize under the treatment condition as opposed to the control condition, while the treatment status itself is held constant. That is,
\deqn{\delta(t) \ = \ E\{Y(t, M(t_1)) - Y(t, M(t_0))\},}{%
      \delta(t) = E[Y(t, M(t1)) - Y(t, M(t0))],}
where \eqn{t, t_1, t_0}{t, t1, t0} are particular values of the treatment \eqn{T} such that \eqn{t_1 \neq t_0}{t1 != t0}, \eqn{M(t)} is the potential mediator, and \eqn{Y(t,m)} is the potential outcome variable. The average direct effect (ADE) is defined similarly as,
\deqn{\zeta(t) \ = \ E\{Y(t_1, M(t)) - Y(t_0, M(t))\},}{%
      \zeta(t) = E[Y(t1, M(t)) - Y(t0, M(t))],}
which represents the expected difference in the potential outcome when the treatment is changed but the mediator is held constant at the value that would realize if the treatment equals \eqn{t}. The two quantities on average add up to the total effect of the treatment on the outcome, \eqn{\tau}. See the references for more details.

When both the mediator model ('model.m') and outcome model ('model.y') are normal linear regressions, the results will be identical to the usual LSEM method by Baron and Kenny (1986).  The function can, however, accommodate other data types including binary, ordered and count outcomes and mediators as well as censored outcomes.  Variables can also be modeled nonparametrically, semiparametrically, or using quantile regression.

If it is desired that inference be made conditional on specific values of the pre-treatment covariates included in the model, the `covariates' argument can be used to set those values as a list or data frame. The list may contain either the entire set or any strict subset of the covariates in the model; in the latter case, the effects will be averaged over the other covariates. The `covariates' argument will be particularly useful when the models contain interactions between the covariates and the treatment and/or mediator (known as ``moderated mediation'').

The prior weights in the mediator and outcome models are taken as sampling weights and the estimated effects will be weighted averages when non-NULL weights are used in fitting 'model.m' and 'model.y'. This will be useful when data does not come from a simple random sample, for example. 

As of version 4.0, the mediator model can be of either 'lm', 'glm' (or `bayesglm'), 'polr' (or `bayespolr'), 'gam', 'rq', `survreg', or `merMod' class, corresponding respectively to the linear regression models, generalized linear models, ordered response models, generalized additive models, quantile regression models, parametric duration models, or multilevel models.. For binary response models, the 'mediator' must be a numeric variable with values 0 or 1 as opposed to a factor. Quasi-likelihood-based inferences are not allowed for the mediator model because the functional form must be exactly specified for the estimation algorithm to work.  The 'binomial' family can only be used for binary response mediators and cannot be used for multiple-trial responses.  This is due to conflicts between how the latter type of models are implemented in \code{\link{glm}} and how 'mediate' is currently written.

For the outcome model, the censored regression model fitted via package \code{VGAM} (of class 'vglm' with 'family@vfamily' equal to "tobit") can be used in addition to the models listed above for the mediator.  The 'mediate' function is not compatible with censored regression models fitted via other packages.  When the quantile regression is used for the outcome model ('rq'), the estimated quantities are quantile causal mediation effects, quantile direct effects and etc., instead of the average effects. If the outcome model is of class 'survreg', the name of the outcome variable must be explicitly supplied in the `outcome' argument. This is due to the fact that 'survreg' objects do not contain that information in an easily extractable form. It should also be noted that for \code{\link{survreg}} models, the \code{\link{Surv}} function must be directly used in the model formula in the call to the survreg function, and that censoring types requiring more than two arguments to Surv (e.g., interval censoring) are not currently supported by 'mediate'.

The quasi-Bayesian approximation (King et al. 2000) cannot be used if 'model.m' is of class 'rq' or 'gam', or if 'model.y' is of class 'gam', 'polr' or 'bayespolr'. In these cases, either an error message is returned or use of the nonparametric bootstrap is forced. Users should note that use of the nonparametric bootstrap often requires significant computing time, especially when 'sims' is set to a large value.

The 'control' argument must be provided when 'gam' is used for the outcome model and user wants to allow ACME and ADE to vary as functions of the treatment (i.e., to relax the "no interaction" assumption). Note that the outcome model must be fitted via package \code{\link{mgcv}} with appropriate formula using \code{\link{s}} constructs (see Imai et al. 2009 in the references). For other model types, the interaction can be allowed by including an interaction term between \eqn{T} and \eqn{M} in the linear predictor of the outcome model. As of version 3.0, the 'INT' argument is deprecated and the existence of the interaction term is automatically detected (except for 'gam' outcome models).

When the treatment variable is continuous or a factor with multiple levels, user must specify the values of \eqn{t_1}{t1} and \eqn{t_0}{t0} using the 'treat.value' and 'control.value' arguments, respectively.  The value of \eqn{t} in the above expressions is set to \eqn{t_0}{t0} for 'd0', 'z0', etc. and to \eqn{t_1}{t1} for 'd1', 'z1', etc.
}


\value{\code{mediate} returns an object of class "\code{mediate}", "\code{mediate.order}" if the outcome model used is 'polr' or 'bayespolr', or "\code{mediate.mer}" if 'lmer' or 'glmer' is used for the outcome or the mediator model, a list that contains the components listed below.  Some of these elements are not available if 'long' is set to 'FALSE' by the user.

The function \code{summary} (i.e., \code{summary.mediate}, \code{summary.mediate.order}, or \code{summary.mediate.mer}) can be used to obtain a table of the results.  The function \code{plot} (i.e., \code{plot.mediate}, \code{plot.mediate.order}, or \code{plot.mediate.mer}) can be used to produce a plot of the estimated average causal mediation, average direct, and total effects along with their confidence intervals. 

\item{d0, d1}{point estimates for average causal mediation effects under the control and treatment conditions.}
\item{d0.ci, d1.ci}{confidence intervals for average causal mediation effects. The confidence level is set at the value specified in 'conf.level'.}
\item{d0.p, d1.p}{two-sided p-values for average causal mediation effects.}
\item{d0.sims, d1.sims}{vectors of length 'sims' containing simulation draws of average causal mediation effects.}
\item{z0, z1}{point estimates for average direct effect under the control and treatment conditions.}
\item{z0.ci, z1.ci}{confidence intervals for average direct effects.}
\item{z0.p, z1.p}{two-sided p-values for average causal direct effects.}
\item{z0.sims, z1.sims}{vectors of length 'sims' containing simulation draws of average direct effects.}
\item{n0, n1}{the "proportions mediated", or the size of the average causal mediation effects relative to the total effect.}
\item{n0.ci, n1.ci}{confidence intervals for the proportions mediated.}
\item{n0.p, n1.p}{two-sided p-values for proportions mediated.}
\item{n0.sims, n1.sims}{vectors of length 'sims' containing simulation draws of the proportions mediated.}
\item{tau.coef}{point estimate for total effect.}
\item{tau.ci}{confidence interval for total effect.}
\item{tau.p}{two-sided p-values for total effect.}
\item{tau.sims}{a vector of length 'sims' containing simulation draws of the total effect.}
\item{d.avg, z.avg, n.avg}{simple averages of d0 and d1, z0 and z1, n0 and n1, respectively, which users may want to use as summary values when those quantities differ.}
\item{d.avg.ci, z.avg.ci, n.avg.ci}{confidence intervals for the above.}
\item{d.avg.p, z.avg.p, n.avg.p}{two-sided p-values for the above.}
\item{d.avg.sims, z.avg.sims, n.avg.sims}{vectors of length 'sims' containing simulation draws of d.avg, z.avg and n.avg, respectively.}
\item{d0.group, d1.group}{group-specific point estimates for average causal mediation effects under the control and treatment conditions.}
\item{d0.ci.group, d1.ci.group}{group-specific confidence intervals for average causal mediation effects. The confidence level is set at the value specified in 'conf.level'.}
\item{d0.p.group, d1.p.group}{group-specific two-sided p-values for average causal mediation effects.}
\item{d0.sims.group, d1.sims.group}{group-specific vectors of length 'sims' containing simulation draws of average causal mediation effects.}
\item{z0.group, z1.group}{group-specific point estimates for average direct effect under the control and treatment conditions.}
\item{z0.ci.group, z1.ci.group}{group-specific confidence intervals for average direct effects.}
\item{z0.p.group, z1.p.group}{group-specific two-sided p-values for average causal direct effects.}
\item{z0.sims.group, z1.sims.group}{group-specific vectors of length 'sims' containing simulation draws of average direct effects.}
\item{n0.group, n1.group}{the group-specific "proportions mediated", or the size of the group-specific average causal mediation effects relative to the total effect.}
\item{n0.ci.group, n1.ci.group}{group-specific confidence intervals for the proportions mediated.}
\item{n0.p.group, n1.p.group}{group-specific two-sided p-values for proportions mediated.}
\item{n0.sims.group, n1.sims.group}{group-specific vectors of length 'sims' containing simulation draws of the proportions mediated.}
\item{tau.coef.group}{group-specific point estimate for total effect.}
\item{tau.ci.group}{group-specific confidence interval for total effect.}
\item{tau.p.group}{group-specific two-sided p-values for total effect.}
\item{tau.sims.group}{a group-specific vector of length 'sims' containing simulation draws of the total effect.}
\item{d.avg.group, z.avg.group, n.avg.group}{group-specific simple averages of d0 and d1, z0 and z1, n0 and n1, respectively, which users may want to use as summary values when those quantities differ.}
\item{d.avg.ci.group, z.avg.ci.group, n.avg.ci.group}{group-specific confidence intervals for the above.}
\item{d.avg.p.group, z.avg.p.group, n.avg.p.group}{group-specific two-sided p-values for the above.}
\item{d.avg.sims.group, z.avg.sims.group, n.avg.sims.group}{group-specific vectors of length 'sims' containing simulation draws of d.avg, z.avg and n.avg, respectively.}
\item{boot}{logical, the 'boot' argument used.}
\item{treat}{a character string indicating the name of the 'treat' variable used.}
\item{mediator}{a character string indicating the name of the 'mediator' variable used.}
\item{INT}{a logical value indicating whether the model specification allows the effects to differ between the treatment and control conditions.}
\item{conf.level}{the confidence level used. }
\item{model.y}{the outcome model used.}
\item{model.m}{the mediator model used.}
\item{group.m}{the name of the mediator group used.}
\item{group.y}{the name of the outcome group used.}
\item{group.name}{the name of the group on which the output is based.}
\item{group.id.m}{the data on the mediator group.}
\item{group.id.y}{the data on the outcome group.}
\item{group.id}{the data on the group on which the output is based.}
\item{control.value}{value of the treatment variable used as the control condition.}
\item{treat.value}{value of the treatment variable used as the treatment condition.}
\item{nobs}{number of observations in the model frame for 'model.m' and 'model.y'. May differ from the numbers in the original models input to 'mediate' if 'dropobs' was 'TRUE'.}
\item{robustSE}{`TRUE' or `FALSE'.}
\item{cluster}{the clusters used.}
}

\references{
Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. (2014). "mediation: R package for Causal Mediation Analysis", Journal of Statistical Software, Vol. 59, No. 5, pp. 1-38. 

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the Black Box of Causality: Learning about Causal Mechanisms from Experimental and Observational Studies, American Political Science Review, Vol. 105, No. 4 (November), pp. 765-789.

Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.} 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}; 
Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}; 
Luke Keele, Penn State University, \email{ljk20@psu.edu}; 
Kosuke Imai, Princeton University, \email{kimai@princeton.edu};
Kentaro Hirose, Princeton University, \email{hirose@princeton.edu}.}
 
\seealso{\code{\link{medsens}}, \code{\link{plot.mediate}}, \code{\link{summary.mediate}}, \code{\link{summary.mediate.mer}}, \code{\link{plot.mediate.mer}}, \code{\link{mediations}}, \code{vcovHC} }

\examples{ 
# Examples with JOBS II Field Experiment

# **For illustration purposes a small number of simulations are used**

data(jobs)

####################################################
# Example 1: Linear Outcome and Mediator Models
####################################################
b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
c <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs)

# Estimation via quasi-Bayesian approximation
contcont <- mediate(b, c, sims=50, treat="treat", mediator="job_seek")
summary(contcont)
plot(contcont)

\dontrun{
# Estimation via nonparametric bootstrap
contcont.boot <- mediate(b, c, boot=TRUE, sims=50, treat="treat", mediator="job_seek")
summary(contcont.boot)
}

# Allowing treatment-mediator interaction
d <- lm(depress2 ~ treat + job_seek + treat:job_seek + econ_hard + sex + age, data=jobs)
contcont.int <- mediate(b, d, sims=50, treat="treat", mediator="job_seek")
summary(contcont.int)

# Allowing ``moderated mediation'' with respect to age
b.int <- lm(job_seek ~ treat*age + econ_hard + sex, data=jobs)
d.int <- lm(depress2 ~ treat*job_seek*age + econ_hard + sex, data=jobs)
contcont.age20 <- mediate(b.int, d.int, sims=50, treat="treat", mediator="job_seek",
			covariates = list(age = 20))
contcont.age70 <- mediate(b.int, d.int, sims=50, treat="treat", mediator="job_seek",
			covariates = list(age = 70))
summary(contcont.age20)
summary(contcont.age70)

# Continuous treatment
jobs$treat_cont <- jobs$treat + rnorm(nrow(jobs))  # (hypothetical) continuous treatment
b.contT <- lm(job_seek ~ treat_cont + econ_hard + sex + age, data=jobs)
c.contT <- lm(depress2 ~ treat_cont + job_seek + econ_hard + sex + age, data=jobs)
contcont.cont <- mediate(b.contT, c.contT, sims=50, 
                    treat="treat_cont", mediator="job_seek",
                    treat.value = 4, control.value = -2)
summary(contcont.cont)

# Categorical treatment 
\dontrun{
b <- lm(job_seek ~ educ + sex, data=jobs)
c <- lm(depress2 ~ educ + job_seek + sex, data=jobs)

# compare two categories of educ --- gradwk and somcol
model.cat <- mediate(b, c, treat="educ", mediator="job_seek", sims=50, 
                     control.value = "gradwk", treat.value = "somcol")
summary(model.cat)
}

######################################################
# Example 2: Binary Outcome and Ordered Mediator
######################################################
\dontrun{
b.ord <- polr(job_disc ~ treat + econ_hard + sex + age, data=jobs,
            method="probit", Hess=TRUE)
d.bin <- glm(work1 ~ treat * job_disc + econ_hard + sex + age, data=jobs,
            family=binomial(link="probit"))
ordbin <- mediate(b.ord, d.bin, sims=50, treat="treat", mediator="job_disc")
summary(ordbin)

# Using heteroskedasticity-consistent standard errors
require(sandwich)
ordbin.rb <- mediate(b.ord, d.bin, sims=50, treat="treat", mediator="job_disc",
            robustSE=TRUE)
summary(ordbin.rb)
}

######################################################
# Example 3: Quantile Causal Mediation Effect
######################################################
require(quantreg)
c.quan <- rq(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs,
            tau = 0.5)  # median
contquan <- mediate(b, c.quan, sims=50, treat="treat", mediator="job_seek")
summary(contquan)

######################################################
# Example 4: GAM Outcome
######################################################
\dontrun{
require(mgcv)
c.gam <- gam(depress2 ~ treat + s(job_seek, bs="cr") + 
            econ_hard + sex + age, data=jobs)
contgam <- mediate(b, c.gam, sims=10, treat="treat", 
                mediator="job_seek", boot=TRUE)
summary(contgam)

# With interaction
d.gam <- gam(depress2 ~ treat + s(job_seek, by = treat) + 
    s(job_seek, by = control) + econ_hard + sex + age, data=jobs)
contgam.int <- mediate(b, d.gam, sims=10, treat="treat", mediator="job_seek",
    control = "control", boot=TRUE)
summary(contgam.int)
}
######################################################
# Example 5: Multilevel Outcome and Mediator Models
######################################################
\dontrun{
require(arm)
 
# educ: mediator group
# occp: outcome group

# Varying intercept for mediator 
model.m <- glmer(job_dich ~ treat + econ_hard + (1 | educ), 
             		     family = binomial(link = "probit"), data = jobs)

# Varying intercept and slope for outcome
model.y <- glmer(work1 ~ treat + job_dich + econ_hard + (1 + treat | occp), 
                family = binomial(link = "probit"), data = jobs)

# Output based on mediator group ("educ")
multilevel <- mediate(model.m, model.y, treat = "treat", 
              mediator = "job_dich", sims=50, group.out="educ")

# Output based on outcome group ("occp")
# multilevel <- mediate(model.m, model.y, treat = "treat", 
              mediator = "job_dich", sims=50) 

# Group-average effects  
summary(multilevel)

# Group-specific effects organized by effect
summary(multilevel, output="byeffect")
# plot(multilevel, group.plots=TRUE)
# See summary.mediate.mer and plot.mediate.mer for detailed explanations 

# Group-specific effects organized by group
summary(multilevel, output="bygroup")
# See summary.mediate.mer for detailed explanations 
}


} 
