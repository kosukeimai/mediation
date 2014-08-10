\name{multimed} 
\alias{multimed} 
\title{Estimation and Sensitivity Analysis for Multiple Causal Mechanisms} 
\description{ 
'multimed' is used for causal mediation analysis when post-treatment 
mediator-outcome confounders, or alternative mediators causally preceding 
the mediator of interest, exist in the hypothesized causal mechanisms.
It estimates the average causal mediation effects (indirect effects) under the
homogeneous interaction assumption based on a varying-coefficient linear
structural equation model. The function also performs sensitivity analysis
with respect to the violation of the homogenous interaction assumption.
The function can be used for both the single experiment design and the
parallel design.
} 

\usage{

multimed(outcome, med.main, med.alt = NULL, treat, 
         covariates = NULL, experiment = NULL,
         data, design = c("single", "parallel"),
         sims = 1000, R2.by = 0.01, conf.level = 0.95)

}

\arguments{
\item{outcome}{ name of the outcome variable in 'data'.}
\item{med.main}{ name of the mediator of interest. Under the parallel design this is the only mediator variable used in the estimation.}
\item{med.alt}{ vector of character strings indicating the names of the post-treatment confounders, i.e., the alternative mediators affecting both the main mediator and outcome. Not needed for the parallel design.}
\item{treat}{ name of the treatment variable in 'data'.}
\item{covariates}{ vector of character strings representing the names of the pre-treatment covariates. Cannot be used for the parallel design.}
\item{experiment}{ name of the binary indicator whether 'med.main' is randomly manipulated or not. Only used under the parallel design.}
\item{data}{ a data frame containing all the above variables.}
\item{design}{ experimental design used. Defaults to 'single'.}
\item{sims}{ number of bootstrap samples used for the calculation of confidence intervals.}
\item{R2.by}{ increment for the "R square tilde" parameter, i.e. the sensitivity parameter representing the proportion of residual outcome variance explained by heterogeneity in treatment-mediator interactions. Must be a numeric value between 0 and 1.}
\item{conf.level}{ level to be used for confidence intervals.}
}

\details{This function implements the framework proposed by Imai and Yamamoto (2012) for the estimation and sensitivity analysis for multiple causal mechanisms. It estimates the average causal mediation effects with respect to the mediator of interest ('med.main'), i.e., the portion of the treatment effect on the outcome that is transmitted through that mediator. Unlike the "standard" causal mediation analysis implemented by \code{\link{mediate}} and \code{\link{medsens}}, this framework allows the existence of post-treatment covariates that confound the relationship between the main mediator and the outcome, or equivalently, alternative mediators ('med.alt') that causally precede the main mediator.
	
When the parallel design was used for the experiment (i.e. when the experiment contained an additional randomly assigned group for which both the treatment and the mediator were randomized), there is no need to specify a particular post-treatment confounder, for any such confounder (observed or unobserved) is allowed to exist by virtue of the design. Similarly, no observed covariates need to be included. The function instead requires an additional variable ('experiment') indicating whether the mediator was randomly manipulated for the unit.

The estimation and sensitivity analysis are both based on a varying-coefficient linear structural equations model, which assumes additivity but allows for an arbitrary degree of heterogeneity in model coefficients across units and thus is substantially more flexible than a traditional SEM framework. For details see Imai and Yamamoto (2012).

The function produces two sets of results. First, point estimates of the average causal mediation effects are calculated, along with their (percentile) bootstrap confidence intervals. These estimates are based on the "homogeneous interaction" assumption, or the assumption that the degree of treatment-mediator interaction is constant across all units. The estimated total treatment effect is also reported.

Second, the bounds on the average causal mediation effects are also estimated and computed for various degrees of interaction heterogeneity (i.e., violation of the identification assumption), which are represented by the values of three alternative sensitivity parameters. These parameters are: (1) sigma, the standard deviation of the (varying) regression coefficient on the interaction term, (2) R square star, the proportion of the residual variance that would be explained by an additional term for interaction heterogeneity, and (3) R square tilde, the proportion of the total variance explained by such a term. The confidence region is also calculated, using the Imbens and Manski (2004) formula with bootstrap standard errors. Further details are given in the above reference.

Note that rows with missing values will be omitted from the calculation of the results. Also note that the treatment variable must be a numeric vector of 1 and 0 and that both mediators and outcome variable must be numeric. The pre-treatment covariates can be of any type that \code{\link{lm}} can handle as predictors.
}

\value{\code{multimed} returns an object of class "\code{multimed}", a list contains the following components. The object can be passed to the \code{summary} and \code{plot} method functions for a summary table and a graphical summary.

\item{sigma}{ values of the sigma sensitivity parameter at which the bounds and confidence intervals are evaluated.}
\item{R2tilde}{ values of the R square tilde parameter.}
\item{R2star}{ values of the R square star parameter.}
\item{d1.lb, d0.lb, d.ave.lb}{ lower bounds on the average causal mediation effects under treatment, control, and the simple average of the two, respectively, corresponding to the values of the sensitivity parameters listed above. Note that the first elements of these vectors equal the point estimates under the homogeneous interaction assumption.}
\item{d1.ub, d0.ub, d.ave.ub}{ upper bounds on the average causal mediation effects.}
\item{d1.ci, d0.ci, d.ave.ci}{ confidence intervals for the average causal mediation effects at different values of the sensitivity parameters.}
\item{tau}{ point estimate of the total treatment effect.}
\item{tau.ci}{ confidence interval for the total treatment effect.}
\item{conf.level}{ confidence level used for the calculation of the confidence intervals.}
}

\references{
Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. (2014). "mediation: R package for Causal Mediation Analysis", Journal of Statistical Software, Vol. 59, No. 5, pp. 1-38. 

Imai, K. and Yamamoto, T. (2012) Identification and Sensitivity Analysis for Multiple Causal Mechanisms: Revisiting Evidence from Framing Experiments, Unpublished manuscript.

}

\author{Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}}

\seealso{\code{\link{plot.multimed}}}

\examples{

# Replicates Figure 3 (right column) of Imai and Yamamoto (2012)
# Note: # of bootstrap samples set low for quick illustration

data(framing)
Xnames <- c("age", "educ", "gender", "income")
res <- multimed("immigr", "emo", "p_harm", "treat", Xnames, 
               data = framing, design = "single", sims = 100)
summary(res)
plot(res, type = "point")
plot(res, type = c("sigma", "R2-total"), tgroup = "average")

# Parallel design example using the simulated data of Imai, Tingley and Yamamoto (2012)

data(boundsdata)
res.para <- multimed(outcome = "out", med.main = "med", treat = "ttt", experiment = "manip",
					 data = boundsdata, design = "parallel", sims = 100)
summary(res.para)
plot(res.para, tg = "av")

}
