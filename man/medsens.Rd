\name{medsens} 
\alias{medsens} 
\title{Sensitivity Analysis for Causal Mediation Effects} 
\description{ 
'medsens' is used to perform sensitivity analysis on the average causal mediation effects and direct effects for violations of the sequential ignorability assumption.  The function takes output from '\code{mediate}' and calculates the true average causal mediation effects and direct effects for different values of the sensitivity parameter representing the degree of the sequantial ignorability violation.
} 
\usage{

medsens(x, rho.by = 0.1, sims = 1000, eps = sqrt(.Machine$double.eps), 
        effect.type = c("indirect", "direct", "both"))

} 

\arguments{ 
\item{x}{ an object of class 'mediate', typically an output from the \code{mediate} function.} 
\item{rho.by}{ a numeric value between 0 and 1 indicating the increment for the sensitivity parameter, rho.}
\item{sims}{ the number of Monte Carlo draws for the calculation of confidence intervals. Only used in cases where either the mediator or outcome variable is binary.}
\item{eps}{ convergence tolerance parameter for the iterative FGLS. Only used when both the mediator and outcome models are linear.}
\item{effect.type}{ a character string indicating which effect(s) to be analyzed.  Default is "indirect". Currently, "direct" and "both" can be used only when both the mediator and outcome models are linear.}
} 

\details{This is the workhorse function for sensitivity analyses for average causal mediation effects. The sensitivity analysis can be used to assess the robustness of the findings from \code{mediate} to the violation of sequantial ignorability, the crucial identification assumption necessary for the estimates to be valid. The analysis proceeds by quantifying the degree of sequantial ignorability violation as the correlation between the error terms of the mediator and outcome models, and then calculating the true values of the average causal mediation effect for given values of this sensitivity parameter, rho. The original findings are deemed sensitive if the true effects are found to vary widely as function of rho.

Users should note that computation time is several minutes for these functions. Setting rho.by to smaller numbers significantly increases computational time, as does increasing eps or the number of simulations.

Function to perform sensitivity analysis on mediation effect for violations of sequential ignorability assumption. This allows for a correlation between the error terms of the outcome model and the mediator model. Sensitivity analysis is possible with 1) continuous mediator and continuous outcome, 2) binary outcome and continuous mediator, and 3) continuous outcome and binary mediator. Output from the function can be passed through summary or plot functions which display estimated mediation effects for given values of rho.
}


\value{\code{medsens} returns an object of class "\code{medsens}", a list containing the following elements. Some of these elements are not available depending on the 'effect.type' argument specified by the user. The output can then be passed to the \code{\link{summary}} (i.e., \code{\link{summary.medsens}}) and \code{\link{plot}} (i.e., \code{\link{plot.medsens}}) functions to produce tabular and graphical summaries of the results.

\item{d0, d1}{vectors of point estimates for average causal mediation effects under the control and treatment conditions for each value of sensitivity parameter rho.}
\item{upper.d0, lower.d0, upper.d1, lower.d1}{vectors of upper and lower confidence limits for average causal mediation effect under the control and treatment conditions for each value of rho.}
\item{z0, z1}{vectors of point estimates for average direct effect under the control and treatment conditions for each value of sensitivity parameter rho.}
\item{upper.z0, lower.z0, upper.z1, lower.z1}{vectors of upper and lower confidence limits for average direct effect under the control and treatment conditions for each value of rho.}
\item{tau}{a vector of point estimates for total effect for each value of rho.  Only present when the outcome model is binary.}
\item{upper.tau, lower.tau}{vectors of upper and lower confidence limits for total effect.  Only present when the outcome model is binary.}
\item{nu}{a vector of point estimates for the proportion mediated for each value of rho.  Only present when the outcome model is binary.}
\item{upper.nu, lower.nu}{vectors of upper and lower confidence limits for the proportion mediated.  Only present when the outcome model is binary.}
\item{rho}{a numeric vector containing the values of sensitivity parameter rho used.}
\item{rho.by}{a numeric value indicating the increment of rho used.}
\item{err.cr.d, err.cr.z}{the values of rho with which the average causal mediation and direct effects are zero.  Vectors of length two if 'INT' is 'TRUE'; numeric values otherwise.}
\item{ind.d0, ind.d1, ind.z0, ind.z1}{vectors of 0s/1s, indicating whether the confidence intervals of d0, d1, z0 and z1 do not cover zero for each value of rho.}
\item{R2star.prod}{a numeric vector containing the values of the products of the two "R square stars", representing the proportions of residual variance in the mediator and outcome explained by the hypothesized unobserved confounder. The values correspond to those of rho.}
\item{R2tilda.prod}{a numeric vector containing the values of the products of the two "R square tildes", representing the proportions of total variance in the mediator and outcome explained by the hypothesized unobserved confounder. The values correspond to those of rho.}
\item{R2star.d.thresh, R2star.z.thresh}{the values of the product of "R square star" for which the average causal mediation and direct effects are zero, respectively.}
\item{R2tilda.d.thresh, R2tilda.z.thresh}{the values of the product of "R square tilde" for which the average causal mediation and direct effects are zero, respectively.}
\item{r.square.y, r.square.m}{the usual R square statistics for the outcome and mediator models.}
\item{INT}{a logical value indicating whether interaction between the treatment and mediator is allowed in the original mediate object.}
\item{effect.type}{ the 'effect.type' argument used.}
\item{type}{a character string indicating the type of the mediator and outcome models used. Currently either "ct" (linear mediator and outcome models), 'bm' (binary mediator and linear outcome models) or 'bo' (linear mediator and binary outcome models).}
}


\references{Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.} 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}, 
Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}, 
Jaquilyn Waddell-Boie, Princeton University, \email{jwaddell@princeton.edu},
Luke Keele, Ohio State University, \email{keele.4@osu.edu}, 
Kosuke Imai, Princeton University, \email{kimai@princeton.edu} }
 
\seealso{\code{\link{mediate}}, \code{\link{summary.medsens}}, \code{\link{plot.medsens}}.}

\examples{ 

#Example with JOBS II Field experiment

#**For illustration purposes simulations set to low number**

data(jobs)

# Fit parametric models
model.m <- lm(job_seek ~ treat + depress1, data=jobs)
model.y <- lm(depress2 ~ treat + job_seek + depress1, data=jobs)

# Pass model objects through mediate function
med.cont <- mediate(model.m, model.y, treat="treat", mediator="job_seek", sims=10)

# Pass mediate output through medsens function
sens.cont <- medsens(med.cont, rho.by=.1, eps=.01, effect.type="both")

# Use summary function to display results
summary(sens.cont)

# Plot true ACMEs and ADEs as function of rho
par.orig <- par(mfrow = c(3,2))
plot(sens.cont, main="JOBS", ylim=c(-.2,.2))

# Plot true ACMEs as function of "R square tildes"
plot(sens.cont, sens.par="R2", r.type=2, sign.prod=1)
par(par.orig)
} 
