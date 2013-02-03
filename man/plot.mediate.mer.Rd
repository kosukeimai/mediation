\name{plot.mediate.mer}
\alias{plot.mediate.mer}

\title{Plotting Indirect, Direct, and Total Effects from Mediation Analysis of Multilevel Models}

\description{
Function to plot group-specific effects derived from causal mediation analysis of multilevel models.
}

\usage{ 
\method{plot}{mediate.mer}(x,  treatment = NULL, xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                             main = NULL, lwd = 1.5, cex = .85,
                             col = "black", group.plot=FALSE, plot.type = 1:5, ...)
}

\arguments{ 
\item{x}{ object of class 'mediate.mer' produced by 'mediate'.}
\item{group.plot}{ group-specific effects if group.plot = "TRUE"; default is group-average effects.}
\item{plot.type}{ select group effects to plot; mediation effect (0) if plot.type = 1; mediation effect (1) if plot.type = 2; direct effect (0) if plot.type = 3; direct effect (1) if plot.type = 4; total effect if plot.type = 5. Combination is possible; for example, mediation effect (1) and direct effect (1) if plot.type = c(2,4). Default is 1 to 5.}
\item{treatment}{ a character string indicating the baseline treatment value of the estimated causal mediation effect and direct effect to plot. Can be either "control", "treated", or "both". If 'NULL' (default), both sets of estimates are plotted if and only if the treatment-mediator interaction term is included in the outcome model. Only relevant if group.plot = "FALSE".}
\item{xlim}{ range of the horizontal axis.}
\item{ylim}{ range of the vertical axis.}
\item{xlab}{ label of the horizontal axis.}
\item{ylab}{ label of the vertical axis.}
\item{main}{ main title.}
\item{lwd}{ width of the horizontal bars for confidence intervals .}
\item{cex}{ size of the dots for point estimates.}
\item{col}{ color of the dots and horizontal bars for the estimates..}
\item{...}{ additional parameters passed to 'plot'.}
}

\author{Kentaro Hirose, Princeton University, \email{hirose@princeton.edu}.}
 
\seealso{\code{\link{mediate}}, \code{\link{summary.mediate.mer}}.}


\examples{ 
# Examples with JOBS II Field Experiment

# **For illustration purposes a small number of simulations are used**

data(jobs)
require(arm)

# educ: mediator group
# occp: outcome group

# Varying intercept for mediator 
model.m <- glmer(job_dich ~ treat + econ_hard + (1 | educ), 
             		     family = binomial(link = "probit"), data = jobs)

# Varying intercept and slope for outcome
model.y <- glmer(work1 ~ treat + job_dich + econ_hard + (1 + econ_hard | occp), 
                family = binomial(link = "probit"), data = jobs)

# Output based on mediator group
multilevel <- mediate(model.m, model.y, treat = "treat", mediator = "job_dich", sims=50, group.m="educ", group.y="occp", group.out="educ")

#plot(multilevel, group.plot=TRUE)     
#plot(multilevel, group.plot=TRUE, plot.type = c(2,4))

