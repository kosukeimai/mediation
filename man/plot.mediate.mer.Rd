\name{plot.mediate.mer}
\alias{plot.mediate.mer}

\title{Plotting Indirect, Direct, and Total Effects from Mediation Analysis of Multilevel Models}

\description{
Function to plot group-specific effects derived from causal mediation analysis of multilevel models.
}

\usage{ 
\method{plot}{mediate.mer}(x, treatment = NULL, group.plots = FALSE,
                             ask = prod(par("mfcol")) < nplots, 
                             xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                             main = NULL, lwd = 1.5, cex = .85,
                             col = "black", ...)
}

\arguments{ 
\item{x}{ object of class 'mediate.mer' produced by 'mediate'.}
\item{treatment}{ a character string indicating the baseline treatment value of the estimated causal mediation effect and direct effect to plot. Can be either "control", "treated", or "both". If 'NULL' (default), both sets of estimates are plotted if and only if they differ.}
\item{group.plots}{ a logical value indicating whether group-specific effects should be plotted in addition to the population-averaged effects.}
\item{ask}{ a logical value. If 'TRUE', the user is asked for input before a new figure is plotted. Default is to ask only if the number of plots on current screen is fewer than necessary.}
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
\dontrun{
data(jobs)
require(arm)

# educ: mediator group
# occp: outcome group

# Varying intercept for mediator 
model.m <- glmer(job_dich ~ treat + econ_hard + (1 | educ), 
             		     family = binomial(link = "probit"), data = jobs)

# Varying intercept and slope for outcome
model.y <- glmer(work1 ~ treat + job_dich + econ_hard + (1 + treat | occp), 
                family = binomial(link = "probit"), data = jobs)

# Output based on mediator group
multilevel <- mediate(model.m, model.y, treat = "treat", 
                      mediator = "job_dich", sims=50, group.out="educ")

#plot(multilevel, group.plots=TRUE)
}
}
