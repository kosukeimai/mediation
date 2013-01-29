\name{plot.mediate.mer}
\alias{plot.mediate.mer}

\title{Plotting Indirect, Direct, and Total Effects from Mediation Analysis of Multilevel Models}

\description{
Function to plot group-specific effects derived from causal mediation analysis of multilevel models.
}

\usage{ 
\method{plot}{mediate.mer}(x, treatment = NULL, xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                           main = NULL, lwd = 1.5, cex = .85,
                           col = "black", group.plot=FALSE,...)
}

\arguments{ 
\item{x}{ object of class 'mediate.mer' produced by 'mediate'.}
\item{group.plot}{ group-specific effects if group.plot = "TRUE"; default is group-average effects.}
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

