\name{plot.mediate} 
\alias{plot.mediate} 
\alias{plot.mediate.order} 
\title{Plotting Indirect, Direct, and Total Effects from Mediation Analysis} 
\description{ 
Function to plot results from \code{mediate}. The vertical axis lists indirect, direct, and total effects and the horizontal axis indicates the respective magnitudes. Most standard options for plot function available.
} 
\usage{
\method{plot}{mediate}(x, treatment = NULL,
    labels = c("ACME","Direct\nEffect","Total\nEffect"), 
    xlim = NULL, ylim = NULL, xlab = "", ylab = "",
    main = NULL, lwd = 1.5, cex = .85,
    col = "black", ...)
}

\arguments{
\item{x}{ object of class \code{mediate} or \code{mediate.order} as produced by \code{mediate}.}
\item{treatment}{ a character string indicating the baseline treatment value of the estimated causal mediation effect and direct effect to plot. Can be either "control", "treated" or "both". If 'NULL' (default), both sets of estimates are plotted if and only if they differ.}
\item{labels}{ a vector of length 3 indicating the labels for the estimated effects.}
\item{xlim}{ range of the horizontal axis.}
\item{ylim}{ range of the vertical axis.}
\item{xlab}{ label of the horizontal axis.}
\item{ylab}{ label of the vertical axis.}
\item{main}{ main title.}
\item{lwd}{ width of the horizontal bars for confidence intervals.}
\item{cex}{ size of the dots for point estimates.}
\item{col}{ color of the dots and horizontal bars for the estimates.}
\item{...}{ additional parameters passed to 'plot'.}
}

\value{\code{mediate} returns an object of class "\code{mediate}". The function \code{summary} is used to obtain a table of the results. The \code{plot} function plots these quantities.
}

\references{Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.} 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu},
Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}.}
 
\seealso{\code{\link{mediate}}, \code{\link{plot}}}

