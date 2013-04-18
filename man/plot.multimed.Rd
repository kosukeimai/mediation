\name{plot.multimed}
\alias{plot.multimed}
\title{Plotting the Results of Causal Mediation Analysis for Multiple Mechanisms}
\description{
  Function to plot results from \code{multimed}. Most standard plotting
  options are available.
}

\usage{

\method{plot}{multimed}(x, type = c("point", "sigma", "R2-residual", "R2-total"),
     tgroup = c("average", "treated", "control"),
     effect.type = c("indirect", "total"),
     ask = prod(par("mfcol")) < nplots,
     xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, main = NULL,
     lwd = par("lwd"), pch = par("pch"), cex = par("cex"), las = par("las"),
     col.eff = "black", col.cbar = "black", col.creg = "gray", ...)

}

\arguments{
\item{x}{ object of class \code{multimed}, typically output from the \code{multimed} funciton.}
\item{type}{ type of plot(s) required. The default is to produce all, i.e., the point estimates of the effects under the homogenous interaction assumpton ("point") and bounds as function of the sigma ("sigma"), R square star ("R2-residual") and R square tilde ("R2-total") parameters.}
\item{tgroup}{ treatment group(s) for which the estimates are produced. The default is to plot all, i.e., the average causal mediation effect when treated ("treated"), control ("control") and the simple average of these two effects ("average").}
\item{effect.type}{ a vector indicating which quantities of interest to plot for the point estimates. Only plotting total effects is not available.}
\item{ask}{  a logical value. If 'TRUE', the user is asked for input before a new figure is plotted.  Default is to ask only if the number of plots on current screen is fewer than necessary.}
\item{xlab}{ label for the x axis. Default labels are used if 'NULL'.}
\item{ylab}{ label for the y axis. Default labels are used if 'NULL'.}
\item{xlim}{ limits of the x axis. If 'NULL' default values are used.}
\item{ylim}{ limits of the y axis. If 'NULL' default values are used.}
\item{main}{ main title for the plot. If 'NULL', default titles are used.}
\item{lwd}{ width of the lines used in graphs. For the "point" plot this is the width of confidence bars. For sensitivity plots this is the width of the lines for the bounds.}
\item{pch}{ plotting points used for the "point" plots.}
\item{cex}{ magnification factor for the plotting points in the "point" plots.}
\item{las}{ style of the y axis labels in the "point" plots.}
\item{col.eff}{ color of the points in the "point" plots and/or the bounds in sensitivity plots.}
\item{col.cbar}{ color of the confidence bars in the "point" plots.}
\item{col.creg}{ color of the confidence regions in sensitivity plots.}
\item{...}{ additional arguments to be passed to plotting functions.}
}

\details{'type' and 'tgroup' can contain multiple character strings, in which case multiple plots are produced. For the use of graphical parameters see \code{\link{plot}} and the links it contains.
}

\references{

Imai, K. and Yamamoto, T. (2012) Identification and Sensitivity Analysis for Multiple Causal Mechanisms: Revisiting Evidence from Framing Experiments, typescript.

}

\author{Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}.}

\seealso{\code{\link{multimed}}, \code{\link{plot}}}

\examples{

# Replicates Figure 3 (right column) of Imai and Yamamoto (2012)
# Note: # of bootstrap samples set low for quick illustration

data(framing)
Xnames <- c("age", "educ", "gender", "income")
res <- multimed("immigr", "emo", "p_harm", "treat", Xnames,
               data = framing, sims = 100)
summary(res)
plot(res, type = "point")
plot(res, type = c("sigma", "R2-total"), tgroup = "average")

}
