\name{plot.mediations}
\alias{plot.mediations}
\alias{plot.mediations.order}
\title{Plotting Indirect, Direct, and Total Effects from Multiple Mediation Analyses}

\description{
Function to plot results from multiple causal mediation analyses conducted via the \code{\link{mediations}} funciton. Output is a series of plots generated via \code{\link{plot.mediate}} for each treatment/mediator/outcome combination specified in the input 'mediations' object.
}



\usage{ 

\method{plot}{mediations}(x, which = names(x),
        ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...)
}



\arguments{ 
\item{x}{ output from the mediations function.}
\item{which}{ subset of names(x), indicating which model combinations to be plotted. Default is to plot all.}
\item{ask}{ logical. If 'TRUE', the user is asked for input before a new figure is plotted.  Default is to ask only if the number of plots on current screen is fewer the number implied by 'which'.}
\item{...}{ arguments passed to the \code{\link{plot.mediate}} function for individual plots.}
}



\value{\code{mediations} returns an object of class \code{mediations}.  The function \code{summary} is used to 
obtain a table of the results. The plot function instead plots these quantities. All additional parameters desired 
for the plotting of an output from \code{mediate} can be passed through.
}


\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu};
Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}.}
 
\seealso{\code{\link{mediations}}, \code{\link{plot.mediate}}, \code{\link{plot}}.}

