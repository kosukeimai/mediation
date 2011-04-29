\name{summary.mediations} 
\alias{summary.mediations}
\alias{summary.mediations.order}
\alias{print.summary.mediations}
\alias{print.summary.mediations.order}

\title{Summarizing Output from Multiple Mediation Analyses} 
\description{ 
The 'summary.mediations' function produces a summary of results from multiple causal analyses conducted via \code{\link{mediations}}.  Output is a series of \code{\link{summary.mediate}} outputs for all the treatment/mediator/outcome combinations used in the input 'mediations' object.
} 
\usage{

  \method{summary}{mediations}(object, ...) 

  \method{print}{summary.mediations}(x, ...)

} 

\arguments{ 
\item{object}{ output from mediations function.} 
\item{x}{ output from summary.mediations function.} 
\item{...}{ additional arguments affecting the summary produced.}
} 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu},
Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}. }
 
\seealso{\code{\link{mediations}}, \code{\link{summary.mediate}}, \code{\link{summary}}.}

