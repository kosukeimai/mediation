\name{summary.multimed}
\alias{summary.multimed}
\alias{print.summary.multimed} 

\title{Summarizing Output from Mediation Analysis with Multiple Mechanisms} 

\description{ 
Function to report results from the \code{\link{multimed}} function.
} 

\usage{

  \method{summary}{multimed}(object, ...)

  \method{print}{summary.multimed}(x, ...)

} 

\arguments{ 
\item{object}{ object of class \code{multimed}, typically output from the \code{multimed} funciton.} 
\item{x}{ output from the summary function.} 
\item{...}{ additional arguments affecting the summary produced.}
} 

\references{
Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. (2014). "mediation: R package for Causal Mediation Analysis", Journal of Statistical Software, Vol. 59, No. 5, pp. 1-38. 

Imai, K. and Yamamoto, T. (2012) Identification and Sensitivity Analysis for Multiple Causal Mechanisms: Revisiting Evidence from Framing Experiments, typescript.

}

\author{Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}.}

\seealso{\code{\link{multimed}}, \code{\link{plot.multimed}}}

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
