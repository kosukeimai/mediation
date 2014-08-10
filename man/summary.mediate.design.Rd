\name{summary.mediate.design}
\alias{summary.mediate.design}
\alias{print.summary.mediate.design}

\title{Summarizing Output from Design Based Mediation Analysis}
\description{
Function to report results from design based mediation analysis. Reported categories differ depending on the design and assumptions used.
}
\usage{

  \method{summary}{mediate.design}(object, ...)

  \method{print}{summary.mediate.design}(x, ...)

}

\arguments{
\item{object}{ object of class \code{mediate.design}, typically output from a function for design-based mediation analysis (such as \code{\link{mediate.sed}}).}
\item{x}{ output from the summary function.}
\item{...}{ additional arguments affecting the summary produced.}
}

\references{
Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. (2014). "mediation: R package for Causal Mediation Analysis", Journal of Statistical Software, Vol. 59, No. 5, pp. 1-38. 

Imai, K., Tingley, D. and Yamamoto, T. (2012) Experimental Designs for Identifying Causal Mechanisms. Journal of the Royal Statistical Society, Series A (Statistics in Society)

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the Black Box of Causality: Learning about Causal Mechanisms from Experimental and Observational Studies, American Political Science Review, Vol. 105, No. 4 (November), pp. 765-789.

Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.

Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer.}


\author{
Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu};
Teppei Yamamoto, Massachusetts Institute of Technology, \email{teppei@mit.edu}.
}

\seealso{\code{\link{mediate}}, \code{\link{plot.mediate}}, \code{\link{summary}}.}

