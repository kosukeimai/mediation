% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medsens.R
\name{plot.medsens}
\alias{plot.medsens}
\title{Plotting Results from Sensitivity Analysis for Causal Mediation Effects}
\usage{
\method{plot}{medsens}(
  x,
  sens.par = c("rho", "R2"),
  r.type = c("residual", "total"),
  sign.prod = c("positive", "negative"),
  pr.plot = FALSE,
  smooth.effect = FALSE,
  smooth.ci = FALSE,
  ask = prod(par("mfcol")) < nplots,
  levels = NULL,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  main = NULL,
  lwd = par("lwd"),
  ...
)
}
\arguments{
\item{x}{'medsens' object, typically output from \code{medsens}.}

\item{sens.par}{a character string indicating the sensitivity parameter to be
used. Default plots effects as functions of "rho". See Details.}

\item{r.type}{type of the R square parameter to be used in "R2" plots. If 
"residual", effects are plotted against the proportions of the residual 
variances that are explained by the unobserved confounder.  If "total", the
proportions of the total variances are used as sensitivity parameters. Only
relevant if 'sens.par' is "R2".}

\item{sign.prod}{a value indicating the direction of hypothesized confounding
in the sensitivity analysis. If "positive", the confounder is assumed to 
affect the mediator and outcome variable in the same direction; if 
"negative" the effects are assumed to be in opposite directions.  Only 
relevant if sens.par is set to "R2".}

\item{pr.plot}{a logical value. If 'TRUE', the "proportions mediated" will be
plotted instead of the average causal mediation effects or direct effects. 
Currently only available if the object 'medsens' is based on the linear 
mediator and binary probit outcome models. Default is 'FALSE'.}

\item{smooth.effect}{a logical value indicating whether the estimated 
mediation effects are smoothed via \code{\link{lowess}} before being 
plotted. Default is 'FALSE'.}

\item{smooth.ci}{a logical value indicating whether the confidence bands are 
smoothed via \code{\link{lowess}} before being plotted. Default is 'FALSE'.}

\item{ask}{a logical value. If 'TRUE', the user is asked for input before a 
new figure is plotted.  Default is to ask only if the number of plots on 
current screen is fewer than necessary.}

\item{levels}{vector of levels at which to draw contour lines. Only relevant 
if 'sens.par' is set to "R2". If 'NULL', default values in 
\code{\link{contour.default}} are used.}

\item{xlab}{label for the x axis. Default labels are used if 'NULL'.}

\item{ylab}{label for the y axis. Default labels are used if 'NULL'.}

\item{xlim}{limits of the x axis. If 'NULL' default values are used.}

\item{ylim}{limits of the y axis. If 'NULL' default values are used.}

\item{main}{main title for the plot. If 'NULL', default titles are used.}

\item{lwd}{width of the lines used in graphs.}

\item{...}{additional arguments to be passed to plotting functions.}
}
\description{
This function is used to plot results from the 'medsens' function. Causal 
average mediation effects (as well as average direct effects and proportions 
mediated for selected models) can be plotted against two alternative 
sensitivity parameters.
}
\details{
The sensitivity analysis for causal mediation effects can be
  conducted in terms of two alternative sensitivity parameters, which both
  quantify the degree of violation of the sequential ignorability assumption.
  The "rho" parameter represents the correlation between the two error terms
  of the (latent) linear models for the mediator and outcome variables. A
  large value of rho indicates the existence of important common unobserved
  predictors for both the mediator and outcome and therefore a high degree of
  sequential ignorability violation, while a value close to zero implies
  there is no such confounders.
  
  The resulting "rho" figures plot the estimated true values of ACME (or ADE,
  proportion mediated) against rho, along with the confidence intervals. When
  rho is zero, sequantial ignorability holds, so the estimated value at that 
  point will be equal to the estimate returned by the \code{\link{mediate}}. 
  The confidence level is determined by the 'conf.level' value of the
  original \code{\link{mediate}} object.
  
  The "R2" parameters represent the proportions of the mediator and outcome 
  variances that are explained by an unobserved pre-treatment confounder, 
  thereby indicating the importance of such a confounder in each model.  When
  'r.type' is "residual", the R2 parameters represent the proportions of the 
  residual variances of the mediator and outcome models that become explained
  by the inclusion of the hypothetical pre-treatment confounder.  These are 
  denoted as "R square stars" in Imai, Keele and Yamamoto (2010) and can also
  be specified as "star" or using a numeric value 1 in \code{medsens.plot}. 
  When 'r.type' is "total", the R2s represent the total mediator and outcome 
  variances the unobserved confounder would explain. This option can also be 
  specified using "tilde" or a numeric value 2.
  
  For both types of the "R2" parameters, 'sign.prod' indicates the
  hypothesized direction in which the unobserved confounder affects the
  mediator and outcome. (The name derives from the fact that this direction
  is mathematically represented by the sign of the product of two regression 
  coefficients.) If "positive" (or a numeric value 1) is given, the
  confounder is assumed to affect the mediator and outcome in the same
  direction.  If "negative" (or a numeric value -1), the effect is assumed to
  be in opposite directions.
  
  The resulting contours in the "R2" plots represent the values of the ACME
  (or ADE) for different combinations of the mediator R2 and outcome R2
  values. When both values are zero (the lower-left corner of the plot), the
  unobserved pre-treatment confounder has no effect on either mediator or
  outcome and therefore sequantial ignorability is satisfied.
}
\section{Warning}{
 The 'smooth.effect' and 'smooth.ci' options should be used 
  with caution since the smoothing could affect substantive implications of 
  the graphical analysis in a significant way.
}

\references{
Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. 
  (2014). "mediation: R package for Causal Mediation Analysis", Journal of 
  Statistical Software, Vol. 59, No. 5, pp. 1-38.
  
  Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal 
  Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 
  309-334.
  
  Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and 
  Sensitivity Analysis for Causal Mediation Effects, Statistical Science, 
  Vol. 25, No. 1 (February), pp. 51-71.
  
  Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation 
  Analysis Using R" in Advances in Social Science Research Using R, ed. H. D.
  Vinod New York: Springer.
}
\seealso{
\code{\link{medsens}}, \code{\link{plot}}, \code{\link{contour}}.
}
\author{
Dustin Tingley, Harvard University, 
  \email{dtingley@gov.harvard.edu}; Teppei Yamamoto, Massachusetts Institute
  of Technology, \email{teppei@mit.edu}; Jaquilyn Waddell-Boie, Princeton 
  University, \email{jwaddell@princeton.edu}; Luke Keele, Penn State 
  University, \email{ljk20@psu.edu}; Kosuke Imai, Princeton University, 
  \email{kimai@princeton.edu}.
}
