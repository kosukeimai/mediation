#' mediation: Parametric and non parametric mediation analysis.
#' @docType package
#' @name mediation
#' 
#' @importFrom MASS mvrnorm gamma.shape polr
#' @importFrom Matrix Diagonal
#' @importMethodsFrom Matrix kronecker
#' @importFrom lpSolve lp
#' @importFrom sandwich vcovHC sandwich estfun
#' @importFrom mvtnorm rmvnorm
#' @importFrom Hmisc wtd.var
#' @importFrom lme4 VarCorr
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics abline axis contour lines par plot plot.default points 
#'   polygon segments text title
#' @importFrom stats approx as.formula coef dnorm family formula getCall glm 
#'   integrate lm lowess model.frame model.matrix model.response model.weights 
#'   na.omit pnorm predict printCoefmat qnorm quantile rgamma sd terms update 
#'   update.formula var vcov weighted.mean binomial delete.response median  
#'   plogis poisson rbinom rlogis rmultinom rnorm rpois rt runif rweibull weights
#' @importFrom utils packageDescription
NULL