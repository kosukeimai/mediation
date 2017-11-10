#' Two-stage Least Squares Estimation of the Average Causal Mediation Effects
#' 
#' Estimate quantities for causal mediation analysis using an instrumental 
#' variable estimator.
#' 
#' @param model.m a fitted model object for mediator, of class \code{lm}.
#' @param model.y a fitted model object for outcome, of class \code{lm}.
#' @param treat a character string indicating the name of the treatment variable
#'   used in the models.
#' @param conf.level level of the returned two-sided confidence intervals. 
#'   Default is to return the 2.5 and 97.5 percentiles of the simulated 
#'   quantities.
#' @param cluster a variable indicating clusters for standard errors. Note that 
#'   this should be a vector of cluster indicators itself, not a character 
#'   string for the name of the variable.
#' @param robustSE a logical value. If 'TRUE', heteroskedasticity-consistent 
#'   standard errors will be used. Default is 'FALSE'.
#' @param ... other arguments passed to vcovCL in the sandwich package:
#'   typically the \code{type} argument for heteroskedasticity-consistent 
#'   standard errors.
#'   
#' @importFrom stats deriv nobs setNames
#'   
#' @return \code{mediate} returns an object of class "\code{mediate}", 
#'   "\code{mediate.tsls}", a list that contains the components listed below.  
#'   
#'   The function \code{summary}  can be used to obtain a table of the results.  
#'   
#'   \item{d1}{point estimate for average causal mediation effects.}
#'   \item{d1.ci}{confidence intervals for average causal mediation 
#'   effect. The confidence level is set at the value specified in 
#'   'conf.level'.}
#'   \item{z0}{point estimates for average direct effect.}
#'   \item{z0.ci}{confidence intervals for average direct effect.}
#'   \item{z0.p}{two-sided p-values for average causal direct effect.}
#'   \item{n0}{the "proportions mediated", or the size of the average causal 
#'   mediation effect relative to the total effect.}
#'   \item{n0.ci}{confidence intervals for the proportion mediated.}
#'   \item{n0.p}{two-sided p-values for proportion mediated.}
#'   \item{tau.coef}{point estimate for total effect.}
#'   \item{tau.ci}{confidence interval for total effect.}
#'   \item{tau.p}{two-sided p-values for total effect.}
#'   \item{boot}{logical, the \code{boot} argument used. Defaults to false, as 
#'   bootstrap is not yet implemented for this function.}
#'   \item{treat}{a character string indicating the name of the 'treat' variable 
#'   used.}
#'   \item{mediator}{a character string indicating the name of the 'mediator' 
#'   variable used.}
#'   \item{INT}{a logical value indicating whether the model specification 
#'   allows the effects to differ between the treatment and control conditions.}
#'   \item{conf.level}{the confidence level used. }
#'   \item{model.y}{the outcome model used.}
#'   \item{model.m}{the mediator model used.}
#'   \item{nobs}{number of observations in the model frame for 'model.m' and 
#'   'model.y'. May differ from the numbers in the original models input to 
#'   'mediate' if 'dropobs' was 'TRUE'.}
#'   \item{cluster}{the clusters used.}
#'   
#' @export
#' 
#' @examples
#' # Generate data
#' set.seed(123)
#' sims <- 10000
#' dat <- data.frame(z = sample(0:1, sims, replace = TRUE), 
#'                   t = sample(0:1, sims, replace = TRUE))
#' dat$m <- rbinom(sims, size = 1, prob = .1 + dat$z * .3 + dat$t * .2)
#' dat$x <- rnorm(sims, mean = 1)
#' dat$y <- 5 + dat$x + 2 * dat$t + 3 * dat$m + rnorm(sims)
#' 
#' model.m <- lm(m ~ t + z, data = dat)
#' model.y <- lm(y ~ t + m + x, data = dat)
#' cluster <- factor(sample(1:3, sims, replace = TRUE))
#' med <- mediate_tsls(model.m, model.y, cluster = NULL, treat = "t")
#' summary(med) 

mediate_tsls <- function(model.m, model.y, treat = "treat.name",
                         conf.level = .95,
                         robustSE = FALSE, cluster = NULL, 
                         ...){
  
  if (!inherits(model.m, "lm") || !inherits(model.y, "lm"))
    stop("Both mediator and outcome models must be of class `lm'.")
  
  m_var <- setdiff(all.vars(formula(model.m)), labels(terms(model.m)))
  y_var <- setdiff(all.vars(formula(model.y)), labels(terms(model.y)))
  t_var <- treat

  if (length(y_var) > 1L || length(m_var) > 1L)
    stop("Left-hand side of model must only have one variable.")
  
  if (!is.null(cluster) & 
      length(unique(c(nobs(model.m), nobs(model.y), length(cluster)))) != 1)
    stop("Length of `cluster' must be equal to number of observations in 
         models.")

  d <- coef(model.y)[m_var] * coef(model.m)[t_var] # mediation effect  
  z <- coef(model.y)[t_var]                        # direct effect
  tau.coef <- d + z                                # total effect
  nu <- d / tau.coef                               # proportion mediated

  # Calculate uncertainty
  if (!is.null(cluster)) {
    vcv_y <- sandwich::vcovCL(model.y, cluster = cluster, ...)
    vcv_m <- sandwich::vcovCL(model.m, cluster = cluster, ...)    
  } else if (robustSE) {
    vcv_y <- sandwich::vcovHC(model.y, ...)
    vcv_m <- sandwich::vcovHC(model.m, ...)    
  } else {
    vcv_y <- vcov(model.y)
    vcv_m <- vcov(model.m)
  }

  se_d <- sqrt(
    coef(model.y)[m_var]^2 * vcv_m[t_var, t_var] +
      coef(model.m)[t_var]^2 * vcv_y[m_var, m_var] +
      vcv_m[t_var, t_var] * vcv_y[m_var, m_var]
  )
  se_z <- sqrt(vcv_y[t_var, t_var])
  se_tau <- sqrt(
    vcv_y[t_var, t_var] + 
      (se_d)^2 + 
      2 * vcv_y[t_var, m_var] * coef(model.m)[t_var]
  )
  
  # Proportion
  delta <- function(f, B, Sigma) {
    ff <- deriv(f, names(B), func = TRUE)
    x <- do.call(ff, as.list(B))
    grad <- as.matrix(attr(x, "gradient"), nr = 1)
    sqrt(grad %*% Sigma %*% t(grad))
  }
  Coefs <- c(coef(model.m)[t_var], coef(model.y)[t_var], coef(model.y)[m_var])
  Coefs <- setNames(Coefs, c("b2", "b3", "gamma"))
  Sigma <- diag(c(vcv_m[t_var, t_var], diag(vcv_y)[c(t_var, m_var)]))
  Sigma[3,2] <- Sigma[2,3] <- vcv_y[t_var, m_var]
  f <- ~b2 * gamma / (b2 * gamma + b3)
  se_n <- as.vector(delta(f, Coefs, Sigma))
  
  qq <- (1 - conf.level) / 2
  qq <- setNames(c(qq, 1 - qq), c("low", "high"))
  d.ci <- d + qnorm(qq) * se_d
  d.p <- pnorm(-d, sd = se_d)
  z.ci <- z + qnorm(qq) * se_z
  z.p <- pnorm(-z, sd = se_z)
  tau.ci <- tau.coef + qnorm(qq) * se_tau
  tau.p <- pnorm(-tau.coef, sd = se_tau)
  n.ci <- nu + qnorm(qq) * se_n
  n.p <- pnorm(-nu, sd = se_n)
  
  # Output
  out <- list(d1 = unname(d), d1.se = se_d, d1.p = d.p, d1.ci = d.ci,
              z0 = unname(z), z0.se = se_z, z0.p = z.p, z0.ci = z.ci,
              tau.coef = unname(tau.coef), tau.ci = tau.ci, tau.p = tau.p,
              n0 = unname(nu), n0.ci = n.ci, n0.p = n.p,
              boot = FALSE, treat = treat, mediator = m_var,
              nobs = nobs(model.y), sims = NA, 
              INT = FALSE, conf.level = conf.level, 
              model.y = model.y, model.m = model.m
  )
  class(out) <- c("mediate", "mediate.tsls")
  return(out)
}