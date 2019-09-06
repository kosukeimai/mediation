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
#' @param boot a logical value. if \code{FALSE} analytic confidence intervals
#'   based on Aroian (1947) will be used; if \code{TRUE} nonparametric
#'   bootstrap will be used. Default is \code{FALSE}.
#' @param sims number of Monte Carlo draws for nonparametric bootstrap.
#' @param est_se estimate standard errors. Primarily for internal use. Default is \code{TRUE}.
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
#'   \item{boot}{logical, the \code{boot} argument used.}
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
#' @references Aroian, L. A. 1947. The probability function of the product of
#'   two normally distributed variables. *Annals of Mathematical Statistics,*
#'   18, 265-271.
#' @export
#' 
#' @examples
#' # Generate data. We use TSLS to address unobserved confounding (n).
#' set.seed(123)
#' sims <- 1000
#' dat <- data.frame(z = sample(0:1, sims, replace = TRUE), 
#'                   t = sample(0:1, sims, replace = TRUE))
#' dat$n <- rnorm(sims, mean = 1)
#' dat$m <- rnorm(sims, mean = dat$z * 0.3 + dat$t * 0.2 + dat$n * 0.7, sd = 0.2)
#' dat$y <- rnorm(sims, mean = 5 + dat$t + dat$m * (-3) + dat$n, sd = 1)
#' model.m <- lm(m ~ t + z, data = dat)
#' model.y <- lm(y ~ t + m, data = dat)
#' cluster <- factor(sample(1:3, sims, replace = TRUE))
#' med <- mediate_tsls(model.m, model.y, cluster = cluster, treat = "t")
#' summary(med) 

mediate_tsls <- function(model.m, model.y, treat = "treat.name",
                         conf.level = .95, 
                         robustSE = FALSE, cluster = NULL,
                         boot = FALSE, sims = 1000, est_se = TRUE,
                         ...){
  
  if (!inherits(model.m, "lm") | !inherits(model.y, "lm"))
    stop("both mediator and outcome models must be of class `lm'.")
  
  m_var <- all.vars(formula(model.m)[[2]])
  y_var <- all.vars(formula(model.y)[[2]])
  t_var <- treat

  if (length(y_var) > 1L || length(m_var) > 1L)
    stop("Left-hand side of model must only have one variable.")
  
  n_y <- nobs(model.y)
  n_m <- nobs(model.m)  
  
  if (n_y != n_m)
    stop("number of observations in both models must be identical.")
  if (!is.null(cluster)) {
    if (NROW(cluster) != n_y) 
      stop("length of `cluster' must be equal to number of observations in models.")
  } else {
    cluster <- seq(n_y)
  }

  # Update y-model using predicted values of mediator
  .dat <- eval(getCall(model.y)$data)
  .dat <- .dat[names(model.m$fitted.values), ]
  # .dat <- model.frame(model.y)
  .dat[[m_var]] <- predict(model.m)
  mod.y <- my_update(model.y, data = .dat)

  # Point estimates
  d <- coef(mod.y)[m_var] * coef(model.m)[t_var]   # mediation effect
  z <- coef(mod.y)[t_var]                          # direct effect
  tau.coef <- d + z                                # total effect
  nu <- d / tau.coef                               # proportion mediated

  if (!est_se) {
    
    se_d <- se_z <- se_tau <- se_n <- NA
    d.ci <- z.ci <- tau.ci <- n.ci <- NA
    d.p <- z.p <- tau.p <- n.p <- NA
    
  } else {

    if (!boot) {
      sims <- NA
      # Analytic CI
      
      if (!is.null(cluster)) {
        vcv_y <- sandwich::vcovCL(mod.y, cluster = cluster, ...)
        vcv_m <- sandwich::vcovCL(model.m, cluster = cluster, ...)    
      } else if (robustSE) {
        vcv_y <- sandwich::vcovHC(mod.y, ...)
        vcv_m <- sandwich::vcovHC(model.m, ...)    
      } else {
        vcv_y <- vcov(mod.y)
        vcv_m <- vcov(model.m)
      }
      
      se_d <- sqrt(
        coef(mod.y)[m_var]^2 * vcv_m[t_var, t_var] +
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
      Coefs <- c(coef(model.m)[t_var], coef(mod.y)[t_var], coef(mod.y)[m_var])
      Coefs <- setNames(Coefs, c("b2", "b3", "gamma"))
      Sigma <- diag(c(vcv_m[t_var, t_var], diag(vcv_y)[c(t_var, m_var)]))
      Sigma[3,2] <- Sigma[2,3] <- vcv_y[t_var, m_var]
      f <- ~b2 * gamma / (b2 * gamma + b3)
      se_n <- as.vector(delta(f, Coefs, Sigma))
      
      # CI and p-values
      qq <- (1 - conf.level) / 2
      qq <- setNames(c(qq, 1 - qq), c("low", "high"))
      d.ci <- d + qnorm(qq) * se_d
      z.ci <- z + qnorm(qq) * se_z
      tau.ci <- tau.coef + qnorm(qq) * se_tau
      n.ci <- nu + qnorm(qq) * se_n
      d.p <- pnorm(-abs(d), sd = se_d)
      z.p <- pnorm(-abs(z), sd = se_z)
      tau.p <- pnorm(-abs(tau.coef), sd = se_tau)
      n.p <- pnorm(-abs(nu), sd = se_n)
      
    } else {

      # clusters
      # taken from sandwich::vcovBS
      cl <- split(seq_along(cluster), cluster)
      
      # matrix for bootstrap estimates
      cf <- matrix(rep.int(0, 4 * sims), ncol = 4,
                   dimnames = list(NULL, c("delta", "zeta", "tau", "nu")))
      
      ## update on bootstrap samples
      for(i in 1:sims) {
        .subset <- unlist(cl[sample(names(cl), length(cl), replace = TRUE)])
        .dat_y <- eval(getCall(model.y)$data)[.subset, ]
        .dat_m <- eval(getCall(model.m)$data)[.subset, ]
        out <- tryCatch({
          up_y <- my_update(model.y, data = .dat_y)
          up_m <- my_update(model.m, data = .dat_m)
          mediate_tsls(up_m, up_y, treat = treat, cluster = NULL, est_se = FALSE)[c("d1", "z0", "tau.coef", "n0")]
        }, error = function(e) {
          setNames(rep(list(NA), 4), c("d1", "z0", "tau.coef", "n0"))
        })
        cf[i, ] <- unlist(out)
      }
  
      se_d <- sd(cf[, "delta"], na.rm = TRUE)
      se_z <- sd(cf[, "zeta"], na.rm = TRUE)
      se_tau <- sd(cf[, "tau"], na.rm = TRUE)
      se_n <- sd(cf[, "nu"], na.rm = TRUE)
      
      qq <- (1 - conf.level) / 2
      qq <- setNames(c(qq, 1 - qq), c("low", "high"))    
      d.ci <- quantile(cf[, "delta"], qq, na.rm = TRUE)
      z.ci <- quantile(cf[, "zeta"], qq, na.rm = TRUE)
      tau.ci <- quantile(cf[, "tau"], qq, na.rm = TRUE)
      n.ci <- quantile(cf[, "nu"], qq, na.rm = TRUE)
      
      d.p <- pval(cf[, "delta"], d)
      z.p <- pval(cf[, "zeta"], z)
      tau.p <- pval(cf[, "tau"], tau.coef)
      n.p <- pval(cf[, "nu"], nu)
      
    }
        
  }

  # Output
  out <- list(d1 = unname(d), d1.se = se_d, d1.p = d.p, d1.ci = d.ci,
              d0 = unname(d), d0.se = se_d, d0.p = d.p, d0.ci = d.ci,
              z1 = unname(z), z1.se = se_z, z1.p = z.p, z1.ci = z.ci,
              z0 = unname(z), z0.se = se_z, z0.p = z.p, z0.ci = z.ci,
              tau.coef = unname(tau.coef), tau.se = se_tau, 
              tau.ci = tau.ci, tau.p = tau.p,
              n0 = unname(nu), n0.se = se_n, n0.ci = n.ci, n0.p = n.p,
              boot = boot, boot.ci.type = "perc",
              treat = treat, mediator = m_var,
              nobs = nobs(model.y), sims = sims, 
              INT = FALSE, conf.level = conf.level, 
              model.y = model.y, model.m = model.m
  )
  class(out) <- c("mediate", "mediate.tsls")
  return(out)
}

# From https://stackoverflow.com/a/13690928/6455166
# This function addresses issues when `update` is called within a function.
my_update <- function(mod, formula = NULL, data = NULL) {
  call <- getCall(mod)
  if (is.null(call)) {
    stop("Model object does not support updating (no call)", call. = FALSE)
  }
  term <- terms(mod)
  if (is.null(term)) {
    stop("Model object does not support updating (no terms)", call. = FALSE)
  }
  
  if (!is.null(data)) call$data <- data
  if (!is.null(formula)) call$formula <- update.formula(call$formula, formula)
  env <- attr(term, ".Environment")

  eval(call, env, parent.frame())
}