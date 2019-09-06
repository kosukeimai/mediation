##################################################################
# Sensitivity Analysis for Multiple Mediators
##################################################################

#' Estimation and Sensitivity Analysis for Multiple Causal Mechanisms
#' 
#' 'multimed' is used for causal mediation analysis when post-treatment 
#' mediator-outcome confounders, or alternative mediators causally preceding the
#' mediator of interest, exist in the hypothesized causal mechanisms. It 
#' estimates the average causal mediation effects (indirect effects) and the 
#' average direct effects under the homogeneous interaction assumption based on 
#' a varying-coefficient linear structural equation model. The function also 
#' performs sensitivity analysis with respect to the violation of the homogenous
#' interaction assumption. The function can be used for both the single
#' experiment design and the parallel design.
#' 
#' @details This function implements the framework proposed by Imai and Yamamoto
#'   (2012) for the estimation and sensitivity analysis for multiple causal
#'   mechanisms. It estimates the average causal mediation effects (indirect
#'   effects) with respect to the mediator of interest ('med.main'), i.e., the
#'   portion of the treatment effect on the outcome that is transmitted through
#'   that mediator, as well as the average direct effects, i.e., the portion of
#'   the treatment effect on the outcome that is not transmitted through the
#'   main mediator. Unlike the "standard" causal mediation analysis implemented
#'   by \code{\link{mediate}} and \code{\link{medsens}}, this framework allows
#'   the existence of post-treatment covariates that confound the relationship 
#'   between the main mediator and the outcome, or equivalently, alternative 
#'   mediators ('med.alt') that causally precede the main mediator.
#'   
#'   When the parallel design was used for the experiment (i.e. when the 
#'   experiment contained an additional randomly assigned group for which both 
#'   the treatment and the mediator were randomized), there is no need to
#'   specify a particular post-treatment confounder, for any such confounder
#'   (observed or unobserved) is allowed to exist by virtue of the design.
#'   Similarly, no observed covariates need to be included. The function instead
#'   requires an additional variable ('experiment') indicating whether the
#'   mediator was randomly manipulated for the unit.
#'   
#'   The estimation and sensitivity analysis are both based on a 
#'   varying-coefficient linear structural equations model, which assumes 
#'   additivity but allows for an arbitrary degree of heterogeneity in model 
#'   coefficients across units and thus is substantially more flexible than a 
#'   traditional SEM framework. For details see Imai and Yamamoto (2012).
#'   
#'   The function produces two sets of results. First, point estimates of the 
#'   average causal mediation effects and the average direct effects are 
#'   calculated, along with their (percentile) bootstrap confidence intervals. 
#'   These estimates are based on the "homogeneous interaction" assumption, or 
#'   the assumption that the degree of treatment-mediator interaction is
#'   constant across all units. The estimated total treatment effect is also
#'   reported.
#'   
#'   Second, the bounds on the average causal mediation effects and the average 
#'   direct effects are also estimated and computed for various degrees of 
#'   interaction heterogeneity (i.e., violation of the identification 
#'   assumption), which are represented by the values of three alternative 
#'   sensitivity parameters. These parameters are: (1) sigma, the standard 
#'   deviation of the (varying) regression coefficient on the interaction term, 
#'   (2) R square star, the proportion of the residual variance that would be 
#'   explained by an additional term for interaction heterogeneity, and (3) R 
#'   square tilde, the proportion of the total variance explained by such a
#'   term. The confidence region is also calculated, using the Imbens and Manski
#'   (2004) formula with bootstrap standard errors. Further details are given in
#'   the above reference.
#'   
#'   Note that rows with missing values will be omitted from the calculation of 
#'   the results. Also note that the treatment variable must be a numeric vector
#'   of 1 and 0 and that both mediators and outcome variable must be numeric.
#'   The pre-treatment covariates can be of any type that \code{\link{lm}} can
#'   handle as predictors.
#'   
#' @param outcome name of the outcome variable in 'data'.
#' @param med.main name of the mediator of interest. Under the parallel design 
#'   this is the only mediator variable used in the estimation.
#' @param med.alt vector of character strings indicating the names of the 
#'   post-treatment confounders, i.e., the alternative mediators affecting both 
#'   the main mediator and outcome. Not needed for the parallel design.
#' @param treat name of the treatment variable in 'data'.
#' @param covariates vector of character strings representing the names of the 
#'   pre-treatment covariates. Cannot be used for the parallel design.
#' @param experiment name of the binary indicator whether 'med.main' is randomly
#'   manipulated or not. Only used under the parallel design.
#' @param data a data frame containing all the above variables.
#' @param design experimental design used. Defaults to 'single'.
#' @param sims number of bootstrap samples used for the calculation of 
#'   confidence intervals.
#' @param R2.by increment for the "R square tilde" parameter, i.e. the 
#'   sensitivity parameter representing the proportion of residual outcome 
#'   variance explained by heterogeneity in treatment-mediator interactions.
#'   Must be a numeric value between 0 and 1.
#' @param conf.level level to be used for confidence intervals.
#' @param weight name of the weights in 'data'.
#' 
#' @return \code{multimed} returns an object of class "\code{multimed}", a list 
#'   contains the following components. The object can be passed to the 
#'   \code{summary} and \code{plot} method functions for a summary table and a 
#'   graphical summary.
#'   
#'   \item{sigma}{ values of the sigma sensitivity parameter at which the bounds 
#'   and confidence intervals are evaluated.}
#'   \item{R2tilde}{ values of the R square tilde parameter.}
#'   \item{R2star}{ values of the R square star parameter.}
#'   \item{d1.lb, d0.lb, d.ave.lb}{ lower bounds on the average causal mediation 
#'   effects under treatment, control, and the simple average of the two, 
#'   respectively, corresponding to the values of the sensitivity parameters 
#'   listed above. Note that the first elements of these vectors equal the point 
#'   estimates under the homogeneous interaction assumption.}
#'   \item{d1.ub, d0.ub, d.ave.ub}{ upper bounds on the average causal mediation 
#'   effects.}
#'   \item{d1.ci, d0.ci, d.ave.ci}{ confidence intervals for the average causal
#'   mediation effects at different values of the sensitivity parameters.}
#'   \item{z1.lb, z0.lb, z.ave.lb}{ lower bounds on the average direct effects 
#'   under treatment, control, and the simple average of the two, respectively, 
#'   corresponding to the values of the sensitivity parameters listed above. 
#'   Note that the first elements of these vectors equal the point estimates 
#'   under the homogeneous interaction assumption.}
#'   \item{z1.ub, z0.ub, z.ave.ub}{ upper bounds on the average direct effects.}
#'   \item{z1.ci, z0.ci, z.ave.ci}{ confidence intervals for the average direct 
#'   effects at different values of the sensitivity parameters.}
#'   \item{tau}{ point estimate of the total treatment effect.}
#'   \item{tau.ci}{ confidence interval for the total treatment effect.}
#'   \item{conf.level}{ confidence level used for the calculation of the 
#'   confidence intervals.}
#'   
#' @author Teppei Yamamoto, Massachusetts Institute of Technology, 
#'   \email{teppei@@mit.edu}
#'   
#' @seealso \code{\link{plot.multimed}}
#' 
#' @references Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. 
#'   (2014). "mediation: R package for Causal Mediation Analysis", Journal of 
#'   Statistical Software, Vol. 59, No. 5, pp. 1-38.
#'   
#'   Imai, K. and Yamamoto, T. (2012) Identification and Sensitivity Analysis
#'   for Multiple Causal Mechanisms: Revisiting Evidence from Framing
#'   Experiments, Unpublished manuscript.
#'   
#' @export
#' @examples
#' \dontrun{
#' # Replicates Figure 3 (right column) of Imai and Yamamoto (2012)
#' # Note: # of bootstrap samples set low for quick illustration
#' 
#' data(framing)
#' Xnames <- c("age", "educ", "gender", "income")
#' res <- multimed("immigr", "emo", "p_harm", "treat", Xnames, 
#'                data = framing, design = "single", sims = 10)
#' summary(res)
#' plot(res, type = "point")
#' plot(res, type = c("sigma", "R2-total"), tgroup = "average")
#' 
#' # Parallel design example using the simulated data of Imai, Tingley and Yamamoto (2012)
#' 
#' data(boundsdata)
#' res.para <- multimed(outcome = "out", med.main = "med", treat = "ttt", experiment = "manip",
#' 					 data = boundsdata, design = "parallel", sims = 10)
#' summary(res.para)
#' plot(res.para, tg = "av")
#' }
multimed <- function(outcome, med.main, med.alt = NULL, treat,
                     covariates = NULL, experiment = NULL,
                     data, design = c("single", "parallel"),
                     sims = 1000, R2.by = 0.01, conf.level = 0.95, weight = NULL){

    if(!is.null(weight) & !is.character(weight)){
        stop("weight must be character string for weights in 'data'")
    } else {
        varnames <- c(outcome, treat, med.main, med.alt, experiment, covariates, weight)
    }
    n.w <- length(med.alt)
    data <- na.omit(data[,varnames])
    design <- match.arg(design)

  if(design == "single"){
    if (is.null(covariates)){
      f.y <- as.formula(paste(outcome, "~", treat, "+", med.main, "+",
                              paste(treat, med.main, sep=":"), "+",
                              paste(med.alt, collapse=" + "), "+",
                              paste(treat, med.alt, sep=":", collapse=" + ")))
      f.ytot <- as.formula(paste(outcome, "~", treat))
      f.m <- as.formula(paste(med.main, "~", treat))
      f.w <- as.list(rep(NA, n.w))
      for(k in 1:n.w){
        f.w[[k]] <- as.formula(paste(med.alt[k], "~", treat))
      }
    } else {
      f.y <- as.formula(paste(outcome, "~", treat, "+", med.main, "+",
                              paste(treat, med.main, sep=":"), "+",
                              paste(med.alt, collapse=" + "), "+",
                              paste(treat, med.alt, sep=":", collapse=" + "), "+",
                              paste(covariates, collapse = " + ")))
      f.ytot <- as.formula(paste(outcome, "~", treat, "+",
                                 paste(covariates, collapse = " + ")))
      f.m <- as.formula(paste(med.main, "~", treat, "+", paste(covariates, collapse = " + ")))
      f.w <- as.list(rep(NA, n.w))
      for(k in 1:n.w){
        f.w[[k]] <- as.formula(paste(med.alt[k], "~", treat, "+", paste(covariates, collapse = " + ")))
      }
    }
  } else if (design == "parallel"){
    if (!is.null(covariates)){
      warning("covariates currently cannot be used for the parallel design; option ignored")
    }
    f.y <- as.formula(paste(outcome, "~", treat, "+", med.main, "+",
                            paste(treat, med.main, sep=":")))
    f.ytot <- as.formula(paste(outcome, "~", treat))
    f.m <- as.formula(paste(med.main, "~", treat))
  } else {
    stop("design unsupported by the multimed procedure")
  }

  # Compute sensitivity parameters
  R2.s <- seq(0, 1, by = R2.by)

  data.1 <- switch(design,
                   single = data,
                   parallel = subset(data, data[[experiment]] == 1))

    if(is.null(weight)){
        ETM2 <- mean(data.1[,treat] * data.1[,med.main]^2)
        model.y <- lm(f.y, data=data.1)
        VY <- var(data.1[,outcome])
    } else {
        ETM2 <- weighted.mean(data.1[,treat] * data.1[,med.main]^2, w = data.1[, weight])
        model.y <- lm(f.y, data=data.1, weights = data.1[, weight])
        VY <- Hmisc::wtd.var(data.1[, outcome], weights = data.1[, weight])
    }
    sigma <- summary(model.y)$sigma * sqrt(R2.s/ETM2)

    R2.t <- ETM2 * sigma^2/VY

  # Bootstrap ACME values

    ACME.1.lo <- ACME.0.lo <- ACME.1.up <- ACME.0.up <-
        ACME.ave.lo <- ACME.ave.up <- matrix(NA, nrow=length(sigma), ncol=sims)
    ADE.1.lo <- ADE.0.lo <- ADE.1.up <- ADE.0.up <-
        ADE.ave.lo <- ADE.ave.up <- matrix(NA, nrow=length(sigma), ncol=sims)
    tau <- rep(NA, length = sims)
  for(b in 1:(sims+1)){
      ## Resample
      data.b <- data[sample(1:nrow(data), nrow(data), replace=TRUE),]
      if(b == sims + 1){  # final iteration is original data
          data.b <- data
      }

    # Subset data for parallel design
    data.b.0 <- switch(design,
                       single = data.b,
                       parallel = subset(data.b, data.b[[experiment]] == 0))
    data.b.1 <- switch(design,
                       single = data.b,
                       parallel = subset(data.b, data.b[[experiment]] == 1))

      ## Fit models
      if(is.null(weight)){
          model.y <- lm(f.y, data=data.b.1)
          model.ytot <- lm(f.ytot, data=data.b.0)
          model.m <- lm(f.m, data=data.b.0)
          if(design == "single"){
              model.w <- as.list(rep(NA, n.w))
              for(k in 1:n.w){
                  model.w[[k]] <- lm(f.w[[k]], data=data.b)
              }
          }
      } else {
          wgt <- data.b[, weight]
          wgt1 <- data.b.1[, weight]
          wgt0 <- data.b.0[, weight]

          model.y <- lm(f.y, weights = wgt1, data=data.b.1)
          model.ytot <- lm(f.ytot, weights = wgt0, data=data.b.0)
          model.m <- lm(f.m, weights = wgt0, data=data.b.0)
          if(design == "single"){
              model.w <- as.list(rep(NA, n.w))
              for(k in 1:n.w){
                  model.w[[k]] <- lm(f.w[[k]], weights = wgt, data=data.b)
              }
          }
      }

    beta3 <- coef(model.y)[treat]
    kappa <- coef(model.y)[paste(treat, med.main, sep=":")]
    if(design == "single"){
      xi3 <- coef(model.y)[med.alt]
      mu3 <- coef(model.y)[paste(treat, med.alt, sep=":")]
    }

    # E(M|T=t)
      mf.m1 <- mf.m0 <- mf.m <- model.frame(model.m)
      mf.m1[,treat] <- 1
      mf.m0[,treat] <- 0
      if(is.null(weight)){
          EM.1 <- mean(predict(model.m, mf.m1))
          EM.0 <- mean(predict(model.m, mf.m0))
      } else { ### model.m from data.b.0 -> wgt0
          EM.1 <- weighted.mean(predict(model.m, mf.m1), w = wgt0)
          EM.0 <- weighted.mean(predict(model.m, mf.m0), w = wgt0)
      }

     # V(M|T=t)
      if(is.null(weight)){
          VM.1 <- sum(model.m$residuals[mf.m[,treat]==1]^2)/(sum(mf.m[,treat]) - length(coef(model.m)))
          VM.0 <- sum(model.m$residuals[mf.m[,treat]==0]^2)/(sum(1-mf.m[,treat]) - length(coef(model.m)))
      } else { ### mf.m from model.m, which is from data.b.0 -> wgt0
          VM.1 <- Hmisc::wtd.var(model.m$residuals[mf.m[,treat]==1], weights = wgt0[mf.m[, treat] == 1])
          VM.0 <- Hmisc::wtd.var(model.m$residuals[mf.m[,treat]==0], weights = wgt0[mf.m[, treat] == 0])
      }

    # E(W|T=t)
    if(design == "single"){
      mf.w1 <- mf.w0 <- as.list(rep(NA, n.w))
      EW.1 <- EW.0 <- rep(NA, n.w)
      for(k in 1:n.w){
        mf.w1[[k]] <- mf.w0[[k]] <- model.frame(model.w[[k]])
        mf.w1[[k]][,treat] <- 1
        mf.w0[[k]][,treat] <- 0
        if(is.null(weight)){
            EW.1[k] <- mean(predict(model.w[[k]], mf.w1[[k]]))
            EW.0[k] <- mean(predict(model.w[[k]], mf.w0[[k]]))
        } else { ### model.w from data.b -> wgt
            EW.1[k] <- weighted.mean(predict(model.w[[k]], mf.w1[[k]]), w = wgt)
            EW.0[k] <- weighted.mean(predict(model.w[[k]], mf.w0[[k]]), w = wgt)
        }
      }
    }

    ## Bounds
    # ACME & ADE
    if(b == sims + 1){
      tau.o <- coef(model.ytot)[treat]
      WXterms <- switch(design,
                        single = (xi3 + mu3) %*% EW.1 - xi3 %*% EW.0,
                        parallel = 0)

      ADE.0.up.o <- beta3 + kappa*EM.0 + sigma*sqrt(VM.0) + as.vector(WXterms)
      ADE.1.up.o <- beta3 + kappa*EM.1 + sigma*sqrt(VM.1) + as.vector(WXterms)
      ADE.0.lo.o <- beta3 + kappa*EM.0 - sigma*sqrt(VM.0) + as.vector(WXterms)
      ADE.1.lo.o <- beta3 + kappa*EM.1 - sigma*sqrt(VM.1) + as.vector(WXterms)
      
      ACME.1.lo.o <- tau.o - ADE.0.up.o
      ACME.0.lo.o <- tau.o - ADE.1.up.o
      ACME.1.up.o <- tau.o - ADE.0.lo.o
      ACME.0.up.o <- tau.o - ADE.1.lo.o
      
      if(is.null(weight)){
          P <- mean(data.b[,treat])
      } else { ### data.b -> wgt
          P <- weighted.mean(data.b[,treat], w = wgt)
      }
      ACME.ave.lo.o <- P * ACME.1.lo.o + (1-P) * ACME.0.lo.o
      ACME.ave.up.o <- P * ACME.1.up.o + (1-P) * ACME.0.up.o

      ADE.ave.lo.o <- P * ADE.1.lo.o + (1-P) * ADE.0.lo.o
      ADE.ave.up.o <- P * ADE.1.up.o + (1-P) * ADE.0.up.o
      
  } else {
      tau[b] <- coef(model.ytot)[treat]
      WXterms <- switch(design,
                        single = (xi3 + mu3) %*% EW.1 - xi3 %*% EW.0,
                        parallel = 0)
      ADE.0.up[,b] <- beta3 + kappa*EM.0 + sigma*sqrt(VM.0) + as.vector(WXterms)
      ADE.1.up[,b] <- beta3 + kappa*EM.1 + sigma*sqrt(VM.1) + as.vector(WXterms)
      ADE.0.lo[,b] <- beta3 + kappa*EM.0 - sigma*sqrt(VM.0) + as.vector(WXterms)
      ADE.1.lo[,b] <- beta3 + kappa*EM.1 - sigma*sqrt(VM.1) + as.vector(WXterms)
      
      ACME.1.lo[,b] <- tau[b] - ADE.0.up[,b]
      ACME.0.lo[,b] <- tau[b] - ADE.1.up[,b]
      ACME.1.up[,b] <- tau[b] - ADE.0.lo[,b]
      ACME.0.up[,b] <- tau[b] - ADE.1.lo[,b]
      
      if(is.null(weight)){
          P <- mean(data.b[,treat])
      } else { ### data.b -> wgt
          P <- weighted.mean(data.b[,treat], w = wgt)
      }
      ACME.ave.lo[,b] <- P * ACME.1.lo[,b] + (1-P) * ACME.0.lo[,b]
      ACME.ave.up[,b] <- P * ACME.1.up[,b] + (1-P) * ACME.0.up[,b]

      ADE.ave.lo[,b] <- P * ADE.1.lo[,b] + (1-P) * ADE.0.lo[,b]
      ADE.ave.up[,b] <- P * ADE.1.up[,b] + (1-P) * ADE.0.up[,b]
    }
  }

    ACME.ave.lo.var <- apply(ACME.ave.lo, 1, var)
    ACME.1.lo.var <- apply(ACME.1.lo, 1, var)
    ACME.0.lo.var <- apply(ACME.0.lo, 1, var)
    ACME.ave.up.var <- apply(ACME.ave.up, 1, var)
    ACME.1.up.var <- apply(ACME.1.up, 1, var)
    ACME.0.up.var <- apply(ACME.0.up, 1, var)

    ADE.ave.lo.var <- apply(ADE.ave.lo, 1, var)
    ADE.1.lo.var <- apply(ADE.1.lo, 1, var)
    ADE.0.lo.var <- apply(ADE.0.lo, 1, var)
    ADE.ave.up.var <- apply(ADE.ave.up, 1, var)
    ADE.1.up.var <- apply(ADE.1.up, 1, var)
    ADE.0.up.var <- apply(ADE.0.up, 1, var)

    ACME.ave.CI <- ACME.1.CI <- ACME.0.CI <- matrix(NA, nrow=2, ncol=length(sigma))
    ADE.ave.CI <- ADE.1.CI <- ADE.0.CI <- matrix(NA, nrow=2, ncol=length(sigma))

    for(i in 1:length(sigma)){
        ACME.ave.CI[,i] <- IMCI(ACME.ave.up.o[i], ACME.ave.lo.o[i],
                                ACME.ave.up.var[i], ACME.ave.lo.var[i], conf.level = conf.level)$ci
        ACME.1.CI[,i] <- IMCI(ACME.1.up.o[i], ACME.1.lo.o[i],
                              ACME.1.up.var[i], ACME.1.lo.var[i], conf.level = conf.level)$ci
        ACME.0.CI[,i] <- IMCI(ACME.0.up.o[i], ACME.0.lo.o[i],
                              ACME.0.up.var[i], ACME.0.lo.var[i], conf.level = conf.level)$ci

        ADE.ave.CI[,i] <- IMCI(ADE.ave.up.o[i], ADE.ave.lo.o[i],
                               ADE.ave.up.var[i], ADE.ave.lo.var[i], conf.level = conf.level)$ci
        ADE.1.CI[,i] <- IMCI(ADE.1.up.o[i], ADE.1.lo.o[i],
                             ADE.1.up.var[i], ADE.1.lo.var[i], conf.level = conf.level)$ci
        ADE.0.CI[,i] <- IMCI(ADE.0.up.o[i], ADE.0.lo.o[i],
                             ADE.0.up.var[i], ADE.0.lo.var[i], conf.level = conf.level)$ci
    }
    tau.CI <- quantile(tau, probs = c((1-conf.level)/2, (1+conf.level)/2), na.rm = TRUE)
    out <- list(sigma = sigma, R2tilde = R2.t, R2star = R2.s, tau = tau.o, tau.ci = tau.CI,
                d1.ci = ACME.1.CI, d0.ci = ACME.0.CI, d.ave.ci = ACME.ave.CI,
                d1.lb = ACME.1.lo.o, d0.lb = ACME.0.lo.o, d.ave.lb = ACME.ave.lo.o,
                d1.ub = ACME.1.up.o, d0.ub = ACME.0.up.o, d.ave.ub = ACME.ave.up.o,
                z1.ci = ADE.1.CI, z0.ci = ADE.0.CI, z.ave.ci = ADE.ave.CI,
                z1.lb = ADE.1.lo.o, z0.lb = ADE.0.lo.o, z.ave.lb = ADE.ave.lo.o,
                z1.ub = ADE.1.up.o, z0.ub = ADE.0.up.o, z.ave.ub = ADE.ave.up.o,
                conf.level = conf.level)
    class(out) <- "multimed"
    out
}


## Calculates Imbens-Manski confidence set for nonparametric bounds
IMCI <- function(upper, lower, var.upper, var.lower,
                 conf.level){
  A <- (upper-lower)/sqrt(max(var.upper,var.lower))
  C <- seq(0,10,by=.001)
  const <- abs(pnorm(C + A) - pnorm(-C) - conf.level)
  Cn <- C[const==min(const)]
  ci <- c(0,0)
  names(ci) <- c("lower","upper")
  ci <- c(lower - Cn*sqrt(var.lower), upper + Cn*sqrt(var.upper))
  if (is.na(ci[1])) ci <- c(NA,NA)
  list(ci=ci, conf.level=conf.level)
}

## Summary

#' Summarizing Output from Mediation Analysis with Multiple Mechanisms
#' 
#' Function to report results from the \code{\link{multimed}} function.
#' 
#' @aliases summary.multimed print.summary.multimed
#' 
#' @param object object of class \code{multimed}, typically output from the 
#'   \code{multimed} funciton.
#' @param x output from the summary function.
#' @param ...  additional arguments affecting the summary produced.
#' 
#' @author Teppei Yamamoto, Massachusetts Institute of Technology, 
#'   \email{teppei@@mit.edu}.
#'   
#' @seealso \code{\link{multimed}}, \code{\link{plot.multimed}}
#' 
#' @references Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. 
#'   (2014). "mediation: R package for Causal Mediation Analysis", Journal of 
#'   Statistical Software, Vol. 59, No. 5, pp. 1-38.
#'   
#'   Imai, K. and Yamamoto, T. (2012) Identification and Sensitivity Analysis
#'   for Multiple Causal Mechanisms: Revisiting Evidence from Framing
#'   Experiments, typescript.
#'   
#' @examples
#' \dontrun{
#' # Replicates Figure 3 (right column) of Imai and Yamamoto (2012)
#' # Note: # of bootstrap samples set low for quick illustration
#' 
#' data(framing)
#' Xnames <- c("age", "educ", "gender", "income")
#' res <- multimed("immigr", "emo", "p_harm", "treat", Xnames, 
#'                data = framing, sims = 10)
#' summary(res)
#' plot(res, type = "point")
#' plot(res, type = c("sigma", "R2-total"), tgroup = "average")
#' }
#' 
#' @export
summary.multimed <- function(object, ...){
  structure(object, class = c("summary.multimed", class(object)))
}

#' @rdname summary.multimed
#' @export
print.summary.multimed <- function(x, ...){
  clp <- 100 * x$conf.level
  cat("\n")
  cat("Causal Mediation Analysis with Confounding by an Alternative Mechanism\n\n")
  cat("Estimates under the Homogeneous Interaction Assumption:\n")

  cmat <- c(x$d1.lb[1], x$d0.lb[1], x$d.ave.lb[1], x$z1.lb[1], x$z0.lb[1], x$z.ave.lb[1], x$tau)
  cmat <- cbind(cmat, rbind(x$d1.ci[,1], x$d0.ci[,1], x$d.ave.ci[,1], x$z1.ci[,1], x$z0.ci[,1], x$z.ave.ci[,1], x$tau.ci))
  colnames(cmat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                          paste(clp, "% CI Upper", sep=""))
  rownames(cmat) <- c("ACME (treated)", "ACME (control)", "ACME (average)", "ADE (treated)", "ADE (control)", "ADE (average)", "Total Effect")
  printCoefmat(cmat[,1:3], digits=3)
  cat("\n")

  cat("Sensitivity Analysis: \n")
  cat("Values of the sensitivity parameters at which ACME first crosses zero:\n")
  ind.d1.b <- sum(sign(x$d1.lb) * sign(x$d1.ub) > 0) + 1
  ind.d1.c <- sum(sign(x$d1.ci[1,]) * sign(x$d1.ci[2,]) > 0) + 1
  ind.d0.b <- sum(sign(x$d0.lb) * sign(x$d0.ub) > 0) + 1
  ind.d0.c <- sum(sign(x$d0.ci[1,]) * sign(x$d0.ci[2,]) > 0) + 1
  ind.d.ave.b <- sum(sign(x$d.ave.lb) * sign(x$d.ave.ub) > 0) + 1
  ind.d.ave.c <- sum(sign(x$d.ave.ci[1,]) * sign(x$d.ave.ci[2,]) > 0) + 1
  smat <- c(x$sigma[ind.d1.b], x$sigma[ind.d1.c],
            x$R2star[ind.d1.b], x$R2star[ind.d1.c], x$R2tilde[ind.d1.b], x$R2tilde[ind.d1.c])
  smat <- rbind(smat, c(x$sigma[ind.d0.b], x$sigma[ind.d0.c],
                        x$R2star[ind.d0.b], x$R2star[ind.d0.c], x$R2tilde[ind.d0.b], x$R2tilde[ind.d0.c]))
  smat <- rbind(smat, c(x$sigma[ind.d.ave.b], x$sigma[ind.d.ave.c],
                        x$R2star[ind.d.ave.b], x$R2star[ind.d.ave.c], x$R2tilde[ind.d.ave.b], x$R2tilde[ind.d.ave.c]))
  colnames(smat) <- c("sigma(bounds)", "sigma(CI)", "R2s(bounds)", "R2s(CI)", "R2t(bounds)", "R2t(CI)")
  rownames(smat) <- c("ACME (treated)", "ACME (control)", "ACME (average)")
  printCoefmat(smat[,1:6], digits=3)

  cat("Values of the sensitivity parameters at which ADE first crosses zero:\n")
  ind.z1.b <- sum(sign(x$z1.lb) * sign(x$z1.ub) > 0) + 1
  ind.z1.c <- sum(sign(x$z1.ci[1,]) * sign(x$z1.ci[2,]) > 0) + 1
  ind.z0.b <- sum(sign(x$z0.lb) * sign(x$z0.ub) > 0) + 1
  ind.z0.c <- sum(sign(x$z0.ci[1,]) * sign(x$z0.ci[2,]) > 0) + 1
  ind.z.ave.b <- sum(sign(x$z.ave.lb) * sign(x$z.ave.ub) > 0) + 1
  ind.z.ave.c <- sum(sign(x$z.ave.ci[1,]) * sign(x$z.ave.ci[2,]) > 0) + 1
  smat <- c(x$sigma[ind.z1.b], x$sigma[ind.z1.c],
            x$R2star[ind.z1.b], x$R2star[ind.z1.c], x$R2tilde[ind.z1.b], x$R2tilde[ind.z1.c])
  smat <- rbind(smat, c(x$sigma[ind.z0.b], x$sigma[ind.z0.c],
                        x$R2star[ind.z0.b], x$R2star[ind.z0.c], x$R2tilde[ind.z0.b], x$R2tilde[ind.z0.c]))
  smat <- rbind(smat, c(x$sigma[ind.z.ave.b], x$sigma[ind.z.ave.c],
                        x$R2star[ind.z.ave.b], x$R2star[ind.z.ave.c], x$R2tilde[ind.z.ave.b], x$R2tilde[ind.z.ave.c]))
  colnames(smat) <- c("sigma(bounds)", "sigma(CI)", "R2s(bounds)", "R2s(CI)", "R2t(bounds)", "R2t(CI)")
  rownames(smat) <- c("ADE (treated)", "ADE (control)", "ADE (average)")
  printCoefmat(smat[,1:6], digits=3)
  cat("\n")

  invisible(x)
}


## Plot

#' Plotting the Results of Causal Mediation Analysis for Multiple Mechanisms
#' 
#' Function to plot results from \code{multimed}. Most standard plotting options
#' are available.
#' 
#' @details 'type' and 'tgroup' can contain multiple character strings, in which
#'   case multiple plots are produced. For the use of graphical parameters see 
#'   \code{\link{plot}} and the links it contains.
#'   
#' @param x object of class \code{multimed}, typically output from the 
#'   \code{multimed} funciton.
#' @param type type of plot(s) required. The default is to produce all, i.e., 
#'   the point estimates of the effects under the homogenous interaction 
#'   assumpton ("point") and bounds as function of the sigma ("sigma"), R square
#'   star ("R2-residual") and R square tilde ("R2-total") parameters.
#' @param tgroup treatment group(s) for which the estimates are produced. The 
#'   default is to plot all, i.e., the average causal mediation effect when 
#'   treated ("treated"), control ("control") and the simple average of these 
#'   two effects ("average").
#' @param effect.type a vector indicating which quantities of interest to plot 
#'   for the point estimates. Only plotting total effects is not available.
#' @param ask a logical value. If 'TRUE', the user is asked for input before a 
#'   new figure is plotted.  Default is to ask only if the number of plots on 
#'   current screen is fewer than necessary.
#' @param xlab label for the x axis. Default labels are used if 'NULL'.
#' @param ylab label for the y axis. Default labels are used if 'NULL'.
#' @param xlim limits of the x axis. If 'NULL' default values are used.
#' @param ylim limits of the y axis. If 'NULL' default values are used.
#' @param main main title for the plot. If 'NULL', default titles are used.
#' @param lwd width of the lines used in graphs. For the "point" plot this is 
#'   the width of confidence bars. For sensitivity plots this is the width of 
#'   the lines for the bounds.
#' @param pch plotting points used for the "point" plots.
#' @param cex magnification factor for the plotting points in the "point" plots.
#' @param las style of the y axis labels in the "point" plots.
#' @param col.eff color of the points in the "point" plots and/or the bounds in 
#'   sensitivity plots.
#' @param col.cbar color of the confidence bars in the "point" plots.
#' @param col.creg color of the confidence regions in sensitivity plots.
#' @param ...  additional arguments to be passed to plotting functions.
#' 
#' @author Teppei Yamamoto, Massachusetts Institute of Technology, 
#'   \email{teppei@@mit.edu}.
#'   
#' @seealso \code{\link{multimed}}, \code{\link{plot}}
#' 
#' @references Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. 
#'   (2014). "mediation: R package for Causal Mediation Analysis", Journal of 
#'   Statistical Software, Vol. 59, No. 5, pp. 1-38.
#'   
#'   Imai, K. and Yamamoto, T. (2012) Identification and Sensitivity Analysis 
#'   for Multiple Causal Mechanisms: Revisiting Evidence from Framing 
#'   Experiments, typescript.
#'   
#' @examples
#' \dontrun{
#' 
#' # Replicates Figure 3 (right column) of Imai and Yamamoto (2012)
#' # Note: # of bootstrap samples set low for quick illustration
#' 
#' data(framing)
#' Xnames <- c("age", "educ", "gender", "income")
#' res <- multimed("immigr", "emo", "p_harm", "treat", Xnames,
#'                data = framing, sims = 10)
#' summary(res)
#' plot(res, type = "point")
#' plot(res, type = c("sigma", "R2-total"), tgroup = "average")
#' 
#' }
#' 
#' @export
plot.multimed <- function(x, type = c("point", "sigma", "R2-residual", "R2-total"),
                          tgroup = c("average", "treated", "control"),
                          effect.type = c("indirect", "direct", "total"),
                          ask = prod(par("mfcol")) < nplots,
                          xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, main = NULL,
                          lwd = par("lwd"), pch = par("pch"), cex = par("cex"), las = par("las"),
                          col.eff = "black", col.cbar = "black", col.creg = "gray", ...){

  type <- match.arg(type, several.ok = TRUE)
  tgroup <- match.arg(tgroup, several.ok = TRUE)
  effect.type <- match.arg(effect.type, several.ok=TRUE)

  IND <- "indirect" %in% effect.type
  DIR <- "direct" %in% effect.type
  TOT <- "total" %in% effect.type

  if(!IND && !DIR && TOT) {
    stop("No support for only plotting total effect.")
  }

  show.point <- "point" %in% type
  nplots <- show.point + (length(type) - show.point) * length(tgroup)
  if(ask){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  eff.up <- eff.lo <- ci.up <- ci.lo <- c()
  d.eff.lo <- d.eff.up <- d.ci.lo <- d.ci.up <- c()
  z.eff.lo <- z.eff.up <- z.ci.lo <- z.ci.up <- c()
  
  if("control" %in% tgroup){
      d.eff.lo <- cbind(eff.lo, x$d0.lb)
      d.eff.up <- cbind(eff.up, x$d0.ub)
      d.ci.lo <- cbind(ci.lo, x$d0.ci[1,])
      d.ci.up <- cbind(ci.up, x$d0.ci[2,])

      z.eff.lo <- cbind(eff.lo, x$z0.lb)
      z.eff.up <- cbind(eff.up, x$z0.ub)
      z.ci.lo <- cbind(ci.lo, x$z0.ci[1,])
      z.ci.up <- cbind(ci.up, x$z0.ci[2,])
  }
  if("treated" %in% tgroup){
      d.eff.lo <- cbind(d.eff.lo, x$d1.lb)
      d.eff.up <- cbind(d.eff.up, x$d1.ub)
      d.ci.lo <- cbind(d.ci.lo, x$d1.ci[1,])
      d.ci.up <- cbind(d.ci.up, x$d1.ci[2,])

      z.eff.lo <- cbind(z.eff.lo, x$z1.lb)
      z.eff.up <- cbind(z.eff.up, x$z1.ub)
      z.ci.lo <- cbind(z.ci.lo, x$z1.ci[1,])
      z.ci.up <- cbind(z.ci.up, x$z1.ci[2,])
  }
  if("average" %in% tgroup){
      d.eff.lo <- cbind(d.eff.lo, x$d.ave.lb)
      d.eff.up <- cbind(d.eff.up, x$d.ave.ub)
      d.ci.lo <- cbind(d.ci.lo, x$d.ave.ci[1,])
      d.ci.up <- cbind(d.ci.up, x$d.ave.ci[2,])

      z.eff.lo <- cbind(z.eff.lo, x$z.ave.lb)
      z.eff.up <- cbind(z.eff.up, x$z.ave.ub)
      z.ci.lo <- cbind(z.ci.lo, x$z.ave.ci[1,])
      z.ci.up <- cbind(z.ci.up, x$z.ave.ci[2,])
  }

  ## 1. Point Estimate under Homogeneous Interaction Assumption
  if(show.point){
    if(is.null(main)){
      ma <- "Point Estimate"
    } else ma <- main
    if(is.null(xlab)){
      xla <- "Average Effects"
    } else xla <- xlab
    if(is.null(ylab)){
        yla <- expression(paste("Total (", bar(tau), ")"))
        if("control" %in% tgroup){
            yla <- c(yla, expression(paste("Control (", bar(zeta)[0], ")")))
        }
        if("treated" %in% tgroup){
            yla <- c(yla, expression(paste("Treated (", bar(zeta)[1], ")")))
        }
        if("average" %in% tgroup){
            yla <- c(yla, expression(paste("Average (", bar(bar(zeta)), ")")))
        }
        if("control" %in% tgroup){
            yla <- c(yla, expression(paste("Control (", bar(delta)[0], ")")))
        }
        if("treated" %in% tgroup){
            yla <- c(yla, expression(paste("Treated (", bar(delta)[1], ")")))
        }
        if("average" %in% tgroup){
            yla <- c(yla, expression(paste("Average (", bar(bar(delta)), ")")))
        }
        if(IND && DIR && !TOT) {
            yla <- yla[2:7]
        }
        if(!IND && !TOT) {
            yla <- yla[2:4]
        }
        if(!DIR && !TOT) {
            yla <- yla[5:7]
        }
        if(!IND && TOT) {
            yla <- yla[c(1, 2:4)]
        }
        if(!DIR && TOT) {
            yla <- yla[c(1, 5:7)]
        }
    } else yla <- ylab

    if(IND && DIR && TOT) {
        eff <- c(x$tau, z.eff.lo[1,], d.eff.lo[1,])
        ci <- cbind(x$tau.ci, rbind(z.ci.lo[1,], z.ci.up[1,]), rbind(d.ci.lo[1,], d.ci.up[1,]))
    } else if(IND && DIR && !TOT) {
        eff <- c(z.eff.lo[1,], d.eff.lo[1,])
        ci <- cbind(rbind(z.ci.lo[1,], z.ci.up[1,]), rbind(d.ci.lo[1,], d.ci.up[1,]))
    } else if(!IND && !TOT) {
        eff <- c(z.eff.lo[1,])
        ci <- rbind(z.ci.lo[1,], z.ci.up[1,])
    } else if(!DIR && !TOT) {
        eff <- c(d.eff.lo[1,])
        ci <- rbind(d.ci.lo[1,], d.ci.up[1,])
    } else if(!IND && TOT) {
        eff <- c(x$tau, z.eff.lo[1,])
        ci <- cbind(x$tau.ci, rbind(z.ci.lo[1,], z.ci.up[1,]))
    } else if(!DIR && TOT) {
        eff <- c(x$tau, d.eff.lo[1,])
        ci <- cbind(x$tau.ci, rbind(d.ci.lo[1,], d.ci.up[1,]))
    }

    if(is.null(xlim)){
      xli <- c(min(ci), max(ci))
    } else xli <- xlim
    if(is.null(ylim)){
      yli <- c(0, length(eff)) + 0.5
    } else yli <- ylim

    plot(0, 0, type = "n", main = ma, xlab = xla, ylab = "",
         xlim = xli, ylim = yli, yaxt = "n", ...)
    for(i in 1:length(eff)){
      segments(ci[1,i], i, ci[2,i], i, lwd = lwd, col = col.cbar)
      points(eff[i], i, pch = pch, cex = cex, col = col.eff)
    }
    abline(v = 0)
    text(y = 1:length(eff), labels = yla, srt = 45, par("usr")[1], xpd = TRUE, pos = 2)
  }

  ## 2. Sensitivity analysis

  if(is.null(ylab)){
    yla <- as.list(rep(NA, 6))
    if("control" %in% tgroup){
      yla[[1]] <- c(expression(paste(bar(delta)[0], "(", sigma, ")")),
                    expression(paste(bar(delta)[0], "(", R^{2}, "*)")),
                    expression(paste(bar(delta)[0], "(", tilde(R)^{2}, ")")))

      yla[[4]] <- c(expression(paste(bar(zeta)[0], "(", sigma, ")")),
                    expression(paste(bar(zeta)[0], "(", R^{2}, "*)")),
                    expression(paste(bar(zeta)[0], "(", tilde(R)^{2}, ")")))
    }
    if("treated" %in% tgroup){
      yla[[2]] <- c(expression(paste(bar(delta)[1], "(", sigma, ")")),
                    expression(paste(bar(delta)[1], "(", R^{2}, "*)")),
                    expression(paste(bar(delta)[1], "(", tilde(R)^{2}, ")")))
      yla[[5]] <- c(expression(paste(bar(zeta)[1], "(", sigma, ")")),
                    expression(paste(bar(zeta)[1], "(", R^{2}, "*)")),
                    expression(paste(bar(zeta)[1], "(", tilde(R)^{2}, ")")))
    }
    if("average" %in% tgroup){
      yla[[3]] <- c(expression(paste(bar(bar(delta)), "(", sigma, ")")),
                    expression(paste(bar(bar(delta)), "(", R^{2}, "*)")),
                    expression(paste(bar(bar(delta)), "(", tilde(R)^{2}, ")")))
      yla[[6]] <- c(expression(paste(bar(bar(zeta)), "(", sigma, ")")),
                    expression(paste(bar(bar(zeta)), "(", R^{2}, "*)")),
                    expression(paste(bar(bar(zeta)), "(", tilde(R)^{2}, ")")))
    }
  } else yla <- ylab

  wh <- c("sigma", "R2-residual", "R2-total") %in% type

  for(j in 1:length(wh)){
    if(!wh[j]) next
    if(is.null(main)){
      ma <- ifelse(j == 1, "Sensitivity with Respect to \n Interaction Heterogeneity",
                   "Sensitivity with Respect to \n Importance of Interaction")
    } else ma <- main
    if(is.null(xlab)){
      xla <- switch(j, expression(sigma),
                    expression(paste(R^{2}, "*")),
                    expression(tilde(R)^{2}))
    } else xla <- xlab
    if(is.null(xlim)){
      if(j == 1) xli <- range(x$sigma) else xli <- c(0,1)
    } else xli <- xlim
    if(is.null(ylim)){
        if(IND && DIR){
            yli <- c(min(c(d.ci.lo, z.ci.lo)), max(c(d.ci.up, z.ci.up)))
            slides <- NULL
            if("control" %in% tgroup){
                slides <- c(1, 4)
            }
            if("treated" %in% tgroup){
                slides <- c(slides, 2, 5)
            }
            if("average" %in% tgroup){
                slides <- c(slides, 3, 6)
            }
            slides <- sort(slides)  
        } else if (!IND && DIR){
            yli <- c(min(z.ci.lo), max(z.ci.up))
            slides <- NULL
            if("control" %in% tgroup){
                slides <- c(4)
            }
            if("treated" %in% tgroup){
                slides <- c(slides, 5)
            }
            if("average" %in% tgroup){
                slides <- c(slides, 6)
            }
            slides <- sort(slides)            
        } else if (IND && !DIR){
            yli <- c(min(d.ci.lo), max(d.ci.up))
            slides <- NULL
            if("control" %in% tgroup){
                slides <- c(1)
            }
            if("treated" %in% tgroup){
                slides <- c(slides, 2)
            }
            if("average" %in% tgroup){
                slides <- c(slides, 3)
            }
            slides <- sort(slides)   
        }
    } else yli <- ylim

    spar <- switch(j, x$sigma, x$R2star, x$R2tilde)

    ci.lo <- cbind(d.ci.lo, z.ci.lo)
    ci.up <- cbind(d.ci.up, z.ci.up)
    eff.lo <- cbind(d.eff.lo, z.eff.lo)
    eff.up <- cbind(d.eff.up, z.eff.up)

    index <- 1:length(slides)
    for(i in slides){
      plot(0, 0, type = "n", main = ma, xlab = xla, ylab = yla[[i]][j],
           xlim = xli, ylim = yli, ...)
      
      ii <- index[i == slides]
      polygon(c(spar, rev(spar)), c(ci.lo[,ii], rev(ci.up[,ii])),
              border = FALSE, col = col.creg)
      lines(spar, eff.lo[,ii], lwd = lwd, col = col.eff)
      lines(spar, eff.up[,ii], lwd = lwd, col = col.eff)
      abline(h = 0)
      abline(h = eff.lo[1,ii], lty = "dashed")
    }
  }
}
