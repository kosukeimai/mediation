test.TMint.mediate <- function(x, conf.level = 0.95){
  if(is.null(x$d0.sims) || is.null(x$d1.sims) || is.null(x$z0.sims) || is.null(x$z1.sims)){
    stop("simulation draws missing; rerun mediate with 'long' set to TRUE")
  }
  if(!x$INT){
    stop("outcome model must include interaction between treatment and mediator")
  }
  d.diff <- x$d1 - x$d0
  d.diff.sims <- x$d1.sims - x$d0.sims  
  cat("\n")
  cat("Null: ACME(1) - ACME(0) = 0 \n")
  cat("Two-sided p-value:", pval(d.diff.sims, d.diff), "\n")
  cat(100*conf.level, "% confidence interval:", 
      quantile(d.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2)), "\n\n")
}


test.modmed.mediate <- function(x, covariates.1, covariates.2, 
                                sims = x$sims, conf.level = x$conf.level){
  
  cl <- getCall(x)
  
  seed <- .Random.seed
  cl$covariates <- covariates.1
  cl$long <- TRUE
  cl$sims <- sims
  out.1 <- eval(cl)
  
  .Random.seed <- seed
  cl$covariates <- covariates.2
  cl$long <- TRUE
  cl$sims <- sims
  out.2 <- eval(cl)
  
  d1.diff <- out.1$d1 - out.2$d1
  d1.diff.sims <- out.1$d1.sims - out.2$d1.sims

  z0.diff <- out.1$z0 - out.2$z0
  z0.diff.sims <- out.1$z0.sims - out.2$z0.sims
  
  if(x$INT){
  
    d0.diff <- out.1$d0 - out.2$d0
    d0.diff.sims <- out.1$d0.sims - out.2$d0.sims

    z1.diff <- out.1$z1 - out.2$z1
    z1.diff.sims <- out.1$z1.sims - out.2$z1.sims
    
  }
  
  # Print results
  
  if(x$INT){
  cat("\n")
  cat("Null: ACME(1|covariates.1) - ACME(1|covariates.2) = 0 \n")
  cat("Two-sided p-value:", pval(d1.diff.sims, d1.diff), "\n")
  cat(100*conf.level, "% confidence interval:", 
      quantile(d1.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2)), "\n\n")
  
  cat("\n")
  cat("Null: ACME(0|covariates.1) - ACME(0|covariates.2) = 0 \n")
  cat("Two-sided p-value:", pval(d0.diff.sims, d0.diff), "\n")
  cat(100*conf.level, "% confidence interval:", 
      quantile(d0.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2)), "\n\n")
  
  cat("\n")
  cat("Null: ADE(1|covariates.1) - ADE(1|covariates.2) = 0 \n")
  cat("Two-sided p-value:", pval(z1.diff.sims, z1.diff), "\n")
  cat(100*conf.level, "% confidence interval:", 
      quantile(z1.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2)), "\n\n")
  
  cat("\n")
  cat("Null: ADE(0|covariates.1) - ADE(0|covariates.2) = 0 \n")
  cat("Two-sided p-value:", pval(z0.diff.sims, z0.diff), "\n")
  cat(100*conf.level, "% confidence interval:", 
      quantile(z0.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2)), "\n\n")

  } else {
  
  cat("\n")
  cat("Null: ACME(covariates.1) - ACME(covariates.2) = 0 \n")
  cat("Two-sided p-value:", pval(d1.diff.sims, d1.diff), "\n")
  cat(100*conf.level, "% confidence interval:", 
      quantile(d1.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2)), "\n\n")
  
  cat("\n")
  cat("Null: ADE(covariates.1) - ADE(covariates.2) = 0 \n")
  cat("Two-sided p-value:", pval(z0.diff.sims, z0.diff), "\n")
  cat(100*conf.level, "% confidence interval:", 
      quantile(z0.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2)), "\n\n")
  
  }
  
}


