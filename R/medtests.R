test.TMint <- function(x, conf.level = 0.95){
  if(is.null(x$d0.sims) || is.null(x$d1.sims) || is.null(x$z0.sims) || is.null(x$z1.sims)){
    stop("simulation draws missing; rerun mediate with 'long' set to TRUE")
  }
  if(!x$INT){
    stop("outcome model must include interaction between treatment and mediator")
  }
  d.diff <- x$d1 - x$d0
  d.diff.sims <- x$d1.sims - x$d0.sims  
  cat("\n")
  cat("Null: Mediation Effect_1 - Mediation Effect_0 = 0 \n")
  cat("Two-sided p-value:", pval(d.diff.sims, d.diff), "\n")
  cat(100*conf.level, "% confidence interval:", 
      quantile(d.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2)), "\n\n")

}

test.modmed <- function(x, moderator.levels, conf.int){
 seed <- .Random.seed 
}

update.mediate <- function(x, covariates){
  cl <- getCall(x)
  cl$covariates <- covariates
  eval(cl)
}