mediate <- function(model.m, model.y, sims = 1000, 
                    boot = FALSE, boot.ci.type = "perc",
                    treat = "treat.name", mediator = "med.name",
                    covariates = NULL, outcome = NULL,
                    control = NULL, conf.level = .95,
                    control.value = 0, treat.value = 1,
                    long = TRUE, dropobs = FALSE,
                    robustSE = FALSE, cluster = NULL, group.out = NULL, ...){
  
  cl <- match.call()
  
  # Warn users who still use INT option
  if(match("INT", names(cl), 0L)){
    warning("'INT' is deprecated - existence of interaction terms is now automatically detected from model formulas")
  }
  
  # Warning for robustSE and cluster used with boot
  if(robustSE && boot){
    warning("'robustSE' is ignored for nonparametric bootstrap")
  }
  
  if(!is.null(cluster) && boot){
    warning("'cluster' is ignored for nonparametric bootstrap")
  }
  
  if(robustSE & !is.null(cluster)){
    stop("choose either `robustSE' or `cluster' option, not both")
  }

  if(boot.ci.type != "bca" & boot.ci.type != "perc"){
      stop("choose either `bca' or `perc' for boot.ci.type")
  }
  
  # Drop observations not common to both mediator and outcome models
  if(dropobs){
    odata.m <- model.frame(model.m)
    odata.y <- model.frame(model.y)
    if(!is.null(cluster)){
      if(is.null(row.names(cluster)) &
           (nrow(odata.m)!=length(cluster) | nrow(odata.y)!=length(cluster))
      ){
        warning("cluster IDs may not correctly match original observations due to missing data")
      }
      odata.y <- merge(odata.y, as.data.frame(cluster), sort=FALSE,
                       by="row.names")
      rownames(odata.y) <- odata.y$Row.names
      odata.y <- odata.y[,-1L]
    }
    newdata <- merge(odata.m, odata.y, sort=FALSE,
                     by=c("row.names", intersect(names(odata.m), names(odata.y))))
    rownames(newdata) <- newdata$Row.names
    newdata <- newdata[,-1L]
    rm(odata.m, odata.y)
    
    call.m <- getCall(model.m)
    call.y <- getCall(model.y)
    
    call.m$data <- call.y$data <- newdata
    if(c("(weights)") %in% names(newdata)){
      call.m$weights <- call.y$weights <- model.weights(newdata)
    }
    model.m <- eval.parent(call.m)
    model.y <- eval.parent(call.y)
    if(!is.null(cluster)){
      cluster <- factor(newdata[, ncol(newdata)])  # factor drops missing levels
    }
  }
  
  # Model type indicators
  isGam.y <- inherits(model.y, "gam")
  isGam.m <- inherits(model.m, "gam")
  isGlm.y <- inherits(model.y, "glm")  # Note gam and bayesglm also inherits "glm"
  isGlm.m <- inherits(model.m, "glm")  # Note gam and bayesglm also inherits "glm"
  isLm.y <- inherits(model.y, "lm")    # Note gam, glm and bayesglm also inherit "lm"
  isLm.m <- inherits(model.m, "lm")    # Note gam, glm and bayesglm also inherit "lm"
  isVglm.y <- inherits(model.y, "vglm")
  isRq.y <- inherits(model.y, "rq")
  isRq.m <- inherits(model.m, "rq")
  isOrdered.y <- inherits(model.y, "polr")  # Note bayespolr also inherits "polr"
  isOrdered.m <- inherits(model.m, "polr")  # Note bayespolr also inherits "polr"
  isSurvreg.y <- inherits(model.y, "survreg")
  isSurvreg.m <- inherits(model.m, "survreg")
  isMer.y <- inherits(model.y, "merMod") # Note lmer and glmer do not inherit "lm" and "glm"
  isMer.m <- inherits(model.m, "merMod") # Note lmer and glmer do not inherit "lm" and "glm"
  
  # Record family and link of model.m if glmer 
  if(isMer.m && getCall(model.m)[[1]] == "glmer"){
    m.family <- as.character(model.m@call$family)
    if(m.family[1] == "binomial" && (is.na(m.family[2]) || m.family[2] == "logit")){
      M.fun <- binomial(link = "logit")
    } else if(m.family[1] == "binomial" && m.family[2] == "probit"){
      M.fun <- binomial(link = "probit")
    } else if(m.family[1] == "binomial" && m.family[2] == "cloglog"){ 
      M.fun <- binomial(link = "cloglog")
    } else if(m.family[1] == "poisson" && (is.na(m.family[2]) || m.family[2] == "log")){
      M.fun <- poisson(link = "log")
    } else if(m.family[1] == "poisson" && m.family[2] == "identity"){
      M.fun <- poisson(link = "identity")
    } else if(m.family[1] == "poisson" && m.family[2] == "sqrt"){
      M.fun <- poisson(link = "sqrt")
    } else {
      stop("glmer family for the mediation model not supported")
    } ### gamma & inverse gaussian excluded (no function to estimate parameters of S4 glmer vs. S3 glm)
  }
  
  # Record family and link of model.y if glmer 
  if(isMer.y && getCall(model.y)[[1]] == "glmer"){
    y.family <- as.character(model.y@call$family)
    if(y.family[1] == "binomial" && (is.na(y.family[2]) || y.family[2] == "logit")){
      Y.fun <- binomial(link = "logit")
    } else if(y.family[1] == "binomial" && y.family[2] == "probit"){
      Y.fun <- binomial(link = "probit")
    } else if(y.family[1] == "binomial" && y.family[2] == "cloglog"){ 
      Y.fun <- binomial(link = "cloglog")
    } else if(y.family[1] == "poisson" && (is.na(y.family[2]) || y.family[2] == "log")){
      Y.fun <- poisson(link = "log")
    } else if(y.family[1] == "poisson" && y.family[2] == "identity"){
      Y.fun <- poisson(link = "identity")
    } else if(y.family[1] == "poisson" && y.family[2] == "sqrt"){
      Y.fun <- poisson(link = "sqrt")
    } else {
      stop("glmer family for the outcome model not supported")
    } ### gamma & inverse gaussian excluded (no function to estimate parameters of S4 glmer vs. S3 glm)
  }
  
  # Record family of model.m if glm
  if(isGlm.m){
    FamilyM <- model.m$family$family
}

  # Record family of model.m if glmer
  if(isMer.m && getCall(model.m)[[1]] == "glmer"){
    FamilyM <- M.fun$family
  }
  
  # Record vfamily of model.y if vglm (currently only tobit)
  if(isVglm.y){
    VfamilyY <- model.y@family@vfamily
  }
  
  # Warning for unused options
  if(!is.null(control) && !isGam.y){
    warning("'control' is only used for GAM outcome models - ignored")
    control <- NULL
  }
  if(!is.null(outcome) && !(isSurvreg.y && boot)){
    warning("'outcome' is only relevant for survival outcome models with bootstrap - ignored")
  }
  
  # Model frames for M and Y models
  m.data <- model.frame(model.m)  # Call.M$data
  y.data <- model.frame(model.y)  # Call.Y$data

  if(!is.null(cluster)){
      row.names(m.data) <- 1:nrow(m.data)
      row.names(y.data) <- 1:nrow(y.data)

      if(!is.null(model.m$weights)){
          m.weights <- as.data.frame(model.m$weights)
          m.name <- as.character(model.m$call$weights)  
          names(m.weights) <- m.name
          m.data <- cbind(m.data, m.weights)
      }

      if(!is.null(model.y$weights)){
          y.weights <- as.data.frame(model.y$weights)
          y.name <- as.character(model.y$call$weights)  
          names(y.weights) <- y.name
          y.data <- cbind(y.data, y.weights)
      }
  }

  # group-level mediator 
  if(isMer.y & !isMer.m){
    m.data <- eval(model.m$call$data, environment(formula(model.m)))  ### add group ID to m.data 
    m.data <- na.omit(m.data)
    y.data <- na.omit(y.data)
  }
  
  # Specify group names
  
  if(isMer.m && isMer.y){
    med.group <- names(model.m@flist)
    out.group <- names(model.y@flist)
    n.med <- length(med.group)
    n.out <- length(out.group)
    if(n.med > 1 || n.out > 1){
      stop("mediate does not support more than two levels per model")
    } else {
      group.m <- med.group
      group.y <- out.group
      if(!is.null(group.out) && !(group.out %in% c(group.m, group.y))){
        warning("group.out does not match group names used in merMod")
      } else if(is.null(group.out)){
        group.out <- group.y
      }
    }
  } else if(!isMer.m && isMer.y){
    out.group <- names(model.y@flist)
    n.out <- length(out.group)
    if(n.out > 1){
      stop("mediate does not support more than two levels per model")
    } else {
      group.m <- NULL
      group.y <- out.group
      group.out <- group.y
    }
  } else if(isMer.m && !isMer.y){
    med.group <- names(model.m@flist)
    n.med <- length(med.group)
    if(n.med > 1){
      stop("mediate does not support more than two levels per model")
    } else {
      group.m <- med.group
      group.y <- NULL
      group.out <- group.m
    }
  } else {
    group.m <- NULL
    group.y <- NULL
    group.out <- NULL
  }
  
  # group data if lmer or glmer 
  if(isMer.m){
    group.id.m <- m.data[,group.m]
  } else {
    group.id.m <- NULL
  } 
  if(isMer.y){
    group.id.y <- y.data[,group.y]	
  } else {
    group.id.y <- NULL
  }
  # group data to be output in summary and plot if lmer or glmer 
  if(isMer.y && isMer.m){
    if(group.out == group.m){
      group.id <- m.data[,group.m]
      group.name <- group.m
    } else {
      group.id <- y.data[,group.y]     
      group.name <- group.y
    }
  } else if(!isMer.y && isMer.m){
    group.id <- m.data[,group.m]   
    group.name <- group.m
  } else if(isMer.y && !isMer.m){   ### group-level mediator
    if(!(group.y %in% names(m.data))){
      stop("specify group-level variable in mediator data")
    } else {
      group.id <- y.data[,group.y]
      group.name <- group.y
      Y.ID<- sort(unique(group.id))
      M.ID <- sort(as.vector(data.matrix(m.data[group.y])))
      if(length(Y.ID) != length(M.ID)){
        stop("groups do not match between mediator and outcome models")
      } else {
        if(FALSE %in% unique(Y.ID == M.ID)){
          stop("groups do not match between mediator and outcome models")
        }
      }
    }
  } else {
    group.id <- NULL
    group.name <- NULL
  }
  
  
  # Numbers of observations and categories
  n.m <- nrow(m.data)
  n.y <- nrow(y.data)
  if(!(isMer.y & !isMer.m)){   ### n.y and n.m are different when group-level mediator is used 
    if(n.m != n.y){
      stop("number of observations do not match between mediator and outcome models")
    } else{
      n <- n.m
    }
    m <- length(sort(unique(model.frame(model.m)[,1])))
  }
  
  # Extracting weights from models
  weights.m <- model.weights(m.data)
  weights.y <- model.weights(y.data)
  
  if(!is.null(weights.m) && isGlm.m && FamilyM == "binomial"){
    message("weights taken as sampling weights, not total number of trials")
  }
  if(!is.null(weights.m) && isMer.m && getCall(model.m)[[1]] == "glmer" && FamilyM == "binomial"){
    message("weights taken as sampling weights, not total number of trials")
  }
  if(is.null(weights.m)){
    weights.m <- rep(1,nrow(m.data))
  }
  if(is.null(weights.y)){
    weights.y <- rep(1,nrow(y.data))
  }
  if(!(isMer.y & !isMer.m)){
    if(!all(weights.m == weights.y)) {
      stop("weights on outcome and mediator models not identical")
    } else {
      weights <- weights.m
    }
  } else{
    weights <- weights.y  ### group-level mediator  
  }

  # Convert character treatment to factor
  if(is.character(m.data[,treat])){
    m.data[,treat] <- factor(m.data[,treat])
  }
  if(is.character(y.data[,treat])){
    y.data[,treat] <- factor(y.data[,treat])
  }
  
  # Convert character mediator to factor
  if(is.character(y.data[,mediator])){
    y.data[,mediator] <- factor(y.data[,mediator])
  }
  
  # Factor treatment indicator
  isFactorT.m <- is.factor(m.data[,treat])
  isFactorT.y <- is.factor(y.data[,treat])
  if(isFactorT.m != isFactorT.y){
    stop("treatment variable types differ in mediator and outcome models")
  } else {
    isFactorT <- isFactorT.y
  }
  
  if(isFactorT){
    t.levels <- levels(y.data[,treat])
    if(treat.value %in% t.levels & control.value %in% t.levels){
      cat.0 <- control.value
      cat.1 <- treat.value
    } else {
      cat.0 <- t.levels[1]
      cat.1 <- t.levels[2]
      warning("treatment and control values do not match factor levels; using ", cat.0, " and ", cat.1, " as control and treatment, respectively")
    }
  } else {
    cat.0 <- control.value
    cat.1 <- treat.value
  }
  
  # Factor mediator indicator
  isFactorM <- is.factor(y.data[,mediator])
  
  if(isFactorM){
    m.levels <- levels(y.data[,mediator])
  }
  
  #####################################
  ## Define functions
  #####################################
  
  indexmax <- function(x){
    ## Return position of largest element in vector x
    order(x)[length(x)]
  }
  
  getvcov <- function(dat, fm, cluster){
    ## Compute cluster robust standard errors
    ## fm is the model object
    cluster <- factor(cluster)  # remove missing levels and NA
    M <- nlevels(cluster)
    N <- sum(!is.na(cluster))
    K <- fm$rank
    dfc <- (M/(M-1))*((N-1)/(N-K))
    uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
    dfc*sandwich(fm, meat. = crossprod(uj)/N)
  }
  
  ############################################################################
  ############################################################################
  ### CASE I: EVERYTHING EXCEPT ORDERED OUTCOME
  ############################################################################
  ############################################################################
  if (!isOrdered.y) {
    
    ########################################################################
    ## Case I-1: Quasi-Bayesian Monte Carlo
    ########################################################################
    if(!boot){
      # Error if gam outcome or quantile mediator
      if(isGam.m | isGam.y | isRq.m){
        stop("'boot' must be 'TRUE' for models used")
      }
      
      # Get mean and variance parameters for mediator simulations
      if(isSurvreg.m && is.null(survival::survreg.distributions[[model.m$dist]]$scale)){
        MModel.coef <- c(coef(model.m), log(model.m$scale))
        scalesim.m <- TRUE
      } else if(isMer.m){
        MModel.fixef <- lme4::fixef(model.m)
        MModel.ranef <- lme4::ranef(model.m)
        scalesim.m <- FALSE    	           	
      } else {
        MModel.coef <- coef(model.m)
        scalesim.m <- FALSE
      }
      
      if(isOrdered.m){
          if(is.null(model.m$Hess)){
              cat("Mediator model object does not contain 'Hessian';")
          }
          k <- length(MModel.coef)
          MModel.var.cov <- vcov(model.m)[(1:k),(1:k)]
      } else if(isSurvreg.m){
          MModel.var.cov <- vcov(model.m)
      } else {
          if(robustSE & !isMer.m){
              MModel.var.cov <- vcovHC(model.m, ...)
          } else if(robustSE & isMer.m){
              MModel.var.cov <- vcov(model.m)
              warning("robustSE does not support mer class: non-robust SEs are computed for model.m")
          } else if(!is.null(cluster)){
              if(nrow(m.data)!=length(cluster)){
                  warning("length of cluster vector differs from # of obs for mediator model")
              }
              dta <- merge(m.data, as.data.frame(cluster), sort=FALSE,
                           by="row.names")
              fm <- update(model.m, data=dta)
              MModel.var.cov <- getvcov(dta, fm, dta[,ncol(dta)])
          } else {
              MModel.var.cov <- vcov(model.m)
          }
      }
      
      # Get mean and variance parameters for outcome simulations
      if(isSurvreg.y && is.null(survival::survreg.distributions[[model.y$dist]]$scale)){
        YModel.coef <- c(coef(model.y), log(model.y$scale))
        scalesim.y <- TRUE  # indicates if survreg scale parameter is simulated
      } else if(isMer.y){
        YModel.fixef <- lme4::fixef(model.y)
        YModel.ranef <- lme4::ranef(model.y)
        scalesim.y <- FALSE    	           	
      } else {
        YModel.coef <- coef(model.y)
        scalesim.y <- FALSE
      }
      
      if(isRq.y){
          YModel.var.cov <- summary(model.y, covariance=TRUE)$cov
      } else if(isSurvreg.y){
          YModel.var.cov <- vcov(model.y)
      } else {
          if(robustSE & !isMer.y){
              YModel.var.cov <- vcovHC(model.y, ...)
          } else if(robustSE & isMer.y){
              YModel.var.cov <- vcov(model.y)
              warning("robustSE does not support mer class: non-robust SEs are computed for model.y")
          } else if(!is.null(cluster)){
              if(nrow(y.data)!=length(cluster)){
                  warning("length of cluster vector differs from # of obs for outcome model")
              }
              dta <- merge(y.data, as.data.frame(cluster), sort=FALSE,
                           by="row.names")
              fm <- update(model.y, data=dta)
              YModel.var.cov <- getvcov(dta, fm, dta[,ncol(dta)])
          } else {
              YModel.var.cov <- vcov(model.y)
          }
      }
      
      # Draw model coefficients from normal
      
      se.ranef.new <- function (object) {
          se.bygroup <- lme4::ranef(object, condVar = TRUE)
          n.groupings <- length(se.bygroup)
          for (m in 1:n.groupings) {
              vars.m <- attr(se.bygroup[[m]], "postVar")
              K <- dim(vars.m)[1]
              J <- dim(vars.m)[3]
              names.full <- dimnames(se.bygroup[[m]])
              se.bygroup[[m]] <- array(NA, c(J, K))
              for (j in 1:J) {
                  se.bygroup[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[, 
                                                                     , j])))
              }
              dimnames(se.bygroup[[m]]) <- list(names.full[[1]], names.full[[2]])
          }
          return(se.bygroup)
      }
      
      if(isMer.m){
          MModel.fixef.vcov <- as.matrix(vcov(model.m))
          MModel.fixef.sim <- rmvnorm(sims,mean=MModel.fixef,sigma=MModel.fixef.vcov)
          Nm.ranef <- ncol(lme4::ranef(model.m)[[1]]) 
          MModel.ranef.sim <- vector("list",Nm.ranef)
          for (d in 1:Nm.ranef){
              MModel.ranef.sim[[d]] <- matrix(rnorm(sims*nrow(lme4::ranef(model.m)[[1]]), mean = lme4::ranef(model.m)[[1]][,d], sd = se.ranef.new(model.m)[[1]][,d]), nrow = sims, byrow = TRUE)
          }
      } else {
          if(sum(is.na(MModel.coef)) > 0){
              stop("NA in model coefficients; rerun models with nonsingular design matrix")
          }
          MModel <- rmvnorm(sims, mean=MModel.coef, sigma=MModel.var.cov)
      }
      
      if(isMer.y){
          YModel.fixef.vcov <- as.matrix(vcov(model.y))
          YModel.fixef.sim <- rmvnorm(sims,mean=YModel.fixef,sigma=YModel.fixef.vcov)
          Ny.ranef <- ncol(lme4::ranef(model.y)[[1]]) 
          YModel.ranef.sim <- vector("list",Ny.ranef)
          for (d in 1:Ny.ranef){
              YModel.ranef.sim[[d]] <- matrix(rnorm(sims*nrow(lme4::ranef(model.y)[[1]]), mean = lme4::ranef(model.y)[[1]][,d], sd = se.ranef.new(model.y)[[1]][,d]), nrow = sims, byrow = TRUE)
          }
      } else {
          if(sum(is.na(YModel.coef)) > 0){
              stop("NA in model coefficients; rerun models with nonsingular design matrix")
          }
          YModel <- rmvnorm(sims, mean=YModel.coef, sigma=YModel.var.cov)
      } 
      
      if(robustSE && (isSurvreg.m | isSurvreg.y)){
        warning("`robustSE' ignored for survival models; fit the model with `robust' option instead\n")
      }
      if(!is.null(cluster) && (isSurvreg.m | isSurvreg.y)){
        warning("`cluster' ignored for survival models; fit the model with 'cluster()' term in the formula\n")
      }
      
      #####################################
      ##  Mediator Predictions
      #####################################
      ### number of observations are different when group-level mediator is used
      if(isMer.y & !isMer.m){
        n <- n.m
      }
      
      pred.data.t <- pred.data.c <- m.data
      
      if(isFactorT){
        pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
        pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
      } else {
        pred.data.t[,treat] <- cat.1
        pred.data.c[,treat] <- cat.0
      }

      if(!is.null(covariates)){
        for(p in 1:length(covariates)){
          vl <- names(covariates[p])
          if(is.factor(pred.data.t[,vl])){
            pred.data.t[,vl] <- pred.data.c[,vl] <- factor(covariates[[p]], levels = levels(m.data[,vl]))
          } else {
            pred.data.t[,vl] <- pred.data.c[,vl] <- covariates[[p]]
          }
        }
      }
      
      mmat.t <- model.matrix(terms(model.m), data=pred.data.t)
      mmat.c <- model.matrix(terms(model.m), data=pred.data.c)
      
      ### Case I-1-a: GLM Mediator
      if(isGlm.m){
        muM1 <- model.m$family$linkinv(tcrossprod(MModel, mmat.t))
        muM0 <- model.m$family$linkinv(tcrossprod(MModel, mmat.c))
        
        if(FamilyM == "poisson"){
          PredictM1 <- matrix(rpois(sims*n, lambda = muM1), nrow = sims)
          PredictM0 <- matrix(rpois(sims*n, lambda = muM0), nrow = sims)
        } else if (FamilyM == "Gamma") {
          shape <- gamma.shape(model.m)$alpha
          PredictM1 <- matrix(rgamma(n*sims, shape = shape,
                                     scale = muM1/shape), nrow = sims)
          PredictM0 <- matrix(rgamma(n*sims, shape = shape,
                                     scale = muM0/shape), nrow = sims)
        } else if (FamilyM == "binomial"){
          PredictM1 <- matrix(rbinom(n*sims, size = 1,
                                     prob = muM1), nrow = sims)
          PredictM0 <- matrix(rbinom(n*sims, size = 1,
                                     prob = muM0), nrow = sims)
        } else if (FamilyM == "gaussian"){
          sigma <- sqrt(summary(model.m)$dispersion)
          error <- rnorm(sims*n, mean=0, sd=sigma)
          PredictM1 <- muM1 + matrix(error, nrow=sims)
          PredictM0 <- muM0 + matrix(error, nrow=sims)
        } else if (FamilyM == "inverse.gaussian"){
          disp <- summary(model.m)$dispersion
          PredictM1 <- matrix(SuppDists::rinvGauss(n*sims, nu = muM1,
                                                   lambda = 1/disp), nrow = sims)
          PredictM0 <- matrix(SuppDists::rinvGauss(n*sims, nu = muM0,
                                                   lambda = 1/disp), nrow = sims)
        } else {
          stop("unsupported glm family")
        }
        
        ### Case I-1-b: Ordered mediator
      } else if(isOrdered.m){
        if(model.m$method=="logistic"){
          linkfn <- plogis
        } else if(model.m$method=="probit") {
          linkfn <- pnorm
        } else {
          stop("unsupported polr method; use 'logistic' or 'probit'")
        }
        
        m.cat <- sort(unique(model.frame(model.m)[,1]))
        lambda <- model.m$zeta
        
        mmat.t <- mmat.t[,-1]
        mmat.c <- mmat.c[,-1]
        
        ystar_m1 <- tcrossprod(MModel, mmat.t)
        ystar_m0 <- tcrossprod(MModel, mmat.c)
        
        PredictM1 <- matrix(,nrow=sims, ncol=n)
        PredictM0 <- matrix(,nrow=sims, ncol=n)
        
        for(i in 1:sims){
          cprobs_m1 <- matrix(NA,n,m)
          cprobs_m0 <- matrix(NA,n,m)
          probs_m1 <- matrix(NA,n,m)
          probs_m0 <- matrix(NA,n,m)
          
          for (j in 1:(m-1)) {  # loop to get category-specific probabilities
            cprobs_m1[,j] <- linkfn(lambda[j]-ystar_m1[i,])
            cprobs_m0[,j] <- linkfn(lambda[j]-ystar_m0[i,])
            # cumulative probabilities
            probs_m1[,m] <- 1-cprobs_m1[,m-1] # top category
            probs_m0[,m] <- 1-cprobs_m0[,m-1] # top category
            probs_m1[,1] <- cprobs_m1[,1]     # bottom category
            probs_m0[,1] <- cprobs_m0[,1]     # bottom category
          }
          
          for (j in 2:(m-1)){  # middle categories
            probs_m1[,j] <- cprobs_m1[,j]-cprobs_m1[,j-1]
            probs_m0[,j] <- cprobs_m0[,j]-cprobs_m0[,j-1]
          }
          
          draws_m1 <- matrix(NA, n, m)
          draws_m0 <- matrix(NA, n, m)
          
          for(ii in 1:n){
            draws_m1[ii,] <- t(rmultinom(1, 1, prob = probs_m1[ii,]))
            draws_m0[ii,] <- t(rmultinom(1, 1, prob = probs_m0[ii,]))
          }
          
          PredictM1[i,] <- apply(draws_m1, 1, indexmax)
          PredictM0[i,] <- apply(draws_m0, 1, indexmax)
        }
        
        ### Case I-1-c: Linear
      } else if(isLm.m){
        sigma <- summary(model.m)$sigma
        error <- rnorm(sims*n, mean=0, sd=sigma)
        muM1 <- tcrossprod(MModel, mmat.t)
        muM0 <- tcrossprod(MModel, mmat.c)
        PredictM1 <- muM1 + matrix(error, nrow=sims)
        PredictM0 <- muM0 + matrix(error, nrow=sims)
        rm(error)
        
        ### Case I-1-d: Survreg
      } else if(isSurvreg.m){
        dd <- survival::survreg.distributions[[model.m$dist]]
        if (is.null(dd$itrans)){
          itrans <- function(x) x
        } else {
          itrans <- dd$itrans
        }
        dname <- dd$dist
        if(is.null(dname)){
          dname <- model.m$dist
        }
        if(scalesim.m){
          scale <- exp(MModel[,ncol(MModel)])
          lpM1 <- tcrossprod(MModel[,1:(ncol(MModel)-1)], mmat.t)
          lpM0 <- tcrossprod(MModel[,1:(ncol(MModel)-1)], mmat.c)
        } else {
          scale <- dd$scale
          lpM1 <- tcrossprod(MModel, mmat.t)
          lpM0 <- tcrossprod(MModel, mmat.c)
        }
        error <- switch(dname,
                        extreme = log(rweibull(sims*n, shape=1, scale=1)),
                        gaussian = rnorm(sims*n),
                        logistic = rlogis(sims*n),
                        t = rt(sims*n, df=dd$parms))
        PredictM1 <- itrans(lpM1 + scale * matrix(error, nrow=sims))
        PredictM0 <- itrans(lpM0 + scale * matrix(error, nrow=sims))
        rm(error)
        
        ### Case I-1-e: Linear Mixed Effect			
      } else if(isMer.m && getCall(model.m)[[1]]=="lmer"){
        M.RANEF1 <- M.RANEF0 <- 0
        for (d in 1:Nm.ranef){
          name <- colnames(lme4::ranef(model.m)[[1]])[d]
          if(name == "(Intercept)"){
            var1 <- var0 <- matrix(1,sims,n) ### RE intercept
          } else if(name == treat){  ### RE slope of treat
            var1 <- matrix(1,sims,n) ### T = 1
            var0 <- matrix(0,sims,n) ### T = 0
          }
          else {
            var1 <- var0 <- matrix(data.matrix(m.data[name]),sims,n,byrow=T) ### RE slope of other variables
          }    
          M.ranef<-matrix(NA,sims,n)
          MModel.ranef.sim.d <- MModel.ranef.sim[[d]]
          Z <- data.frame(MModel.ranef.sim.d)
          if(is.factor(group.id.m)){
            colnames(Z) <- levels(group.id.m)
            for (i in 1:n){
              M.ranef[,i]<-Z[,group.id.m[i]==levels(group.id.m)] 
            }
          } else {
            colnames(Z) <- sort(unique(group.id.m)) 
            for (i in 1:n){
              M.ranef[,i]<-Z[,group.id.m[i]==sort(unique(group.id.m))]
            }
          }
          M.RANEF1 <- M.ranef*var1 + M.RANEF1   # sum of (random effects*corresponding covarites)  
          M.RANEF0 <- M.ranef*var0 + M.RANEF0
        }
        sigma <- attr(lme4::VarCorr(model.m), "sc")        
        error <- rnorm(sims*n, mean=0, sd=sigma)
        muM1 <- tcrossprod(MModel.fixef.sim, mmat.t) + M.RANEF1
        muM0 <- tcrossprod(MModel.fixef.sim, mmat.c) + M.RANEF0
        PredictM1 <- muM1 + matrix(error, nrow=sims)
        PredictM0 <- muM0 + matrix(error, nrow=sims)
        rm(error)          	
        
        ### Case I-1-f: Generalized Linear Mixed Effect                  	
      } else if(isMer.m && getCall(model.m)[[1]]=="glmer"){
        M.RANEF1 <-M.RANEF0 <- 0 ### 1=RE for M(1); 0=RE for M(0)
        for (d in 1:Nm.ranef){
          name <- colnames(lme4::ranef(model.m)[[1]])[d]
          if(name == "(Intercept)"){
            var1 <- var0 <- matrix(1,sims,n) ### RE intercept
          } else if(name == treat){ ### RE slope of treat
            var1 <- matrix(1,sims,n) ### T = 1
            var0 <- matrix(0,sims,n) ### T = 0
          } else {
            var1 <- var0 <- matrix(data.matrix(m.data[name]),sims,n,byrow=T) ### RE slope of other variables
          }    
          M.ranef<-matrix(NA,sims,n)
          MModel.ranef.sim.d <- MModel.ranef.sim[[d]]
          Z <- data.frame(MModel.ranef.sim.d)
          if(is.factor(group.id.m)){
            colnames(Z) <- levels(group.id.m)
            for (i in 1:n){
              M.ranef[,i]<-Z[,group.id.m[i]==levels(group.id.m)] 
            }
          } else {
            colnames(Z) <- sort(unique(group.id.m)) 
            for (i in 1:n){
              M.ranef[,i]<-Z[,group.id.m[i]==sort(unique(group.id.m))]
            }
          }
          M.RANEF1 <- M.ranef*var1 + M.RANEF1   # sum of (random effects*corresponding covarites)  
          M.RANEF0 <- M.ranef*var0 + M.RANEF0
        }       	
        muM1 <- M.fun$linkinv(tcrossprod(MModel.fixef.sim,mmat.t) + M.RANEF1)
        muM0 <- M.fun$linkinv(tcrossprod(MModel.fixef.sim,mmat.c) + M.RANEF0)
        FamilyM <- M.fun$family
        if(FamilyM == "poisson"){
          PredictM1 <- matrix(rpois(sims*n, lambda = muM1), nrow = sims)
          PredictM0 <- matrix(rpois(sims*n, lambda = muM0), nrow = sims)
        } else if (FamilyM == "binomial"){
          PredictM1 <- matrix(rbinom(n*sims, size = 1,
                                     prob = muM1), nrow = sims)
          PredictM0 <- matrix(rbinom(n*sims, size = 1,
                                     prob = muM0), nrow = sims)
        } 
      } else {
        stop("mediator model is not yet implemented")
      }

      ### group-level mediator : J -> NJ
      if(isMer.y & !isMer.m){
        J <- nrow(m.data)
        group.id.m <- as.vector(data.matrix(m.data[group.y]))
        v1 <- v0 <- matrix(NA, sims, length(group.id.y))
        num.m <- 1:J
        num.y <- 1:length(group.id.y)
        for (j in 1:J){
          id.y <- unique(group.id.y)[j]
          NUM.M <- num.m[group.id.m == id.y]
          NUM.Y <- num.y[group.id.y == id.y]
          v1[, NUM.Y] <- PredictM1[, NUM.M]
          v0[, NUM.Y] <- PredictM0[, NUM.M]
        }
        PredictM1 <- v1
        PredictM0 <- v0
      }
      
      rm(mmat.t, mmat.c)
      
      #####################################
      ##  Outcome Predictions
      #####################################
      ### number of observations are different when group-level mediator is used
      if(isMer.y & !isMer.m){
        n <- n.y
      }

      effects.tmp <- array(NA, dim = c(n, sims, 4))

      if(isMer.y){
        Y.RANEF1 <- Y.RANEF2 <- Y.RANEF3 <- Y.RANEF4 <- 0
        ### 1=RE for Y(1,M(1)); 2=RE for Y(1,M(0)); 3=RE for Y(0,M(1)); 4=RE for Y(0,M(0))
        for (d in 1:Ny.ranef){
          name <- colnames(lme4::ranef(model.y)[[1]])[d]
          if(name == "(Intercept)"){
            var1 <- var2 <- var3 <- var4 <- matrix(1,sims,n)
          } else if(name == treat){
            var1 <- matrix(1,sims,n)
            var2 <- matrix(1,sims,n)
            var3 <- matrix(0,sims,n)
            var4 <- matrix(0,sims,n)
          } else if(name == mediator){
            var1 <- PredictM1
            var2 <- PredictM0
            var3 <- PredictM1
            var4 <- PredictM0
          } else {
            var1 <- var2 <- var3 <- var4 <- matrix(data.matrix(y.data[name]),sims,n,byrow=T)
          }    
          Y.ranef<-matrix(NA,sims,n)
          YModel.ranef.sim.d <- YModel.ranef.sim[[d]]
          Z <- data.frame(YModel.ranef.sim.d)
          if(is.factor(group.id.y)){
            colnames(Z) <- levels(group.id.y)
            for (i in 1:n){
              Y.ranef[,i]<-Z[,group.id.y[i]==levels(group.id.y)] 
            }
          } else {
            colnames(Z) <- sort(unique(group.id.y)) 
            for (i in 1:n){
              Y.ranef[,i]<-Z[,group.id.y[i]==sort(unique(group.id.y))]
            }
          }
          Y.RANEF1 <- Y.ranef*var1 + Y.RANEF1   # sum of (random effects*corresponding covarites)
          Y.RANEF2 <- Y.ranef*var2 + Y.RANEF2
          Y.RANEF3 <- Y.ranef*var3 + Y.RANEF3
          Y.RANEF4 <- Y.ranef*var4 + Y.RANEF4
        }       	
      }

      for(e in 1:4){
        tt <- switch(e, c(1,1,1,0), c(0,0,1,0), c(1,0,1,1), c(1,0,0,0))
        Pr1 <- matrix(, nrow=n, ncol=sims)
        Pr0 <- matrix(, nrow=n, ncol=sims)
        
        for(j in 1:sims){
          pred.data.t <- pred.data.c <- y.data
          
          if(!is.null(covariates)){
            for(p in 1:length(covariates)){
              vl <- names(covariates[p])
              if(is.factor(pred.data.t[,vl])){
                pred.data.t[,vl] <- pred.data.c[,vl] <- factor(covariates[[p]], levels = levels(y.data[,vl]))
              } else {
                pred.data.t[,vl] <- pred.data.c[,vl] <- covariates[[p]]
              }
            }
          }
          
          # Set treatment values
          cat.t <- ifelse(tt[1], cat.1, cat.0)
          cat.c <- ifelse(tt[2], cat.1, cat.0)
          cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
          cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
          if(isFactorT){
            pred.data.t[,treat] <- factor(cat.t, levels = t.levels)
            pred.data.c[,treat] <- factor(cat.c, levels = t.levels)
            if(!is.null(control)){
              pred.data.t[,control] <- factor(cat.t.ctrl, levels = t.levels)
              pred.data.c[,control] <- factor(cat.c.ctrl, levels = t.levels)
            }
          } else {
            pred.data.t[,treat] <- cat.t
            pred.data.c[,treat] <- cat.c
            if(!is.null(control)){
              pred.data.t[,control] <- cat.t.ctrl
              pred.data.c[,control] <- cat.c.ctrl
            }
          }
          
          # Set mediator values
          PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
          PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
          if(isFactorM) {
            pred.data.t[,mediator] <- factor(PredictMt, levels=1:m, labels=m.levels)
            pred.data.c[,mediator] <- factor(PredictMc, levels=1:m, labels=m.levels)
          } else {
            pred.data.t[,mediator] <- PredictMt
            pred.data.c[,mediator] <- PredictMc
          }
          
          ymat.t <- model.matrix(terms(model.y), data=pred.data.t)
          ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
          
          if(isVglm.y){
            if(VfamilyY=="tobit") {
              Pr1.tmp <- ymat.t %*% YModel[j,-2]
              Pr0.tmp <- ymat.c %*% YModel[j,-2]
              Pr1[,j] <- pmin(pmax(Pr1.tmp, model.y@misc$Lower), model.y@misc$Upper)
              Pr0[,j] <- pmin(pmax(Pr0.tmp, model.y@misc$Lower), model.y@misc$Upper)
            } else {
              stop("outcome model is in unsupported vglm family")
            }
          } else if(scalesim.y){
            Pr1[,j] <- t(as.matrix(YModel[j,1:(ncol(YModel)-1)])) %*% t(ymat.t)
            Pr0[,j] <- t(as.matrix(YModel[j,1:(ncol(YModel)-1)])) %*% t(ymat.c)
          } else if(isMer.y){ 
            if(e == 1){             ### mediation(1)
              Y.RANEF.A <- Y.RANEF1
              Y.RANEF.B <- Y.RANEF2
            } else if(e == 2){      ### mediation(0)
              Y.RANEF.A <- Y.RANEF3
              Y.RANEF.B <- Y.RANEF4
            } else if(e == 3){      ### direct(1)
              Y.RANEF.A <- Y.RANEF1
              Y.RANEF.B <- Y.RANEF3
            } else {                ### direct(0)
              Y.RANEF.A <- Y.RANEF2
              Y.RANEF.B <- Y.RANEF4
            }
            Pr1[,j] <- t(as.matrix(YModel.fixef.sim[j,])) %*% t(ymat.t) + Y.RANEF.A[j,]
            Pr0[,j] <- t(as.matrix(YModel.fixef.sim[j,])) %*% t(ymat.c) + Y.RANEF.B[j,]                  
          } else {
            Pr1[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.t)
            Pr0[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.c)
          }
          
          rm(ymat.t, ymat.c, pred.data.t, pred.data.c)
        }
        
        if(isGlm.y){
          Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
          Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
        } else if(isSurvreg.y){
          dd <- survival::survreg.distributions[[model.y$dist]]
          if (is.null(dd$itrans)){
            itrans <- function(x) x
          } else {
            itrans <- dd$itrans
          }
          Pr1 <- apply(Pr1, 2, itrans)
          Pr0 <- apply(Pr0, 2, itrans)
        } else if(isMer.y && getCall(model.y)[[1]] == "glmer"){
          Pr1 <- apply(Pr1, 2, Y.fun$linkinv)
          Pr0 <- apply(Pr0, 2, Y.fun$linkinv)                	
        }
        
        effects.tmp[,,e] <- Pr1 - Pr0 ### e=1:mediation(1); e=2:mediation(0); e=3:direct(1); e=4:direct(0)
        rm(Pr1, Pr0)
      }
      
      if(!isMer.m && !isMer.y){
        rm(PredictM1, PredictM0, YModel, MModel)
      } else if(!isMer.m && isMer.y){
        rm(PredictM1, PredictM0, YModel.ranef.sim)
      } else {
        rm(PredictM1, PredictM0, MModel.ranef.sim)
      }
      
      et1<-effects.tmp[,,1] ### mediation effect (1)
      et2<-effects.tmp[,,2] ### mediation effect (0)
      et3<-effects.tmp[,,3] ### direct effect (1)
      et4<-effects.tmp[,,4] ### direct effect (0)
      
      delta.1 <- t(as.matrix(apply(et1, 2, weighted.mean, w=weights)))
      delta.0 <- t(as.matrix(apply(et2, 2, weighted.mean, w=weights)))
      zeta.1 <- t(as.matrix(apply(et3, 2, weighted.mean, w=weights)))
      zeta.0 <- t(as.matrix(apply(et4, 2, weighted.mean, w=weights)))
      rm(effects.tmp)
      
      tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
      nu.0 <- delta.0/tau
      nu.1 <- delta.1/tau
      delta.avg <- (delta.1 + delta.0)/2
      zeta.avg <- (zeta.1 + zeta.0)/2
      nu.avg <- (nu.1 + nu.0)/2
      
      d0 <- mean(delta.0)			# mediation effect
      d1 <- mean(delta.1)
      z1 <- mean(zeta.1)			# direct effect
      z0 <- mean(zeta.0)
      tau.coef <- mean(tau)	  	        # total effect
      n0 <- median(nu.0)
      n1 <- median(nu.1)
      d.avg <- (d0 + d1)/2
      z.avg <- (z0 + z1)/2
      n.avg <- (n0 + n1)/2
      
      if(isMer.y | isMer.m){
        if(!is.null(group.m) && group.name == group.m){
          G<-length(unique(group.id.m))
          delta.1.group<-matrix(NA,G,sims)
          delta.0.group<-matrix(NA,G,sims)
          zeta.1.group<-matrix(NA,G,sims)
          zeta.0.group<-matrix(NA,G,sims)
          for (g in 1:G){0
                         delta.1.group[g,] <- t(apply(matrix(et1[group.id.m==unique(group.id.m)[g],], ncol=sims), 2, weighted.mean, w=weights[group.id.m==unique(group.id.m)[g]]))
                         delta.0.group[g,] <- t(apply(matrix(et2[group.id.m==unique(group.id.m)[g],], ncol=sims), 2, weighted.mean, w=weights[group.id.m==unique(group.id.m)[g]]))
                         zeta.1.group[g,] <- t(apply(matrix(et3[group.id.m==unique(group.id.m)[g],], ncol=sims), 2, weighted.mean, w=weights[group.id.m==unique(group.id.m)[g]]))
                         zeta.0.group[g,] <- t(apply(matrix(et4[group.id.m==unique(group.id.m)[g],], ncol=sims), 2, weighted.mean, w=weights[group.id.m==unique(group.id.m)[g]]))
          }
        } else {
          G<-length(unique(group.id.y))
          delta.1.group<-matrix(NA,G,sims)
          delta.0.group<-matrix(NA,G,sims)
          zeta.1.group<-matrix(NA,G,sims)
          zeta.0.group<-matrix(NA,G,sims)
          for (g in 1:G){
            delta.1.group[g,] <- t(apply(matrix(et1[group.id.y==unique(group.id.y)[g],], ncol=sims), 2, weighted.mean, w=weights[group.id.y==unique(group.id.y)[g]]))
            delta.0.group[g,] <- t(apply(matrix(et2[group.id.y==unique(group.id.y)[g],], ncol=sims), 2, weighted.mean, w=weights[group.id.y==unique(group.id.y)[g]]))
            zeta.1.group[g,] <- t(apply(matrix(et3[group.id.y==unique(group.id.y)[g],], ncol=sims), 2, weighted.mean, w=weights[group.id.y==unique(group.id.y)[g]]))
            zeta.0.group[g,] <- t(apply(matrix(et4[group.id.y==unique(group.id.y)[g],], ncol=sims), 2, weighted.mean, w=weights[group.id.y==unique(group.id.y)[g]]))
          }
        } 
        tau.group <- (zeta.1.group + delta.0.group + zeta.0.group + delta.1.group)/2
        nu.0.group <- delta.0.group/tau.group
        nu.1.group <- delta.1.group/tau.group
        delta.avg.group <- (delta.1.group + delta.0.group)/2
        zeta.avg.group <- (zeta.1.group + zeta.0.group)/2
        nu.avg.group <- (nu.1.group + nu.0.group)/2
        
        d0.group <- apply(delta.0.group,1,mean)			# mediation effect
        d1.group <- apply(delta.1.group,1,mean)
        z1.group <- apply(zeta.1.group,1,mean)			# direct effect
        z0.group <- apply(zeta.0.group,1,mean)
        tau.coef.group <- apply(tau.group,1,mean)		# total effect
        n0.group <- apply(nu.0.group,1,median)
        n1.group <- apply(nu.1.group,1,median)
        d.avg.group <- (d0.group + d1.group)/2
        z.avg.group <- (z0.group + z1.group)/2
        n.avg.group <- (n0.group + n1.group)/2
      } 
      ########################################################################
      ## Case I-2: Nonparametric Bootstrap
      ########################################################################
    } else {
      # Error if lmer or glmer 
      if(isMer.m | isMer.y){
        stop("'boot' must be 'FALSE' for models used")
      }
      
      Call.M <- getCall(model.m)
      Call.Y <- getCall(model.y)
      
      # Storage
      delta.1 <- matrix(NA, sims, 1)
      delta.0 <- matrix(NA, sims, 1)
      zeta.1 <- matrix(NA, sims, 1)
      zeta.0 <- matrix(NA, sims, 1)
      tau <- matrix(NA, sims, 1)
      
      # Bootstrap loop begins
      for(b in 1:(sims+1)){
        index <- sample(1:n, n, replace = TRUE)
        
        if(b == sims+1){  # in the final run, use the original data
          index <- 1:n
        }
        
        if(isSurvreg.m){
          if(ncol(model.m$y) > 2){
            stop("unsupported censoring type")
          }
          mname <- names(m.data)[1]
          if(substr(mname, 1, 4) != "Surv"){
            stop("refit the survival model with `Surv' used directly in model formula")
          }
          nc <- nchar(mediator)
          eventname <- substr(mname, 5 + nc + 3, nchar(mname) - 1)
          if(nchar(eventname) == 0){
            m.data.tmp <- data.frame(m.data,
                                     as.numeric(m.data[,1L][,1L]))
            names(m.data.tmp)[c(1L, ncol(m.data)+1)] <- c(mname, mediator)
          } else {
            m.data.tmp <- data.frame(m.data,
                                     as.numeric(m.data[,1L][,1L]),
                                     as.numeric(model.m$y[,2]))
            names(m.data.tmp)[c(1L, ncol(m.data)+(1:2))] <- c(mname, mediator, eventname)
          }
          Call.M$data <- m.data.tmp[index,]
        } else {
          Call.M$data <- m.data[index,]
        }
        
        if(isSurvreg.y){
          if(ncol(model.y$y) > 2){
            stop("unsupported censoring type")
          }
          yname <- names(y.data)[1]
          if(substr(yname, 1, 4) != "Surv"){
            stop("refit the survival model with `Surv' used directly in model formula")
          }
          if(is.null(outcome)){
            stop("`outcome' must be supplied for survreg outcome with boot")
          }
          nc <- nchar(outcome)
          eventname <- substr(yname, 5 + nc + 3, nchar(yname) - 1)
          if(nchar(eventname) == 0){
            y.data.tmp <- data.frame(y.data,
                                     as.numeric(y.data[,1L][,1L]))
            names(y.data.tmp)[c(1L, ncol(y.data)+1)] <- c(yname, outcome)
          } else {
            y.data.tmp <- data.frame(y.data,
                                     as.numeric(y.data[,1L][,1L]),
                                     as.numeric(model.y$y[,2]))
            names(y.data.tmp)[c(1L, ncol(y.data)+(1:2))] <- c(yname, outcome, eventname)
          }
          Call.Y$data <- y.data.tmp[index,]
        } else {
          Call.Y$data <- y.data[index,]
        }
        
        Call.M$weights <- m.data[index,"(weights)"]
        Call.Y$weights  <- y.data[index,"(weights)"]
        
        if(isOrdered.m && length(unique(y.data[index,mediator]))!=m){
          stop("insufficient variation on mediator")
        }
        
        # Refit Models with Resampled Data
        new.fit.M <- eval.parent(Call.M)
        new.fit.Y <- eval.parent(Call.Y)
        
        #####################################
        #  Mediator Predictions
        #####################################
        pred.data.t <- pred.data.c <- m.data
        
        if(isFactorT){
          pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
          pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
        } else {
          pred.data.t[,treat] <- cat.1
          pred.data.c[,treat] <- cat.0
        }
        
        if(!is.null(covariates)){
          for(p in 1:length(covariates)){
            vl <- names(covariates[p])
            if(is.factor(pred.data.t[,vl])){
              pred.data.t[,vl] <- pred.data.c[,vl] <- factor(covariates[[p]], levels = levels(m.data[,vl]))
            } else {
              pred.data.t[,vl] <- pred.data.c[,vl] <- covariates[[p]]
            }
          }
        }
        
        ### Case I-2-a: GLM Mediator (including GAMs)
        if(isGlm.m){
          
          muM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
          muM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)
          
          if(FamilyM == "poisson"){
            PredictM1 <- rpois(n, lambda = muM1)
            PredictM0 <- rpois(n, lambda = muM0)
          } else if (FamilyM == "Gamma") {
            shape <- gamma.shape(new.fit.M)$alpha
            PredictM1 <- rgamma(n, shape = shape, scale = muM1/shape)
            PredictM0 <- rgamma(n, shape = shape, scale = muM0/shape)
          } else if (FamilyM == "binomial"){
            PredictM1 <- rbinom(n, size = 1, prob = muM1)
            PredictM0 <- rbinom(n, size = 1, prob = muM0)
          } else if (FamilyM == "gaussian"){
            sigma <- sqrt(summary(new.fit.M)$dispersion)
            error <- rnorm(n, mean=0, sd=sigma)
            PredictM1 <- muM1 + error
            PredictM0 <- muM0 + error
          } else if (FamilyM == "inverse.gaussian"){
            disp <- summary(new.fit.M)$dispersion
            PredictM1 <- SuppDists::rinvGauss(n, nu = muM1, lambda = 1/disp)
            PredictM0 <- SuppDists::rinvGauss(n, nu = muM0, lambda = 1/disp)
          } else {
            stop("unsupported glm family")
          }
          
          ### Case I-2-b: Ordered Mediator
        } else if(isOrdered.m) {
          probs_m1 <- predict(new.fit.M, newdata=pred.data.t, type="probs")
          probs_m0 <- predict(new.fit.M, newdata=pred.data.c, type="probs")
          draws_m1 <- matrix(NA, n, m)
          draws_m0 <- matrix(NA, n, m)
          for(ii in 1:n){
            draws_m1[ii,] <- t(rmultinom(1, 1, prob = probs_m1[ii,]))
            draws_m0[ii,] <- t(rmultinom(1, 1, prob = probs_m0[ii,]))
          }
          PredictM1 <- apply(draws_m1, 1, indexmax)
          PredictM0 <- apply(draws_m0, 1, indexmax)
          
          ### Case I-2-c: Quantile Regression for Mediator
        } else if(isRq.m){
          # Use inverse transform sampling to predict M
          call.new <- new.fit.M$call
          call.new$tau <- runif(n)
          newfits <- eval.parent(call.new)
          tt <- delete.response(terms(new.fit.M))
          m.t <- model.frame(tt, pred.data.t, xlev = new.fit.M$xlevels)
          m.c <- model.frame(tt, pred.data.c, xlev = new.fit.M$xlevels)
          X.t <- model.matrix(tt, m.t, contrasts = new.fit.M$contrasts)
          X.c <- model.matrix(tt, m.c, contrasts = new.fit.M$contrasts)
          rm(tt, m.t, m.c)
          PredictM1 <- rowSums(X.t * t(newfits$coefficients))
          PredictM0 <- rowSums(X.c * t(newfits$coefficients))
          rm(newfits, X.t, X.c)
          
          ### Case I-2-d: Linear
        } else if(isLm.m){
          sigma <- summary(new.fit.M)$sigma
          error <- rnorm(n, mean=0, sd=sigma)
          PredictM1 <- predict(new.fit.M, type="response",
                               newdata=pred.data.t) + error
          PredictM0 <- predict(new.fit.M, type="response",
                               newdata=pred.data.c) + error
          rm(error)
          
          ### Case I-2-e: Survreg
        } else if(isSurvreg.m){
          dd <- survival::survreg.distributions[[new.fit.M$dist]]
          if (is.null(dd$itrans)){
            itrans <- function(x) x
          } else {
            itrans <- dd$itrans
          }
          dname <- dd$dist
          if(is.null(dname)){
            dname <- new.fit.M$dist
          }
          scale <- new.fit.M$scale
          lpM1 <- predict(new.fit.M, newdata=pred.data.t, type="linear")
          lpM0 <- predict(new.fit.M, newdata=pred.data.c, type="linear")
          error <- switch(dname,
                          extreme = log(rweibull(n, shape=1, scale=1)),
                          gaussian = rnorm(n),
                          logistic = rlogis(n),
                          t = rt(n, df=dd$parms))
          PredictM1 <- as.numeric(itrans(lpM1 + scale * error))
          PredictM0 <- as.numeric(itrans(lpM0 + scale * error))
          rm(error)
          
        } else {
          stop("mediator model is not yet implemented")
        }
        
        #####################################
        #  Outcome Predictions
        #####################################
        effects.tmp <- matrix(NA, nrow = n, ncol = 4)
        for(e in 1:4){
          tt <- switch(e, c(1,1,1,0), c(0,0,1,0), c(1,0,1,1), c(1,0,0,0))
          pred.data.t <- pred.data.c <- y.data
          
          if(!is.null(covariates)){
            for(p in 1:length(covariates)){
              vl <- names(covariates[p])
              if(is.factor(pred.data.t[,vl])){
                pred.data.t[,vl] <- pred.data.c[,vl] <- factor(covariates[[p]], levels = levels(y.data[,vl]))
              } else {
                pred.data.t[,vl] <- pred.data.c[,vl] <- covariates[[p]]
              }
            }
          }
          
          # Set treatment values
          cat.t <- ifelse(tt[1], cat.1, cat.0)
          cat.c <- ifelse(tt[2], cat.1, cat.0)
          cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
          cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
          if(isFactorT){
            pred.data.t[,treat] <- factor(cat.t, levels = t.levels)
            pred.data.c[,treat] <- factor(cat.c, levels = t.levels)
            if(!is.null(control)){
              pred.data.t[,control] <- factor(cat.t.ctrl, levels = t.levels)
              pred.data.c[,control] <- factor(cat.c.ctrl, levels = t.levels)
            }
          } else {
            pred.data.t[,treat] <- cat.t
            pred.data.c[,treat] <- cat.c
            if(!is.null(control)){
              pred.data.t[,control] <- cat.t.ctrl
              pred.data.c[,control] <- cat.c.ctrl
            }
          }
          
          # Set mediator values
          PredictM1.tmp <- PredictM1
          PredictM0.tmp <- PredictM0
          PredictMt <- PredictM1 * tt[3] + PredictM0 * (1 - tt[3])
          PredictMc <- PredictM1 * tt[4] + PredictM0 * (1 - tt[4])
          if(isFactorM) {
            pred.data.t[,mediator] <- factor(PredictMt, levels=1:m, labels=m.levels)
            pred.data.c[,mediator] <- factor(PredictMc, levels=1:m, labels=m.levels)
          } else {
            pred.data.t[,mediator] <- PredictMt
            pred.data.c[,mediator] <- PredictMc
          }
          
          if(isRq.y){
            pr.1 <- predict(new.fit.Y, type="response",
                            newdata=pred.data.t, interval="none")
            pr.0 <- predict(new.fit.Y, type="response",
                            newdata=pred.data.c, interval="none")
          } else {
            pr.1 <- predict(new.fit.Y, type="response",
                            newdata=pred.data.t)
            pr.0 <- predict(new.fit.Y, type="response",
                            newdata=pred.data.c)
          }
          pr.mat <- as.matrix(cbind(pr.1, pr.0))
          effects.tmp[,e] <- pr.mat[,1] - pr.mat[,2]
          
          rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
        }
        
        # Compute all QoIs
        if(b == sims+1){
          d1 <- weighted.mean(effects.tmp[,1], weights)
          d0 <- weighted.mean(effects.tmp[,2], weights)
          z1 <- weighted.mean(effects.tmp[,3], weights)
          z0 <- weighted.mean(effects.tmp[,4], weights)
        } else {
          delta.1[b] <- weighted.mean(effects.tmp[,1], weights)
          delta.0[b] <- weighted.mean(effects.tmp[,2], weights)
          zeta.1[b] <- weighted.mean(effects.tmp[,3], weights)
          zeta.0[b] <- weighted.mean(effects.tmp[,4], weights)
        }
      }  # bootstrap loop ends
      
      tau.coef <- (d1 + d0 + z1 + z0)/2
      n0 <- d0/tau.coef
      n1 <- d1/tau.coef
      d.avg <- (d1 + d0)/2
      z.avg <- (z1 + z0)/2
      n.avg <- (n0 + n1)/2
      
      tau <- (delta.1 + delta.0 + zeta.1 + zeta.0)/2
      nu.0 <- delta.0/tau
      nu.1 <- delta.1/tau
      delta.avg <- (delta.0 + delta.1)/2
      zeta.avg <- (zeta.0 + zeta.1)/2
      nu.avg <- (nu.0 + nu.1)/2
      
    }  # nonpara boot branch ends
    
    ########################################################################
    ## Compute Outputs and Put Them Together
    ########################################################################
    
    low <- (1 - conf.level)/2
    high <- 1 - low

    if (boot & boot.ci.type == "bca"){
        BC.CI <- function(theta){
            z.inv <- length(theta[theta < mean(theta)])/sims
            z <- qnorm(z.inv)
            U <- (sims - 1) * (mean(theta) - theta)
            top <- sum(U^3)
            under <- 6 * (sum(U^2))^{3/2}
            a <- top / under
            lower.inv <-  pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
            lower2 <- lower <- quantile(theta, lower.inv)
            upper.inv <-  pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
            upper2 <- upper <- quantile(theta, upper.inv)
            return(c(lower, upper))      
        }
        d0.ci <- BC.CI(delta.0)
        d1.ci <- BC.CI(delta.1)
        tau.ci <- BC.CI(tau)
        z1.ci <- BC.CI(zeta.1)
        z0.ci <- BC.CI(zeta.0)
        n0.ci <- BC.CI(nu.0)
        n1.ci <- BC.CI(nu.1)
        d.avg.ci <- BC.CI(delta.avg)
        z.avg.ci <- BC.CI(zeta.avg)
        n.avg.ci <- BC.CI(nu.avg)
    } else {
        d0.ci <- quantile(delta.0, c(low,high), na.rm=TRUE)
        d1.ci <- quantile(delta.1, c(low,high), na.rm=TRUE)
        tau.ci <- quantile(tau, c(low,high), na.rm=TRUE)
        z1.ci <- quantile(zeta.1, c(low,high), na.rm=TRUE)
        z0.ci <- quantile(zeta.0, c(low,high), na.rm=TRUE)
        n0.ci <- quantile(nu.0, c(low,high), na.rm=TRUE)
        n1.ci <- quantile(nu.1, c(low,high), na.rm=TRUE)
        d.avg.ci <- quantile(delta.avg, c(low,high), na.rm=TRUE)
        z.avg.ci <- quantile(zeta.avg, c(low,high), na.rm=TRUE)
        n.avg.ci <- quantile(nu.avg, c(low,high), na.rm=TRUE)
    }
    
    # p-values
    d0.p <- pval(delta.0, d0)
    d1.p <- pval(delta.1, d1)
    d.avg.p <- pval(delta.avg, d.avg)
    z0.p <- pval(zeta.0, z0)
    z1.p <- pval(zeta.1, z1)
    z.avg.p <- pval(zeta.avg, z.avg)        
    n0.p <- pval(nu.0, n0)
    n1.p <- pval(nu.1, n1)
    n.avg.p <- pval(nu.avg, n.avg)
    tau.p <- pval(tau, tau.coef)
    
    if(isMer.y | isMer.m){
      QUANT<-function(object){
        z<-quantile(object,c(low,high),na.rm=TRUE)
        return(z)
      }
      d0.ci.group <- t(apply(delta.0.group,1,QUANT))
      d1.ci.group <- t(apply(delta.1.group,1,QUANT))
      tau.ci.group <- t(apply(tau.group,1,QUANT))
      z1.ci.group <- t(apply(zeta.1.group,1,QUANT))
      z0.ci.group <- t(apply(zeta.0.group,1,QUANT))
      n0.ci.group <- t(apply(nu.0.group,1,QUANT))
      n1.ci.group <- t(apply(nu.1.group,1,QUANT))
      d.avg.ci.group <- t(apply(delta.avg.group,1,QUANT))
      z.avg.ci.group <- t(apply(zeta.avg.group,1,QUANT))
      n.avg.ci.group <- t(apply(nu.avg.group,1,QUANT))
      
      d0.p.group<-rep(NA,G)
      d1.p.group<-rep(NA,G)
      d.avg.p.group<-rep(NA,G)
      z0.p.group<-rep(NA,G)
      z1.p.group<-rep(NA,G)
      z.avg.p.group<-rep(NA,G)
      n0.p.group<-rep(NA,G)
      n1.p.group<-rep(NA,G)
      n.avg.p.group<-rep(NA,G)
      tau.p.group<-rep(NA,G)
      for (g in 1:G){
        d0.p.group[g] <- pval(delta.0.group[g,], d0.group[g])
        d1.p.group[g] <- pval(delta.1.group[g,], d1.group[g])
        d.avg.p.group[g] <- pval(delta.avg.group[g,], d.avg.group[g])
        z0.p.group[g] <- pval(zeta.0.group[g,], z0.group[g])
        z1.p.group[g] <- pval(zeta.1.group[g,], z1.group[g])
        z.avg.p.group[g] <- pval(zeta.avg.group[g,], z.avg.group[g])        
        n0.p.group[g] <- pval(nu.0.group[g,], n0.group[g])
        n1.p.group[g] <- pval(nu.1.group[g,], n1.group[g])
        n.avg.p.group[g] <- pval(nu.avg.group[g,], n.avg.group[g])
        tau.p.group[g] <- pval(tau.group[g,], tau.coef.group[g])
      }
    }
    
    # Detect whether models include T-M interaction
    INT <- paste(treat,mediator,sep=":") %in% attr(terms(model.y),"term.labels") |
      paste(mediator,treat,sep=":") %in% attr(terms(model.y),"term.labels")
    if(!INT & isGam.y){
      INT <- !isTRUE(all.equal(d0, d1))  # if gam, determine empirically
    }
    
    if(long && !isMer.y && !isMer.m) {
      out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                  d0.p=d0.p, d1.p=d1.p,
                  d0.sims=delta.0, d1.sims=delta.1,
                  z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
                  z0.p=z0.p, z1.p=z1.p,
                  z0.sims=zeta.0, z1.sims=zeta.1,
                  n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
                  n0.p=n0.p, n1.p=n1.p,
                  n0.sims=nu.0, n1.sims=nu.1,
                  tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                  tau.sims=tau,
                  d.avg=d.avg, d.avg.p=d.avg.p, d.avg.ci=d.avg.ci, d.avg.sims=delta.avg,
                  z.avg=z.avg, z.avg.p=z.avg.p, z.avg.ci=z.avg.ci, z.avg.sims=zeta.avg,
                  n.avg=n.avg, n.avg.p=n.avg.p, n.avg.ci=n.avg.ci, n.avg.sims=nu.avg,
                  boot=boot, boot.ci.type=boot.ci.type,
                  treat=treat, mediator=mediator,
                  covariates=covariates,
                  INT=INT, conf.level=conf.level,
                  model.y=model.y, model.m=model.m,
                  control.value=control.value, treat.value=treat.value,
                  nobs=n, sims=sims, call=cl,
                  robustSE = robustSE, cluster = cluster)
      class(out) <- "mediate"
    } 
    if(!long && !isMer.y && !isMer.m){
      out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                  d0.p=d0.p, d1.p=d1.p,
                  z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
                  z0.p=z0.p, z1.p=z1.p,
                  n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
                  n0.p=n0.p, n1.p=n1.p,
                  tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                  d.avg=d.avg, d.avg.p=d.avg.p, d.avg.ci=d.avg.ci,
                  z.avg=z.avg, z.avg.p=z.avg.p, z.avg.ci=z.avg.ci,
                  n.avg=n.avg, n.avg.p=n.avg.p, n.avg.ci=n.avg.ci,
                  boot=boot, boot.ci.type=boot.ci.type,
                  treat=treat, mediator=mediator,
                  covariates=covariates,
                  INT=INT, conf.level=conf.level,
                  model.y=model.y, model.m=model.m,
                  control.value=control.value, treat.value=treat.value,
                  nobs=n, sims=sims, call=cl,
                  robustSE = robustSE, cluster = cluster)
      class(out) <- "mediate"
    }
    if(long && (isMer.y || isMer.m)) {
      out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                  d0.p=d0.p, d1.p=d1.p,
                  d0.sims=delta.0, d1.sims=delta.1,
                  z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
                  z0.p=z0.p, z1.p=z1.p,
                  z0.sims=zeta.0, z1.sims=zeta.1,
                  n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
                  n0.p=n0.p, n1.p=n1.p,
                  n0.sims=nu.0, n1.sims=nu.1,
                  tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                  tau.sims=tau,
                  d.avg=d.avg, d.avg.p=d.avg.p, d.avg.ci=d.avg.ci, d.avg.sims=delta.avg,
                  z.avg=z.avg, z.avg.p=z.avg.p, z.avg.ci=z.avg.ci, z.avg.sims=zeta.avg,
                  n.avg=n.avg, n.avg.p=n.avg.p, n.avg.ci=n.avg.ci, n.avg.sims=nu.avg,
                  d0.group=d0.group, d1.group=d1.group, d0.ci.group=d0.ci.group, d1.ci.group=d1.ci.group,
                  d0.p.group=d0.p.group, d1.p.group=d1.p.group,
                  d0.sims.group=delta.0.group, d1.sims.group=delta.1.group,
                  z0.group=z0.group, z1.group=z1.group, z0.ci.group=z0.ci.group, z1.ci.group=z1.ci.group,
                  z0.p.group=z0.p.group, z1.p.group=z1.p.group,
                  z0.sims.group=zeta.0.group, z1.sims.group=zeta.1.group,
                  n0.group=n0.group, n1.group=n1.group, n0.ci.group=n0.ci.group, n1.ci.group=n1.ci.group,
                  n0.p.group=n0.p.group, n1.p.group=n1.p.group,
                  n0.sims.group=nu.0.group, n1.sims.group=nu.1.group,
                  tau.coef.group=tau.coef.group, tau.ci.group=tau.ci.group, tau.p.group=tau.p.group,
                  tau.sims.group=tau.group,
                  d.avg.group=d.avg.group, d.avg.p.group=d.avg.p.group, d.avg.ci.group=d.avg.ci.group, d.avg.sims.group=delta.avg.group,
                  z.avg.group=z.avg.group, z.avg.p.group=z.avg.p.group, z.avg.ci.group=z.avg.ci.group, z.avg.sims.group=zeta.avg.group,
                  n.avg.group=n.avg.group, n.avg.p.group=n.avg.p.group, n.avg.ci.group=n.avg.ci.group, n.avg.sims.group=nu.avg.group,
                  boot=boot, boot.ci.type=boot.ci.type,
                  treat=treat, mediator=mediator,
                  covariates=covariates,
                  INT=INT, conf.level=conf.level,
                  model.y=model.y, model.m=model.m,
                  control.value=control.value, treat.value=treat.value,
                  nobs=n, sims=sims, call=cl,
                  group.m=group.m,group.y=group.y,group.name=group.name,
                  group.id.m=group.id.m,group.id.y=group.id.y,group.id=group.id,
                  robustSE = robustSE, cluster = cluster)  
      class(out) <- "mediate.mer"
    }
    if(!long && (isMer.y || isMer.m)){
      out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                  d0.p=d0.p, d1.p=d1.p,
                  z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
                  z0.p=z0.p, z1.p=z1.p,
                  n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
                  n0.p=n0.p, n1.p=n1.p,
                  tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                  d.avg=d.avg, d.avg.p=d.avg.p, d.avg.ci=d.avg.ci,
                  z.avg=z.avg, z.avg.p=z.avg.p, z.avg.ci=z.avg.ci,
                  n.avg=n.avg, n.avg.p=n.avg.p, n.avg.ci=n.avg.ci,
                  d0.group=d0.group, d1.group=d1.group, d0.ci.group=d0.ci.group, d1.ci.group=d1.ci.group,
                  d0.p.group=d0.p.group, d1.p.group=d1.p.group,
                  z0.group=z0.group, z1.group=z1.group, z0.ci.group=z0.ci.group, z1.ci.group=z1.ci.group,
                  z0.p.group=z0.p.group, z1.p.group=z1.p.group,
                  n0.group=n0.group, n1.group=n1.group, n0.ci.group=n0.ci.group, n1.ci.group=n1.ci.group,
                  n0.p.group=n0.p.group, n1.p.group=n1.p.group,
                  tau.coef.group=tau.coef.group, tau.ci.group=tau.ci.group, tau.p.group=tau.p.group,
                  d.avg.group=d.avg.group, d.avg.p.group=d.avg.p.group, d.avg.ci.group=d.avg.ci.group,
                  z.avg.group=z.avg.group, z.avg.p.group=z.avg.p.group, z.avg.ci.group=z.avg.ci.group,
                  n.avg.group=n.avg.group, n.avg.p.group=n.avg.p.group, n.avg.ci.group=n.avg.ci.group,
                  boot=boot, boot.ci.type=boot.ci.type,
                  treat=treat, mediator=mediator,
                  covariates=covariates,
                  INT=INT, conf.level=conf.level,
                  model.y=model.y, model.m=model.m,
                  control.value=control.value, treat.value=treat.value,
                  nobs=n, sims=sims, call=cl,
                  group.m=group.m,group.y=group.y,group.name=group.name,
                  group.id.m=group.id.m,group.id.y=group.id.y,group.id=group.id,
                  robustSE = robustSE, cluster = cluster)
      class(out) <- "mediate.mer"
    }
    
    out
    
    ############################################################################
    ############################################################################
    ### CASE II: ORDERED OUTCOME
    ############################################################################
    ############################################################################
  } else {
    if(boot != TRUE){
      warning("ordered outcome model can only be used with nonparametric bootstrap - option forced")
      boot <- TRUE
    }
    
    if(isMer.m){
      stop("merMod class is not supported for ordered outcome")
    }
    
    n.ycat <- length(unique(model.response(y.data)))
    
    # Storage
    delta.1 <- matrix(NA, sims, n.ycat)
    delta.0 <- matrix(NA, sims, n.ycat)
    zeta.1 <- matrix(NA, sims, n.ycat)
    zeta.0 <- matrix(NA, sims, n.ycat)
    tau <- matrix(NA, sims, n.ycat)
    
    # Bootstrap loop begins
    for(b in 1:(sims+1)){
      
      # Resampling Step
      index <- sample(1:n, n, replace = TRUE)
      if(b == sims + 1){  # use original data for the last iteration
        index <- 1:n
      }
      
      Call.M <- model.m$call
      Call.Y <- model.y$call
      
      if(isSurvreg.m){
        if(ncol(model.m$y) > 2){
          stop("unsupported censoring type")
        }
        mname <- names(m.data)[1]
        if(substr(mname, 1, 4) != "Surv"){
          stop("refit the survival model with `Surv' used directly in model formula")
        }
        nc <- nchar(mediator)
        eventname <- substr(mname, 5 + nc + 3, nchar(mname) - 1)
        if(nchar(eventname) == 0){
          m.data.tmp <- data.frame(m.data,
                                   as.numeric(m.data[,1L][,1L]))
          names(m.data.tmp)[c(1L, ncol(m.data)+1)] <- c(mname, mediator)
        } else {
          m.data.tmp <- data.frame(m.data,
                                   as.numeric(m.data[,1L][,1L]),
                                   as.numeric(model.m$y[,2]))
          names(m.data.tmp)[c(1L, ncol(m.data)+(1:2))] <- c(mname, mediator, eventname)
        }
        Call.M$data <- m.data.tmp[index,]
      } else {
        Call.M$data <- m.data[index,]
      }
      
      Call.Y$data <- y.data[index,]
      Call.M$weights <- m.data[index,"(weights)"]
      Call.Y$weights  <- y.data[index,"(weights)"]
      new.fit.M <- eval.parent(Call.M)
      new.fit.Y <- eval.parent(Call.Y)
      
      if(isOrdered.m && length(unique(y.data[index,mediator]))!=m){
        # Modify the coefficients when mediator has empty cells
        coefnames.y <- names(model.y$coefficients)
        coefnames.new.y <- names(new.fit.Y$coefficients)
        new.fit.Y.coef <- rep(0, length(coefnames.y))
        names(new.fit.Y.coef) <- coefnames.y
        new.fit.Y.coef[coefnames.new.y] <- new.fit.Y$coefficients
        new.fit.Y$coefficients <- new.fit.Y.coef
      }
      
      #####################################
      # Mediator Predictions
      #####################################
      pred.data.t <- pred.data.c <- m.data
      
      if(isFactorT){
        pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
        pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
      } else {
        pred.data.t[,treat] <- cat.1
        pred.data.c[,treat] <- cat.0
      }
      
      if(!is.null(covariates)){
        for(p in 1:length(covariates)){
          vl <- names(covariates[p])
          if(is.factor(pred.data.t[,vl])){
            pred.data.t[,vl] <- pred.data.c[,vl] <- factor(covariates[[p]], levels = levels(m.data[,vl]))
          } else {
            pred.data.t[,vl] <- pred.data.c[,vl] <- covariates[[p]]
          }
        }
      }
      
      ### Case II-a: GLM Mediator (including GAMs)
      if(isGlm.m){
        
        muM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
        muM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)
        
        if(FamilyM == "poisson"){
          PredictM1 <- rpois(n, lambda = muM1)
          PredictM0 <- rpois(n, lambda = muM0)
        } else if (FamilyM == "Gamma") {
          shape <- gamma.shape(model.m)$alpha
          PredictM1 <- rgamma(n, shape = shape, scale = muM1/shape)
          PredictM0 <- rgamma(n, shape = shape, scale = muM0/shape)
        } else if (FamilyM == "binomial"){
          PredictM1 <- rbinom(n, size = 1, prob = muM1)
          PredictM0 <- rbinom(n, size = 1, prob = muM0)
        } else if (FamilyM == "gaussian"){
          sigma <- sqrt(summary(model.m)$dispersion)
          error <- rnorm(n, mean=0, sd=sigma)
          PredictM1 <- muM1 + error
          PredictM0 <- muM0 + error
        } else if (FamilyM == "inverse.gaussian"){
          disp <- summary(model.m)$dispersion
          PredictM1 <- SuppDists::rinvGauss(n, nu = muM1, lambda = 1/disp)
          PredictM0 <- SuppDists::rinvGauss(n, nu = muM0, lambda = 1/disp)
        } else {
          stop("unsupported glm family")
        }
        
        ### Case II-b: Ordered Mediator
      } else if(isOrdered.m) {
        probs_m1 <- predict(new.fit.M, type="probs", newdata=pred.data.t)
        probs_m0 <- predict(new.fit.M, type="probs", newdata=pred.data.c)
        draws_m1 <- matrix(NA, n, m)
        draws_m0 <- matrix(NA, n, m)
        
        for(ii in 1:n){
          draws_m1[ii,] <- t(rmultinom(1, 1, prob = probs_m1[ii,]))
          draws_m0[ii,] <- t(rmultinom(1, 1, prob = probs_m0[ii,]))
        }
        PredictM1 <- apply(draws_m1, 1, indexmax)
        PredictM0 <- apply(draws_m0, 1, indexmax)
        
        ### Case II-c: Quantile Regression for Mediator
      } else if(isRq.m){
        # Use inverse transform sampling to predict M
        call.new <- new.fit.M$call
        call.new$tau <- runif(n)
        newfits <- eval.parent(call.new)
        tt <- delete.response(terms(new.fit.M))
        m.t <- model.frame(tt, pred.data.t, xlev = new.fit.M$xlevels)
        m.c <- model.frame(tt, pred.data.c, xlev = new.fit.M$xlevels)
        X.t <- model.matrix(tt, m.t, contrasts = new.fit.M$contrasts)
        X.c <- model.matrix(tt, m.c, contrasts = new.fit.M$contrasts)
        rm(tt, m.t, m.c)
        PredictM1 <- rowSums(X.t * t(newfits$coefficients))
        PredictM0 <- rowSums(X.c * t(newfits$coefficients))
        rm(newfits, X.t, X.c)
        
        ### Case II-d: Linear
      } else if(isLm.m){
        sigma <- summary(new.fit.M)$sigma
        error <- rnorm(n, mean=0, sd=sigma)
        PredictM1 <- predict(new.fit.M, type="response",
                             newdata=pred.data.t) + error
        PredictM0 <- predict(new.fit.M, type="response",
                             newdata=pred.data.c) + error
        rm(error)
        
        ### Case I-2-e: Survreg
      } else if(isSurvreg.m){
        dd <- survival::survreg.distributions[[new.fit.M$dist]]
        if (is.null(dd$itrans)){
          itrans <- function(x) x
        } else {
          itrans <- dd$itrans
        }
        dname <- dd$dist
        if(is.null(dname)){
          dname <- new.fit.M$dist
        }
        scale <- new.fit.M$scale
        lpM1 <- predict(new.fit.M, newdata=pred.data.t, type="linear")
        lpM0 <- predict(new.fit.M, newdata=pred.data.c, type="linear")
        error <- switch(dname,
                        extreme = log(rweibull(n, shape=1, scale=1)),
                        gaussian = rnorm(n),
                        logistic = rlogis(n),
                        t = rt(n, df=dd$parms))
        PredictM1 <- as.numeric(itrans(lpM1 + scale * error))
        PredictM0 <- as.numeric(itrans(lpM0 + scale * error))
        rm(error)
        
      } else {
        stop("mediator model is not yet implemented")
      }
      
      #####################################
      #  Outcome Predictions
      #####################################
      effects.tmp <- array(NA, dim = c(n, n.ycat, 4))
      for(e in 1:4){
        tt <- switch(e, c(1,1,1,0), c(0,0,1,0), c(1,0,1,1), c(1,0,0,0))
        pred.data.t <- pred.data.c <- y.data
        
        if(!is.null(covariates)){
          for(p in 1:length(covariates)){
            vl <- names(covariates[p])
            if(is.factor(pred.data.t[,vl])){
              pred.data.t[,vl] <- pred.data.c[,vl] <- factor(covariates[[p]], levels = levels(y.data[,vl]))
            } else {
              pred.data.t[,vl] <- pred.data.c[,vl] <- covariates[[p]]
            }
          }
        }
        
        # Set treatment values
        cat.t <- ifelse(tt[1], cat.1, cat.0)
        cat.c <- ifelse(tt[2], cat.1, cat.0)
        cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
        cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
        if(isFactorT){
          pred.data.t[,treat] <- factor(cat.t, levels = t.levels)
          pred.data.c[,treat] <- factor(cat.c, levels = t.levels)
          if(!is.null(control)){
            pred.data.t[,control] <- factor(cat.t.ctrl, levels = t.levels)
            pred.data.c[,control] <- factor(cat.c.ctrl, levels = t.levels)
          }
        } else {
          pred.data.t[,treat] <- cat.t
          pred.data.c[,treat] <- cat.c
          if(!is.null(control)){
            pred.data.t[,control] <- cat.t.ctrl
            pred.data.c[,control] <- cat.c.ctrl
          }
        }
        
        # Set mediator values
        PredictM1.tmp <- PredictM1
        PredictM0.tmp <- PredictM0
        PredictMt <- PredictM1 * tt[3] + PredictM0 * (1 - tt[3])
        PredictMc <- PredictM1 * tt[4] + PredictM0 * (1 - tt[4])
        if(isFactorM) {
          pred.data.t[,mediator] <- factor(PredictMt, levels=1:m, labels=m.levels)
          pred.data.c[,mediator] <- factor(PredictMc, levels=1:m, labels=m.levels)
        } else {
          pred.data.t[,mediator] <- PredictMt
          pred.data.c[,mediator] <- PredictMc
        }
        probs_p1 <- predict(new.fit.Y, newdata=pred.data.t, type="probs")
        probs_p0 <- predict(new.fit.Y, newdata=pred.data.c, type="probs")
        effects.tmp[,,e] <- probs_p1 - probs_p0
        rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
      }
      
      # Compute all QoIs
      if(b == sims+1){
        d1 <- apply(effects.tmp[,,1], 2, weighted.mean, w=weights)
        d0 <- apply(effects.tmp[,,2], 2, weighted.mean, w=weights)
        z1 <- apply(effects.tmp[,,3], 2, weighted.mean, w=weights)
        z0 <- apply(effects.tmp[,,4], 2, weighted.mean, w=weights)
      } else {
        delta.1[b,] <- apply(effects.tmp[,,1], 2, weighted.mean, w=weights)
        delta.0[b,] <- apply(effects.tmp[,,2], 2, weighted.mean, w=weights)
        zeta.1[b,] <- apply(effects.tmp[,,3], 2, weighted.mean, w=weights)
        zeta.0[b,] <- apply(effects.tmp[,,4], 2, weighted.mean, w=weights)
      }
      
    }  # Bootstrap loop ends
    
    tau.coef <- (d1 + d0 + z1 + z0)/2
    tau <- (zeta.1 + zeta.0 + delta.0 + delta.1)/2
    
    ########################################################################
    ## Compute Outputs and Put Them Together
    ########################################################################
    low <- (1 - conf.level)/2
    high <- 1 - low

    if(boot.ci.type == "bca"){
        BC.CI <- function(theta){
            z.inv <- length(theta[theta < mean(theta)])/sims
            z <- qnorm(z.inv)
            U <- (sims - 1) * (mean(theta) - theta)
            top <- sum(U^3)
            under <- 6 * (sum(U^2))^{3/2}
            a <- top / under
            lower.inv <-  pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
            lower2 <- lower <- quantile(theta, lower.inv)
            upper.inv <-  pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
            upper2 <- upper <- quantile(theta, upper.inv)
            return(c(lower, upper))      
        }
        d0.ci <- BC.CI(delta.0)
        d1.ci <- BC.CI(delta.1)
        tau.ci <- BC.CI(tau)
        z1.ci <- BC.CI(zeta.1)
        z0.ci <- BC.CI(zeta.0)
    } else {
        CI <- function(theta){
            return(quantile(theta, c(low, high), na.rm = TRUE))
        }
        d0.ci <- apply(delta.0, 2, CI)
        d1.ci <- apply(delta.1, 2, CI)
        tau.ci <- apply(tau, 2, CI)
        z1.ci <- apply(zeta.1, 2, CI)
        z0.ci <- apply(zeta.0, 2, CI)
    }
    
    # p-values
    d0.p <- d1.p <- z0.p <- z1.p <- tau.p <- rep(NA, n.ycat)
    for(i in 1:n.ycat){
      d0.p[i] <- pval(delta.0[,i], d0[i])
      d1.p[i] <- pval(delta.1[,i], d1[i])
      z0.p[i] <- pval(zeta.0[,i], z0[i])
      z1.p[i] <- pval(zeta.1[,i], z1[i])
      tau.p[i] <- pval(tau[,i], tau.coef[i])
    }
    
    # Detect whether models include T-M interaction
    INT <- paste(treat,mediator,sep=":") %in% attr(model.y$terms,"term.labels") |
      paste(mediator,treat,sep=":") %in% attr(model.y$terms,"term.labels")
    
    if(long) {
      out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                  d0.p=d0.p, d1.p=d1.p,
                  d0.sims=delta.0, d1.sims=delta.1,
                  tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                  z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
                  z0.p=z0.p, z1.p=z1.p,
                  z1.sims=zeta.1, z0.sims=zeta.0, tau.sims=tau,
                  boot=boot, boot.ci.type=boot.ci.type,
                  treat=treat, mediator=mediator,
                  covariates=covariates,
                  INT=INT, conf.level=conf.level,
                  model.y=model.y, model.m=model.m,
                  control.value=control.value, treat.value=treat.value, 
                  nobs=n, sims=sims, call=cl,
                  robustSE = robustSE, cluster = cluster)
    } else {
      out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                  d0.p=d0.p, d1.p=d1.p,
                  tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                  z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
                  z0.p=z0.p, z1.p=z1.p,
                  boot=boot, boot.ci.type=boot.ci.type,
                  treat=treat, mediator=mediator,
                  covariates=covariates,
                  INT=INT, conf.level=conf.level,
                  model.y=model.y, model.m=model.m,
                  control.value=control.value, treat.value=treat.value, 
                  nobs=n, sims=sims, call=cl,
                  robustSE = robustSE, cluster = cluster)
    }
    class(out) <- "mediate.order"
    out
  }
}

##################################################################
summary.mediate <- function(object, ...){
  structure(object, class = c("summary.mediate", class(object)))
}

print.summary.mediate <- function(x, ...){
  clp <- 100 * x$conf.level
  cat("\n")
  cat("Causal Mediation Analysis \n\n")
  if(x$boot){
    if(x$boot.ci.type == "perc"){
      cat("Nonparametric Bootstrap Confidence Intervals with the Percentile Method\n\n")
    } else if(x$boot.ci.type == "bca"){
      cat("Nonparametric Bootstrap Confidence Intervals with the BCa Method\n\n") 
    }
  } else {
    cat("Quasi-Bayesian Confidence Intervals\n\n")
  }
  
  if(!is.null(x$covariates)){
    cat("(Inference Conditional on the Covariate Values Specified in `covariates')\n\n")
  }
  
  isLinear.y <- (	(class(x$model.y)[1] %in% c("lm", "rq")) ||
                    (inherits(x$model.y, "glm") &&
                       x$model.y$family$family == "gaussian" &&
                       x$model.y$family$link == "identity") ||
                    (inherits(x$model.y, "survreg") &&
                       x$model.y$dist == "gaussian") )
  
  printone <- !x$INT && isLinear.y
  
  if (printone){
    # Print only one set of values if lmY/quanY/linear gamY without interaction
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    rownames(smat) <- c("ACME", "ADE",
                        "Total Effect", "Prop. Mediated")
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    rownames(smat) <- c("ACME (control)", "ACME (treated)",
                        "ADE (control)", "ADE (treated)",
                        "Total Effect",
                        "Prop. Mediated (control)",
                        "Prop. Mediated (treated)",
                        "ACME (average)",
                        "ADE (average)",
                        "Prop. Mediated (average)")
  }
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                      paste(clp, "% CI Upper", sep=""), "p-value")
  printCoefmat(smat, digits=3)
  cat("\n")
  cat("Sample Size Used:", x$nobs,"\n\n")
  cat("\n")
  cat("Simulations:", x$sims,"\n\n")
  invisible(x)
}

#########################################################################
summary.mediate.mer <- function(object, output=c("default","byeffect","bygroup"),...){
  output <- match.arg(output)
  switch(output,
         default = structure(object, class = c("summary.mediate.mer", class(object))),
         byeffect = structure(object, class = c("summary.mediate.mer.2", class(object))),   
         bygroup = structure(object, class = c("summary.mediate.mer.3", class(object))))
}

#########################################################################
print.summary.mediate.mer <- function(x,...){  
  clp <- 100 * x$conf.level
  cat("\n")
  cat("Causal Mediation Analysis \n\n")

  if(x$boot){
    if(x$boot.ci.type == "perc"){
      cat("Nonparametric Bootstrap Confidence Intervals with the Percentile Method\n\n")
    } else if(x$boot.ci.type == "bca"){
      cat("Nonparametric Bootstrap Confidence Intervals with the BCa Method\n\n") 
    }
  } else {
    cat("Quasi-Bayesian Confidence Intervals\n\n")
  }
  
  if(!is.null(x$covariates)){
    cat("(Inference Conditional on the Covariate Values Specified in `covariates')\n\n")
  }
  
  cat("Mediator Groups:", x$group.m,"\n\n")
  cat("Outcome Groups:", x$group.y,"\n\n")
  cat("Output Based on Overall Averages Across Groups","\n\n")
  
  isLinear.y <- (	(class(x$model.y)[1] %in% c("lm", "rq")) || # lm or quantile
                    (inherits(x$model.y, "glm") &&
                       x$model.y$family$family == "gaussian" &&
                       x$model.y$family$link == "identity") ||      # glm normal
                    (inherits(x$model.y, "survreg") &&
                       x$model.y$dist == "gaussian") ||             # surv normal
                    (inherits(x$model.y, "merMod") &&
                       x$model.y@call[[1]] == "lmer") )          # lmer
  
  printone <- !x$INT && isLinear.y
  
  if(printone){
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    rownames(smat) <- c("ACME", "ADE",
                        "Total Effect", "Prop. Mediated")
  }else{
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    rownames(smat) <- c("ACME (control)", "ACME (treated)",
                        "ADE (control)", "ADE (treated)",
                        "Total Effect",
                        "Prop. Mediated (control)",
                        "Prop. Mediated (treated)",
                        "ACME (average)",
                        "ADE (average)",
                        "Prop. Mediated (average)")
  }
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                      paste(clp, "% CI Upper", sep=""), "p-value")
  printCoefmat(smat, digits=3)
  cat("\n")
  cat("Sample Size Used:", x$nobs,"\n\n")
  cat("\n")
  cat("Simulations:", x$sims,"\n\n")
  invisible(x)
}
#########################################################################
print.summary.mediate.mer.2 <- function(x,...){
  clp <- 100 * x$conf.level
  cat("\n")
  cat("Causal Mediation Analysis \n\n")
  
  if(x$boot){
    if(x$boot.ci.type == "perc"){
      cat("Nonparametric Bootstrap Confidence Intervals with the Percentile Method\n\n")
    } else if(x$boot.ci.type == "bca"){
      cat("Nonparametric Bootstrap Confidence Intervals with the BCa Method\n\n") 
    }
  } else {
    cat("Quasi-Bayesian Confidence Intervals\n\n")
  }
  
  if(!is.null(x$covariates)){
    cat("(Inference Conditional on the Covariate Values Specified in `covariates')\n\n")
  }
  
  cat("Mediator Groups:", x$group.m,"\n\n")
  cat("Outcome Groups:", x$group.y,"\n\n")
  cat("Output Based on", x$group.name,"\n\n")
  
  if(is.factor(x$group.id)){
    gname <- levels(x$group.id)
  } else {
    gname <- sort(unique(x$group.id))
  }
  
  isLinear.y <- (	(class(x$model.y)[1] %in% c("lm", "rq")) || # lm or quantile
                    (inherits(x$model.y, "glm") &&
                       x$model.y$family$family == "gaussian" &&
                       x$model.y$family$link == "identity") ||      # glm normal
                    (inherits(x$model.y, "survreg") &&
                       x$model.y$dist == "gaussian") ||             # surv normal
                    (inherits(x$model.y, "merMod") &&
                       x$model.y@call[[1]] == "lmer") )          # lmer
  
  printone <- !x$INT && isLinear.y
  
  if(printone){
    cat("ACME","\n\n")
    smat<-cbind(x$d0.group,x$d0.ci.group,x$d0.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")
    cat("ADE","\n\n")
    smat<-cbind(x$z0.group, x$z0.ci.group, x$z0.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")
    cat("Total Effect","\n\n")
    smat<-cbind(x$tau.coef.group, x$tau.ci.group, x$tau.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n") 
    cat("Prop. Mediated","\n\n")
    smat<-cbind(x$n0.group, x$n0.ci.group, x$n0.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")  
    
    cat("\n")
    cat("Sample Size Used:", x$nobs,"\n\n")
    cat("\n")
    cat("Simulations:", x$sims,"\n\n")
    invisible(x) 
  }else{
    cat("ACME (control)","\n\n")
    smat<-cbind(x$d0.group,x$d0.ci.group,x$d0.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")
    cat("ACME (treated)","\n\n")
    smat<-cbind(x$d1.group, x$d1.ci.group, x$d1.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")
    cat("ADE (control)","\n\n")
    smat<-cbind(x$z0.group, x$z0.ci.group, x$z0.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")
    cat("ADE (treated)","\n\n")
    smat<-cbind(x$z1.group, x$z1.ci.group, x$z1.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")    
    cat("Total Effect","\n\n")
    smat<-cbind(x$tau.coef.group, x$tau.ci.group, x$tau.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n") 
    cat("Prop. Mediated (control)","\n\n")
    smat<-cbind(x$n0.group, x$n0.ci.group, x$n0.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")  
    cat("Prop. Mediated (treated)","\n\n")
    smat<-cbind(x$n1.group, x$n1.ci.group, x$n1.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n") 
    cat("ACME (average)","\n\n")
    smat<-cbind(x$d.avg.group, x$d.avg.ci.group, x$d.avg.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n") 
    cat("ADE (average)","\n\n")
    smat<-cbind(x$z.avg.group, x$z.avg.ci.group, x$z.avg.p.group)
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")
    cat("Prop. Mediated (average)","\n\n")
    smat<-cbind(x$n.avg.group, x$n.avg.ci.group, x$n.avg.p.group)   
    rownames(smat) <- gname
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                        paste(clp, "% CI Upper", sep=""), "p-value") 
    printCoefmat(smat, digits=3)
    cat("\n")                          
    
    cat("\n")
    cat("Sample Size Used:", x$nobs,"\n\n")
    cat("\n")
    cat("Simulations:", x$sims,"\n\n")
    invisible(x)
  }
}
#########################################################################
print.summary.mediate.mer.3 <- function(x,...){
  clp <- 100 * x$conf.level
  cat("\n")
  cat("Causal Mediation Analysis \n\n")

  if(x$boot){
    if(x$boot.ci.type == "perc"){
      cat("Nonparametric Bootstrap Confidence Intervals with the Percentile Method\n\n")
    } else if(x$boot.ci.type == "bca"){
      cat("Nonparametric Bootstrap Confidence Intervals with the BCa Method\n\n") 
    }
  } else {
    cat("Quasi-Bayesian Confidence Intervals\n\n")
  }
  
  if(!is.null(x$covariates)){
    cat("(Inference Conditional on the Covariate Values Specified in `covariates')\n\n")
  }
  
  cat("Mediator Groups:", x$group.m,"\n\n")
  cat("Outcome Groups:", x$group.y,"\n\n")
  cat("Output Based on", x$group.name,"\n\n")
  
  if(is.factor(x$group.id)){
    gname <- levels(x$group.id)
  } else {
    gname <- sort(unique(x$group.id))
  }
  G<-length(gname)
  isLinear.y <- (	(class(x$model.y)[1] %in% c("lm", "rq")) || # lm or quantile
                    (inherits(x$model.y, "glm") &&
                       x$model.y$family$family == "gaussian" &&
                       x$model.y$family$link == "identity") ||      # glm normal
                    (inherits(x$model.y, "survreg") &&
                       x$model.y$dist == "gaussian") ||             # surv normal
                    (inherits(x$model.y, "merMod") &&
                       x$model.y@call[[1]] == "lmer") )          # lmer
  
  printone <- !x$INT && isLinear.y
  
  if(printone){
    for (g in 1:G){  
      cat("\n")
      cat("Group:",gname[g],"\n") 
      smat <- c(x$d0.group[g], x$d0.ci.group[g,], x$d0.p.group[g])
      smat <- rbind(smat, c(x$z0.group[g], x$z0.ci.group[g,], x$z0.p.group[g]))
      smat <- rbind(smat, c(x$tau.coef.group[g], x$tau.ci.group[g,], x$tau.p.group[g]))
      smat <- rbind(smat, c(x$n0.group[g], x$n0.ci.group[g,], x$n0.p.group[g]))
      rownames(smat) <- c("ACME", 
                          "ADE", 
                          "Total Effect",
                          "Prop. Mediated") 
      colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                          paste(clp, "% CI Upper", sep=""), "p-value")             
      printCoefmat(smat, digits=3)
    }
    cat("\n")
    cat("Sample Size Used:", x$nobs,"\n\n")
    cat("\n")
    cat("Simulations:", x$sims,"\n\n")
    invisible(x)    
  }else{
    for (g in 1:G){  
      cat("\n")
      cat("Group:",gname[g],"\n") 
      smat <- c(x$d0.group[g], x$d0.ci.group[g,], x$d0.p.group[g])
      smat <- rbind(smat, c(x$d1.group[g], x$d1.ci.group[g,], x$d1.p.group[g]))
      smat <- rbind(smat, c(x$z0.group[g], x$z0.ci.group[g,], x$z0.p.group[g]))
      smat <- rbind(smat, c(x$z1.group[g], x$z1.ci.group[g,], x$z1.p.group[g]))
      smat <- rbind(smat, c(x$tau.coef.group[g], x$tau.ci.group[g,], x$tau.p.group[g]))
      smat <- rbind(smat, c(x$n0.group[g], x$n0.ci.group[g,], x$n0.p.group[g]))
      smat <- rbind(smat, c(x$n1.group[g], x$n1.ci.group[g,], x$n1.p.group[g]))
      smat <- rbind(smat, c(x$d.avg.group[g], x$d.avg.ci.group[g,], x$d.avg.p.group[g]))
      smat <- rbind(smat, c(x$z.avg.group[g], x$z.avg.ci.group[g,], x$z.avg.p.group[g]))
      smat <- rbind(smat, c(x$n.avg.group[g], x$n.avg.ci.group[g,], x$n.avg.p.group[g]))   
      rownames(smat) <- c("ACME (control)", "ACME (treated)",
                          "ADE (control)", "ADE (treated)",
                          "Total Effect",
                          "Prop. Mediated (control)",
                          "Prop. Mediated (treated)",
                          "ACME (average)",
                          "ADE (average)",
                          "Prop. Mediated (average)") 
      colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                          paste(clp, "% CI Upper", sep=""), "p-value")             
      printCoefmat(smat, digits=3)
    }
    cat("\n")
    cat("Sample Size Used:", x$nobs,"\n\n")
    cat("\n")
    cat("Simulations:", x$sims,"\n\n")
    invisible(x)
  }
}
#########################################################################
summary.mediate.order <- function(object, ...){
  structure(object, class = c("summary.mediate.order", class(object)))
}



print.summary.mediate.order <- function(x, ...){
  tab.d0 <- rbind(x$d0, x$d0.ci, x$d0.p)
  tab.d1 <- rbind(x$d1, x$d1.ci, x$d1.p)
  tab.z0 <- rbind(x$z0, x$z0.ci, x$z0.p)
  tab.z1 <- rbind(x$z1, x$z1.ci, x$z1.p)
  tab.tau <- rbind(x$tau.coef, x$tau.ci, x$tau.p)
  
  # Outcome Table Labels
  y.lab <- sort(unique(levels(model.frame(x$model.y)[,1])))
  out.names <- c()
  for(i in 1:length(y.lab)){
    out.names.tmp <- paste("Pr(Y=",y.lab[i],")",sep="")
    out.names <- c(out.names, out.names.tmp)
  }
  
  # Label Tables
  rownames(tab.d0)[1] <- "ACME (control) "
  rownames(tab.d0)[4] <- "p-value "
  colnames(tab.d0) <- out.names
  rownames(tab.d1)[1] <- "ACME (treated) "
  rownames(tab.d1)[4] <- "p-value "
  colnames(tab.d1) <- out.names
  rownames(tab.z0)[1] <- "ADE (control)  "
  rownames(tab.z0)[4] <- "p-value "
  colnames(tab.z0) <- out.names
  rownames(tab.z1)[1] <- "ADE (treated)  "
  rownames(tab.z1)[4] <- "p-value "
  colnames(tab.z1) <- out.names
  rownames(tab.tau)[1] <- "Total Effect  "
  rownames(tab.tau)[4] <- "p-value "
  colnames(tab.tau) <- out.names
  
  cat("\n")
  cat("Causal Mediation Analysis \n\n")

  if(x$boot.ci.type == "perc"){
    cat("Nonparametric Bootstrap Confidence Intervals with the Percentile Method\n\n")
  } else if(x$boot.ci.type == "bca"){
    cat("Nonparametric Bootstrap Confidence Intervals with the BCa Method\n\n") 
  }
  
  if(!is.null(x$covariates)){
    cat("(Inference Conditional on the Covariate Values Specified in `covariates')\n\n")
  }
  print(tab.d0, digits=3)
  cat("\n")
  print(tab.d1, digits=3)
  cat("\n")
  print(tab.z0, digits=3)
  cat("\n")
  print(tab.z1, digits=3)
  cat("\n")
  print(tab.tau, digits=3)
  cat("\n\n")
  cat("Sample Size Used:", x$nobs,"\n\n")
  cat("\n\n")
  cat("Simulations:", x$sims,"\n\n")
  invisible(x)
}
#########################################################################
plot.process <- function(model) {
  coef.vec.1 <- c(model$d1, model$z1)
  lower.vec.1 <- c(model$d1.ci[1], model$z1.ci[1])
  upper.vec.1 <- c(model$d1.ci[2], model$z1.ci[2])
  tau.vec <- c(model$tau.coef,model$tau.ci[1],model$tau.ci[2])
  range.1 <- range(model$d1.ci[1], model$z1.ci[1],model$tau.ci[1],
                   model$d1.ci[2], model$z1.ci[2],model$tau.ci[2])
  
  coef.vec.0 <- c(model$d0, model$z0)
  lower.vec.0 <- c(model$d0.ci[1], model$z0.ci[1])
  upper.vec.0 <- c(model$d0.ci[2], model$z0.ci[2])
  range.0 <- range(model$d0.ci[1], model$z0.ci[1],model$tau.ci[1],
                   model$d0.ci[2], model$z0.ci[2],model$tau.ci[2])
  
  return(list(coef.vec.1=coef.vec.1, lower.vec.1=lower.vec.1,
              upper.vec.1=upper.vec.1, coef.vec.0=coef.vec.0,
              lower.vec.0=lower.vec.0, upper.vec.0=upper.vec.0, tau.vec=tau.vec,
              range.1=range.1, range.0=range.0))
}

#########################################################################
plot.mediate <- function(x, treatment = NULL, labels = NULL,
                         effect.type = c("indirect","direct","total"),
                         xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                         main = NULL, lwd = 1.5, cex = .85,
                         col = "black", ...){
  # Determine which graph to plot
  isLinear.y <- (	(class(x$model.y)[1] %in% c("lm", "rq")) ||
                    (inherits(x$model.y, "glm") &&
                       x$model.y$family$family == "gaussian" &&
                       x$model.y$family$link == "identity") ||
                    (inherits(x$model.y, "survreg") &&
                       x$model.y$dist == "gaussian") )
  
  printone <- !x$INT && isLinear.y
  
  effect.type <- match.arg(effect.type, several.ok=TRUE)
  IND <- "indirect" %in% effect.type
  DIR <- "direct" %in% effect.type
  TOT <- "total" %in% effect.type
  
  if(is.null(treatment)){
    if(printone){
      treatment <- 1
    } else {
      treatment <- c(0,1)
    }
  } else {
    treatment <- switch(treatment,
                        control = 0,
                        treated = 1,
                        both = c(0,1))
  }
  
  param <- plot.process(x)
  
  y.axis <- (IND + DIR + TOT):1
  
  # Set xlim
  if(is.null(xlim)){
    if(length(treatment) > 1) {
      xlim <- range(param$range.1, param$range.0) * 1.2
    } else if (treatment == 1){
      xlim <- param$range.1 * 1.2
    } else {
      xlim <- param$range.0 * 1.2
    }
  }
  
  # Set ylim
  if(is.null(ylim)){
    ylim <- c(min(y.axis) - 0.5, max(y.axis) + 0.5)
  }
  
  # Create blank plot first
  plot(rep(0,IND+DIR+TOT), y.axis, type = "n", xlab = xlab, ylab = ylab,
       yaxt = "n", xlim = xlim, ylim = ylim, main = main, ...)
  
  # Set offset values depending on number of bars to plot
  if(length(treatment) == 1){
    adj <- 0
  } else {
    adj <- 0.05
  }
  
  if(1 %in% treatment){
    if(IND && DIR) {
      points(param$coef.vec.1, y.axis[1:2] + adj, type = "p", pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1, y.axis[1:2] + adj, param$upper.vec.1, y.axis[1:2] + adj,
               lwd = lwd, col = col)
    }
    if(IND && !DIR) {
      points(param$coef.vec.1[1], y.axis[1] + adj, type = "p", pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1[1], y.axis[1] + adj, param$upper.vec.1[1], y.axis[1] + adj,
               lwd = lwd, col = col)
    }
    if(!IND && DIR) {
      points(param$coef.vec.1[2], y.axis[1] + adj, type = "p", pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1[2], y.axis[1] + adj, param$upper.vec.1[2], y.axis[1] + adj,
               lwd = lwd, col = col)
    }    
  }
  if(0 %in% treatment) {
    if(IND && DIR) {
      points(param$coef.vec.0, y.axis[1:2] - adj, type = "p", pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0, y.axis[1:2] - adj, param$upper.vec.0, y.axis[1:2] - adj,
               lwd = lwd, lty = 3, col = col)
    }
    if(IND && !DIR) {
      points(param$coef.vec.0[1], y.axis[1] - adj, type = "p", pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0[1], y.axis[1] - adj, param$upper.vec.0[1], y.axis[1] - adj,
               lwd = lwd, lty = 3, col = col)
    }
    if(!IND && DIR) {
      points(param$coef.vec.0[2], y.axis[1] - adj, type = "p", pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0[2], y.axis[1] - adj, param$upper.vec.0[2], y.axis[1] - adj,
               lwd = lwd, lty = 3, col = col)
    }
  }
  if (TOT) {
    points(param$tau.vec[1], 1 , type = "p", pch = 19, cex = cex, col = col)
    segments(param$tau.vec[2], 1 , param$tau.vec[3], 1 ,
             lwd = lwd, col = col)
  }
  
  if(is.null(labels)){
    labels <- c("ACME","ADE","Total\nEffect")[c(IND,DIR,TOT)]
  }
  axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE, ...)
  abline(v = 0, lty = 2)
}

#########################################################################
plot.process.mer <- function(model) {
  coef.vec.1 <- c(model$d1, model$z1)
  lower.vec.1 <- c(model$d1.ci[1], model$z1.ci[1])
  upper.vec.1 <- c(model$d1.ci[2], model$z1.ci[2])
  tau.vec <- c(model$tau.coef,model$tau.ci[1],model$tau.ci[2])
  range.1 <- range(model$d1.ci[1], model$z1.ci[1],model$tau.ci[1],
                   model$d1.ci[2], model$z1.ci[2],model$tau.ci[2])
  
  coef.vec.0 <- c(model$d0, model$z0)
  lower.vec.0 <- c(model$d0.ci[1], model$z0.ci[1])
  upper.vec.0 <- c(model$d0.ci[2], model$z0.ci[2])
  range.0 <- range(model$d0.ci[1], model$z0.ci[1],model$tau.ci[1],
                   model$d0.ci[2], model$z0.ci[2],model$tau.ci[2])
  
  coef.vec.0.group <- cbind(model$d0.group, model$z0.group)
  lower.vec.0.group <- cbind(model$d0.ci.group[,1], model$z0.ci.group[,1])
  upper.vec.0.group <- cbind(model$d0.ci.group[,2], model$z0.ci.group[,2])
  
  coef.vec.1.group <- cbind(model$d1.group, model$z1.group)
  lower.vec.1.group <- cbind(model$d1.ci.group[,1], model$z1.ci.group[,1])
  upper.vec.1.group <- cbind(model$d1.ci.group[,2], model$z1.ci.group[,2])
  
  tau.vec.group <- cbind(model$tau.coef.group, model$tau.ci.group[,1], model$tau.ci.group[,2])
  
  return(list(coef.vec.1=coef.vec.1, lower.vec.1=lower.vec.1,
              upper.vec.1=upper.vec.1, coef.vec.0=coef.vec.0,
              lower.vec.0=lower.vec.0, upper.vec.0=upper.vec.0, tau.vec=tau.vec,
              range.1=range.1, range.0=range.0,coef.vec.1.group=coef.vec.1.group, lower.vec.1.group=lower.vec.1.group,
              upper.vec.1.group=upper.vec.1.group, coef.vec.0.group=coef.vec.0.group,
              lower.vec.0.group=lower.vec.0.group, upper.vec.0.group=upper.vec.0.group, tau.vec.group=tau.vec.group))
}

#########################################################################
plot.mediate.mer <- function(x, treatment = NULL, group.plots = FALSE,
                             ask = prod(par("mfcol")) < nplots, 
                             xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                             main = NULL, lwd = 1.5, cex = .85,
                             col = "black", ...){
  
  param <- plot.process.mer(x)
  
  isLinear.y <- (	(class(x$model.y)[1] %in% c("lm", "rq")) || # lm or quantile
                    (inherits(x$model.y, "glm") &&
                       x$model.y$family$family == "gaussian" &&
                       x$model.y$family$link == "identity") ||      # glm normal
                    (inherits(x$model.y, "survreg") &&
                       x$model.y$dist == "gaussian") ||             # surv normal
                    (inherits(x$model.y, "merMod") &&
                       x$model.y@call[[1]] == "lmer") )          # lmer
  
  printone <- !x$INT && isLinear.y
  
  if(is.null(treatment)){
    if(printone){
      treatment <- 1
    } else {
      treatment <- c(0,1)
    }
  } else {
    treatment <- switch(treatment,
                        control = 0,
                        treated = 1,
                        both = c(0,1))
  }
  
  nplots <- 1 + group.plots * (2 * length(treatment) + 1)  # n of plots necessary
  
  if(ask){
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  # 1. Summary plot for population effects
  
  labels = c("ACME","ADE","Total\nEffect")
  y.axis <- c(length(param$coef.vec.1):.5)
  y.axis <- y.axis + 1
  
  if(is.null(xlim)){
    if(length(treatment) > 1) {
      xlim <- range(param$range.1, param$range.0) * 1.2
    } else if (treatment == 1){
      xlim <- param$range.1 * 1.2
    } else {
      xlim <- param$range.0 * 1.2
    }
  }
  
  if(is.null(ylim)){
    ylim <- c(min(y.axis) -1- 0.5, max(y.axis) + 0.5)
  }
  
  plot(param$coef.vec.1, y.axis, type = "n", xlab = xlab, ylab = ylab,
       yaxt = "n", xlim = xlim, ylim = ylim, main = main, ...)
  
  if(length(treatment) == 1){
    adj <- 0
  } else {
    adj <- 0.05
  }
  
  if(1 %in% treatment){
    points(param$coef.vec.1, y.axis + adj, type = "p", pch = 19, cex = cex, col = col)
    segments(param$lower.vec.1, y.axis + adj, param$upper.vec.1, y.axis + adj,
             lwd = lwd, col = col)
    points(param$tau.vec[1], 1, type = "p", pch = 19, cex = cex, col = col)
    segments(param$tau.vec[2], 1 , param$tau.vec[3], 1 ,
             lwd = lwd, col = col)
  }
  if(0 %in% treatment) {
    points(param$coef.vec.0, y.axis - adj, type = "p", pch = 1, cex = cex, col = col)
    segments(param$lower.vec.0, y.axis - adj, param$upper.vec.0, y.axis - adj,
             lwd = lwd, lty = 3, col = col)
  }
  y.axis.new <- c(3,2,1)
  axis(2, at = y.axis.new, labels = labels, las = 1, tick = TRUE, ...)
  abline(v = 0, lty = 2)
  
  # 2. Group effects    
  if(group.plots){
    
    #########################################################################     
    if(is.factor(x$group.id)){
      labels <- levels(x$group.id)
    } else {
      labels = sort(unique(x$group.id))
    }
    #########################################################################
    if(0 %in% treatment){
      y.axis <- c(length(param$coef.vec.0.group[,1]):.5)
      y.axis <- y.axis + 1
      
      xlim <- c(param$lower.vec.0.group[,1], param$coef.vec.0.group[,1], param$upper.vec.0.group[,1]) 
      MIN <- ifelse(min(xlim)>0, min(xlim)*0.9, min(xlim)*1.1)
      MAX <- ifelse(max(xlim)>0, max(xlim)*1.1, max(xlim)*0.9)
      xlim <- c(MIN, MAX)
      
      ylim <- c(min(y.axis), max(y.axis))
      
      plot(param$coef.vec.0.group[,1], y.axis, type = "n", xlab = xlab, ylab = ylab,
           yaxt = "n", xlim = xlim, ylim = ylim, main = "ACME (control)", ...)
      
      points(param$coef.vec.0.group[,1], y.axis, type = "p", pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0.group[,1], y.axis, param$upper.vec.0.group[,1], y.axis,
               lwd = lwd, col = col)
      
      axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE, ...)
      abline(v = 0, lty = 2)
    }
    ######################################################################### 
    if(1 %in% treatment){
      y.axis <- c(length(param$coef.vec.1.group[,1]):.5)
      y.axis <- y.axis + 1
      
      xlim <- c(param$lower.vec.1.group[,1], param$coef.vec.1.group[,1], param$upper.vec.1.group[,1])
      MIN <- ifelse(min(xlim)>0, min(xlim)*0.9, min(xlim)*1.1)
      MAX <- ifelse(max(xlim)>0, max(xlim)*1.1, max(xlim)*0.9)
      xlim <- c(MIN, MAX)
      
      ylim <- c(min(y.axis), max(y.axis))
      
      plot(param$coef.vec.1.group[,1], y.axis, type = "n", xlab = xlab, ylab = ylab,
           yaxt = "n", xlim = xlim, ylim = ylim, main = ifelse(printone, "ACME", "ACME (treated)"), ...)
      
      points(param$coef.vec.1.group[,1], y.axis, type = "p", pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1.group[,1], y.axis, param$upper.vec.1.group[,1], y.axis,
               lwd = lwd, col = col)
      
      axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE, ...)
      abline(v = 0, lty = 2)
    }
    ######################################################################### 
    if(0 %in% treatment){
      y.axis <- c(length(param$coef.vec.0.group[,2]):.5)
      y.axis <- y.axis + 1
      
      xlim <- c(param$lower.vec.0.group[,2], param$coef.vec.0.group[,2], param$upper.vec.0.group[,2])
      MIN <- ifelse(min(xlim)>0, min(xlim)*0.9, min(xlim)*1.1)
      MAX <- ifelse(max(xlim)>0, max(xlim)*1.1, max(xlim)*0.9)
      xlim <- c(MIN, MAX)
      
      ylim <- c(min(y.axis), max(y.axis))
      
      plot(param$coef.vec.0.group[,2], y.axis, type = "n", xlab = xlab, ylab = ylab,
           yaxt = "n", xlim = xlim, ylim = ylim, main = "ADE (control)", ...)
      
      points(param$coef.vec.0.group[,2], y.axis, type = "p", pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0.group[,2], y.axis, param$upper.vec.0.group[,2], y.axis,
               lwd = lwd, lty = 3, col = col)
      
      axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE, ...)
      abline(v = 0, lty = 2)
    }    
    ######################################################################### 
    if(1 %in% treatment){
      y.axis <- c(length(param$coef.vec.1.group[,2]):.5)
      y.axis <- y.axis + 1
      
      xlim <- c(param$lower.vec.1.group[,2], param$coef.vec.1.group[,2], param$upper.vec.1.group[,2])
      MIN <- ifelse(min(xlim)>0, min(xlim)*0.9, min(xlim)*1.1)
      MAX <- ifelse(max(xlim)>0, max(xlim)*1.1, max(xlim)*0.9)
      xlim <- c(MIN, MAX)
      
      ylim <- c(min(y.axis), max(y.axis))
      
      plot(param$coef.vec.1.group[,2], y.axis, type = "n", xlab = xlab, ylab = ylab,
           yaxt = "n", xlim = xlim, ylim = ylim, main = ifelse(printone, "ADE", "ADE (treated)"), ...)
      
      points(param$coef.vec.1.group[,2], y.axis, type = "p", pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1.group[,2], y.axis, param$upper.vec.1.group[,2], y.axis,
               lwd = lwd, lty = 3, col = col)
      
      axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE, ...)
      abline(v = 0, lty = 2)
    }    
    ######################################################################### 
    y.axis <- c(length(param$tau.vec.group[,1]):.5)
    y.axis <- y.axis + 1
    
    xlim <- c(param$tau.vec.group[,1], param$tau.vec.group[,2], param$tau.vec.group[,3])
    MIN <- ifelse(min(xlim)>0, min(xlim)*0.9, min(xlim)*1.1)
    MAX <- ifelse(max(xlim)>0, max(xlim)*1.1, max(xlim)*0.9)
    xlim <- c(MIN, MAX)
    
    ylim <- c(min(y.axis), max(y.axis))
    
    plot(param$tau.vec.group[,1], y.axis, type = "n", xlab = xlab, ylab = ylab,
         yaxt = "n", xlim = xlim, ylim = ylim, main = "Total\nEffect", ...)
    
    points(param$tau.vec.group[,1], y.axis , type = "p", pch = 19, cex = cex, col = col)
    segments(param$tau.vec.group[,2], y.axis , param$tau.vec.group[,3], y.axis ,
             lwd = lwd, col = col)
    
    axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE, ...)
    abline(v = 0, lty = 2)
    
  }
}

#########################################################################
plot.process.order <- function(model){
  length <- length(model$d1)
  coef.vec.1 <- lower.vec.1 <- upper.vec.1 <-
    coef.vec.0 <- lower.vec.0 <- upper.vec.0 <- matrix(NA,ncol=2,nrow=length)
  tau.vec<-matrix(NA,ncol=3,nrow=length)
  for(j in 1:length){
    coef.vec.1[j,] <- c(model$d1[j], model$z1[j])
    lower.vec.1[j,] <- c(model$d1.ci[1,j], model$z1.ci[1,j])
    upper.vec.1[j,] <- c(model$d1.ci[2,j], model$z1.ci[2,j])
    
    coef.vec.0[j,] <- c(model$d0[j], model$z0[j])
    lower.vec.0[j,] <- c(model$d0.ci[1,j], model$z0.ci[1,j])
    upper.vec.0[j,] <- c(model$d0.ci[2,j], model$z0.ci[2,j])
    
    tau.vec[j,] <- c(model$tau.coef[j], model$tau.ci[1,j], model$tau.ci[2,j])
    
  }
  
  range.1 <- range(model$d1.ci[1,], model$z1.ci[1,],model$tau.ci[1,],
                   model$d1.ci[2,], model$z1.ci[2,],model$tau.ci[2,])
  range.0 <- range(model$d0.ci[1,], model$z0.ci[1,],model$tau.ci[1,],
                   model$d0.ci[2,], model$z0.ci[2,],model$tau.ci[2,])
  
  return(list(coef.vec.1=coef.vec.1, lower.vec.1=lower.vec.1,
              upper.vec.1=upper.vec.1, coef.vec.0=coef.vec.0,
              lower.vec.0=lower.vec.0, upper.vec.0=upper.vec.0,
              tau.vec=tau.vec,
              range.1=range.1, range.0=range.0, length=length))
}

#########################################################################
plot.mediate.order <- function(x, treatment = NULL,
                               labels = c("ACME","ADE","Total\nEffect"),
                               xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                               main = NULL, lwd = 1.5, cex = .85,
                               col = "black", ...){
  # Determine which graph to plot
  if(is.null(treatment)){
    if(x$INT){
      treatment <- c(0,1)
    } else {
      treatment <- 1
    }
  } else {
    treatment <- switch(treatment,
                        control = 0,
                        treated = 1,
                        both = c(0,1))
  }
  
  param <- plot.process.order(x)
  y.axis <- c(ncol(param$coef.vec.1):.5)
  y.axis <- y.axis + 1
  # create indicator for y.axis, descending so labels go from top to bottom
  
  # Set xlim
  if(is.null(xlim)){
    if(length(treatment) > 1) {
      xlim <- range(param$range.1, param$range.0) * 1.2
    } else if (treatment == 1){
      xlim <- param$range.1 * 1.2
    } else {
      xlim <- param$range.0 * 1.2
    }
  }
  
  # Set ylim
  if(is.null(ylim)){
    ylim <- c(min(y.axis) - 1 - 0.5, max(y.axis) + 0.5)
  }
  
  # Plot
  plot(param$coef.vec.1[1,], y.axis, type = "n", xlab = xlab, ylab = ylab,
       yaxt = "n", xlim = xlim, ylim = ylim, main = main, ...)
  
  # Set offset values depending on number of bars to plot
  if(length(treatment) == 1){
    adj <- 0
  } else {
    adj <- 0.05
  }
  
  if(1 %in% treatment){
    adj.1 <- adj * nrow(param$coef.vec.1)
    for(z in 1:nrow(param$coef.vec.1)){
      points(param$coef.vec.1[z,], y.axis + adj.1,
             type = "p", pch = 19, cex = cex, col = col)
      segments(param$lower.vec.1[z,], y.axis + adj.1,
               param$upper.vec.1[z,], y.axis + adj.1,
               lwd = lwd, col = col)
      points(param$tau.vec[z,1], 1 + adj.1 ,
             type = "p", pch = 19, cex = cex, col = col)
      segments(param$tau.vec[z,2], 1 + adj.1 ,
               param$tau.vec[z,3], 1 + adj.1 ,
               lwd = lwd, col = col)
      adj.1 <- adj.1 - 0.05
    }
    
  }
  if(0 %in% treatment) {
    adj.0 <- adj
    for(z in 1:nrow(param$coef.vec.0)){
      points(param$coef.vec.0[z,], y.axis - adj.0,
             type = "p", pch = 1, cex = cex, col = col)
      segments(param$lower.vec.0[z,], y.axis - adj.0,
               param$upper.vec.0[z,], y.axis - adj.0,
               lwd = lwd, lty = 3, col = col)
      adj.0 <- adj.0 + 0.05
    }
  }
  if (treatment[1]==0 & length(treatment)==1){
    print("test")
    adj.1 <- adj * nrow(param$coef.vec.1)
    for(z in 1:nrow(param$tau.vec)){
      points(param$tau.vec[z,1], 1 + adj.1 ,
             type = "p", pch = 19, cex = cex, col = col)
      segments(param$tau.vec[z,2], 1 + adj.1 ,
               param$tau.vec[z,3], 1 +adj.1 ,
               lwd = lwd, col = col)
      adj.1 <- adj.1 - 0.05
    }
  }
  
  y.axis.new <- c(3,2,1)
  axis(2, at = y.axis.new, labels = labels, las = 1, tick = TRUE, ...)
  abline(v = 0, lty = 2)
}

pval <- function(x, xhat){
  ## Compute p-values
  if (xhat == 0) out <- 1
  else {
    out <- 2 * min(sum(x > 0), sum(x < 0)) / length(x)
  }
  return(min(out, 1))
}
