med.boot <- function(y.data, index, m.data) {
          
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
  d1 <- weighted.mean(effects.tmp[,1], weights)
  d0 <- weighted.mean(effects.tmp[,2], weights)
  z1 <- weighted.mean(effects.tmp[,3], weights)
  z0 <- weighted.mean(effects.tmp[,4], weights)
  
  c(d1 = d1, d0 = d0, z1 = z1, z0 = z0)
}