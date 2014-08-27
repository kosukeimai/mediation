ivmediate <- function(model.t, model.m, model.y, ci = TRUE, sims = 1000, boot = TRUE,
                      enc = "enc.name", treat = "treat.name", mediator = "med.name",
                      conf.level = .95, long = TRUE, dropobs = FALSE,
                      multicore = FALSE, mc.cores = getOption("mc.cores", 2L)){
  
  # Drop observations not common to all three models
  if(dropobs){
    odata.t <- model.frame(model.t)
    odata.m <- model.frame(model.m)
    odata.y <- model.frame(model.y)
    odata.tm <- merge(odata.t, odata.m, sort=FALSE,
                      by=c("row.names", intersect(names(odata.t), names(odata.m))))
    rownames(odata.tm) <- odata.tm$Row.names
    newdata <- merge(odata.tm[,-1L], odata.y, sort=FALSE,
                     by=c("row.names", intersect(names(odata.tm), names(odata.y))))
    rownames(newdata) <- newdata$Row.names
    newdata <- newdata[,-1L]
    rm(odata.t, odata.m, odata.y, odata.tm)
    
    call.t <- getCall(model.t)
    call.m <- getCall(model.m)
    call.y <- getCall(model.y)
    
    call.t$data <- call.m$data <- call.y$data <- newdata
    if(c("(weights)") %in% names(newdata)){
      call.t$weights <- call.m$weights <- call.y$weights <- model.weights(newdata)
    }
    model.t <- eval.parent(call.t)
    model.m <- eval.parent(call.m)
    model.y <- eval.parent(call.y)
  }
  n <- nrow(model.frame(model.y))
  
  # Get point estimates
  out.point <- ivmediate.fit(NULL, model.t, model.m, model.y, enc, treat, mediator)
  
  
  # Get confidence intervals
  
  dc1.ci <- dc0.ci <- zc1.ci <- zc0.ci <- tauc.ci <- 
    dc1.sims <- dc0.sims <- zc1.sims <- zc0.sims <- tauc.sims <- NULL
  
  if(ci){
    if(boot){
      if(multicore){
        out.sim.full <- parallel::mclapply(1:sims, ivmediate.fit.b, 
                                 model.t = model.t, model.m = model.m, model.y = model.y,
                                 enc = enc, treat = treat, mediator = mediator,
                                 mc.cores = mc.cores)
      } else {
        out.sim.full <- lapply(1:sims, ivmediate.fit.b, 
                               model.t = model.t, model.m = model.m, model.y = model.y,
                               enc = enc, treat = treat, mediator = mediator)
      }
    } else {
      if(multicore){
        out.sim.full <- parallel::mclapply(1:sims, ivmediate.fit, 
                                 model.t = model.t, model.m = model.m, model.y = model.y,
                                 enc = enc, treat = treat, mediator = mediator, sim = TRUE,
                                 mc.cores = mc.cores)
      } else {
        out.sim.full <- lapply(1:sims, ivmediate.fit,
                               model.t = model.t, model.m = model.m, model.y = model.y,
                               enc = enc, treat = treat, mediator = mediator, sim = TRUE)
      }
    }
    
    out.sim <- lapply(out.sim.full, function(x) if (is.character(x)) rep(NA, length(x)) else x)
    out.sim <- simplify2array(out.sim)
    
    out.ci <- array(dim=c(2,nrow(out.sim),length(conf.level)))
    for(i in 1:length(conf.level)){
      prob <- (1-conf.level[i])/2
      out.ci[,,i] <- apply(out.sim, 1, quantile, probs = c(prob, 1-prob), na.rm=TRUE)
    }
    dimnames(out.ci) <- list(c("lower","upper"), NULL,
                             paste(100*conf.level, rep("%",length(conf.level)), sep=""))
    dc1.ci <- as.matrix(out.ci[,1,])
    dc0.ci <- as.matrix(out.ci[,2,])
    zc1.ci <- as.matrix(out.ci[,3,])
    zc0.ci <- as.matrix(out.ci[,4,])
    tauc.ci <- as.matrix(out.ci[,5,])
    
    if(long){
      dc1.sims <- out.sim[1,]
      dc0.sims <- out.sim[2,]
      zc1.sims <- out.sim[3,]
      zc0.sims <- out.sim[4,]
      tauc.sims <- out.sim[5,]
    }
  }  
  
  # Wrap things up
  
  out <- list(dc1 = out.point[1], dc0 = out.point[2], 
              zc1 = out.point[3], zc0 = out.point[4],
              tauc = out.point[5],
              dc1.issue = out.point[6], dc0.issue = out.point[7], 
              zc1.issue = out.point[8], zc0.issue = out.point[9],
              dc1.inf = out.point[10], dc0.inf = out.point[11],
              zc1.inf = out.point[12], zc0.inf = out.point[13],
              dc1.ci = dc1.ci, dc0.ci = dc0.ci, 
              zc1.ci = zc1.ci, zc0.ci = zc0.ci,
              tauc.ci = tauc.ci,
              dc1.sims = dc1.sims, dc0.sims = dc0.sims, 
              zc1.sims = zc1.sims, zc0.sims = zc0.sims,
              tauc.sims = tauc.sims,
              boot = boot, enc = enc, treat = treat, mediator = mediator,
              conf.level = conf.level,
              nobs = n, sims = sims
  )
  class(out) <- "ivmediate"
  out
}


ivmediate.fit.b <- function(x, model.t, model.m, model.y, 
                            enc, treat, mediator){
  
  data <- model.frame(model.y)    
  sb <- sample(1:nrow(data), nrow(data), replace=TRUE)
  Data.b <- data[sb,]
  
  model.t.b <- update(model.t, data = Data.b)
  model.m.b <- update(model.m, data = Data.b)
  model.y.b <- update(model.y, data = Data.b)
  
  out <- ivmediate.fit(NULL, model.t.b, model.m.b, model.y.b, enc, treat, mediator)
  return(out)  
}


ivmediate.fit <- function(x, model.t, model.m, model.y, 
                          enc, treat, mediator, sim = FALSE){
  
  # Model type indicators
  isGlm.t <- inherits(model.t, "glm")  # Note gam and bayesglm also inherits "glm"
  isGlm.m <- inherits(model.m, "glm")  # Note gam and bayesglm also inherits "glm"
  isGlm.y <- inherits(model.y, "glm")  # Note gam and bayesglm also inherits "glm"
  isLm.t <- inherits(model.t, "lm")    # Note gam, glm and bayesglm also inherit "lm"
  isLm.m <- inherits(model.m, "lm")    # Note gam, glm and bayesglm also inherit "lm"
  isLm.y <- inherits(model.y, "lm")    # Note gam, glm and bayesglm also inherit "lm"
  #  isVglm.y <- inherits(model.y, "vglm")
  #  isRq.y <- inherits(model.y, "rq")
  #  isRq.m <- inherits(model.m, "rq")
  #  isOrdered.y <- inherits(model.y, "polr")  # Note bayespolr also inherits "polr"
  #  isOrdered.m <- inherits(model.m, "polr")  # Note bayespolr also inherits "polr"
  #  isSurvreg.y <- inherits(model.y, "survreg")
  #  isSurvreg.m <- inherits(model.m, "survreg")
  #  isMer.y <- inherits(model.y, "merMod") # Note lmer and glmer do not inherit "lm" and "glm"
  #  isMer.m <- inherits(model.m, "merMod") # Note lmer and glmer do not inherit "lm" and "glm"
  
  # Extract or simulate coefs
  name.coef.t <- names(coef(model.t))
  name.coef.m <- names(coef(model.m))	
  name.coef.y <- names(coef(model.y))	
  coef.t <- as.matrix(coef(model.t))
  coef.m <- as.matrix(coef(model.m))
  coef.y <- as.matrix(coef(model.y))
  if(sim){
    vcov.t <- as.matrix(vcov(model.t))
    vcov.m <- as.matrix(vcov(model.m))
    vcov.y <- as.matrix(vcov(model.y))
    coef.t <- t(as.matrix(rmvnorm(1, coef.t, vcov.t)))
    coef.m <- t(as.matrix(rmvnorm(1, coef.m, vcov.m)))
    coef.y <- t(as.matrix(rmvnorm(1, coef.y, vcov.y)))
  }
  
  mf.T.1 <- mf.T.0 <- mf.T <- model.frame(model.t)
  mf.T.1[,enc] <- 1
  mf.T.0[,enc] <- 0
  
  Q11 <- tcrossprod(model.matrix(terms(model.t), mf.T.1), t(coef.t))
  Q01 <- 1 - Q11
  Q10 <- tcrossprod(model.matrix(terms(model.t), mf.T.0), t(coef.t))
  Q00 <- 1 - Q10
  
  if(isGlm.t){
    Q11 <- model.t$family$linkinv(Q11)
    Q10 <- model.t$family$linkinv(Q10)
    Q01 <- model.t$family$linkinv(Q01)
    Q00 <- model.t$family$linkinv(Q00)
  }
  
  mf.M.11 <- mf.M.10 <- mf.M.01 <- mf.M.00 <- mf.M <- model.frame(model.m)
  mf.M.11[,treat] <- mf.M.11[,enc] <- mf.M.10[,treat] <- mf.M.01[,enc] <- 1
  mf.M.10[,enc] <- mf.M.01[,treat] <- mf.M.00[,treat] <- mf.M.00[,enc] <- 0
  mu.M.11 <- tcrossprod(model.matrix(terms(model.m), mf.M.11), t(coef.m))
  mu.M.10 <- tcrossprod(model.matrix(terms(model.m), mf.M.10), t(coef.m))
  mu.M.01 <- tcrossprod(model.matrix(terms(model.m), mf.M.01), t(coef.m))
  mu.M.00 <- tcrossprod(model.matrix(terms(model.m), mf.M.00), t(coef.m))
  
  if(isGlm.m){
    if(family(model.m)$family == "binomial"){
      g11 <- function(m,i) (1-m) + family(model.m)$linkinv(mu.M.11[i])*(-1)^(1-m)
      g10 <- function(m,i) (1-m) + family(model.m)$linkinv(mu.M.10[i])*(-1)^(1-m)
      g01 <- function(m,i) (1-m) + family(model.m)$linkinv(mu.M.01[i])*(-1)^(1-m)
      g00 <- function(m,i) (1-m) + family(model.m)$linkinv(mu.M.00[i])*(-1)^(1-m)
    } else stop("mediator model is not of supported GLM family")
  } else if(isLm.m){
    g11 <- function(m,i) dnorm(m, mu.M.11[i], summary(model.m)$sigma)  # density of M given T=Z=1
    g10 <- function(m,i) dnorm(m, mu.M.10[i], summary(model.m)$sigma)  # density of M given T=1,Z=0
    g01 <- function(m,i) dnorm(m, mu.M.01[i], summary(model.m)$sigma)  # density of M given T=0,Z=1
    g00 <- function(m,i) dnorm(m, mu.M.00[i], summary(model.m)$sigma)  # density of M given T=Z=0
  } else stop("mediator model not supported")
  
  mf.Y.11 <- mf.Y.10 <- mf.Y.01 <- mf.Y.00 <- mf.Y <- model.frame(model.y)
  n <- nrow(mf.Y)
  mf.Y.11[,treat] <- mf.Y.11[,enc] <- mf.Y.10[,treat] <- mf.Y.01[,enc] <- 1
  mf.Y.10[,enc] <- mf.Y.01[,treat] <- mf.Y.00[,treat] <- mf.Y.00[,enc] <- 0
  
  labs.y <- attr(terms(model.y), "term.labels")
  labs.y.ind <- sapply(strsplit(labs.y, ":"), function(x) mediator %in% x)
  labs.y.a <- labs.y[!labs.y.ind]
  labs.y.b <- sapply(strsplit(labs.y[labs.y.ind], ":"), function(x) paste(x[-match(mediator,x)],collapse=":"))
  form.y.a <- formula(c("~", paste(labs.y.a, collapse="+")))
  form.y.b <- formula(c("~", paste(labs.y.b, collapse="+")))  
  mm.Y.11.a <- model.matrix(form.y.a, mf.Y.11)
  mm.Y.11.b <- model.matrix(form.y.b, mf.Y.11)
  mm.Y.10.a <- model.matrix(form.y.a, mf.Y.10)
  mm.Y.10.b <- model.matrix(form.y.b, mf.Y.10)
  mm.Y.01.a <- model.matrix(form.y.a, mf.Y.01)
  mm.Y.01.b <- model.matrix(form.y.b, mf.Y.01)
  mm.Y.00.a <- model.matrix(form.y.a, mf.Y.00)
  mm.Y.00.b <- model.matrix(form.y.b, mf.Y.00)

  mu.Y.ind <- sapply(strsplit(name.coef.y , ":"), function(x) mediator %in% x)
  mu.Y.11.a <- tcrossprod(mm.Y.11.a, t(coef.y[!mu.Y.ind]))
  mu.Y.11.b <- tcrossprod(mm.Y.11.b, t(coef.y[mu.Y.ind]))
  mu.Y.10.a <- tcrossprod(mm.Y.10.a, t(coef.y[!mu.Y.ind]))
  mu.Y.10.b <- tcrossprod(mm.Y.10.b, t(coef.y[mu.Y.ind]))
  mu.Y.01.a <- tcrossprod(mm.Y.01.a, t(coef.y[!mu.Y.ind]))
  mu.Y.01.b <- tcrossprod(mm.Y.01.b, t(coef.y[mu.Y.ind]))
  mu.Y.00.a <- tcrossprod(mm.Y.00.a, t(coef.y[!mu.Y.ind]))
  mu.Y.00.b <- tcrossprod(mm.Y.00.b, t(coef.y[mu.Y.ind]))
  
  if(isGlm.y){
    if(family(model.y)$family == "binomial"){
      S11 <- function(m,i) family(model.y)$linkinv(mu.Y.11.a[i] + m * mu.Y.11.b[i])  # mean of Y given M=m, T=1 and Z=1
      S10 <- function(m,i) family(model.y)$linkinv(mu.Y.10.a[i] + m * mu.Y.10.b[i])  # mean of Y given M=m, T=1 and Z=0
      S01 <- function(m,i) family(model.y)$linkinv(mu.Y.01.a[i] + m * mu.Y.01.b[i])  # mean of Y given M=m, T=0 and Z=1
      S00 <- function(m,i) family(model.y)$linkinv(mu.Y.00.a[i] + m * mu.Y.00.b[i])  # mean of Y given M=m, T=0 and Z=0
    } else stop("outcome model is not of supported GLM family")
  } else if (isLm.y){
    S11 <- function(m,i) mu.Y.11.a[i] + m * mu.Y.11.b[i]  # mean of Y given M=m, T=1 and Z=1
    S10 <- function(m,i) mu.Y.10.a[i] + m * mu.Y.10.b[i]  # mean of Y given M=m, T=1 and Z=0
    S01 <- function(m,i) mu.Y.01.a[i] + m * mu.Y.01.b[i]  # mean of Y given M=m, T=0 and Z=1
    S00 <- function(m,i) mu.Y.00.a[i] + m * mu.Y.00.b[i]  # mean of Y given M=m, T=0 and Z=0
  } else stop("outcome model not supported")
  
  LACME1.integrand <- function(m,i){
    val <- ((Q11[i]*g11(m,i)*S11(m,i) - Q10[i]*g10(m,i)*S10(m,i))/(Q11[i]*g11(m,i) - Q10[i]*g10(m,i))) *
      (Q11[i]*g11(m,i) - Q10[i]*g10(m,i) - Q00[i]*g00(m,i) + Q01[i]*g01(m,i))/(Q11[i] - Q10[i])
    val <- ifelse(is.nan(val), 0, val)
    val <- ifelse(val == Inf, .Machine$double.xmax, val)
    ifelse(val == -Inf, -.Machine$double.xmax, val)
  }
  
  LACME0.integrand <- function(m,i){
    val <- ((Q00[i]*g00(m,i)*S00(m,i) - Q01[i]*g01(m,i)*S01(m,i))/(Q00[i]*g00(m,i) - Q01[i]*g01(m,i))) *
      (Q11[i]*g11(m,i) - Q10[i]*g10(m,i) - Q00[i]*g00(m,i) + Q01[i]*g01(m,i))/(Q11[i] - Q10[i])
    val <- ifelse(is.nan(val), 0, val)
    val <- ifelse(val == Inf, .Machine$double.xmax, val)
    ifelse(val == -Inf, -.Machine$double.xmax, val)
  }
  
  LANDE1.integrand <- function(m,i){
    val <- ((Q11[i]*g11(m,i)*S11(m,i) - Q10[i]*g10(m,i)*S10(m,i))/(Q11[i]*g11(m,i) - Q10[i]*g10(m,i)) -
              (Q00[i]*g00(m,i)*S00(m,i) - Q01[i]*g01(m,i)*S01(m,i))/(Q00[i]*g00(m,i) - Q01[i]*g01(m,i))) *
      (Q11[i]*g11(m,i) - Q10[i]*g10(m,i))/(Q11[i] - Q10[i])
    val <- ifelse(is.nan(val), 0, val)
    val <- ifelse(val == Inf, .Machine$double.xmax, val)
    ifelse(val == -Inf, -.Machine$double.xmax, val)
  }
  
  LANDE0.integrand <- function(m,i){
    val <- ((Q11[i]*g11(m,i)*S11(m,i) - Q10[i]*g10(m,i)*S10(m,i))/(Q11[i]*g11(m,i) - Q10[i]*g10(m,i)) -
              (Q00[i]*g00(m,i)*S00(m,i) - Q01[i]*g01(m,i)*S01(m,i))/(Q00[i]*g00(m,i) - Q01[i]*g01(m,i))) *
      (Q00[i]*g00(m,i) - Q01[i]*g01(m,i))/(Q11[i] - Q10[i])
    val <- ifelse(is.nan(val), 0, val)
    val <- ifelse(val == Inf, .Machine$double.xmax, val)
    ifelse(val == -Inf, -.Machine$double.xmax, val)
  }
  
  LACME1.issue <- LACME0.issue <- LANDE1.issue <- LANDE0.issue <- 0
  LACME1.inf <- LACME0.inf <- LANDE1.inf <- LANDE0.inf <- 0
  
  if(isGlm.m){
    if(family(model.m)$family == "binomial"){
      
      LACME1 <- mean(LACME1.integrand(1,1:n) + LACME1.integrand(0,1:n))
      LACME0 <- mean(LACME0.integrand(1,1:n) + LACME0.integrand(0,1:n))
      LANDE1 <- mean(LANDE1.integrand(1,1:n) + LANDE1.integrand(0,1:n))
      LANDE0 <- mean(LANDE0.integrand(1,1:n) + LANDE0.integrand(0,1:n))
      
    } else stop("mediator model is not of supported GLM family")
  } else if(isLm.m){
    
    varnames <- c(names(mf.T),names(mf.M),names(mf.Y)[-1L])
    ind.mtz <- varnames %in% c(mediator,treat,enc)
    
    if(!sum(!ind.mtz)){  # no covariates
      
      LACME1 <- integrate(LACME1.integrand, -Inf, Inf, i=1, stop.on.error=FALSE)$value
      LACME0 <- integrate(LACME0.integrand, -Inf, Inf, i=1, stop.on.error=FALSE)$value
      LANDE1 <- integrate(LANDE1.integrand, -Inf, Inf, i=1, stop.on.error=FALSE)$value
      LANDE0 <- integrate(LANDE0.integrand, -Inf, Inf, i=1, stop.on.error=FALSE)$value
      
    } else {  # with covariates
      
      LACME1.hat <- LACME0.hat <- LANDE1.hat <- LANDE0.hat <-
        LACME1.hat.issue <- LACME0.hat.issue <- LANDE1.hat.issue <- LANDE0.hat.issue <- rep(NA, n)
      for(i in 1:n){
        
        LACME1.out <- integrate(LACME1.integrand, -Inf, Inf, i=i, stop.on.error=F)
        LACME0.out <- integrate(LACME0.integrand, -Inf, Inf, i=i, stop.on.error=F)
        LANDE1.out <- integrate(LANDE1.integrand, -Inf, Inf, i=i, stop.on.error=F)
        LANDE0.out <- integrate(LANDE0.integrand, -Inf, Inf, i=i, stop.on.error=F)
        
        LACME1.hat[i] <- LACME1.out$value
        LACME0.hat[i] <- LACME0.out$value
        LANDE1.hat[i] <- LANDE1.out$value
        LANDE0.hat[i] <- LANDE0.out$value
        
        LACME1.hat.issue[i] <- LACME1.out$message
        LACME0.hat.issue[i] <- LACME0.out$message
        LANDE1.hat.issue[i] <- LANDE1.out$message
        LANDE0.hat.issue[i] <- LANDE0.out$message
        
        #        cat(i, "\n")
      }
      
      LACME1.finite <- is.finite(LACME1.hat)
      LACME0.finite <- is.finite(LACME0.hat)
      LANDE1.finite <- is.finite(LANDE1.hat)
      LANDE0.finite <- is.finite(LANDE0.hat)
      
      LACME1 <- sum(LACME1.hat[LACME1.finite])/sum(LACME1.finite)
      LACME0 <- sum(LACME0.hat[LACME0.finite])/sum(LACME0.finite)
      LANDE1 <- sum(LANDE1.hat[LANDE1.finite])/sum(LANDE1.finite)
      LANDE0 <- sum(LANDE0.hat[LANDE0.finite])/sum(LANDE0.finite)
      
      LACME1.issue <- sum(LACME1.hat.issue != "OK")
      LACME0.issue <- sum(LACME0.hat.issue != "OK")
      LANDE1.issue <- sum(LANDE1.hat.issue != "OK")
      LANDE0.issue <- sum(LANDE0.hat.issue != "OK")
      
      LACME1.inf <- sum(!LACME1.finite)
      LACME0.inf <- sum(!LACME0.finite)
      LANDE1.inf <- sum(!LANDE1.finite)
      LANDE0.inf <- sum(!LANDE0.finite)
      
    }
  } else stop("mediator model not supported")
  
  LATE <- LACME1 + LANDE0
  
  out <- c(LACME1, LACME0, LANDE1, LANDE0, LATE,
           LACME1.issue, LACME0.issue, LANDE1.issue, LANDE0.issue,
           LACME1.inf, LACME0.inf, LANDE1.inf, LANDE0.inf)
  names(out) <- c("dc1", "dc0", "zc1", "zc0", "tauc",
                  "dc1.issue", "dc0.issue", "zc1.issue", "zc0.issue",
                  "dc1.inf", "dc0.inf", "zc1.inf", "zc0.inf")
  return(out)
  
}


summary.ivmediate <- function(object, conf.level = object$conf.level[1], ...){
  if(!is.null(object$dc1.ci)){
      cl.ind <- match(conf.level, object$conf.level)
      if(is.na(cl.ind)){
        stop("conf.level must be one of the levels contained in the ivmediate output")
      }
      object$conf.level <- object$conf.level[cl.ind]
      object$dc1.ci <- object$dc1.ci[,cl.ind]
      object$dc0.ci <- object$dc0.ci[,cl.ind]
      object$zc1.ci <- object$zc1.ci[,cl.ind]
      object$zc0.ci <- object$zc0.ci[,cl.ind]
      object$tauc.ci <- object$tauc.ci[,cl.ind]
  }
  structure(object, class = c("summary.ivmediate", class(object)))
}


print.summary.ivmediate <- function(x, ...){
  cat("\n")
  cat("Causal Mediation Analysis with Treatment Noncompliance\n\n")
  if(is.null(x$dc1.ci)){
    smat <- as.matrix(c(x$dc0, x$dc1, x$zc0, x$zc1, x$tauc))
    rownames(smat) <- c("LACME (control)", "LACME (treated)",
                        "LANDE (control)", "LANDE (treated)",
                        "LATE")
    colnames(smat) <- "Estimate"
  } else {
    clp <- 100 * x$conf.level
    if(x$boot){
        cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
    } else {
        cat("Confidence Intervals Based on Quasi-Bayesian Monte Carlo\n\n")
    }
    smat <- c(x$dc0, x$dc0.ci)
    smat <- rbind(smat, c(x$dc1, x$dc1.ci))
    smat <- rbind(smat, c(x$zc0, x$zc0.ci))
    smat <- rbind(smat, c(x$zc1, x$zc1.ci))
    smat <- rbind(smat, c(x$tauc, x$tauc.ci))
    rownames(smat) <- c("LACME (control)", "LACME (treated)",
                        "LANDE (control)", "LANDE (treated)",
                        "LATE")
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                      paste(clp, "% CI Upper", sep=""))
  }
  
  printCoefmat(smat)
  cat("\n")
  cat("Sample Size Used:", x$nobs,"\n\n")
  
  if(!is.null(x$dc1.ci)){
      cat("\n")
      cat("Simulations:", x$sims,"\n\n")
  }
  
  invisible(x)
}


plot.ivmediate <- function(x, treatment = NULL, labels = NULL,
                         effect.type = c("indirect","direct","total"),
                         conf.level = x$conf.level[1],
                         xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                         main = NULL, lwd = 1.5, cex = .85,
                         col = "black", ...){
  if(is.null(x$dc1.ci)){
    stop("plot cannot be used for an ivmediate object without confidence intervals")
  }

  cl.ind <- match(conf.level, x$conf.level)
  if(is.na(cl.ind)){
    stop("conf.level must be one of the levels contained in the ivmediate output")
  }
  
  effect.type <- match.arg(effect.type, several.ok=TRUE)
  IND <- "indirect" %in% effect.type
  DIR <- "direct" %in% effect.type
  TOT <- "total" %in% effect.type
  
  if(is.null(treatment)){
    treatment <- c(0,1)
  } else {
    treatment <- switch(treatment,
                        control = 0,
                        treated = 1,
                        both = c(0,1))
  }
  
  range.1 <- range(x$dc1.ci, x$zc1.ci, x$tauc.ci)
  range.0 <- range(x$dc0.ci, x$zc0.ci, x$tauc.ci)

  y.axis <- (IND + DIR + TOT):1
  
  # Set xlim
  if(is.null(xlim)){
    if(length(treatment) > 1) {
      xlim <- range(range.1, range.0) * 1.2
    } else if (treatment == 1){
      xlim <- range.1 * 1.2
    } else {
      xlim <- range.0 * 1.2
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
      points(c(x$dc1,x$zc1), y.axis[1:2] + adj, type = "p", pch = 19, cex = cex, col = col)
      segments(c(x$dc1.ci[1,cl.ind],x$zc1.ci[1,cl.ind]), y.axis[1:2] + adj, 
                c(x$dc1.ci[2,cl.ind],x$zc1.ci[2,cl.ind]), y.axis[1:2] + adj,
                lwd = lwd, col = col)
    }
    if(IND && !DIR) {
      points(x$dc1, y.axis[1] + adj, type = "p", pch = 19, cex = cex, col = col)
      segments(x$dc1.ci[1,cl.ind], y.axis[1] + adj, 
                x$dc1.ci[2,cl.ind], y.axis[1] + adj,
                lwd = lwd, col = col)
    }
    if(!IND && DIR) {
      points(x$zc1, y.axis[1] + adj, type = "p", pch = 19, cex = cex, col = col)
      segments(x$zc1.ci[1,cl.ind], y.axis[1] + adj, 
                x$zc1.ci[2,cl.ind], y.axis[1] + adj,
                lwd = lwd, col = col)
    }    
  }
  
  if(0 %in% treatment) {
    if(IND && DIR) {
      points(c(x$dc0,x$zc0), y.axis[1:2] - adj, type = "p", pch = 19, cex = cex, col = col)
      segments(c(x$dc0.ci[1,cl.ind],x$zc0.ci[1,cl.ind]), y.axis[1:2] - adj, 
                c(x$dc0.ci[2,cl.ind],x$zc0.ci[2,cl.ind]), y.axis[1:2] - adj,
                lwd = lwd, lty = 3, col = col)
    }
    if(IND && !DIR) {
      points(x$dc0, y.axis[1] - adj, type = "p", pch = 19, cex = cex, col = col)
      segments(x$dc0.ci[1,cl.ind], y.axis[1] - adj, 
                x$dc0.ci[2,cl.ind], y.axis[1] - adj,
                lwd = lwd, lty = 3, col = col)
    }
    if(!IND && DIR) {
      points(x$zc0, y.axis[1] - adj, type = "p", pch = 19, cex = cex, col = col)
      segments(x$zc0.ci[1,cl.ind], y.axis[1] - adj, 
                x$zc0.ci[2,cl.ind], y.axis[1] - adj,
                lwd = lwd, lty = 3, col = col)
    }    
  }
  
  if (TOT) {
    points(x$tauc, 1 , type = "p", pch = 19, cex = cex, col = col)
    segments(x$tauc.ci[1,cl.ind], 1, x$tauc.ci[2,cl.ind], 1,
             lwd = lwd, col = col)
  }
  
  if(is.null(labels)){
    labels <- c("LACME","LANDE","LATE")[c(IND,DIR,TOT)]
  }
  axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE, ...)
  abline(v = 0, lty = 2)
}


