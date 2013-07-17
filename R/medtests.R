test.TMint <- function(x, ...){
  UseMethod("test.TMint")
}

test.TMint.default <- function(x, ...){
  stop("currently no test.TMint method exists for the input object.")
}

test.TMint.mediate <- function(x, conf.level = x$conf.level, ...){
  if(is.null(x$d0.sims) || is.null(x$d1.sims) || is.null(x$z0.sims) || is.null(x$z1.sims)){
    stop("simulation draws missing; rerun mediate with 'long' set to TRUE")
  }
  if(!x$INT){
    stop("outcome model must include interaction between treatment and mediator")
  }
  d.diff <- x$d1 - x$d0
  d.diff.sims <- x$d1.sims - x$d0.sims  
  
  # Format results in the htest format
  
  pv <- pval(d.diff.sims, d.diff)
  ci <- quantile(d.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2))
  
  null.value <- 0
  names(d.diff) <- names(null.value) <- "ACME(1) - ACME(0)"
  attr(ci, "conf.level") <- conf.level
  
  res <- list(statistic = d.diff, p.value = pv, conf.int = ci,
               null.value = null.value, alternative = "two.sided",
               method = "Test of ACME(1) - ACME(0) = 0",
               data.name = paste("estimates from", deparse(substitute(x))))
               
  class(res) <- "htest"
  return(res)
}

test.TMint.mediate.order <- function(x, conf.level = x$conf.level, ...){
  if(is.null(x$d0.sims) || is.null(x$d1.sims) || is.null(x$z0.sims) || is.null(x$z1.sims)){
    stop("simulation draws missing; rerun mediate with 'long' set to TRUE")
  }
  if(!x$INT){
    stop("outcome model must include interaction between treatment and mediator")
  }
  d.diff <- x$d1 - x$d0
  d.diff.sims <- x$d1.sims - x$d0.sims  
  
  # Format results in the htest format
  
  pv <- pval(d.diff.sims, d.diff)
  ci <- quantile(d.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2))
  
  null.value <- 0
  names(d.diff) <- names(null.value) <- "ACME(1) - ACME(0)"
  attr(ci, "conf.level") <- conf.level
  
  res <- list(statistic = d.diff, p.value = pv, conf.int = ci,
               null.value = null.value, alternative = "two.sided",
               method = "Test of ACME(1) - ACME(0) = 0",
               data.name = paste("estimates from", deparse(substitute(x))))
               
  class(res) <- "htest"
  return(res)
}



test.modmed <- function(object, ...){
  UseMethod("test.modmed")
}

test.modmed.default <- function(object, ...){
  stop("currently no test.modmed method exists for the input object.")
}


test.modmed.mediate <- function(object, covariates.1, covariates.2,
                                sims = object$sims, conf.level = object$conf.level, ...){
  
  cl <- getCall(object)
  cl$long <- TRUE
  cl$sims <- sims
  
  seed <- .Random.seed
  cl$covariates <- covariates.1
  out.1 <- eval(cl)
  
  .Random.seed <- seed
  cl$covariates <- covariates.2
  out.2 <- eval(cl)
  
  d1.diff <- out.1$d1 - out.2$d1
  d1.diff.sims <- out.1$d1.sims - out.2$d1.sims

  z0.diff <- out.1$z0 - out.2$z0
  z0.diff.sims <- out.1$z0.sims - out.2$z0.sims
  
  if(object$INT){
  
    d0.diff <- out.1$d0 - out.2$d0
    d0.diff.sims <- out.1$d0.sims - out.2$d0.sims

    z1.diff <- out.1$z1 - out.2$z1
    z1.diff.sims <- out.1$z1.sims - out.2$z1.sims
    
  }
  
  # Format results
  
  null.value <- 0
  
  if(object$INT){
    
    pv <- pval(d1.diff.sims, d1.diff)
    ci <- quantile(d1.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2))
    names(d1.diff) <- names(null.value) <- "ACME(1|covariates.1) - ACME(1|covariates.2)"
    attr(ci, "conf.level") <- conf.level
    res.d1 <- list(statistic = d1.diff, p.value = pv, conf.int = ci,
                 null.value = null.value, alternative = "two.sided",
                 method = "Test of ACME(1|covariates.1) - ACME(1|covariates.2) = 0",
                 data.name = paste("estimates from", deparse(substitute(object))))
  
    pv <- pval(d0.diff.sims, d0.diff)
    ci <- quantile(d0.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2))
    names(d0.diff) <- names(null.value) <- "ACME(0|covariates.1) - ACME(0|covariates.2)"
    attr(ci, "conf.level") <- conf.level
    res.d0 <- list(statistic = d0.diff, p.value = pv, conf.int = ci,
                 null.value = null.value, alternative = "two.sided",
                 method = "Test of ACME(0|covariates.1) - ACME(0|covariates.2) = 0",
                 data.name = paste("estimates from", deparse(substitute(object))))

    pv <- pval(z1.diff.sims, z1.diff)
    ci <- quantile(z1.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2))
    names(z1.diff) <- names(null.value) <- "ADE(1|covariates.1) - ADE(1|covariates.2)"
    attr(ci, "conf.level") <- conf.level
    res.z1 <- list(statistic = z1.diff, p.value = pv, conf.int = ci,
                 null.value = null.value, alternative = "two.sided",
                 method = "Test of ADE(1|covariates.1) - ADE(1|covariates.2) = 0",
                 data.name = paste("estimates from", deparse(substitute(object))))
  
    pv <- pval(z0.diff.sims, z0.diff)
    ci <- quantile(z0.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2))
    names(z0.diff) <- names(null.value) <- "ADE(0|covariates.1) - ADE(0|covariates.2)"
    attr(ci, "conf.level") <- conf.level
    res.z0 <- list(statistic = z0.diff, p.value = pv, conf.int = ci,
                 null.value = null.value, alternative = "two.sided",
                 method = "Test of ADE(0|covariates.1) - ADE(0|covariates.2) = 0",
                 data.name = paste("estimates from", deparse(substitute(object))))
    
    class(res.d1) <- class(res.d0) <- class(res.z1) <- class(res.z0) <- "htest"
    res <- list(res.d1, res.d0, res.z1, res.z0)

  } else {
    
    pv <- pval(d1.diff.sims, d1.diff)
    ci <- quantile(d1.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2))
    names(d1.diff) <- names(null.value) <- "ACME(covariates.1) - ACME(covariates.2)"
    attr(ci, "conf.level") <- conf.level
    res.d1 <- list(statistic = d1.diff, p.value = pv, conf.int = ci,
                 null.value = null.value, alternative = "two.sided",
                 method = "Test of ACME(covariates.1) - ACME(covariates.2) = 0",
                 data.name = paste("estimates from", deparse(substitute(object))))

    pv <- pval(z0.diff.sims, z0.diff)
    ci <- quantile(z0.diff.sims, c((1 - conf.level)/2, (1 + conf.level)/2))
    names(z0.diff) <- names(null.value) <- "ADE(covariates.1) - ADE(covariates.2)"
    attr(ci, "conf.level") <- conf.level
    res.z0 <- list(statistic = z0.diff, p.value = pv, conf.int = ci,
                 null.value = null.value, alternative = "two.sided",
                 method = "Test of ADE(covariates.1) - ADE(covariates.2) = 0",
                 data.name = paste("estimates from", deparse(substitute(object))))
  
    class(res.d1) <- class(res.z0) <- "htest"
    res <- list(res.d1, res.z0)
    
  }
  
  class(res) <- "test.modmed.mediate"
  return(res)
  
}

print.test.modmed.mediate <- function(x, ...){
  for(i in 1:length(x)){
    print(x[[i]], ...)
  }
}


##############################Order Functions##################################

test.TMint.mediate.order <- function(x, conf.level = x$conf.level, ...){
  if(is.null(x$d0.sims) || is.null(x$d1.sims) || is.null(x$z0.sims) || is.null(x$z1.sims)){
    stop("simulation draws missing; rerun mediate with 'long' set to TRUE")
  }
  if(!x$INT){
    stop("outcome model must include interaction between treatment and mediator")
  }
  d.diff <- x$d1 - x$d0
  d.diff.sims <- x$d1.sims - x$d0.sims  
  
  # Format results in the htest.order format
  # p-values
    y.lab <- sort(unique(levels(model.frame(x$model.y)[,1])))
    n.ycat <- length(y.lab)
    int.p <- rep(NA, n.ycat)
    ci.lo <- rep(NA, n.ycat)
    ci.up <- rep(NA, n.ycat)
    for(i in 1:n.ycat){
     int.p[i] <- pval(d.diff.sims[,i], d.diff[i])
     ci.lo[i] <- quantile(d.diff.sims[,i], (1 - conf.level)/2)
     ci.up[i] <- quantile(d.diff.sims[,i], (1 + conf.level)/2)
    }
    
    tab.1 <- rbind(d.diff, ci.lo, ci.up, int.p)
    rownames(tab.1)[1] <- "ACME(1) - ACME(0)"
    rownames(tab.1)[2] <- "2.5%"
    rownames(tab.1)[3] <- "97.5%"
    rownames(tab.1)[4] <- "p-value"


  null.value <- 0
  res <- list(statistic = d.diff, p.value = int.p, conf.int = c(ci.lo, ci.up),
               null.value = null.value, alternative = "two.sided")  
  class(res) <- "htest.order"
  print(tab.1, digits=3)
}


test.modmed.mediate.order <- function(object, covariates.1, covariates.2,
                                sims = object$sims, conf.level = object$conf.level, ...){
  
  cl <- getCall(object)
  cl$long <- TRUE
  cl$sims <- sims
  
  seed <- .Random.seed
  cl$covariates <- covariates.1
  out.1 <- eval(cl)
  
  .Random.seed <- seed
  cl$covariates <- covariates.2
  out.2 <- eval(cl)
  
  d1.diff <- out.1$d1 - out.2$d1
  d1.diff.sims <- out.1$d1.sims - out.2$d1.sims

  z0.diff <- out.1$z0 - out.2$z0
  z0.diff.sims <- out.1$z0.sims - out.2$z0.sims
  
  if(object$INT){
  
    d0.diff <- out.1$d0 - out.2$d0
    d0.diff.sims <- out.1$d0.sims - out.2$d0.sims

    z1.diff <- out.1$z1 - out.2$z1
    z1.diff.sims <- out.1$z1.sims - out.2$z1.sims
    
  }
  
  # Format results	
  	y.lab <- sort(unique(levels(model.frame(object$model.y)[,1])))
    n.ycat <- length(y.lab)
    out.names <- c()
    for(i in 1:length(y.lab)){
    out.names.tmp <- paste("Pr(Y=",y.lab[i],")",sep="")
    out.names <- c(out.names, out.names.tmp)
  }
    
  null.value <- 0
  
  if(object$INT){

    d1.int.p <- d1.ci.lo <- d1.ci.up <- d0.ci.lo <- d0.ci.up <- d0.int.p <- z1.int.p <- z0.int.p <- z1.ci.lo <- z1.ci.up <- z0.ci.lo <- z0.ci.up <- rep(NA, n.ycat)
    
    for(i in 1:n.ycat){
     d1.int.p[i] <- pval(d1.diff.sims[,i], d1.diff[i])
     d1.ci.lo[i] <- quantile(d1.diff.sims[,i], (1 - conf.level)/2)
     d1.ci.up[i] <- quantile(d1.diff.sims[,i], (1 + conf.level)/2)
     
     d0.int.p[i] <- pval(d0.diff.sims[,i], d0.diff[i])
     d0.ci.lo[i] <- quantile(d0.diff.sims[,i], (1 - conf.level)/2)
     d0.ci.up[i] <- quantile(d0.diff.sims[,i], (1 + conf.level)/2)
     
     z1.int.p[i] <- pval(z1.diff.sims[,i], z1.diff[i])
     z1.ci.lo[i] <- quantile(z1.diff.sims[,i], (1 - conf.level)/2)
     z1.ci.up[i] <- quantile(z1.diff.sims[,i], (1 + conf.level)/2)
     
     z0.int.p[i] <- pval(z0.diff.sims[,i], d0.diff[i])
     z0.ci.lo[i] <- quantile(z0.diff.sims[,i], (1 - conf.level)/2)
     z0.ci.up[i] <- quantile(z0.diff.sims[,i], (1 + conf.level)/2)
     
    }
    
    tab.1 <- rbind(d1.diff, d1.ci.lo, d1.ci.up, d1.int.p)
    rownames(tab.1)[1] <- "ACME(1) - ACME(0)"
    rownames(tab.1)[2] <- "2.5%"
    rownames(tab.1)[3] <- "97.5%"
    rownames(tab.1)[4] <- "p-value"
    colnames(tab.1) <-  out.names
    tab.2 <- rbind(d0.diff, d0.ci.lo, d0.ci.up, d0.int.p)
    rownames(tab.2)[1] <- "ACME(1) - ACME(0)"
    rownames(tab.2)[2] <- "2.5%"
    rownames(tab.2)[3] <- "97.5%"
    rownames(tab.2)[4] <- "p-value"
    colnames(tab.2) <-  out.names
    tab.3 <- rbind(z1.diff, z1.ci.lo, z1.ci.up, z1.int.p)
    rownames(tab.3)[1] <- "ACME(1) - ACME(0)"
    rownames(tab.3)[2] <- "2.5%"
    rownames(tab.3)[3] <- "97.5%"
    rownames(tab.3)[4] <- "p-value"
    colnames(tab.3) <-  out.names
    tab.4 <- rbind(z0.diff, z0.ci.lo, z0.ci.up, z0.int.p)
    rownames(tab.4)[1] <- "ACME(1) - ACME(0)"
    rownames(tab.4)[2] <- "2.5%"
    rownames(tab.4)[3] <- "97.5%"
    rownames(tab.4)[4] <- "p-value"
    colnames(tab.4) <-  out.names
        
  cat("Test of ACME(1|covariates.1) - ACME(1|covariates.2) = 0 \n")
  print(tab.1, digits=3)
  cat("\n")
  cat("Test of ACME(0|covariates.1) - ACME(0|covariates.2) = 0 \n")
  print(tab.2, digits=3)
  cat("\n")
  cat("Test of ADE(1|covariates.1) - ADE(1|covariates.2) = 0 \n")
  print(tab.3, digits=3)
  cat("\n")
  cat("Test of ADE(0|covariates.1) - ADE(0|covariates.2) = 0 \n")
  print(tab.4, digits=3)
  
        res <- list(statistic = c(d1.diff, d0.diff, z1.diff, z0.diff), p.value = c(d1.int.p,d0.int.p, z1.int.p, z0.int.p), conf.int = c(d1.ci.lo, d1.ci.up, d0.ci.lo, d0.ci.up, z1.ci.lo, z1.ci.up, z0.ci.lo, z0.ci.up), null.value = null.value)  


  } else {
  	
  	y.lab <- sort(unique(levels(model.frame(object$model.y)[,1])))
  	y.lab <- sort(unique(levels(model.frame(contord$model.y)[,1])))
    n.ycat <- length(y.lab)
    d1.int.p <- d1.ci.lo <- d1.ci.up <- d0.int.p <- d1.int.p <- z1.int.p <- z0.int.p <- z1.ci.lo <- z1.ci.up <- z0.ci.lo <- z0.ci.up <- rep(NA, n.ycat)
    
    for(i in 1:n.ycat){
     d1.int.p[i] <- pval(d1.diff.sims[,i], d1.diff[i])
     d1.ci.lo[i] <- quantile(d1.diff.sims[,i], (1 - conf.level)/2)
     d1.ci.up[i] <- quantile(d1.diff.sims[,i], (1 + conf.level)/2)
     
     d0.int.p[i] <- pval(d0.diff.sims[,i], d0.diff[i])
     d0.ci.lo[i] <- quantile(d0.diff.sims[,i], (1 - conf.level)/2)
     d0.ci.up[i] <- quantile(d0.diff.sims[,i], (1 + conf.level)/2)
     
     z1.int.p[i] <- pval(z1.diff.sims[,i], z1.diff[i])
     z1.ci.lo[i] <- quantile(z1.diff.sims[,i], (1 - conf.level)/2)
     z1.ci.up[i] <- quantile(z1.diff.sims[,i], (1 + conf.level)/2)
     
     z0.int.p[i] <- pval(z0.diff.sims[,i], d0.diff[i])
     z0.ci.lo[i] <- quantile(z0.diff.sims[,i], (1 - conf.level)/2)
     z0.ci.up[i] <- quantile(z0.diff.sims[,i], (1 + conf.level)/2)
     
    }
    
    tab.1 <- rbind(d1.diff, d1.ci.lo, d1.ci.up, d1.int.p)
    rownames(tab.1)[1] <- "ACME(1) - ACME(0)"
    rownames(tab.1)[2] <- "2.5%"
    rownames(tab.1)[3] <- "97.5%"
    rownames(tab.1)[4] <- "p-value"
    colnames(tab.1) <-  out.names
    tab.2 <- rbind(z0.diff, z0.ci.lo, z0.ci.up, z0.int.p)
    rownames(tab.2)[1] <- "ACME(1) - ACME(0)"
    rownames(tab.2)[2] <- "2.5%"
    rownames(tab.2)[3] <- "97.5%"
    rownames(tab.2)[4] <- "p-value"
    colnames(tab.2) <-  out.names
    cat("Test of ACME(covariates.1) - ACME(covariates.2) = 0 \n")
    print(tab.1, digits=3)
    cat("\n")
    cat("Test of ADE(covariates.1) - ADE(covariates.2) = 0\n")
    print(tab.2, digits=3)
  
      res <- list(statistic = c(d1.diff,z0.diff), p.value = c(d1.int.p,z0.int.p), conf.int = c(d1.ci.lo, d1.ci.up, z0.ci.lo, z0.ci.up),
               null.value = null.value)   
  }
  class(res) <- "test.modmed.mediate.order"
  return(res)
}


