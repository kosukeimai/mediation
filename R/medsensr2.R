medsensr2 <- function(model.m, model.y, T="treat.name", M="med.name")
{
  ## Step 0: Setting Variable labels
  ## Uppercase letters (e.g. T) = labels in the input matrix
  ## Uppercase letters + ".out" (e.g. T.out) = labels in the regression output
  y.t.data <- model.frame(model.y)
  Y <- colnames(y.t.data)[1]
  
  if(is.factor(y.t.data[,T])==TRUE){
    cat.c <- levels(y.t.data[,T])[1] 
    cat.t <- levels(y.t.data[,T])[2]
    T.out <- paste(T,cat.t, sep="") 
  } else {
    cat.c <- NULL
    cat.t <- NULL
    T.out <- paste(T,cat.t, sep="")
  }
  
  if(is.factor(y.t.data[,M])==TRUE){
    cat.m0 <- levels(y.t.data[,M])[1] 
    cat.m1 <- levels(y.t.data[,M])[2]
    M.out <- paste(M,cat.m1, sep="") 
  } else {
    cat.m0 <- NULL
    cat.m1 <- NULL
    M.out <- paste(M,cat.m1, sep="")
  }
  
                                        #Functions for Calculating Rho Values
  rho.pos <- function(a,b){
    rho.sq <- a*b
    sqrt(rho.sq)
  }
  
  rho.neg <- function(a,b){
    rho.sq <- a*b
    sqrt(rho.sq) *-1
  }
  
  rho.pos.tilde <- function(a,b){
    rho.sq <- a*b/((1 - r.sq.m)*(1 - r.sq.y))
    sqrt(rho.sq)
  }
  
  rho.neg.tilde <- function(a,b){
    rho.sq <- a*b/((1 - r.sq.m)*(1 - r.sq.y))
    sqrt(rho.sq) *-1
  }
#########################################################
  ## CASE 1: Continuous Outcome + Continuous Mediator
#########################################################
  if(class(model.y)[1]=="lm" & class(model.m)[1]=="lm") {
                                        #R^2 Values
    r.sq.m <- summary(model.m)$adj.r.squared
    r.sq.y <- summary(model.y)$adj.r.squared
    r.sq.m.star <- seq(0.01,.98,.01)
    r.sq.y.star <- seq(0.01,.98,.01)
    r.tilde.m <- (1 - r.sq.m)*r.sq.m.star
    r.tilde.y <- (1 - r.sq.y)*r.sq.y.star   
    sigma.2 <- var(model.m$resid)
    e.3.star <- model.y$coef[M]*model.m$resid + model.y$resid 
    sigma.3.star <- var(e.3.star)
    sigma.23.star <- cov(model.m$resid, e.3.star) 
    
    acme.lm <- function(rho){
      model.m$coef[T.out]*(sigma.23.star/sigma.2) - rho/sqrt(sigma.2) * sqrt(1/(1-rho^2)*(sigma.3.star - sigma.23.star/sigma.2))
    }
                                        #Matrix of Rho Values
    rho.p <- outer(r.sq.m.star, r.sq.y.star, rho.pos)
    rho.n <- outer(r.sq.m.star, r.sq.y.star, rho.neg)
    rho.til.p <- outer(r.tilde.m, r.tilde.y, rho.pos.tilde)
    rho.til.n <- outer(r.tilde.m, r.tilde.y, rho.neg.tilde)
                                        #Compute ACME
    r.p <- sapply(rho.p, acme.lm)
    r.p <- matrix(r.p, nrow(rho.p), ncol(rho.p))
    r.n <- sapply(rho.n, acme.lm)
    r.n <- matrix(r.n, nrow(rho.n), ncol(rho.n))
    r.til.p <- sapply(rho.til.p, acme.lm)
    r.til.p <- matrix(r.til.p, nrow(rho.til.p), ncol(rho.til.p))
    r.til.n <- sapply(rho.til.n, acme.lm)
    r.til.n <- matrix(r.til.n, nrow(rho.til.n), ncol(rho.til.n))
    
    type <- "ct"
    ## END OF CASE 1: Continuous Outcome + Continuous Mediator    
  } else if(class(model.y)[1]=="lm" & class(model.m)[1]=="glm") {
#########################################################
    ## CASE 2: Continuous Outcome + Binary Mediator
########################################################
    ##R^2 Values
    Mmodel.coef <- model.m$coef
    m.mat <- model.matrix(model.m)
    fitted<- m.mat %*% Mmodel.coef
    var.mstar <- var(fitted)
    r.sq.m <-var.mstar/(1+var.mstar)
                                        #r.sq.m <- 1 - (model.m$deviance/model.m$null.deviance)
    
    r.sq.y <- summary(model.y)$adj.r.squared
    r.sq.m.star <- seq(0.01, 1, length.out=25)
    r.sq.y.star <- seq(0.01, 1, length.out=25)
    r.tilde.m <- (1 - r.sq.m)*r.sq.m.star
    r.tilde.y <- (1 - r.sq.y)*r.sq.y.star
    
    ## Variable values (LABEL.value)
    Y.value <- y.t.data[,1]
    y.k <- length(model.y$coef)
    
                                        # Step 1: ACME Preliminaries
                                        # Step 1-1: Med model coefficients
    Mmodel.coef <- model.m$coef
    
                                        # Step 1-2: lambda_0 and lambda_1 values; 
    m.mat <- model.matrix(model.m)
    m.mat.1 <- model.matrix(model.m)
    m.mat.1[,T.out] <- 1 # M-model matrix with t=1
    m.mat.0 <- model.matrix(model.m)
    m.mat.0[,T.out] <- 0 # M-model matrix with t=0
    mu <- m.mat %*% Mmodel.coef # E(M|T,X)
    mu.1 <- m.mat.1 %*% Mmodel.coef # E(M|T=1,X)
    mu.0 <- m.mat.0 %*% Mmodel.coef # E(M|T=0,X)
    lambda11 <- dnorm(-mu.1) / pnorm(mu.1) #lambda for m=1,t=1
    lambda10 <- dnorm(-mu.0) / pnorm(mu.0) #lambda for m=1,t=0
    lambda01 <- -dnorm(-mu.1) / pnorm(-mu.1) #lambda for m=0,t=1
    lambda00 <- -dnorm(-mu.0) / pnorm(-mu.0) #lambda for m=0,t=0
    
                                        # Step 1-3: Define lambda function
    lambda <- function(mmodel, mcoef) {
      mu <- model.matrix(mmodel) %*% mcoef
      m <- mmodel$y #this is M
      return((m*dnorm(-mu)-(1-m)*dnorm(-mu))/(m*pnorm(mu)+(1-m)*pnorm(-mu)))
    }
    
    acme <- function(rho){
      ## Step 2-1: Obtain the initial Y model with the correction term
      adj <- lambda(model.m, Mmodel.coef) * rho # the adjustment term
      w <- 1 - rho^2*lambda(model.m, Mmodel.coef)*(lambda(model.m, Mmodel.coef) + mu)
      y.t.data.adj <- data.frame(y.t.data, w, adj)
      model.y.adj <- update(model.y, as.formula(paste(". ~ . + adj")), weights=w, data=y.t.data.adj)
      sigma.3 <- summary(model.y.adj)$sigma
      
      ## Step 2-2: Update the Y model via Iterative FGLS
      eps <- .1
      sigma.dif <- 1
      while(abs(sigma.dif) > eps){
        Y.star <- Y.value - sigma.3 * adj
        y.t.data.star <- data.frame(Y.star, y.t.data.adj)
        model.y.update <- update(model.y, as.formula(paste("Y.star ~ .")), weights=w, data=y.t.data.star)
        sigma.3.temp <- summary(model.y.update)$sigma
        sigma.dif <- sigma.3.temp - sigma.3
        sigma.3 <- sigma.3.temp
      }
      
      ## Step 2-3: Adjusted Y model parameters
      Ymodel.coef <- model.y.update$coef
      ## Step 2-4: ACMEs; means are over observations
      mean((Ymodel.coef[M.out]) * (pnorm(mu.1) - pnorm(mu.0)))
    }
                                        #Matrix of Rho Values
    rho.p <- outer(r.sq.m.star, r.sq.y.star, rho.pos)
    rho.n <- outer(r.sq.m.star, r.sq.y.star, rho.neg)
    rho.til.p <- outer(r.tilde.m, r.tilde.y, rho.pos.tilde)
    rho.til.n <- outer(r.tilde.m, r.tilde.y, rho.neg.tilde)
                                        #ACMEs 
    r.p <- sapply(rho.p, acme)
    r.p <- matrix(r.p, nrow(rho.p), ncol(rho.p))
    r.n <- sapply(rho.n, acme)
    r.n <- matrix(r.n, nrow(rho.n), ncol(rho.n))
    r.til.p <- sapply(rho.til.p, acme)
    r.til.p <- matrix(r.til.p, nrow(rho.til.p), ncol(rho.til.p))
    r.til.n <- sapply(rho.til.n, acme)
    r.til.n <- matrix(r.til.n, nrow(rho.til.n), ncol(rho.til.n))
    type <- "bm"
    
    ## END OF CASE 2: Continuous Outcome + Binary Mediator
  } else if(class(model.y)[1]=="glm" & class(model.m)[1]=="lm") {
#########################################################
    ## CASE 3: Binary Outcome + Continuous Mediator
#########################################################
    
                                        #R^2 Values
    r.sq.m <- summary(model.m)$adj.r.squared
    
    Ymodel.coef <- model.y$coef
    y.mat <- model.matrix(model.y)
    fitted<- y.mat %*% Ymodel.coef
    var.ystar <- var(fitted)
    r.sq.y <-var.ystar/(1+var.ystar)
                                        #r.sq.y <- 1 - (model.y$deviance/model.y$null.deviance)
    
    r.sq.m.star <- seq(.01,1,.01)
    r.sq.y.star <- seq(.01,1,.01)
    r.tilde.m <- (1 - r.sq.m)*r.sq.m.star
    r.tilde.y <- (1 - r.sq.y)*r.sq.y.star
    
                                        # Step 1: Obtain Model Parameters
    Mmodel.coef <- model.m$coef
    m.k <- length(model.m$coef)
    if(is.factor(y.t.data[,T])==TRUE){
      beta2.temp <- Mmodel.coef[T.out] 
    } else {
      beta2.temp <- Mmodel.coef[T]    
    }
    sigma.2 <- summary(model.m)$sigma
    sig2.shape <- model.m$df/2
    sig2.invscale <- (model.m$df/2) * sigma.2^2
    Ymodel.coef <- model.y$coef
    y.k <- length(Ymodel.coef)
    gamma.tilde <- Ymodel.coef[M.out]
    
                                        # Step 2: Compute ACME via the procedure in IKT
    ## Step 2-1: Estimate Error Correlation from inconsistent estimate of Y on M
    rho12 <- (sigma.2 * gamma.tilde) / (1 + sqrt(sigma.2^2*gamma.tilde^2))
    
    ## Step 2-2: Calculate alpha_1, beta_1 and xi_1
    YTmodel.coef <- Ymodel.coef[!names(Ymodel.coef)%in%M.out] * sqrt(1-rho12^2) %x% t(rep(1,y.k-1)) + Mmodel.coef * (rho12/sigma.2) %x% t(rep(1,y.k-1))
    
    ## Step 2-3: Calculate Gamma
    ## Data matrices for the Y model less M
    y.mat.1 <- model.matrix(model.y)[,!colnames(model.matrix(model.y))%in%M.out]
    y.mat.1[,T.out] <- 1
    y.mat.0 <- model.matrix(model.y)[,!colnames(model.matrix(model.y))%in%M.out]
    y.mat.0[,T.out] <- 0
    
    acme <- function(rho){
      gamma.temp <- (-rho + rho12*sqrt((1-rho^2)/(1-rho12^2)))/sigma.2
      d0.temp <- mean( pnorm(y.mat.0 %*% t(YTmodel.coef) + 
                             gamma.temp*beta2.temp/sqrt(gamma.temp^2*sigma.2^2+2*gamma.temp*rho*sigma.2+1)) -
                      pnorm(y.mat.0 %*% t(YTmodel.coef)) )
      d1.temp <- mean( pnorm(y.mat.1 %*% t(YTmodel.coef)) - pnorm(y.mat.1 %*% t(YTmodel.coef) - 
                                                                  gamma.temp*beta2.temp/sqrt(gamma.temp^2*sigma.2^2+2*gamma.temp*rho*sigma.2+1)) )
      (d0.temp + d1.temp)/(2)
    }
    
                                        #Rho Values    
    rho.p <- outer(r.sq.m.star, r.sq.y.star, rho.pos)
    rho.n <- outer(r.sq.m.star, r.sq.y.star, rho.neg)
    rho.til.p <- outer(r.tilde.m, r.tilde.y, rho.pos.tilde)
    rho.til.n <- outer(r.tilde.m, r.tilde.y, rho.neg.tilde)
                                        #Calculate ACMEs
    r.p <- sapply(rho.p, acme)
    r.p <- matrix(r.p, nrow(rho.p), ncol(rho.p))
    r.n <- sapply(rho.n, acme)
    r.n <- matrix(r.n, nrow(rho.n), ncol(rho.n))
    r.til.p <- sapply(rho.til.p, acme)
    r.til.p <- matrix(r.til.p, nrow(rho.til.p), ncol(rho.til.p))
    r.til.n <- sapply(rho.til.n, acme)
    r.til.n <- matrix(r.til.n, nrow(rho.til.n), ncol(rho.til.n))
    type <- "bo"
    ## END OF CASE 3: Binary Outcome + Continuous Mediator    
  }

        #Extract R^2 Values Where ACME=0
  n.r <- nrow(r.til.p)
  n.c <- ncol(r.til.p)
  r.c <- n.r*n.c
  
                                        #Vectorize Matrices
  r.til.p.vec <- matrix(r.til.p, r.c, 1)
  out.1 <- cbind(round(r.tilde.m, 4), round(r.tilde.y, 4), r.til.p.vec)
  ii <- which(abs(r.til.p.vec-0)==min(abs(r.til.p.vec-0)))
  out.til.p <- out.1[ii,1:2]
  
  r.til.n.vec <- matrix(r.til.n, r.c, 1)
  out.2 <- cbind(round(r.tilde.m, 4), round(r.tilde.y, 4), r.til.n.vec)
  ii <- which(abs(r.til.n.vec-0)==min(abs(r.til.n.vec-0)))
  out.til.n <- out.2[ii,1:2]
  
  r.n.vec <- matrix(r.n, r.c, 1)
  out.3 <- cbind(round(r.sq.m.star, 4), round(r.sq.y.star, 4), r.n.vec)
  ii <- which(abs(r.n.vec-0)==min(abs(r.n.vec-0)))
  out.r.n <- out.3[ii,1:2]
  
  r.p.vec <- matrix(r.p, r.c, 1)
  out.4 <- cbind(round(r.sq.m.star, 4), round(r.sq.y.star, 4), r.p.vec)
  ii <- which(abs(r.p.vec-0)==min(abs(r.p.vec-0)))
  out.r.p <- out.4[ii,1:2]
  
  out <- list(r.p = r.p, r.n = r.n, r.til.p = r.til.p, r.til.n = r.til.n, out.til.p=out.til.p, out.til.n=out.til.n, out.r.n=out.r.n,  out.r.p=out.r.p, type=type, r.sq.m.star=r.sq.m.star, r.sq.y.star=r.sq.y.star, r.tilde.m=r.tilde.m, r.tilde.y=r.tilde.y)
    class(out) <- "medsensr2"
    out
## END OF SENSITIVITY FUNCTION
}

print.medsensr2 <- function(x, ...){
    print(unlist(x[1:5]))
    invisible(x)
    }

summary.medsensr2 <- function(object, ...)
    structure(object, class = c("summary.medsensr2", class(object)))
 
print.summary.medsensr2 <- function(x, ...){
    cat("\nMediation Sensitivity Analysis\n")
    
    if(length(x$out.r.p)==2){
    x$out.r.p <- as.matrix(x$out.r.p)
    x$out.r.p <- t(x$out.r.p)
    }
    if(length(x$out.r.n)==2){
    x$out.r.n <- as.matrix(x$out.r.n)
    x$out.r.n <- t(x$out.r.n)
    }
    if(length(x$out.til.p)==2){
    x$out.til.p <- as.matrix(x$out.til.p)
    x$out.til.p <- t(x$out.til.p)
    }
    if(length(x$out.til.n)==2){
    x$out.til.n <- as.matrix(x$out.til.n)
    x$out.til.n <- t(x$out.til.n)
    }
    
    colnames(x$out.r.p) <-  c("M R-squared Star", "Y R-squared Star")
    rownames(x$out.r.p) <- NULL
    colnames(x$out.r.n) <-  c("M R-squared Star", "Y R-squared Star")
    rownames(x$out.r.n) <- NULL
          
    colnames(x$out.til.p) <-  c("M R-squared Tilde", "Y R-squared Tilde")
    rownames(x$out.til.p) <- NULL
    colnames(x$out.til.n) <-  c("M R-squared Tilde", "Y R-squared Tilde")
    rownames(x$out.til.n) <- NULL
      
    cat("Sign of Unobserved Confounder is Positive\n")
    print(x$out.til.p)
    print(x$out.r.p)
    
    cat("Sign of Unobserved Confounder is Negative\n")
    print(x$out.til.n)
    print(x$out.r.n)
    invisible(x)
         
}


plot.medsensr2 <- function(x, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL,... ){
    #if(r.type== 1){
        
        if(x$type=="bo"){
            lev.1 <- c(0,.1)
            lev.2 <- c(.1,-.1)
            } else if(x$type=="bm"){
            lev.1 <- c(-.1,.1)
            lev.2 <- c(-.15,.1)
            } else {
            lev.1 = c(-.5,1.7)
            lev.2 = c(0,-3)
                    }
        
        par(mfrow=c(2,2))
        par(mar=c(3.55,7,4,0))  #(b,l,t,r)
        #Column 1 Sgn -1
        contour(x$r.sq.m.star, x$r.sq.y.star, x$r.n, ylim=c(0,1), xlim=c(0,1), levels=pretty(lev.1, 25), asp=1)
        title(xlab=expression(paste(R[M]^{2},"*")), line=2.5, cex.lab=.9)
        title(ylab=expression(paste(R[Y]^2,"*")), line=2.5, cex.lab=.9 )
        title(main=expression(paste("sgn", (lambda[2]*lambda[3])==-1)))
        axis(2,at=seq(0,1,by=.1))
        axis(1,at=seq(0,1,by=.1))
        mtext("Proportion of unexplained variance \n explained by an unobserved confounder", side=2, line=4, cex=.9)

        #Column 2 Sgn 1
        par(mar=c(3.55,5,4,2)) 
        contour(x$r.sq.m.star, x$r.sq.y.star, x$r.p, ylim=c(0,1), xlim=c(0,1), levels = pretty(lev.2, 25), asp=1)
        title(xlab=expression(paste(R[M]^2,"*")), line=2.5, cex.lab=.9)
        title(ylab=expression(paste(R[Y]^2,"*")), line=2.5, cex.lab=.9 )
        title(main=expression(paste("sgn", (lambda[2]*lambda[3])==1)))
        axis(2,at=seq(0,1,by=.1))
        axis(1,at=seq(0,1,by=.1))

        #} else if(r.type==2) {
        #par(mfrow=c(1,2))
        #Lower Left
        par(mar=c(5,7,2,0))
        contour(x$r.tilde.m, x$r.tilde.y, x$r.til.n, ylim=c(0,1), xlim=c(0,1), levels=pretty(lev.1, 25), asp=1)
        title(xlab=expression(paste(tilde(R)[M]^{2})), line=2.5, cex.lab=.9)
        title(ylab=expression(paste(tilde(R)[Y]^2)), line=2.5, cex.lab=.9 )
        #title(main=expression(paste("sgn", (lambda[2]*lambda[3])==-1)))
        mtext("Proportion of original variance \n explained by an unobserved confounder", side=2, line=4, cex=.9)
        axis(2,at=seq(0,1,by=.1))
        axis(1,at=seq(0,1.1,by=.1))

        #Lower Right
        par(mar=c(5,5,2,2)) 
        contour(x$r.tilde.m, x$r.tilde.y, x$r.til.p, ylim=c(0,1), xlim=c(0,1), levels = pretty(lev.2, 25), asp=1)
        title(xlab=expression(paste(tilde(R)[M]^2)), line=2.5, cex.lab=.9)
        title(ylab=expression(paste(tilde(R)[Y]^2)), line=2.5, cex.lab=.9 )
        #title(main=expression(paste("sgn", (lambda[2]*lambda[3])==1)))
        axis(2,at=seq(0,1,by=.1))
        axis(1,at=seq(0,2,by=.1))
            
        #}
    
        }
