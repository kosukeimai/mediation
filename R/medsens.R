medsens <- function(model.m, model.y, T="treat.name", M="med.name", INT=FALSE, DETAIL=FALSE, sims=1000, eps=.Machine$double.eps)
{

    #########################################################
    ## Setting Up Sensitivity Parameter
    #########################################################
    if(DETAIL==TRUE){
        rho <- round(seq(-.99, .99, by = .01),2)    
        } else {
        rho <- round(seq(-.90, .90, by = .1),2)
            }
    
    #########################################################
    ## CASE 1: Continuous Outcome + Continuous Mediator
    #########################################################
    if(class(model.y)[1]=="lm" & class(model.m)[1]=="lm") {

        d0 <- matrix(NA, length(rho), 1)
        d1 <- matrix(NA, length(rho), 1)
        d0.var <- matrix(NA, length(rho), 1)
        d1.var <- matrix(NA, length(rho), 1)
        
        y.t.data <- model.frame(model.y)    
        if(is.factor(y.t.data[,paste(T)])==TRUE){
            cat.c <- levels(y.t.data[,T])[1] 
            cat.t <- levels(y.t.data[,T])[2]
            T.cat <- paste(T,cat.t, sep="") 
            } else {
            cat.c <- NULL
            cat.t <- NULL
            T.cat <- paste(T,cat.t, sep="")
            }
            
        if(INT==TRUE){
        int.lab <- paste(T.cat,M, sep=":")
        t.m <- paste(T,M, sep=":")
            }
                
        #Estimate Error Correlation
        if(INT==TRUE){
            mod.y <- update(model.y,as.formula(paste(". ~ . -", t.m, "-", M)))
            } else {
            mod.y <- update(model.y,as.formula(paste(". ~ . -", M)))
                }
        err.cr <- cor(model.m$resid, mod.y$resid)
                    
        
        for(i in 1:length(rho)){
        
        e.cor <- rho[i]
        
        b.dif <- 1
        #eps <- .001 #.Machine$double.eps
        
        #Stacked Equations
        m.mat <- model.matrix(model.m)
        y.mat <- model.matrix(model.y)
        m.k <- ncol(m.mat)
        m.n <- nrow(m.mat)
        y.k <- ncol(y.mat)
        y.n <- nrow(y.mat)
        n <- y.n
        m.zero <- matrix(0, m.n, y.k)
        y.zero <- matrix(0, y.n, m.k)
        X.1 <- cbind(m.mat, m.zero)
        X.2 <- cbind(y.zero, y.mat)
        X <- rbind(X.1, X.2)
        
        m.frame <- model.frame(model.m)
        y.frame <- model.frame(model.y)
        Y.c <- rbind(as.matrix(m.frame[,1]), as.matrix(y.frame[,1]))
        
        #Estimates of OLS Start Values
        inxx <- solve(crossprod(X))
        b.ols <- inxx %*% crossprod(X,Y.c)
        b.tmp <- b.ols
        
        while(abs(b.dif) > eps){
        
        e.hat <- as.matrix(Y.c - (X %*% b.tmp))
        
        e.1 <- e.hat[1:n]
        e.2 <- e.hat[(n+1): (2*n)]
        
        sd.1 <- sd(e.1)
        sd.2 <- sd(e.2)
        
        omega <- matrix(NA, 2,2)
        
        omega[1,1] <- crossprod(e.1)/(n-1)
        omega[2,2] <- crossprod(e.2)/(n-1) 
        omega[2,1] <- e.cor*sd.1*sd.2
        omega[1,2] <- e.cor*sd.1*sd.2
        
        I <- diag(1,n)
        omega.i <- solve(omega)
        v.i <-  kronecker(omega.i, I)
        
        X.sur <- crossprod(X, v.i) %*% X 
        b.sur <- solve(X.sur) %*% crossprod(X, v.i) %*% Y.c
        
        #Variance-Covariance Matrix
        v.cov <- crossprod(X, v.i) %*% X
        v.cov <- solve(v.cov)
        
        b.old <- b.tmp
        b.dif <- sum((b.sur - b.old)^2)
        b.tmp <- b.sur
        
        }
        
        #Name Elements - Extract Quantities
        m.names <- names(model.m$coef)
        y.names <- names(model.y$coef)
        b.names <- c(m.names, y.names)
        row.names(b.sur) <- b.names
        m.coefs <- as.matrix(b.sur[1:m.k])
        y.coefs <- as.matrix(b.sur[(m.k+1):(m.k+y.k)])
        row.names(m.coefs) <- m.names
        row.names(y.coefs) <- y.names
        rownames(v.cov) <- b.names
        colnames(v.cov) <- b.names
        v.m <- v.cov[1:m.k,1:m.k]
        v.y <- v.cov[(m.k+1):(m.k+y.k),(m.k+1):(m.k+y.k)]
        
        #Save Estimates
        if(INT==TRUE){
            d0[i,] <- m.coefs[paste(T.cat),]*y.coefs[paste(M),] 
            d1[i,] <- m.coefs[paste(T.cat),]*(y.coefs[paste(M),] + y.coefs[paste(int.lab),])
            } else {
                d0[i,] <- m.coefs[paste(T.cat),]*y.coefs[paste(M),]
                d1[i,] <- m.coefs[paste(T.cat),]*y.coefs[paste(M),] 
                }
        
        #Save Variance Estimates
        if(INT==TRUE){
            d0.var[i,] <- (y.coefs[paste(M),] + 0*y.coefs[paste(int.lab),])^2*v.m[T.cat,T.cat] + m.coefs[paste(T.cat),]^2*(v.y[M,M] + 0*v.y[int.lab, int.lab] + 0*2*v.y[M, int.lab])
            d1.var[i,] <- (y.coefs[paste(M),] + y.coefs[paste(int.lab),])^2*v.m[T.cat,T.cat] + m.coefs[paste(T.cat),]^2*(v.y[M,M] + v.y[int.lab, int.lab] + 2*v.y[M, int.lab])
            } else {
            d0.var[i,] <- (m.coefs[paste(T.cat),]^2*v.y[M,M]) + (y.coefs[paste(M),]^2*v.m[T.cat,T.cat])
            d1.var[i,] <- (m.coefs[paste(T.cat),]^2*v.y[M,M]) + (y.coefs[paste(M),]^2*v.m[T.cat,T.cat])
                }
                
        rm(b.sur, m.coefs, y.coefs, v.cov, v.m, v.y)
        
        }
        
        if(INT==FALSE){
        upper.d0 <- d0 + qnorm(0.975) * sqrt(d0.var)
        lower.d0 <- d0 - qnorm(0.975) * sqrt(d0.var)
        upper.d1 <- NULL
        lower.d1 <- NULL
        ind.d0 <- as.numeric(lower.d0 < 0 & upper.d0 > 0)
        ind.d1 <- NULL  
            } else {
        upper.d0 <- d0 + qnorm(0.975) * sqrt(d0.var)
        lower.d0 <- d0 - qnorm(0.975) * sqrt(d0.var)
        upper.d1 <- d1 + qnorm(0.975) * sqrt(d1.var)
        lower.d1 <- d1 - qnorm(0.975) * sqrt(d1.var)    
        ind.d0 <- as.numeric(lower.d0 < 0 & upper.d0 > 0)
        ind.d1 <- as.numeric(lower.d1 < 0 & upper.d1 > 0)
                }
        if(INT==TRUE){
        ii <- which(abs(d0-0)==min(abs(d0-0)))
        kk <- which(abs(d1-0)==min(abs(d1-0)))
		err.cr.1 <- rho[ii]
		err.cr.2 <- rho[kk]
		err.cr <- c(err.cr.1, err.cr.2)
		} 
        type <- "ct"
        out <- list(rho = rho, err.cr=err.cr, d0=d0, d1=d1, upper.d0=upper.d0, lower.d0=lower.d0, upper.d1=upper.d1, lower.d1=lower.d1, ind.d0=ind.d0, ind.d1=ind.d1, INT=INT, DETAIL=DETAIL, sims=sims,tau=NULL, upper.tau=NULL, lower.tau=NULL, nu=NULL, upper.nu=NULL, lower.nu=NULL, type=type)
        class(out) <- "medsens"
        out

    ## END OF CASE 1: Continuous Outcome + Continuous Mediator    
    } else

    #########################################################
    ## CASE 2: Continuous Outcome + Binary Mediator
    #########################################################
    if(class(model.y)[1]=="lm" & class(model.m)[1]=="glm") {

        # Step 0: Setting Variable labels
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
        
        if(INT==TRUE){
        TM <- paste(T,M, sep=":")
        TM.out <- paste(T.out,M.out, sep=":")
            }
        
        ## Variable values (LABEL.value)
        Y.value <- y.t.data[,1]
        y.k <- length(model.y$coef)
        
        # Step 1: Pre-loop computations
        # Step 1-1: Bootstrap M model parameters
        Mmodel.coef <- model.m$coef
        Mmodel.var.cov <- vcov(model.m)
        Mmodel.coef.boot <- mvrnorm(sims, mu=Mmodel.coef, Sigma=Mmodel.var.cov) # bootstrap M-model parameters
        
        # Step 1-2: Bootstrap lambda_0 and lambda_1; lambdas are (n x sims) matrix
        m.mat <- model.matrix(model.m)
        m.mat.1 <- model.matrix(model.m)
        m.mat.1[,T.out] <- 1 # M-model matrix with t=1
        m.mat.0 <- model.matrix(model.m)
        m.mat.0[,T.out] <- 0 # M-model matrix with t=0
        mu.boot <- m.mat %*% t(Mmodel.coef.boot) # E(M|T,X)
        mu.1.boot <- m.mat.1 %*% t(Mmodel.coef.boot) # E(M|T=1,X)
        mu.0.boot <- m.mat.0 %*% t(Mmodel.coef.boot) # E(M|T=0,X)
        lambda11 <- dnorm(-mu.1.boot) / pnorm(mu.1.boot) #lambda for m=1,t=1
        lambda10 <- dnorm(-mu.0.boot) / pnorm(mu.0.boot) #lambda for m=1,t=0
        lambda01 <- -dnorm(-mu.1.boot) / pnorm(-mu.1.boot) #lambda for m=0,t=1
        lambda00 <- -dnorm(-mu.0.boot) / pnorm(-mu.0.boot) #lambda for m=0,t=0
        
        # Step 1-3: Define lambda function
        lambda <- function(mmodel, mcoef) {
            mu <- model.matrix(mmodel) %*% mcoef
            m <- mmodel$y #this is M
            return((m*dnorm(-mu)-(1-m)*dnorm(-mu))/(m*pnorm(mu)+(1-m)*pnorm(-mu)))
        }

        # Step 2: Rho loop
        ## Step 2-0: Initialize containers
        d0 <- d1 <- matrix(NA, length(rho), 1)
        upper.d0 <- upper.d1 <- lower.d0 <- lower.d1 <- matrix(NA, length(rho), 1)
        ind.d0 <- ind.d1 <- matrix(NA, length(rho), 1)
        Ymodel.coef.boot <- matrix(NA, sims, y.k)
        colnames(Ymodel.coef.boot) <- names(model.y$coef)
        sigma.3.boot <- rep(NA, sims)
        d0.boot <- d1.boot <- rep(NA, sims)
        ## START OF RHO LOOP
        for(i in 1:length(rho)){
            
            ## START OF BOOTSTRAP LOOP
            for(k in 1:sims){
            ## Step 2-1: Obtain the initial Y model with the correction term
            adj <- lambda(model.m, Mmodel.coef.boot[k,]) * rho[i] # the adjustment term
            w <- 1 - rho[i]^2*lambda(model.m, Mmodel.coef.boot[k,])*(lambda(model.m, Mmodel.coef.boot[k,]) + mu.boot[,k])
            y.t.data.adj <- data.frame(y.t.data, w, adj)
            model.y.adj <- update(model.y, as.formula(paste(". ~ . + adj")), weights=w, data=y.t.data.adj)
            sigma.3 <- summary(model.y.adj)$sigma
            
            ## Step 2-2: Update the Y model via Iterative FGLS
            #eps <- .Machine$double.eps
            sigma.dif <- 1
            while(abs(sigma.dif) > eps){
                Y.star <- Y.value - sigma.3 * adj
                y.t.data.star <- data.frame(Y.star, y.t.data.adj)
                model.y.update <- update(model.y, as.formula(paste("Y.star ~ .")), weights=w, data=y.t.data.star)
                sigma.3.temp <- summary(model.y.update)$sigma
                sigma.dif <- sigma.3.temp - sigma.3
                sigma.3 <- sigma.3.temp
            }
            
            ## Step 2-3: Bootstrap Y model parameters
            Ymodel.coef <- model.y.update$coef
            Ymodel.var.cov <- vcov(model.y.update)
            Ymodel.coef.boot[k,] <- mvrnorm(1, mu=Ymodel.coef, Sigma=Ymodel.var.cov) #draw one bootstrap sample of Y-model parameters for each k
            #sig3.shape <- model.y.update$df/2
            #sig3.invscale <- (model.y.update$df/2) * sigma.3^2
            #sigma.3.boot[k] <- sqrt(1 / rgamma(1, shape = sig3.shape, scale = 1/sig3.invscale)) #one sample of sigma.3 via inverse-gamma posterior
            
            ## Step 2-4: Bootstrap ACMEs; means are over observations
            d0.boot[k] <- mean( (Ymodel.coef.boot[k,M.out]) * (pnorm(mu.1.boot[,k]) - pnorm(mu.0.boot[,k])) )
            if(INT==TRUE){
                d1.boot[k] <- mean( (Ymodel.coef.boot[k,M.out] + Ymodel.coef.boot[k,TM.out]) * 
                    (pnorm(mu.1.boot[,k]) - pnorm(mu.0.boot[,k])) )
                } else {
                d1.boot[k] <- mean( (Ymodel.coef.boot[k,M.out]) * 
                    (pnorm(mu.1.boot[,k]) - pnorm(mu.0.boot[,k])) )
                }

            ## END OF BOOTSTAP LOOP
            }
            
        ## Step 2-5: Compute Outputs
        d0[i] <- mean(d0.boot)
        d1[i] <- mean(d1.boot)
        upper.d0[i] <- quantile(d0.boot, 0.975)
        upper.d1[i] <- quantile(d1.boot, 0.975)
        lower.d0[i] <- quantile(d0.boot, 0.025)
        lower.d1[i] <- quantile(d1.boot, 0.025)
        ind.d0[i] <- as.numeric(lower.d0[i] < 0 & upper.d0[i] > 0)
        ind.d1[i] <- as.numeric(lower.d1[i] < 0 & upper.d1[i] > 0)
        
        ## END OF RHO LOOP
        }
        if(INT==TRUE){
        ii <- which(abs(d0-0)==min(abs(d0-0)))
        kk <- which(abs(d1-0)==min(abs(d1-0)))
		err.cr.1 <- rho[ii]
		err.cr.2 <- rho[kk]
		err.cr <- c(err.cr.1, err.cr.2)
		} else {
		ii <- which(abs(d0-0)==min(abs(d0-0)))
		err.cr <- rho[ii]	
			}
        type <- "bm"
        # Step 3: Output
        out <- list(rho = rho, d0=d0, d1=d1, upper.d0=upper.d0, lower.d0=lower.d0, upper.d1=upper.d1, lower.d1=lower.d1, ind.d0=ind.d0, ind.d1=ind.d1, INT=INT, DETAIL=DETAIL, sims=sims,tau=NULL, upper.tau=NULL, lower.tau=NULL, nu=NULL, upper.nu=NULL, lower.nu=NULL,type=type, err.cr=err.cr)
        class(out) <- "medsens"
        out

    ## END OF CASE 2: Continuous Outcome + Binary Mediator
    } else
    
    #########################################################
    ## CASE 3: Binary Outcome + Continuous Mediator
    #########################################################
    if(class(model.y)[1]=="glm" & class(model.m)[1]=="lm") {
        	if(INT==TRUE){
		stop("Sensitivity Analysis Not Available Binary Mediator With Interactions \n")
		}
        # Step 0: Setting Variable labels
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
        
        if(INT==TRUE){
        TM <- paste(T,M, sep=":")
        TM.out <- paste(T.out,M.out, sep=":")
            }
        
        # Step 1: Obtain Model Parameters
        ## Step 1-1: Bootstrap M model parameters
        Mmodel.coef <- model.m$coef
        m.k <- length(model.m$coef)
        Mmodel.var.cov <- vcov(model.m)
        Mmodel.coef.boot <- mvrnorm(sims, mu=Mmodel.coef, Sigma=Mmodel.var.cov)
        if(is.factor(y.t.data[,T])==TRUE){
         beta2.boot <- Mmodel.coef.boot[,T.out] 
            } else {
          beta2.boot <- Mmodel.coef.boot[,T]    
                }

        sigma.2 <- summary(model.m)$sigma
        sig2.shape <- model.m$df/2
        sig2.invscale <- (model.m$df/2) * sigma.2^2
        sigma.2.boot <- sqrt(1 / rgamma(sims, shape = sig2.shape, scale = 1/sig2.invscale))
        
        ## Step 1-2: Bootstrap Y model parameters
        Ymodel.coef <- model.y$coef
        Ymodel.var.cov <- vcov(model.y)
        y.k <- length(Ymodel.coef)
        Ymodel.coef.boot <- mvrnorm(sims, mu=Ymodel.coef, Sigma=Ymodel.var.cov)
        colnames(Ymodel.coef.boot) <- names(Ymodel.coef)
        gamma.tilde <- Ymodel.coef.boot[,M.out]

        # Step 2: Compute ACME via the procedure in IKT
        ## Step 2-1: Estimate Error Correlation from inconsistent estimate of Y on M
        rho12.boot <- (sigma.2.boot * gamma.tilde) / (1 + sqrt(sigma.2.boot^2*gamma.tilde^2))
        
        ## Step 2-2: Calculate alpha_1, beta_1 and xi_1
        YTmodel.coef.boot <- Ymodel.coef.boot[,!colnames(Ymodel.coef.boot)%in%M.out] * sqrt(1-rho12.boot^2) %x% t(rep(1,y.k-1)) + Mmodel.coef.boot * (rho12.boot/sigma.2.boot) %x% t(rep(1,y.k-1))
        
        ## Step 2-3: Calculate Gamma
        ## Data matrices for the Y model less M
        y.mat.1 <- model.matrix(model.y)[,!colnames(model.matrix(model.y))%in%M.out]
        y.mat.1[,T.out] <- 1
        y.mat.0 <- model.matrix(model.y)[,!colnames(model.matrix(model.y))%in%M.out]
        y.mat.0[,T.out] <- 0
        
        ## Initialize objects before the rho loop
        d0 <- d1 <- rep(NA, length(rho))
        upper.d0 <- upper.d1 <- lower.d0 <- lower.d1 <- rep(NA, length(rho))
        ind.d0 <- ind.d1 <- rep(NA, length(rho))
        tau <- nu <- rep(NA, length(rho))
        upper.tau <- upper.nu <- lower.tau <- lower.nu <- rep(NA, length(rho))
        d0.boot <- d1.boot <- rep(NA, sims)
        tau.boot <- nu.boot <- rep(NA, sims)
        
        ## START OF RHO LOOP
        for(i in 1:length(rho)){
            gamma.boot <- (-rho[i] + rho12.boot*sqrt((1-rho[i]^2)/(1-rho12.boot^2)))/sigma.2.boot
            for(k in 1:sims){
            d0.boot[k] <- mean( pnorm(y.mat.0 %*% YTmodel.coef.boot[k,] + 
                gamma.boot[k]*beta2.boot[k]/sqrt(gamma.boot[k]^2*sigma.2.boot[k]^2+2*gamma.boot[k]*rho[i]*sigma.2.boot[k]+1)) -
                pnorm(y.mat.0 %*% YTmodel.coef.boot[k,]) )
            d1.boot[k] <- mean( pnorm(y.mat.1 %*% YTmodel.coef.boot[k,]) - pnorm(y.mat.1 %*% YTmodel.coef.boot[k,] - 
                gamma.boot[k]*beta2.boot[k]/sqrt(gamma.boot[k]^2*sigma.2.boot[k]^2+2*gamma.boot[k]*rho[i]*sigma.2.boot[k]+1)) )
            tau.boot[k] <- mean( pnorm(y.mat.1 %*% YTmodel.coef.boot[k,]) - pnorm(y.mat.0 %*% YTmodel.coef.boot[k,]) )
            nu.boot[k] <- (d0.boot[k] + d1.boot[k])/(2*tau.boot[k])
            }
            
        ## Step 2-4: Compute Outputs
        d0[i] <- mean(d0.boot) # ACME(t=0)
        d1[i] <- mean(d1.boot) # ACME(t=1)
        upper.d0[i] <- quantile(d0.boot, 0.975)
        upper.d1[i] <- quantile(d1.boot, 0.975)
        lower.d0[i] <- quantile(d0.boot, 0.025)
        lower.d1[i] <- quantile(d1.boot, 0.025)
        tau[i] <- mean(tau.boot) # ATE
        nu[i] <- mean(nu.boot) # Proportion Mediated
        upper.tau[i] <- quantile(tau.boot, 0.975)
        upper.nu[i] <- quantile(nu.boot, 0.975)
        lower.tau[i] <- quantile(tau.boot, 0.025)
        lower.nu[i] <- quantile(nu.boot, 0.025)
        ind.d0[i] <- as.numeric(lower.d0[i] < 0 & upper.d0[i] > 0)
        ind.d1[i] <- as.numeric(lower.d1[i] < 0 & upper.d1[i] > 0)
        
        ## END OF RHO LOOP
        }
        if(INT==TRUE){
        ii <- which(abs(d0-0)==min(abs(d0-0)))
        kk <- which(abs(d1-0)==min(abs(d1-0)))
		err.cr.1 <- rho[ii]
		err.cr.2 <- rho[kk]
		err.cr <- c(err.cr.1, err.cr.2)
		} else {
		ii <- which(abs(d0-0)==min(abs(d0-0)))
		err.cr <- rho[ii]	
			}
        type <- "bo"
        ## Step 3: Output
        err.cr <- mean(rho12.boot) # Rho_12 estimate
        out <- list(rho = rho, err.cr=err.cr, d0=d0, d1=d1, upper.d0=upper.d0, lower.d0=lower.d0, 
            upper.d1=upper.d1, lower.d1=lower.d1, ind.d0=ind.d0, ind.d1=ind.d1, 
            tau=tau, upper.tau=upper.tau, lower.tau=lower.tau, nu=nu, upper.nu=upper.nu, lower.nu=lower.nu,
            INT=INT, DETAIL=DETAIL, sims=sims, type=type)
        class(out) <- "medsens"
        out
    ## END OF CASE 3: Binary Outcome + Continuous Mediator    
    }

## END OF SENSITIVITY FUNCTION
}



print.medsens <- function(x, ...){
    print(unlist(x[1:16]))
    invisible(x)
    }

summary.medsens <- function(object, ...)
    structure(object, class = c("summary.medsens", class(object)))
 
print.summary.medsens <- function(x, ...){
 if(x$type=="ct"){
        if(x$INT==FALSE){
            tab <- cbind(x$rho, round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            if(sum(x$ind.d0)==1){
    tab <- as.matrix(tab[x$ind.d0==1, -5])
    tab <- t(tab)
    } else {
    tab <- tab[x$ind.d0==1, -5] 
        }
            colnames(tab) <-  c("Rho", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region\n\n")
            print(tab)
            cat("\nRho at which Delta = 0:", round(x$err.cr, 4), "\n\n")
            invisible(x)    
                } else { #Interaction Tables Start Here
            tab.d0 <- cbind(x$rho, round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            if(sum(x$ind.d0)==1){
            tab.d0 <- as.matrix(tab.d0[x$ind.d0==1, -5])
            tab.d0 <- t(tab.d0)
            } else {
            tab.d0 <- tab.d0[x$ind.d0==1, -5] 
            }
            colnames(tab.d0) <-  c("Rho", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d0) <- NULL
            tab.d1 <- cbind(x$rho, round(x$d1,4), round(x$lower.d1,4), round(x$upper.d1, 4), x$ind.d1)
            if(sum(x$ind.d1)==1){
            tab.d1 <- as.matrix(tab.d1[x$ind.d1==1, -5])
            tab.d1 <- t(tab.d1)
            } else {
            tab.d1 <- tab.d1[x$ind.d1==1, -5] 
            }
            colnames(tab.d1) <-  c("Rho", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d1) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region: d0\n\n")
            print(tab.d0)
            cat("\nRho at which Delta_0 = 0:", round(x$err.cr[1], 4), "\n\n")
            cat("\nSensitivity Region: d1\n\n")
            print(tab.d1)
            cat("\nRho at which Delta_1 = 0:", round(x$err.cr[2], 4), "\n\n")
            invisible(x)        
            }
    } else if(x$type=="bm") {
        if(x$INT==FALSE){
            tab <- cbind(x$rho, round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            if(sum(x$ind.d0)==1){
            tab <- as.matrix(tab[x$ind.d0==1, -5])
            tab <- t(tab)
            } else {
            tab <- tab[x$ind.d0==1, -5] 
            }
            colnames(tab) <-  c("Rho", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region\n\n")
            print(tab)
            cat("\nRho at which Delta = 0:", round(x$err.cr, 4), "\n\n")
            invisible(x)    
                } else {
            tab.d0 <- cbind(x$rho, round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            if(sum(x$ind.d0)==1){
            tab.d0 <- as.matrix(tab.d0[x$ind.d0==1, -5])
            tab.d0 <- t(tab.d0)
            } else {
            tab.d0 <- tab.d0[x$ind.d0==1, -5] 
            }
            colnames(tab.d0) <-  c("Rho","Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d0) <- NULL
            tab.d1 <- cbind(x$rho, round(x$d1,4), round(x$lower.d1,4), round(x$upper.d1, 4), x$ind.d1)
                if(sum(x$ind.d1)==1){
                tab.d1 <- as.matrix(tab.d1[x$ind.d1==1, -5])
                tab.d1 <- t(tab.d1)
                } else {
            tab.d1 <- tab.d1[x$ind.d1==1, -5] 
            }
            colnames(tab.d1) <-  c("Rho","Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d1) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region: d0\n\n")
            print(tab.d0)
            cat("\nRho at which Delta_0 = 0:", round(x$err.cr[1], 4), "\n\n")
            cat("\nSensitivity Region: d1\n\n")
            print(tab.d1)
            cat("\nRho at which Delta_1 = 0:", round(x$err.cr[2], 4), "\n\n")
            invisible(x)        
            }
        } else if(x$type=="bo") {
        if(x$INT==FALSE){
            tab <- cbind(x$rho, round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            if(sum(x$ind.d0)==1){
            tab <- as.matrix(tab[x$ind.d0==1, -5])
            tab <- t(tab)
            } else {
            tab <- tab[x$ind.d0==1, -5] 
            }
            colnames(tab) <-  c("Rho", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region\n\n")
            print(tab)
            cat("\nRho at which Delta = 0:", round(x$err.cr, 4), "\n\n")
            invisible(x)    
                } else {
            tab.d0 <- cbind(x$rho, round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            if(sum(x$ind.d0)==1){
            tab <- as.matrix(tab[x$ind.d0==1, -5])
            tab <- t(tab)
            } else {
            tab <- tab[x$ind.d0==1, -5] 
            }
            colnames(tab.d0) <-  c("Rho","Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d0) <- NULL
            tab.d1 <- cbind(x$rho, round(x$d1,4), round(x$lower.d1,4), round(x$upper.d1, 4), x$ind.d1)
            if(sum(x$ind.d1)==1){
            tab <- as.matrix(tab[x$ind.d1==1, -5])
            tab <- t(tab)
            } else {
            tab <- tab[x$ind.d1==1, -5] 
            }
            colnames(tab.d1) <-  c("Rho","Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d1) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region: d0\n\n")
            print(tab.d0)
            cat("\nRho at which Delta_0 = 0:", round(x$err.cr[1], 4), "\n\n")
            cat("\nSensitivity Region: d1\n\n")
            print(tab.d1)
            cat("\nRho at which Delta_1 = 0:", round(x$err.cr[2], 4), "\n\n")
            invisible(x)        
            }
     }
        
}

plot.medsens <- function(x, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL, pr.plot=FALSE,...){
    if(pr.plot==TRUE){
        plot(x$rho, x$nu, type="n", xlab="", ylab = "", main=main, xlim=xlim, ylim=ylim)
        polygon(x=c(x$rho, rev(x$rho)), y=c(x$lower.nu, rev(x$upper.nu)), border=FALSE, col=8, lty=2)
        lines(x$rho, x$d0, lty=1)
        abline(h=0)
        abline(v=0)
        title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
        title(ylab = expression(paste("Proportion Mediated: ", bar(delta))), cex.lab=.9)
        } else {
   if(x$INT==FALSE){
        plot(x$rho, x$d0, type="n", xlab="", ylab = "", main=main, xlim=xlim, ylim=ylim)
        polygon(x=c(x$rho, rev(x$rho)), y=c(x$lower.d0, rev(x$upper.d0)), border=FALSE, col=8, lty=2)
        lines(x$rho, x$d0, lty=1)
        abline(h=0)
        abline(v=0)
        title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
        title(ylab = expression(paste("Average Mediation Effect: ", bar(delta)(t))), cex.lab=.9)
        if(x$DETAIL==TRUE){
        abline(h=x$d0[96], lty=2)
            } else {
        abline(h=x$d0[10], lty=2)       
                }
        } else {
        par(mfrow=c(1,2))
        plot(x$rho, x$d0, type="n", xlab="", ylab = "", ylim = c(-.2,.2))
        polygon(x=c(x$rho, rev(x$rho)), y=c(x$lower.d0, rev(x$upper.d0)), border=FALSE, col=8, lty=2)
        lines(x$rho, x$d0, lty=1)
        abline(h=0)
        if(x$DETAIL==TRUE){
        abline(h=x$d0[96], lty=2)
            } else {
        abline(h=x$d0[10], lty=2)       
                }
        abline(v=0)
        title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
        title(ylab = expression(paste("Average Mediation Effect: ", bar(delta[0]))), cex.lab=.9)

        #Delta_1
        plot(x$rho, x$d1, type="n", xlab="", ylab = "", ylim = c(-.2,.2))
        polygon(x=c(x$rho, rev(x$rho)), y=c(x$lower.d1, rev(x$upper.d1)), border=FALSE, col=8, lty=2)
        lines(x$rho, x$d1, lty=1)
        abline(h=0)
                if(x$DETAIL==TRUE){
        abline(h=x$d1[96], lty=2)
            } else {
        abline(h=x$d1[10], lty=2)       
                }
        abline(v=0)
        title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
        title(ylab = expression(paste("Average Mediation Effect: ", bar(delta[1]))), cex.lab=.9)
            }
        }
    }
