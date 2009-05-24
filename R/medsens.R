medsens <- function(z, ...){ 
    UseMethod("medsens", z)
    }

medsens.default <- function(z, model.y, T="treat.name", M="med.name", INT=FALSE, DETAIL=TRUE, nboot=1000)
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
    if(class(model.y)=="lm" & class(model.m)=="lm") {

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
        eps <- .001 #.Machine$double.eps
        
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
        
        if(INT==TRUE){
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
        
        out <- list(rho = rho, err.cr=err.cr, d0=d0, d1=d1, upper.d0=upper.d0, lower.d0=lower.d0, upper.d1=upper.d1, lower.d1=lower.d1, ind.d0=ind.d0, ind.d1=ind.d1, INT=INT, DETAIL=DETAIL, nboot=nboot)
        class(out) <- "sens.c"
        out

    ## END OF CASE 1: Continuous Outcome + Continuous Mediator    
    } else

    #########################################################
    ## CASE 2: Continuous Outcome + Binary Mediator
    #########################################################
    if(class(model.y)=="lm" & class(model.m)=="glm") {

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
        Mmodel.coef.boot <- mvrnorm(nboot, mu=Mmodel.coef, Sigma=Mmodel.var.cov)
        
        # Step 1-2: Bootstrap lambda_0 and lambda_1; lambdas are (n x nboot) matrix
        m.mat.1 <- model.matrix(model.m)
        m.mat.1[,T.out] <- 1
        m.mat.0 <- model.matrix(model.m)
        m.mat.0[,T.out] <- 0
        mu.1.boot <- m.mat.1 %*% t(Mmodel.coef.boot) 
        mu.0.boot <- m.mat.0 %*% t(Mmodel.coef.boot) 
        lambda11 <- apply(-mu.1,2,dnorm) / apply(mu.1,2,pnorm) #m=1,t=1
        lambda10 <- apply(-mu.0,2,dnorm) / apply(mu.0,2,pnorm) #m=1,t=0
        lambda01 <- -apply(-mu.1,2,dnorm) / apply(-mu.1,2,pnorm) #m=0,t=1
        lambda00 <- -apply(-mu.0,2,dnorm) / apply(-mu.0,2,pnorm) #m=0,t=0
        
        # Step 1-3: Define lambda functions
        lambda <- function(mmodel, mcoef) {
            mu <- model.matrix(mmodel) %*% mcoef
            m <- mmodel$y
            return((m*dnorm(-mu)-(1-m)*dnorm(-mu))/(m*pnorm(mu)+(1-m)*pnorm(-mu)))
        }

        # Step 2: Rho loop
        ## Step 2-0: Initialize containers
        d0 <- d1 <- matrix(NA, length(rho), 1)
        upper.d0 <- upper.d1 <- lower.d0 <- lower.d1 <- matrix(NA, length(rho), 1)
        ind.d0 <- ind.d1 <- matrix(NA, length(rho), 1)
        Ymodel.coef.boot <- matrix(NA, nboot, y.k)
        sigma.3.boot <- rep(NA, nboot)
        d0.boot <- d1.boot <- rep(NA, nboot)
        
        ## START OF RHO LOOP
        for(i in 1:length(rho)){
        
            ## START OF BOOTSTRAP LOOP
            for(k in 1:nboot){
            
            ## Step 2-1: Obtain the initial Y model with the correction term
            adj <- lambda(mmodel, Mmodel.coef.boot[k]) * rho[i]
            model.y.adj <- update(model.y, as.formula(paste(". ~ . + adj")))
            sigma.3 <- summary(model.y.adj)$sigma
            
            ## Step 2-2: Update the Y model via Iterative OLS
            eps <- .001
            sigma.dif <- 1
            while(abs(sigma.dif) > eps){
                Y.star <- Y.value - sigma.3 * adj
                model.y.update <- update(model.y, as.formula(paste("Y.star ~ .")))
                sigma.3.temp <- summary(model.y.update)$sigma
                sigma.dif <- sigma.3.temp - sigma.3
                sigma.3 <- sigma.3.temp
            }
            
            ## Step 2-3: Bootstrap Y model parameters
            Ymodel.coef <- y.update$coef
            Ymodel.var.cov <- vcov(y.update)
            Ymodel.coef.boot[k,] <- mvrnorm(1, mu=Ymodel.coef, Sigma=Ymodel.var.cov)
            sig3.shape <- y.update$df/2
            sig3.invscale <- (y.update$df/2) * sigma.3^2
            sigma.3.boot[k] <- sqrt(1 / rgamma(1, shape = sig3.shape, scale = 1/sig3.invscale))
            
            ## Step 2-4: Bootstrap ACMEs
            d0.boot[k] <- mean( (Ymodel.coef.boot[k,M.out] + rho[i]*sigma.3.boot[k]*(lambda10[,k] - lambda00[,k])) * 
                (dnorm(mu.1.boot[,k]) - dnorm(mu.0.boot[,k])) )
            if(INT==TRUE){
                d1.boot[k] <- mean( (Ymodel.coef.boot[k,M.out] + Ymodel.coef.boot[k,TM.out] + rho[i]*sigma.3.boot[k]*(lambda11[,k] - lambda01[,k])) *
                    (dnorm(mu.1,boot[,k]) - dnorm(mu.0.boot[,k])) )
                } else {
                d1.boot[k] <- mean( (Ymodel.coef.boot[k,M.out] + rho[i]*sigma.3.boot[k]*(lambda11[,k] - lambda01[,k])) *
                    (dnorm(mu.1,boot[,k]) - dnorm(mu.0.boot[,k])) )
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
        
        # Step 3: Output
        err.cr <- NULL
        out <- list(rho = rho, err.cr=err.cr, d0=d0, d1=d1, upper.d0=upper.d0, lower.d0=lower.d0, upper.d1=upper.d1, lower.d1=lower.d1, ind.d0=ind.d0, ind.d1=ind.d1, INT=INT, DETAIL=DETAIL, nboot=nboot)
        class(out) <- "sens.bm"
        out
        
    ## END OF CASE 2: Continuous Outcome + Binary Mediator    
    }
    
    #########################################################
    ## CASE 3: Binary Outcome + Continuous Mediator
    #########################################################
    if(class(model.y)=="glm" & class(model.m)=="lm") {
        nboot<-1000 #Luke, for the dichotomous outcome case the nboot needs to be set. I do it within the function here but this be user specified so it can match the number of nboot used with mediate.R

        Mmodel.coef <- model.m$coef
        m.k <- length(model.m$coef)
        Mmodel.var.cov <- vcov(model.m)
        mdraws <- mvrnorm(nboot, mu=Mmodel.coef, Sigma=Mmodel.var.cov)
        alpha2.est <- mdraws[,1]
        beta2.est <- mdraws[,2]
        xi2.est <- mdraws[,3:m.k]
        sigma <- rep(summary(model.m)$sigma, nboot)
        sigma.sq <- sigma^2
        
        #General xi.2 Quantity
        m.data <- model.frame(model.m)
        mmat <- model.matrix(model.m, m.data)
        mmat[,1] <- 0
        mmat[,2] <- 0
        xi2.X <-  mdraws %*% t(mmat)

        Ymodel.coef <- model.y$coef
        Ymodel.var.cov <- vcov(model.y)
        y.k <- length(Ymodel.coef)
        ydraws <- mvrnorm(nboot, mu=Ymodel.coef, Sigma=Ymodel.var.cov)
        alpha3.tilde <- ydraws[,1]
        beta3.tilde <- ydraws[,T] #Luke-will it know to use T and M?
        gamma.tilde <- ydraws[,M]
        xi3.tilde <- ydraws[,4:y.k] #Dimension by extra covariates or use predict or what i use in binary code.
            
        ###############################################
        # Calculations Start Here
        #Step 2 - Estimate Error Correlation from inconsistent estimate of Y on M
        rho12 <- (sigma*gamma.tilde)/(sqrt(sigma.sq*gamma.tilde^2+1))
        
        #Step 2 - Calculate Alpha_1
        alpha.1 <- alpha3.tilde*sqrt(1 - rho12^2) + (alpha2.est*rho12)/sigma
        
        #Step 3 - Calculate Beta_1
        beta.1 <- beta3.tilde*sqrt(1 - rho12^2) + (beta2.est*rho12)/sigma
        
        #calculate xi1
        xi.1 <- (xi3.tilde*sqrt(1 - rho12^2)) + ((xi2.est*rho12)/sigma)
        
            #Luke-this appears to require that the first variable in the model for the outcome is T. We want to avoid this.                      

        y.data <- model.frame(model.y)
        ymat <- model.matrix(model.y, y.data)
        ymat[,1] <- 0
        ymat[,2] <- 0
        xi3.X <- ydraws %*% t(ymat)
        
        tau <- matrix(, nrow=nboot, ncol=length(rho))
        #d0 <- matrix(, nrow=nboot, ncol=length(rho))
        #d1 <- matrix(, nrow=nboot, ncol=length(rho))
        d0.tmp <- matrix(, nrow=nboot, ncol=length(rho))
        d1.tmp <- matrix(, nrow=nboot, ncol=length(rho))                    
        d.sum <- matrix(, nrow=nboot, ncol=length(rho))

        #d0 <- matrix(NA, length(rho), 1)
        #d1 <- matrix(NA, length(rho), 1)
        #tau<-matrix(NA, length(rho), 1)
    
    
        #Luke -why do we have to set i here?
        i <- 1
            #Teppei-I think this is wrong...we need to loop over nboot too..right?
            #Also, I think this messes up the way we are using the covars, bc xi2.X is quite what we want. Is the fix to also do xi2.X[k]?
    
        for(i in 1:length(rho)){
            #Step 4 - Calculate Gamma
            gamma <- (-rho[i] + rho12*sqrt((1-rho[i]^2)/(1-rho12^2)))/sigma
            
            #Step 5 - Alpha_3
            alpha3 <- alpha.1*sqrt(gamma^2*sigma.sq + 2*gamma*rho[i]*sigma + 1) - gamma*alpha2.est
            
            #Step 6 - Beta_3
            beta3 <- beta.1*sqrt(gamma^2*sigma.sq + 2*gamma*rho[i]*sigma + 1) - gamma*beta2.est
            
            t <- 0
            d0.covartemp <- pnorm((alpha3 + beta3*t + xi3.X + gamma*(alpha2.est + beta2.est + xi2.X))/sqrt(gamma^2*sigma.sq+2*gamma*rho[i]*sigma+1)) -pnorm((alpha3 + beta3*t + xi3.X + gamma*(alpha2.est + xi2.X))/sqrt(gamma^2*sigma.sq+2*gamma*rho[i]*sigma+1))
            
            t <- 1
            d1.covartemp <- pnorm((alpha3 + beta3*t + xi3.X + gamma*(alpha2.est + beta2.est + xi2.X))/sqrt(gamma^2*sigma.sq+2*gamma*rho[i]*sigma+1)) -pnorm((alpha3 + beta3*t + xi3.X + gamma*(alpha2.est + xi2.X))/sqrt(gamma^2*sigma.sq+2*gamma*rho[i]*sigma+1))
            
            #d0[,i] <- d0.tmp <- apply(d0.covartemp, 1, mean)
            #d1[,i] <- d1.tmp <- apply(d1.covartemp, 1, mean)
            d0.tmp[,i] <- apply(d0.covartemp, 1, mean)
            d1.tmp[,i] <- apply(d1.covartemp, 1, mean)
        }
                            
    #This calculates various quantities of interest that are required for the summary and plot functions
        tau <- pnorm(alpha.1 + beta.1) - pnorm(alpha.1)
        
        #d0.ci <- apply(d0, 2, quantile, probs=c(0.025, 0.975),na.rm=TRUE)
        #d1.ci <- apply(d1, 2, quantile, probs=c(0.025, 0.975),na.rm=TRUE)
        #d.avg <- (d0 + d1)/2
        #d.sum <- d.avg / tau       
        #pr.ci <- apply(d.sum, 2, quantile, probs=c(0.025, 0.975),na.rm=TRUE)

        d0.ci <- apply(d0.tmp, 2, quantile, probs=c(0.025, 0.975),na.rm=TRUE)
        d1.ci <- apply(d1.tmp, 2, quantile, probs=c(0.025, 0.975),na.rm=TRUE)
        d.avg <- (d0.tmp + d1.tmp)/2
        d.sum <- d.avg / tau       
        pr.ci <- apply(d.sum, 2, quantile, probs=c(0.025, 0.975),na.rm=TRUE)

        

        pr.med <- apply(d.sum, 2, mean)
        pr.ci <- apply(d.sum, 2, quantile, probs=c(0.025, 0.975))

        #d0 <- apply(d0.tmp,2,mean)
        #d1 <- apply(d1.tmp,2,mean) 

        d0 <- apply(d0.tmp,2,mean)
        d1 <- apply(d1.tmp,2,mean)

        if(INT==TRUE){
            upper.d0 <- d0.ci[2,]
            lower.d0 <- d0.ci[1,]
            upper.d1 <- NULL
            lower.d1 <- NULL
            ind.d0 <- as.numeric(lower.d0 < 0 & upper.d0 > 0)
            ind.d1 <- NULL
            upper.pr<-pr.ci[2,]
            lower.pr<-pr.ci[1,]
            } else {
                upper.d0 <-  d0.ci[2,]
                lower.d0 <-  d0.ci[1,]
                upper.d1 <-  d1.ci[2,]
                lower.d1 <-  d1.ci[1,]   
                ind.d0 <- as.numeric(lower.d0 < 0 & upper.d0 > 0)
                ind.d1 <- as.numeric(lower.d1 < 0 & upper.d1 > 0)
                upper.pr<-pr.ci[2,]
                lower.pr<-pr.ci[1,]
                }

        #Luke -I have to create err.cr to get the summary to work. Can't figure out why.            
        err.cr <- matrix(1, length(rho), 1)


        out <- list(rho = rho, err.cr=err.cr, d0=d0, d1=d1, upper.d0=upper.d0, lower.d0=lower.d0, upper.d1=upper.d1, lower.d1=lower.d1, ind.d0=ind.d0, ind.d1=ind.d1, INT=INT, DETAIL=DETAIL)
        #out <- list(rho = rho, d0=d0, d1=d1, upper.d0=upper.d0, lower.d0=lower.d0, upper.d1=upper.d1, lower.d1=lower.d1, ind.d0=ind.d0, ind.d1=ind.d1, INT=INT, DETAIL=DETAIL)

        class(out) <- "sens.c"
        out
    ## END OF CASE 3: Binary Outcome + Continuous Mediator    
    }

## END OF SENSITIVITY FUNCTION
}



print.sens.c <- function(x, ...){
    print(unlist(x[1:16]))
    invisible(x)
    }

summary.sens.c <- function(object)
    structure(object, class = c("sum.sens.c", class(object)))
 
print.sum.sens.c <- function(x, ...){
     if(class(model.y)=="lm") {
        if(x$INT==FALSE){
            tab <- cbind(x$rho, round(x$err.cr,4), round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            tab <- tab[x$ind.d0==1, -6]
            colnames(tab) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region\n\n")
            print(tab)
            invisible(x)    
                } else {
            tab.d0 <- cbind(x$rho, round(x$err.cr,4), round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            tab.d0 <- tab.d0[x$ind==1, -6]
            colnames(tab.d0) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d0) <- NULL
            tab.d1 <- cbind(x$rho, round(x$err.cr,4), round(x$d1,4), round(x$lower.d1,4), round(x$upper.d1, 4), x$ind.d1)
            tab.d1 <- tab[x$ind.d1==1, -6]
            colnames(tab.d1) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d1) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region: d0\n\n")
            print(tab.d0)
            cat("\nSensitivity Region: d1\n\n")
            print(tab.d1)
            invisible(x)        
            }
      } else
        
     if(class(model.y)=="glm") {
#Luke -for some reason I can only get the summary command to work by creating a vector for err.cr in the binary function above, and then include it in the printout below
#Not sure why, can you try to fix this?
        if(x$INT==FALSE){
            tab <- cbind(x$rho, round(x$err.cr,4), round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            tab <- tab[x$ind.d0==1, -6]
            colnames(tab) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region\n\n")
            print(tab)
            invisible(x)    
                } else {
            tab.d0 <- cbind(x$rho, round(x$err.cr,4), round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
            tab.d0 <- tab.d0[x$ind==1, -6]
            colnames(tab.d0) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d0) <- NULL
            tab.d1 <- cbind(x$rho, round(x$err.cr,4), round(x$d1,4), round(x$lower.d1,4), round(x$upper.d1, 4), x$ind.d1)
            tab.d1 <- tab[x$ind.d1==1, -6]
            colnames(tab.d1) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
            rownames(tab.d1) <- NULL
            cat("\nMediation Sensitivity Analysis\n")
            cat("\nSensitivity Region: d0\n\n")
            print(tab.d0)
            cat("\nSensitivity Region: d1\n\n")
            print(tab.d1)
            invisible(x)        
            }
     }
        
}

plot.sens.c <- function(x, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL,...){
    if(x$INT==FALSE){
        plot(x$rho, x$d0, type="n", xlab="", ylab = "", main=main, xlim=xlim, ylim=ylim)
        polygon(x=c(x$rho, rev(x$rho)), y=c(x$lower.d0, rev(x$upper.d0)), border=FALSE, col=8, lty=2)
        lines(x$rho, x$d0, lty=1)
        abline(h=0)
        abline(v=0)
        title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
        title(ylab = expression(paste("Average Mediation Effect: ", bar(delta))), cex.lab=.9)
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
        lines(rho, d1, lty=1)
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
