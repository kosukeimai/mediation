mediate <- function(model.m, model.y, sims=1000, boot=FALSE, 
                    treat="treat.name", mediator="med.name", 
                    control=NULL, conf.level=.95,
                    control.value=0, treat.value=1, 
                    long=TRUE, dropobs=FALSE, robustSE = FALSE, ...){
    
    # Warn users who still use INT option
    if(match("INT", names(match.call()), 0L)){
        warning("'INT' is deprecated - existence of interaction term
        is now automatically detected")
    }
    
    # Warning for robustSE used with boot
    if(robustSE && boot){
        warning("'robustSE' is ignored for nonparametric bootstrap")
    }
    
    
    # Model type indicators
    isGam.y <- class(model.y)[1] == "gam"
    isGam.m <- class(model.m)[1] == "gam"
    isVglm.y <- class(model.y)[1] == "vglm"
    isRq.y <- class(model.y)[1] == "rq"
    isRq.m <- class(model.m)[1] == "rq"
    isOrdered.y <- class(model.y)[1] == "polr"
    isOrdered.m <- class(model.m)[1] == "polr"
    
    # Drop observations not common to both mediator and outcome models
    if(dropobs){
        odata.m <- model.frame(model.m)
        odata.y <- model.frame(model.y)
        newdata <- merge(odata.m, odata.y, sort=FALSE,
                    by=c("row.names", intersect(names(odata.m), names(odata.y))))
        rownames(newdata) <- newdata$Row.names
        newdata <- newdata[,-1L]
        rm(odata.m, odata.y)
                
        if(isS4(model.m)){
            call.m <- model.m@call
        } else {
            call.m <- model.m$call
        }
        if(isS4(model.y)){
            call.y <- model.y@call
        } else {
            call.y <- model.y$call
        }
        call.m$data <- call.y$data <- newdata
        if(c("(weights)") %in% names(newdata)){
            call.m$weights <- call.y$weights <- model.weights(newdata) 
        }
        model.m <- eval.parent(call.m)
        model.y <- eval.parent(call.y)
    }
            
    # Record class of model.m as "ClassM"
    if(isGam.m){
        ClassM <- class(model.m)[2]
    } else {
        ClassM <- class(model.m)[1]
    }
    
    # Record family of model.m
    if(ClassM=="glm" || isGam.m){
        FamilyM <- model.m$family$family
    }
    
    # Record class of model.y as "ClassY"
    if(isGam.y){
        ClassY <- class(model.y)[2]
    } else {
        ClassY <- class(model.y)[1]
        if(isVglm.y){
            vfamily <- model.y@family@vfamily
            # Indicates which vglm model (currently always tobit)
        }
    }
    
    # Warning for control option in non-GAM outcome models
    if(!is.null(control) && !isGam.y){
        warning("'control' is only used for GAM outcome models - ignored")
        control <- NULL
    }
    
    # Model frames for M and Y models
    m.data <- model.frame(model.m)  # Call.M$data
    y.data <- model.frame(model.y)  # Call.Y$data
    
    # Numbers of observations/variables/categories
    n.m <- length(m.data[,1])
    n <- n.y <- length(y.data[,1])
    if(n.m != n.y){
        stop("number of observations do not match between mediator and outcome models")
    }
    m <- length(sort(unique(model.frame(model.m)[,1])))
    m.min <- as.numeric(sort(model.frame(model.m)[,1])[1])
    
    # Extracting weights from models
    weights.m <- model.weights(m.data)
    weights.y <- model.weights(y.data)
    
    if(!is.null(weights.m) && ClassM == "glm" && FamilyM == "binomial"){
        message("weights taken as sampling weights, not total number of trials")
    }
    
    if(is.null(weights.m)){
        weights.m <- rep(1,nrow(m.data))
    }
    if(is.null(weights.y)){
        weights.y <- rep(1,nrow(y.data))
    }
    
    if(!all(weights.m == weights.y)) {
        stop("weights on outcome and mediator models not identical")
    } else {
        weights <- weights.m
    }
    
    # Factor treatment indicator
    isFactorT.m <- is.factor(m.data[,treat])
    isFactorT.y <- is.factor(y.data[,treat])
    
    if(isFactorT.y){
        # TODO: Allow non-binary factor treatment
        t.levels <- levels(y.data[,treat])
        if(sum(c(treat.value, control.value) %in% t.levels)){
            cat.0 <- control.value
            cat.1 <- treat.value
        } else {
            cat.0 <- t.levels[1]
            cat.1 <- t.levels[2]
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
    
    predictY.dataprep <- function(T.t,T.c,M.t,M.c){
        ## Prepare model matrix for outcome predictions; arguments 0 or 1
        ## Objects that must exist in surrounding environment:
        ##    y.data, cat.1, cat.0,
        ##    t.levels, m.levels, isFactorT.y, isFactorM,
        ##    PredictM1, PredictM0, treat, mediator, control,
        ##    boot, j (if !boot)
        
        pred.data.t <- pred.data.c <- y.data
        
        # Set treatment values
        cat.t <- ifelse(T.t, cat.1, cat.0)
        cat.c <- ifelse(T.c, cat.1, cat.0)
        cat.t.ctrl <- ifelse(T.t, cat.0, cat.1)
        cat.c.ctrl <- ifelse(T.c, cat.0, cat.1)
        if(isFactorT.y){
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
        if(boot){
            PredictM1.tmp <- PredictM1
            PredictM0.tmp <- PredictM0
        } else {
            PredictM1.tmp <- PredictM1[j,]
            PredictM0.tmp <- PredictM0[j,]
        }
        if(M.t){
            PredictMt <- PredictM1.tmp
        } else {
            PredictMt <- PredictM0.tmp
        }
        if(M.c){
            PredictMc <- PredictM1.tmp
        } else {
            PredictMc <- PredictM0.tmp
        }
        if(isFactorM) {
            pred.data.t[,mediator] <- factor(PredictMt, levels=1:m, labels=m.levels)
            pred.data.c[,mediator] <- factor(PredictMc, levels=1:m, labels=m.levels)
        } else {
            pred.data.t[,mediator] <- PredictMt
            pred.data.c[,mediator] <- PredictMc
        }
        return(list(t=pred.data.t, c=pred.data.c))
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
        if(boot == FALSE){
            # Error if gam outcome or quantile mediator
            if(isGam.m | isGam.y | isRq.m){
                stop("'boot' must be 'TRUE' for models used")
            }
            
            # Get mean and variance parameters for simulations
            MModel.coef <- model.m$coef
            if(ClassM=="polr"){  # TRUE if model.m is ordered
                if(is.null(model.m$Hess)){
                    cat("mediator model object does not contain 'Hessian';")
                }
                k <- length(model.m$coef)
                MModel.var.cov <- vcov(model.m)[(1:k),(1:k)]
            } else {
                if(robustSE){
                    MModel.var.cov <- vcovHC(model.m, ...)
                } else {
                    MModel.var.cov <- vcov(model.m)
                }
            }
            
            if(isS4(model.y)){
                TMmodel.coef <- model.y@coefficients
            } else {
                TMmodel.coef <- model.y$coef
            }
            if(isRq.y){
                TMmodel.var.cov <- summary(model.y, covariance=TRUE)$cov
            } else{
                if(robustSE){
                    TMmodel.var.cov <- vcovHC(model.y, ...)
                } else {
                    TMmodel.var.cov <- vcov(model.y)
                }
            }
            
            # Draw model coefficients from normal
            MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov)
            TMmodel <- mvrnorm(sims, mu=TMmodel.coef, Sigma=TMmodel.var.cov)
            
            #####################################
            #  Mediator Predictions
            #####################################
            pred.data.t <- m.data
            pred.data.c <- m.data
            
            if(isFactorT.m){
                pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
            } else {
                pred.data.t[,treat] <- cat.1
                pred.data.c[,treat] <- cat.0
            }
            
            mmat.t <- model.matrix(terms(model.m), data=pred.data.t)
            mmat.c <- model.matrix(terms(model.m), data=pred.data.c)
            
            ### Case I-1-a: GLM Mediator
            if(ClassM == "glm"){
                
                muM1 <- model.m$family$linkinv(MModel %*% t(mmat.t))
                muM0 <- model.m$family$linkinv(MModel %*% t(mmat.c))
                
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
            } else if(ClassM=="polr"){
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
                
                ystar_m1 <- MModel %*% t(mmat.t) 
                ystar_m0 <- MModel %*% t(mmat.c)
                
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
            } else if(ClassM=="lm"){
                sigma <- summary(model.m)$sigma
                error <- rnorm(sims*n, mean=0, sd=sigma)
                muM1 <- MModel %*% t(mmat.t)
                muM0 <- MModel %*% t(mmat.c)
                PredictM1 <- muM1 + matrix(error, nrow=sims)
                PredictM0 <- muM0 + matrix(error, nrow=sims)
                rm(error)
                
            } else {
                stop("mediator model is not yet implemented")
            }
            rm(mmat.t, mmat.c)
            
            #####################################
            ##  Outcome Predictions
            #####################################
            # ACME(1)
            Pr1 <- matrix(, nrow=n, ncol=sims)
            Pr0 <- matrix(, nrow=n, ncol=sims)
            
            for(j in 1:sims){
                pred.data <- predictY.dataprep(1,1,1,0)
                pred.data.t <- pred.data$t
                pred.data.c <- pred.data$c
                
                ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
                ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
                
                if(isVglm.y){
                    if(vfamily=="tobit") {
                        Pr1.tmp <- ymat.t %*% TMmodel[j,-2]
                        Pr0.tmp <- ymat.c %*% TMmodel[j,-2]
                        Pr1[,j] <- pmin(pmax(Pr1.tmp, model.y@misc$Lower), model.y@misc$Upper)
                        Pr0[,j] <- pmin(pmax(Pr0.tmp, model.y@misc$Lower), model.y@misc$Upper)
                    } else {
                        stop("outcome model is in unsupported vglm family")
                    }
                } else {
                    Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
                    Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
                }
                
                rm(ymat.t, ymat.c, pred.data.t, pred.data.c, pred.data)
            }
            
            if(ClassY=="glm"){
                Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                delta.1.tmp <- Pr1 - Pr0
            } else {
                delta.1.tmp <- Pr1 - Pr0
            }
            rm(Pr1, Pr0)
                    
            # ACME(0)
            Pr1 <- matrix(, nrow=n, ncol=sims)
            Pr0 <- matrix(, nrow=n, ncol=sims)
            
            for(j in 1:sims){
                pred.data <- predictY.dataprep(0,0,1,0)
                pred.data.t <- pred.data$t
                pred.data.c <- pred.data$c
                
                ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
                ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
                
                if(isVglm.y){
                    if(vfamily=="tobit") {
                        Pr1.tmp <- ymat.t %*% TMmodel[j,-2]
                        Pr0.tmp <- ymat.c %*% TMmodel[j,-2]
                        Pr1[,j] <- pmin(pmax(Pr1.tmp, model.y@misc$Lower), model.y@misc$Upper)
                        Pr0[,j] <- pmin(pmax(Pr0.tmp, model.y@misc$Lower), model.y@misc$Upper)
                    } else {
                        stop("outcome model is in unsupported vglm family")
                    }
                } else {
                        Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
                }
                rm(ymat.t, ymat.c, pred.data.t, pred.data.c, pred.data)
            }
            
            if(ClassY=="glm"){
                Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                delta.0.tmp <- Pr1 - Pr0
            } else {
                delta.0.tmp <- Pr1 - Pr0
            }
            rm(Pr1, Pr0)
            
            # DE(1)
            Pr1 <- matrix(, nrow=n, ncol=sims)
            Pr0 <- matrix(, nrow=n, ncol=sims)
            
            for(j in 1:sims){
                pred.data <- predictY.dataprep(1,0,1,1)
                pred.data.t <- pred.data$t
                pred.data.c <- pred.data$c
                
                ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
                ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
                
                if(isVglm.y){
                    if(vfamily=="tobit") {
                        Pr1.tmp <- ymat.t %*% TMmodel[j,-2]
                        Pr0.tmp <- ymat.c %*% TMmodel[j,-2]
                        Pr1[,j] <- pmin(pmax(Pr1.tmp, model.y@misc$Lower), model.y@misc$Upper)
                        Pr0[,j] <- pmin(pmax(Pr0.tmp, model.y@misc$Lower), model.y@misc$Upper)
                    } else {
                        stop("outcome model is in unsupported vglm family")
                    }
                } else {
                    Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
                    Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
                }
                rm(ymat.t, ymat.c, pred.data.t, pred.data.c, pred.data)
            }    
            
            if(ClassY=="glm"){
                Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                zeta.1.tmp <- Pr1 - Pr0
            } else {
                zeta.1.tmp <- Pr1 - Pr0
            }
            rm(Pr1, Pr0)
            
            # DE(0)
            Pr1 <- matrix(, nrow=n, ncol=sims)
            Pr0 <- matrix(, nrow=n, ncol=sims)
            
            for(j in 1:sims){
                pred.data <- predictY.dataprep(1,0,0,0)
                pred.data.t <- pred.data$t
                pred.data.c <- pred.data$c
                
                ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
                ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
                
                if(isVglm.y){
                    if(vfamily=="tobit") {
                        Pr1.tmp <- ymat.t %*% TMmodel[j,-2]
                        Pr0.tmp <- ymat.c %*% TMmodel[j,-2]
                        Pr1[,j] <- pmin(pmax(Pr1.tmp, model.y@misc$Lower), model.y@misc$Upper)
                        Pr0[,j] <- pmin(pmax(Pr0.tmp, model.y@misc$Lower), model.y@misc$Upper)
                    } else {
                        stop("outcome model is in unsupported vglm family")
                    }
                } else {        
                    Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
                    Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
                }
                rm(ymat.t, ymat.c, pred.data.t, pred.data.c, pred.data)
            }
            
            if(ClassY=="glm"){
                Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                zeta.0.tmp <- Pr1 - Pr0
            } else {
                zeta.0.tmp <- Pr1 - Pr0
            }
            
            rm(Pr1, Pr0, PredictM1, PredictM0, TMmodel, MModel) 
            zeta.1 <- t(as.matrix(apply(zeta.1.tmp, 2, weighted.mean, w=weights)))
            rm(zeta.1.tmp)
            zeta.0 <- t(as.matrix(apply(zeta.0.tmp, 2, weighted.mean, w=weights)))
            rm(zeta.0.tmp)
            delta.1 <- t(as.matrix(apply(delta.1.tmp, 2, weighted.mean, w=weights)))
            rm(delta.1.tmp)
            delta.0 <- t(as.matrix(apply(delta.0.tmp, 2, weighted.mean, w=weights)))
            rm(delta.0.tmp)
            tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
            
            
        ########################################################################
        ## Case I-2: Nonparametric Bootstrap
        ########################################################################
        } else {
            
            if(isS4(model.m)){
                Call.M <- model.m@call
            } else {
                Call.M <- model.m$call
            }
            if(isS4(model.y)){
                Call.Y.t <- model.y@call
            } else {
                Call.Y.t <- model.y$call
            }
            
            # Storage
            delta.1 <- matrix(NA, sims, 1)
            delta.0 <- matrix(NA, sims, 1)
            zeta.1 <- matrix(NA, sims, 1)
            zeta.0 <- matrix(NA, sims, 1)
            tau <- matrix(NA, sims, 1)
            
            
            
            # Bootstrap loop begins
            for(b in 1:sims){
                index <- sample(1:n, n, repl=TRUE)
                Call.M$data <- m.data[index,]
                Call.Y.t$data <- y.data[index,]
                Call.M$weights <- m.data[index,"(weights)"]
                Call.Y.t$weights  <- y.data[index,"(weights)"]
                
                # TODO: Verify if below is unnecessary with all outcome model types.
                if(ClassM=="polr" && length(unique(y.data[index,mediator]))!=m){
                        stop("insufficient variation on mediator")
                }
                
                # Refit Models with Resampled Data
                new.fit.M <- model.m <- eval.parent(Call.M)
                new.fit.t <- model.y <- eval.parent(Call.Y.t)
                
                #####################################
                #  Mediator Predictions
                #####################################
                pred.data.t <- m.data
                pred.data.c <- m.data
                
                if(isFactorT.m){
                    pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                } else {
                    pred.data.t[,treat] <- cat.1
                    pred.data.c[,treat] <- cat.0
                }
                
                ### Case I-2-a: GLM Mediator
                if(ClassM=="glm" || isGam.m){
                
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
            
                ### Case I-2-b: Ordered Mediator
                } else if(ClassM=="polr") {
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
                } else if(ClassM=="lm"){
                    sigma <- summary(new.fit.M)$sigma
                    error <- rnorm(n, mean=0, sd=sigma)
                    PredictM1 <- predict(new.fit.M, type="response", 
                                          newdata=pred.data.t) + error
                    PredictM0 <- predict(new.fit.M, type="response", 
                                          newdata=pred.data.c) + error
                    rm(error)
                    
                } else {
                    stop("mediator model is not yet implemented")
                }
               
                #####################################
                #  Outcome Predictions
                #####################################
                # ACME(1)
                pred.data <- predictY.dataprep(1,1,1,0)
                pred.data.t <- pred.data$t
                pred.data.c <- pred.data$c
                
                if(isRq.y){
                    pr.1 <- predict(new.fit.t, type="response",
                                         newdata=pred.data.t, interval="none")
                    pr.0 <- predict(new.fit.t, type="response",
                                         newdata=pred.data.c, interval="none")
                } else {
                    pr.1 <- predict(new.fit.t, type="response",
                                         newdata=pred.data.t)
                    pr.0 <- predict(new.fit.t, type="response",
                                         newdata=pred.data.c)
                }
                pr.mat <- as.matrix(cbind(pr.1, pr.0))
                delta.1.tmp <- pr.mat[,1] - pr.mat[,2]
                
                rm(pred.data.t, pred.data.c, pred.data, pr.1, pr.0, pr.mat)
                
                
                # ACME(0)
                pred.data <- predictY.dataprep(0,0,1,0)
                pred.data.t <- pred.data$t
                pred.data.c <- pred.data$c
                
                if(isRq.y){
                    pr.1 <- predict(new.fit.t, type="response",
                                        newdata=pred.data.t, interval="none")
                    pr.0 <- predict(new.fit.t, type="response",
                                        newdata=pred.data.c, interval="none");
                } else {
                    pr.1 <- predict(new.fit.t, type="response",
                                        newdata=pred.data.t)
                    pr.0 <- predict(new.fit.t, type="response",
                                        newdata=pred.data.c)
                }
                pr.mat <- as.matrix(cbind(pr.1, pr.0))
                delta.0.tmp <-pr.mat[,1] - pr.mat[,2]
                
                rm(pred.data.t, pred.data.c, pred.data, pr.1, pr.0, pr.mat)
               
                
                # DE(1)
                pred.data <- predictY.dataprep(1,0,1,1)
                pred.data.t <- pred.data$t
                pred.data.c <- pred.data$c
                
                if(isRq.y){
                    pr.1 <- predict(new.fit.t, type="response",
                                        newdata=pred.data.t, interval="none")
                    pr.0 <- predict(new.fit.t, type="response",
                                        newdata=pred.data.c, interval="none")
                } else {
                    pr.1 <- predict(new.fit.t, type="response",
                                        newdata=pred.data.t)
                    pr.0 <- predict(new.fit.t, type="response",
                                        newdata=pred.data.c)
                }
                pr.mat <- as.matrix(cbind(pr.1, pr.0))
                zeta.1.tmp <- pr.mat[,1] - pr.mat[,2]
                
                rm(pred.data.t, pred.data.c, pred.data, pr.1, pr.0, pr.mat)
                
                
                # DE(0)
                pred.data <- predictY.dataprep(1,0,0,0)
                pred.data.t <- pred.data$t
                pred.data.c <- pred.data$c
                
                if(isRq.y){
                    pr.1 <- predict(new.fit.t, type="response",
                                            newdata=pred.data.t, interval="none")
                    pr.0 <- predict(new.fit.t, type="response",
                                            newdata=pred.data.c, interval="none")
                } else {
                    pr.1 <- predict(new.fit.t, type="response",
                                            newdata=pred.data.t)
                    pr.0 <- predict(new.fit.t, type="response",
                                            newdata=pred.data.c)
                }
                pr.mat <- as.matrix(cbind(pr.1, pr.0))
                zeta.0.tmp <-pr.mat[,1] - pr.mat[,2]
                
                rm(pred.data.t, pred.data.c, pred.data, pr.1, pr.0, pr.mat)
                
                # Compute all QoIs
                zeta.1[b] <- weighted.mean(zeta.1.tmp, weights)
                zeta.0[b] <- weighted.mean(zeta.0.tmp, weights)
                delta.1[b] <- weighted.mean(delta.1.tmp, weights)
                delta.0[b] <- weighted.mean(delta.0.tmp, weights)
                tau[b] <- (zeta.1[b] + zeta.0[b] + delta.0[b] + delta.1[b])/2
                
            }  # bootstrap loop ends
        }  # nonpara boot branch ends
        
        ########################################################################
        ## Compute Outputs and Put Them Together
        ########################################################################
        d0 <- mean(delta.0)
        d1 <- mean(delta.1)
        low <- (1 - conf.level)/2
        high <- 1 - low
        d0.ci <- quantile(delta.0,c(low,high), na.rm=TRUE)
        d1.ci <- quantile(delta.1,c(low,high), na.rm=TRUE)
        tau.coef <- mean(tau)
        tau.ci <- quantile(tau,c(low,high), na.rm=TRUE)
        z1 <- mean(zeta.1)
        z1.ci <- quantile(zeta.1,c(low,high), na.rm=TRUE)
        z0 <- mean(zeta.0)
        z0.ci <- quantile(zeta.0,c(low,high), na.rm=TRUE)
        nu.0 <- delta.0/tau
        n0 <- median(nu.0)
        n0.ci <- quantile(nu.0, c(low,high), na.rm=TRUE)
        nu.1 <- delta.1/tau
        n1 <- median(nu.1)
        n1.ci <- quantile(nu.1, c(low,high), na.rm=TRUE)
        
        delta.avg <- (delta.0 + delta.1)/2
        zeta.avg <- (zeta.0 + zeta.1)/2
        d.avg <- mean(delta.avg)
        z.avg <- mean(zeta.avg)
        d.avg.ci <- quantile(delta.avg, c(low,high), na.rm=TRUE)
        z.avg.ci <- quantile(zeta.avg, c(low,high), na.rm=TRUE)
        nu.avg <- delta.avg/tau
        n.avg <- median(nu.avg)
        n.avg.ci <- quantile(nu.avg, c(low,high), na.rm=TRUE)
        
        # Detect whether models include T-M interaction
        if(!isS4(model.y)) {
            INT <- paste(treat,mediator,sep=":") %in% attr(model.y$terms,"term.labels") | 
                 paste(mediator,treat,sep=":") %in% attr(model.y$terms,"term.labels") 
        } else {
            INT <- paste(treat,mediator,sep=":") %in% attr(model.y@terms,"term.labels") |
                 paste(mediator,treat,sep=":") %in% attr(model.y@terms,"term.labels") 
        }
        if(!INT & isGam.y){
            INT <- !isTRUE(all.equal(d0, d1))  # if gam, determine empirically
        }
        
        if(long) {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                        d0.sims=delta.0, d1.sims=delta.1,
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci, 
                        z0.sims=zeta.0, z1.sims=zeta.1, 
                        n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
                        n0.sims=nu.0, n1.sims=nu.1,
                        tau.coef=tau.coef, tau.ci=tau.ci, 
                        tau.sims=tau,
                        d.avg=d.avg, d.avg.ci=d.avg.ci, d.avg.sims=delta.avg,
                        z.avg=z.avg, z.avg.ci=z.avg.ci, z.avg.sims=zeta.avg,
                        n.avg=n.avg, n.avg.ci=n.avg.ci, n.avg.sims=nu.avg,
                        boot=boot, treat=treat, mediator=mediator, 
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m, 
                        control.value=control.value, treat.value=treat.value,
                        nobs=n)
        } else {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci, 
                        n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
                        tau.coef=tau.coef, tau.ci=tau.ci, 
                        d.avg=d.avg, d.avg.ci=d.avg.ci,
                        z.avg=z.avg, z.avg.ci=z.avg.ci,
                        n.avg=n.avg, n.avg.ci=n.avg.ci,
                        boot=boot, treat=treat, mediator=mediator, 
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m, 
                        control.value=control.value, treat.value=treat.value,
                        nobs=n)
        }
        class(out) <- "mediate"
        out
        
    ############################################################################
    ############################################################################
    ### CASE II: ORDERED OUTCOME
    ############################################################################
    ############################################################################
    } else {
        if(boot != TRUE){
            warning("ordered outcome model can only be used with 
                     nonparametric bootstrap - option forced")
            boot <- TRUE
        }
        
        n.ycat <- length(unique(model.response(y.data)))
        
        # Storage - Now Dynamic
        delta.1 <- matrix(NA, sims, n.ycat)
        delta.0 <- matrix(NA, sims, n.ycat)
        zeta.1 <- matrix(NA, sims, n.ycat)
        zeta.0 <- matrix(NA, sims, n.ycat)
        tau <- matrix(NA, sims, n.ycat)
        
        # Bootstrap loop begins
        for(b in 1:sims){
            
            # Resampling Step
            # Pull off data.star.
            index <- sample(1:n, n, repl=TRUE)
            call.m <- model.m$call
            call.y <- model.y$call
            call.m$data <- m.data[index,]
            call.y$data <- y.data[index,]
            call.m$weights <- m.data[index,"(weights)"]
            call.y$weights  <- y.data[index,"(weights)"]
            new.fit.M <- eval.parent(call.m)
            new.fit.t <- eval.parent(call.y)
            
            if(ClassM=="polr" && length(unique(y.data[index,mediator]))!=m){
                # Modify the coefficients when mediator has empty cells
                coefnames.y <- names(model.y$coefficients)
                coefnames.new.y <- names(new.fit.t$coefficients)
                new.fit.t.coef <- rep(0, length(coefnames.y))
                names(new.fit.t.coef) <- coefnames.y
                new.fit.t.coef[coefnames.new.y] <- new.fit.t$coefficients
                new.fit.t$coefficients <- new.fit.t.coef
            }
            
            #####################################
            # Mediator Predictions
            #####################################
            pred.data.t <- m.data
            pred.data.t[,treat] <- cat.1
            pred.data.c <- m.data
            pred.data.c[,treat] <- cat.0
            
            if(isFactorT.m){
                pred.data.t[,treat] <- as.factor(pred.data.t[,treat])
                pred.data.c[,treat] <- as.factor(pred.data.c[,treat])
            } else { 
                pred.data.t[,treat] <- as.numeric(pred.data.t[,treat])
                pred.data.c[,treat] <- as.numeric(pred.data.c[,treat])
            } 
            
            ### Case II-a: GLM Mediator
            if(ClassM=="glm" || isGam.m){
                
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
            } else if(ClassM=="polr") {
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
            } else if(ClassM=="lm"){
                sigma <- summary(new.fit.M)$sigma
                error <- rnorm(n, mean=0, sd=sigma)
                PredictM1 <- predict(new.fit.M, type="response",
                                      newdata=pred.data.t) + error
                PredictM0 <- predict(new.fit.M, type="response",
                                      newdata=pred.data.c) + error
                rm(error)
            } else {
                stop("mediator model is not yet implemented")
            }
            
            #####################################
            #  Outcome Predictions
            #####################################
            # ACME(1)
            pred.data <- predictY.dataprep(1,1,1,0)
            pred.data.t <- pred.data$t
            pred.data.c <- pred.data$c
            
            probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
            probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
            delta.1.tmp <-  probs_p1 - probs_p0
            rm(pred.data.t, pred.data.c, pred.data, probs_p1, probs_p0)
            
            # ACME(0)
            pred.data <- predictY.dataprep(0,0,1,0)
            pred.data.t <- pred.data$t
            pred.data.c <- pred.data$c
            
            probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
            probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
            delta.0.tmp <- probs_p1 - probs_p0
            rm(pred.data.t, pred.data.c, pred.data, probs_p1, probs_p0)
            
            # DE(1)
            pred.data <- predictY.dataprep(1,0,1,1)
            pred.data.t <- pred.data$t
            pred.data.c <- pred.data$c
            
            probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
            probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
            zeta.1.tmp <- probs_p1 - probs_p0
            rm(pred.data.t, pred.data.c, pred.data, probs_p1, probs_p0)
            
            # DE(0)
            pred.data <- predictY.dataprep(1,0,0,0)
            pred.data.t <- pred.data$t
            pred.data.c <- pred.data$c
            
            probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
            probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
            zeta.0.tmp <- probs_p1 - probs_p0
            rm(pred.data.t, pred.data.c, pred.data, probs_p1, probs_p0)
            
            # Compute all QoIs
            zeta.1[b,] <- apply(zeta.1.tmp, 2, weighted.mean, w=weights)
            zeta.0[b,] <- apply(zeta.0.tmp, 2, weighted.mean, w=weights)
            delta.1[b,] <- apply(delta.1.tmp, 2, weighted.mean, w=weights)
            delta.0[b,] <- apply(delta.0.tmp, 2, weighted.mean, w=weights)
            tau[b,] <- (zeta.1[b,] + zeta.0[b,] + delta.0[b,] + delta.1[b,])/2
            
        }  # Bootstrap loop ends
            
        ########################################################################
        ## Compute Outputs and Put Them Together
        ########################################################################
        d0 <- apply(delta.0, 2, mean)
        d1 <- apply(delta.1, 2, mean)
        low <- (1 - conf.level)/2
        high <- 1 - low
        d0.ci <- apply(delta.0, 2, quantile, c(low,high))
        d1.ci <- apply(delta.1, 2, quantile, c(low,high))
        tau.coef <- apply(tau, 2, mean)
        tau.ci <- apply(tau, 2, quantile, c(low,high))
        z1 <- apply(zeta.1, 2, mean)
        z1.ci <- apply(zeta.1,2, quantile, c(low,high))
        z0 <- apply(zeta.0,2, mean)
        z0.ci <- apply(zeta.0,2, quantile, c(low,high))
        
        # Detect whether models include T-M interaction
        INT <- paste(treat,mediator,sep=":") %in% attr(model.y$terms,"term.labels") | 
             paste(mediator,treat,sep=":") %in% attr(model.y$terms,"term.labels") 
        
        if(long) {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                        d0.sims=delta.0, d1.sims=delta.1,
                        tau.coef=tau.coef, tau.ci=tau.ci, 
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci, 
                        z1.sims=zeta.1, z0.sims=zeta.0, tau.sims=tau,
                        boot=boot, treat=treat, mediator=mediator, 
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m, 
                        control.value=control.value, treat.value=treat.value, nobs=n)
        } else {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                        tau.coef=tau.coef, tau.ci=tau.ci, 
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,  
                        boot=boot, treat=treat, mediator=mediator, 
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m, 
                        control.value=control.value, treat.value=treat.value, nobs=n)
        }
        class(out) <- "mediate.order"
        out
    }
}



summary.mediate <- function(object, ...){
    structure(object, class = c("summary.mediate", class(object)))
}



print.summary.mediate <- function(x, ...){
    clp <- 100 * x$conf.level
    cat("\n Causal Mediation Analysis \n\n")
    if(x$boot==TRUE){
        cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
    } else {
        cat("Quasi-Bayesian Confidence Intervals\n\n")
    }
    
    printone <- x$INT == FALSE && (class(x$model.y)[1] %in% c("lm", "rq") ||
        (inherits(x$model.y, "glm") && x$model.y$family$family == "gaussian"
         && x$model.y$family$link == "identity"))
    
    if (printone){
        # Print only one set of values if lmY/quanY/linear gamY without interaction
        cat("Mediation Effect: ", format(x$d1, digits=4), clp, "% CI ", 
                format(x$d1.ci, digits=4), "\n")
        cat("Direct Effect: ", format(x$z0, digits=4), clp, "% CI ", 
                format(x$z0.ci, digits=4), "\n")
        cat("Total Effect: ", format(x$tau.coef, digits=4), clp, "% CI ", 
                format(x$tau.ci, digits=4), "\n")
        cat("Proportion of Total Effect via Mediation: ", format(x$n0, digits=4), 
                clp, "% CI ", format(x$n0.ci, digits=4),"\n\n")
        cat("Sample Size Used:", x$nobs,"\n\n")        
    } else {
        cat("Mediation Effect_0: ", format(x$d0, digits=4), clp, "% CI ", 
                format(x$d0.ci, digits=4), "\n")
        cat("Mediation Effect_1: ", format(x$d1, digits=4), clp, "% CI ", 
                format(x$d1.ci, digits=4), "\n")
        cat("Direct Effect_0: ", format(x$z0, digits=4), clp, "% CI ", 
                format(x$z0.ci, digits=4), "\n")
        cat("Direct Effect_1: ", format(x$z1, digits=4), clp, "% CI ", 
                format(x$z1.ci, digits=4), "\n")
        cat("Total Effect: ", format(x$tau.coef, digits=4), clp, "% CI ", 
                format(x$tau.ci, digits=4), "\n")
        cat("Proportion of Total Effect via Mediation_0: ", format(x$n0, digits=4), 
                clp, "% CI ", format(x$n0.ci, digits=4),"\n")
        cat("Proportion of Total Effect via Mediation_1: ", format(x$n1, digits=4), 
                clp, "% CI ", format(x$n1.ci, digits=4),"\n\n")
        cat("Mediation Effect (Average): ", format(x$d.avg, digits=4), clp, "% CI ", 
                format(x$d.avg.ci, digits=4), "\n")
        cat("Direct Effect (Average): ", format(x$z.avg, digits=4), clp, "% CI ", 
                format(x$z.avg.ci, digits=4), "\n")
        cat("Proportion of Total Effect via Mediation (Average): ", format(x$n.avg, digits=4), 
                clp, "% CI ", format(x$n.avg.ci, digits=4),"\n\n")
        cat("Sample Size Used:", x$nobs,"\n\n") 
    } 
    invisible(x)
}



summary.mediate.order <- function(object, ...){
    structure(object, class = c("summary.mediate.order", class(object)))
}



print.summary.mediate.order <- function(x, ...){
    tab.d0 <- rbind(x$d0, x$d0.ci)
    tab.d1 <- rbind(x$d1, x$d1.ci)
    tab.z0 <- rbind(x$z0, x$z0.ci)
    tab.z1 <- rbind(x$z1, x$z1.ci)
    tab.tau <- rbind(x$tau.coef, x$tau.ci)
    
    # Outcome Table Labels
    y.lab <- sort(unique(levels(model.frame(x$model.y)[,1])))
    out.names <- c()
    for(i in 1:length(y.lab)){
        out.names.tmp <- paste("Pr(Y=",y.lab[i],")",sep="")
        out.names <- c(out.names, out.names.tmp)
    }
    
    # Label Tables
    rownames(tab.d0)[1] <- "Mediation Effect_0: "
    colnames(tab.d0) <- out.names
    rownames(tab.d1)[1] <- "Mediation Effect_1: "
    colnames(tab.d1) <- out.names
    rownames(tab.z0)[1] <- "Direct Effect_0: "
    colnames(tab.z0) <- out.names
    rownames(tab.z1)[1] <- "Direct Effect_1: "
    colnames(tab.z1) <- out.names
    rownames(tab.tau)[1] <- "Total Effect: "
    colnames(tab.tau) <- out.names
    
    cat("\n Causal Mediation Analysis \n\n")
    cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
    print(tab.d0, digits=4)
    cat("\n")
    print(tab.d1, digits=4)
    cat("\n")
    print(tab.z0, digits=4)
    cat("\n")
    print(tab.z1, digits=4)
    cat("\n")
    print(tab.tau, digits=4)
    cat("\n\n")
    cat("Sample Size Used:", x$nobs,"\n\n") 
    invisible(x)
}



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



plot.mediate <- function(x, treatment = NULL,
                        labels = c("ACME","Direct\nEffect","Total\nEffect"), 
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
    
    param <- plot.process(x)
    y.axis <- c(length(param$coef.vec.1):.5)
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
        ylim <- c(min(y.axis) -1- 0.5, max(y.axis) + 0.5)
    }
    
    # Plot
    plot(param$coef.vec.1, y.axis, type = "n", xlab = xlab, ylab = ylab,
            yaxt = "n", xlim = xlim, ylim = ylim, main = main, ...)
    
    # Set offset values depending on number of bars to plot
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
    if(treatment[1]==0 & length(treatment)==1) {
        points(param$tau.vec[1], 1 , type = "p", pch = 19, cex = cex, col = col)
        segments(param$tau.vec[2], 1 , param$tau.vec[3], 1 , 
                lwd = lwd, col = col) 
    }
    y.axis.new <- c(3,2,1)
    axis(2, at = y.axis.new, labels = labels, las = 1, tick = TRUE, ...)
    abline(v = 0, lty = 2)
}



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



plot.mediate.order <- function(x, treatment = NULL,
                        labels = c("ACME","Direct\nEffect","Total\nEffect"), 
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
