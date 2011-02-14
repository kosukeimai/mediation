mediate <- function(model.m, model.y, sims=1000, boot=FALSE, treat="treat.name",
                    mediator="med.name", control=NULL, conf.level=.95, 
                    control.value=0, treat.value=1, long=TRUE, INT=NULL){
    
    # Warn users who still use INT option
    if(!is.null(INT)){
        warning("Argument INT no longer necessary; existence of interaction term
        is now automatically detected.")
    }
    
    # Model type indicators
    isGam.y <- class(model.y)[1] == "gam"
    isGam.m <- class(model.m)[1] == "gam"
    isVglm.y <- class(model.y)[1] == "vglm"
    isRq.y <- class(model.y)[1] == "rq"
    isOrdered.y <- class(model.y)[1] == "polr"
    
    # Record class of model.m as "ClassM"
    if(isGam.m){
        ClassM <- class(model.m)[2]
    } else {
        ClassM <- class(model.m)[1]
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
        
    # Detect whether models include T-M interaction
    if(isGam.y){
        INT <- !is.null(control)
    } else if(!isS4(model.y)) {
        INT <- paste(treat,mediator,sep=":") %in% attr(model.y$terms,"term.labels") | 
             paste(mediator,treat,sep=":") %in% attr(model.y$terms,"term.labels") 
    } else {
        INT <- paste(treat,mediator,sep=":") %in% attr(model.y@terms,"term.labels") |
             paste(mediator,treat,sep=":") %in% attr(model.y@terms,"term.labels") 
    }
    
    # Warning for control option in non-GAM outcome models
    if(!is.null(control) && !isGam.y){
        warning("argument control irrelevant except for GAM outcome models, 
        therefore ignored")
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
    
    # Define functions
    
    indexmax <- function(x){
        ## Return position of largest element in vector x
        order(x)[length(x)]
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
            # Error if gam
            if(isGam.y){
                stop("boot must be TRUE for gam models")
            }
            
            # Get mean and variance parameters for simulations
            MModel.coef <- model.m$coef
            if(ClassM=="polr"){  # TRUE if model.m is ordered
                if(is.null(model.m$Hess)){
                    cat("\nMediator model fitted without Hess = TRUE;")
                }
                k <- length(model.m$coef)
                MModel.var.cov <- vcov(model.m)[(1:k),(1:k)]
            } else {
                MModel.var.cov <- vcov(model.m)
            }
            if(isS4(model.y)){
                TMmodel.coef <- model.y@coefficients
            } else {
                TMmodel.coef <- model.y$coef
            }
            TMmodel.var.cov <- vcov(model.y)
            
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
            if(ClassM=="glm"){
                PredictM1 <- model.m$family$linkinv(MModel %*% t(mmat.t))
                PredictM0 <- model.m$family$linkinv(MModel %*% t(mmat.c))
            
            ### Case I-1-b: Ordered mediator
            } else if(ClassM=="polr"){
                if(model.m$method=="logit"){
                    linkfn <- plogis
                } else {
                    linkfn <- pnorm
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
            
            ### Case I-1-c: Other 
            } else {
                sigma <- summary(model.m)$sigma
                error <- rnorm(n, mean=0, sd=sigma)
                PredictM1 <- MModel %*% t(mmat.t)
                PredictM0 <- MModel %*% t(mmat.c)
                PredictM1 <- PredictM1 + error
                PredictM0 <- PredictM0 + error
                rm(error)
            }
            rm(mmat.t, mmat.c)
            
            #####################################
            ##  Outcome Predictions
            #####################################
            # ACME(1) Predictions Data
            Pr1 <- matrix(, nrow=n, ncol=sims)
            Pr0 <- matrix(, nrow=n, ncol=sims)
            
            for(j in 1:sims){
                pred.data.t <- y.data
                pred.data.c <- y.data
                if(isFactorT.y){
                    pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.1, levels = t.levels)
                } else {
                    pred.data.t[,treat] <- cat.1
                    pred.data.c[,treat] <- cat.1
                }
                
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictM1[j,], labels = m.levels)
                    pred.data.c[,mediator] <- factor(PredictM0[j,], labels = m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictM1[j,]
                    pred.data.c[,mediator] <- PredictM0[j,]
                }
                if(!isVglm.y){
                    ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
                    ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
                } else {
                    ymat.t <- model.matrix(model.y@terms, data=pred.data.t) 
                        # There seems to be a problem here in how things work
                    ymat.c <- model.matrix(model.y@terms, data=pred.data.c)        
                }
               
                #print(as.matrix(TMmodel[j,]))  # the second row of this is a Log(scale) 
                                                # parameter from tobit.
                #print(t(ymat.t))  # this only has 3 rows, hence we must take this 
                                   # second row out in order to create the predictions 
                                   # we want.
                #print(dim(ymat.t))
                #print(as.matrix(TMmodel[j,]))
                #print(dim(as.matrix(TMmodel[j,])))
              
                # ACME(1) Predictions
                if(isVglm.y){
                    if(vfamily=="tobit") {
                        # TODO: This part should be made into a function
                        Pr1[,j] <- t(as.matrix(TMmodel[j,-2])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(TMmodel[j,-2])) %*% t(ymat.c)
                        model.y.upper <- model.y@misc$Upper
                        model.y.lower <- model.y@misc$Lower
                        censored <- (Pr1[,j] > model.y.upper)
                        Pr1[censored,j] <- model.y.upper
                        censored <- (Pr1[,j] < model.y.lower)
                        Pr1[censored,j] <- model.y.lower
                        censored <- (Pr0[,j] > model.y.upper)
                        Pr0[censored,j] <- model.y.upper
                        censored <- (Pr0[,j] < model.y.lower)
                        Pr0[censored,j] <- model.y.lower
                    } else {
                        stop("outcome model is in unsupported vglm family")
                    }
                } else {
                    Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
                    Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
                }
                
                rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
            }
            
            if(ClassY=="glm"){
                Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                delta.1.tmp <- Pr1 - Pr0
            } else {
                delta.1.tmp <- Pr1 - Pr0
            }
            rm(Pr1, Pr0)
                    
            # ACME(0) Predictions Data
            Pr1 <- matrix(, nrow=n, ncol=sims)
            Pr0 <- matrix(, nrow=n, ncol=sims)
            
            for(j in 1:sims){
                pred.data.t <- y.data
                pred.data.c <- y.data
                if(isFactorT.y){
                    pred.data.t[,treat] <- factor(cat.0, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                } else{
                    pred.data.t[,treat] <- cat.0
                    pred.data.c[,treat] <- cat.0
                }
                
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictM1[j,], labels = m.levels)
                    pred.data.c[,mediator] <- factor(PredictM0[j,], labels = m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictM1[j,]
                    pred.data.c[,mediator] <- PredictM0[j,]
                }
                
                ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
                ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
                
                # ACME(0) Predictions
                if(isVglm.y){
                    if(vfamily=="tobit") {
                        # TODO: This part should be made into a function
                        Pr1[,j] <- t(as.matrix(TMmodel[j,-2])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(TMmodel[j,-2])) %*% t(ymat.c)
                        model.y.upper <- model.y@misc$Upper
                        model.y.lower <- model.y@misc$Lower
                        censored <- (Pr1[,j] > model.y.upper)
                        Pr1[censored,j] <- model.y.upper
                        censored <- (Pr1[,j] < model.y.lower)
                        Pr1[censored,j] <- model.y.lower
                        censored <- (Pr0[,j] > model.y.upper)
                        Pr0[censored,j] <- model.y.upper
                        censored <- (Pr0[,j] < model.y.lower)
                        Pr0[censored,j] <- model.y.lower
                    } else {
                        stop("outcome model is in unsupported vglm family")
                    }
                } else {
                        Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
                }
                rm(ymat.t, ymat.c)
            }
            
            if(ClassY=="glm"){
                Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                delta.0.tmp <- Pr1 - Pr0
            } else {
                delta.0.tmp <- Pr1 - Pr0
            }
            rm(Pr1, Pr0)
            
            # DE(1) Predictions Data
            Pr1 <- matrix(, nrow=n, ncol=sims)
            Pr0 <- matrix(, nrow=n, ncol=sims)
            
            for(j in 1:sims){
                pred.data.t <- y.data
                pred.data.c <- y.data
                if(isFactorT.y){
                    pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                } else {
                    pred.data.t[,treat] <- cat.1
                    pred.data.c[,treat] <- cat.0
                }
                
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictM1[j,], labels = m.levels)
                    pred.data.c[,mediator] <- factor(PredictM1[j,], labels = m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictM1[j,]
                    pred.data.c[,mediator] <- PredictM1[j,]
                }
                
                ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
                ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
                
                # DE(1) Predictions
                if(isVglm.y){
                    if(vfamily=="tobit") {
                        # TODO: This part should be made into a function
                        Pr1[,j] <- t(as.matrix(TMmodel[j,-2])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(TMmodel[j,-2])) %*% t(ymat.c)
                        model.y.upper <- model.y@misc$Upper
                        model.y.lower <- model.y@misc$Lower
                        censored <- (Pr1[,j] > model.y.upper)
                        Pr1[censored,j] <- model.y.upper
                        censored <- (Pr1[,j] < model.y.lower)
                        Pr1[censored,j] <- model.y.lower
                        censored <- (Pr0[,j] > model.y.upper)
                        Pr0[censored,j] <- model.y.upper
                        censored <- (Pr0[,j] < model.y.lower)
                        Pr0[censored,j] <- model.y.lower
                    } else {
                        stop("outcome model is in unsupported vglm family")
                    }
                } else {
                    Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
                    Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
                }
                rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
            }    
            
            if(ClassY=="glm"){
                Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                zeta.1.tmp <- Pr1 - Pr0
            } else {
                zeta.1.tmp <- Pr1 - Pr0
            }
            rm(Pr1, Pr0)
            
            # DE(0) Predictions Data
            Pr1 <- matrix(, nrow=n, ncol=sims)
            Pr0 <- matrix(, nrow=n, ncol=sims)
            
            for(j in 1:sims){
                pred.data.t <- y.data
                pred.data.c <- y.data
                if(isFactorT.y){
                    pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                } else {
                    pred.data.t[,treat] <- cat.1
                    pred.data.c[,treat] <- cat.0
                }
                
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictM0[j,], labels = m.levels)
                    pred.data.c[,mediator] <- factor(PredictM0[j,], labels = m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictM0[j,]
                    pred.data.c[,mediator] <- PredictM0[j,]
                }
                
                ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
                ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
                
                # DE(0) Predictions
                if(isVglm.y){
                    if(vfamily=="tobit") {
                        # TODO: This part should be made into a function
                        Pr1[,j] <- t(as.matrix(TMmodel[j,-2])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(TMmodel[j,-2])) %*% t(ymat.c)
                        model.y.upper <- model.y@misc$Upper
                        model.y.lower <- model.y@misc$Lower
                        censored <- (Pr1[,j] > model.y.upper)
                        Pr1[censored,j] <- model.y.upper
                        censored <- (Pr1[,j] < model.y.lower)
                        Pr1[censored,j] <- model.y.lower
                        censored <- (Pr0[,j] > model.y.upper)
                        Pr0[censored,j] <- model.y.upper
                        censored <- (Pr0[,j] < model.y.lower)
                        Pr0[censored,j] <- model.y.lower
                    } else {
                        stop("outcome model is in unsupported vglm family")
                    }
                } else {        
                    Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
                    Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
                }
                rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
            }
            
            if(ClassY=="glm"){
                Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                zeta.0.tmp <- Pr1 - Pr0
            } else {
                zeta.0.tmp <- Pr1 - Pr0
            }
            
            rm(Pr1, Pr0, PredictM1, PredictM0, TMmodel, MModel) 
            zeta.1 <- t(as.matrix(apply(zeta.1.tmp, 2, mean)))
            rm(zeta.1.tmp)
            zeta.0 <- t(as.matrix(apply(zeta.0.tmp, 2, mean)))
            rm(zeta.0.tmp)
            delta.1 <- t(as.matrix(apply(delta.1.tmp, 2, mean)))
            rm(delta.1.tmp)
            delta.0 <- t(as.matrix(apply(delta.0.tmp, 2, mean)))
            rm(delta.0.tmp)
            tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
            
            
        ########################################################################
        ## Case I-2: Nonparametric Bootstrap
        ########################################################################
        } else {
            
            Call.M <- model.m$call
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
                
                if(ClassM=="polr" && length(unique(y.data[index,mediator]))!=m){
                        stop("Insufficient Variation on Mediator")
                }
                
               #print("here tobit error")  # After X iterations, where X can vary, 
                                           # we get the following: 
                                           # Error in lm.fit(X_vlm, z_vlm, ...) :   
                                           #     NA/NaN/Inf in foreign function call (arg 4)
                                           # Printing the call objects doesn't 
                                           # indicate an issue.
                                           # print(eval.parent(Call.M))
                                           # print(eval.parent(Call.Y.t))
                                       
                # Refit Models with Resampled Data
                new.fit.M <- eval.parent(Call.M)
                new.fit.t <- eval.parent(Call.Y.t)
                
                #####################################
                #  Mediator Predictions
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
                
                ### Case I-2-a: GLM Mediator
                if(ClassM=="glm"){
                    PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
                    PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)
                    
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
                    
                ### Case I-2-c: Other
                } else {
                    if(isGam.m){
                        sigma <- summary(new.fit.M)$scale
                    } else {
                        # print(summary(new.fit.M))  # when of class rq (quantile regression)
                                                     # there is not such sigma variable
                        # print(names(summary(new.fit.M)))  # this produces  
                                                            #   "call"
                                                            #   "terms"
                                                            #   "coefficients"
                                                            #   "rdf"
                                                            #   "tau"
                        sigma <- summary(new.fit.M)$sigma
                    }
                    if(class(model.m)[1]!="rq"){
                        error <- rnorm(n, mean=0, sd=sigma)
                        PredictM1 <- predict(new.fit.M, type="response", 
                                              newdata=pred.data.t) + error
                        PredictM0 <- predict(new.fit.M, type="response", 
                                              newdata=pred.data.c) + error
                        rm(error)
                    } else {
                        # For quantile regression we are unable to introduce error 
                        # into the prediction. Hence sample quantiles from 0 1 and
                        # then predict.
                        tau.list <- runif(n)
                        PredictM1 <- PredictM0 <- vector(length=n)
                        for(zz in 1:n){
                            PredictM1[zz] <- predict(update(new.fit.M, tau=tau.list[zz]),
                                                      type="response", newdata=pred.data.t)
                            PredictM0[zz] <- predict(update(new.fit.M, tau=tau.list[zz]),
                                                      type="response", newdata=pred.data.c)
                        }
                    }
                }
               
                #####################################
                #  Outcome Predictions
                #####################################
                # ACME(1) Predictions Data
                pred.data.t <- y.data
                pred.data.c <- y.data
                    
                if(isFactorT.y){
                    pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.1, levels = t.levels)
                    if(!is.null(control)){
                        pred.data.t[,control] <- factor(cat.0, levels = t.levels)
                        pred.data.c[,control] <- factor(cat.0, levels = t.levels)
                    }
                } else {
                    pred.data.t[,treat] <- cat.1
                    pred.data.c[,treat] <- cat.1
                    if(!is.null(control)){
                        pred.data.t[,control] <- cat.0
                        pred.data.c[,control] <- cat.0
                    }
                }
                
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictM1, labels = m.levels)
                    pred.data.c[,mediator] <- factor(PredictM0, labels = m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictM1
                    pred.data.c[,mediator] <- PredictM0
                }
                
                # ACME(1) Predictions
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
                
                rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
                
                
                # ACME(0) Predictions Data
                pred.data.t <- y.data
                pred.data.c <- y.data
                    
                if(isFactorT.y){
                    pred.data.t[,treat] <- factor(cat.0, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                    if(!is.null(control)){
                        pred.data.t[,control] <- factor(cat.1, levels = t.levels)
                        pred.data.c[,control] <- factor(cat.1, levels = t.levels)
                    }
                } else {
                    pred.data.t[,treat] <- cat.0
                    pred.data.c[,treat] <- cat.0
                    if(!is.null(control)){
                        pred.data.t[,control] <- cat.1
                        pred.data.c[,control] <- cat.1
                    }
                }
                
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictM1, labels = m.levels)
                    pred.data.c[,mediator] <- factor(PredictM0, labels = m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictM1
                    pred.data.c[,mediator] <- PredictM0
                }
                
                # ACME(0) Predictions
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
                
                rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
               
                
                # DE(1) Predictions Data
                pred.data.t <- y.data
                pred.data.c <- y.data
                if(isFactorT.y){
                    pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                    if(!is.null(control)){
                        pred.data.t[,control] <- factor(cat.0, levels = t.levels)
                        pred.data.c[,control] <- factor(cat.1, levels = t.levels)
                    }
                } else {
                    pred.data.t[,treat] <- cat.1
                    pred.data.c[,treat] <- cat.0
                    if(!is.null(control)){
                        pred.data.t[,control] <- cat.0
                        pred.data.c[,control] <- cat.1
                    }
                }
                
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictM1, labels = m.levels)
                    pred.data.c[,mediator] <- factor(PredictM1, labels = m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictM1
                    pred.data.c[,mediator] <- PredictM1
                }
                
                # DE(1) Predictions
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
                
                rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
                
                
                # DE(0) Predictions Data
                pred.data.t <- y.data
                pred.data.c <- y.data
                if(isFactorT.y){
                    pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                    if(!is.null(control)){
                        pred.data.t[,control] <- factor(cat.0, levels = t.levels)
                        pred.data.c[,control] <- factor(cat.1, levels = t.levels)
                    }
                } else {
                    pred.data.t[,treat] <- cat.1
                    pred.data.c[,treat] <- cat.0
                    if(!is.null(control)){
                        pred.data.t[,control] <- cat.0
                        pred.data.c[,control] <- cat.1
                    }
                }
                
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictM0, labels = m.levels)
                    pred.data.c[,mediator] <- factor(PredictM0, labels = m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictM0
                    pred.data.c[,mediator] <- PredictM0
                }
                
                # DE(0) Predictions
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
                
                rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
                
                
                # Compute all QoIs
                zeta.1[b] <- mean(zeta.1.tmp)
                zeta.0[b] <- mean(zeta.0.tmp)
                delta.1[b] <- mean(delta.1.tmp)
                delta.0[b] <- mean(delta.0.tmp)
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
        
        avg.delta <- (d0 + d1)/2
        pct.dist <- avg.delta/tau
        pct.coef <- median(pct.dist)
        pct.ci <- quantile(pct.dist,c(low,high), na.rm=TRUE)
        
        if(long) {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                        d0.sims=delta.0, d1.sims=delta.1,
                        pct.coef=pct.coef, pct.ci=pct.ci,
                        tau.coef=tau.coef, tau.ci=tau.ci, 
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci, 
                        z1.sims=zeta.1, z0.sims=zeta.0, tau.sims=tau,
                        boot=boot, treat=treat, mediator=mediator, 
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m, 
                        control.value=control.value, treat.value=treat.value)
        } else {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                        pct.coef=pct.coef, pct.ci=pct.ci,
                        tau.coef=tau.coef, tau.ci=tau.ci, 
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,  
                        boot=boot, treat=treat, mediator=mediator, 
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m, 
                        control.value=control.value, treat.value=treat.value)    
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
            warning("Ordered outcome model can only be used with 
                     nonparametric bootstrap, option forced")
            boot <- TRUE
        }
        
        # Storage - Now Dynamic
        delta.1 <- matrix(NA, sims, m)
        delta.0 <- matrix(NA, sims, m)
        zeta.1 <- matrix(NA, sims, m)
        zeta.0 <- matrix(NA, sims, m)
        tau <- matrix(NA, sims, m)
        
        # Bootstrap loop begins
        for(b in 1:sims){
            if(ClassM=="polr" && length(unique(y.data[index,mediator]))!=m){
                stop("Insufficient Variation on Mediator")
            }
            
            # Resampling Step
            # Pull off data.star.
            index <- sample(1:n, n, repl=TRUE)
            data.star.m <- m.data[index,]
            data.star.y <- y.data[index,]
            new.fit.M <- update(model.m, data=data.star.m)
            new.fit.t <- update(model.y, data=data.star.y)
            
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
            if(ClassM=="glm"){
                PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
                PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)
            
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
#                max.pr_m1 <- apply(probs_m1, 1, max)
#                max.pr_m0 <- apply(probs_m0, 1, max)
#                cat.m1 <- predict(new.fit.M, newdata=pred.data.t)
#                cat.m0 <- predict(new.fit.M, newdata=pred.data.c)
#                draws_m1 <- round(runif(n, m.min, m), 0)
#                draws_m0 <- round(runif(n, m.min, m), 0)
#                # where is .75 from?
#                PredictM1 <- ifelse(max.pr_m1 > runif(n, 0, .75), draws_m1, cat.m1)
#                PredictM0 <- ifelse(max.pr_m0 > runif(n, 0, .75), draws_m0, cat.m0)
                PredictM1 <- apply(draws_m1, 1, indexmax)
                PredictM0 <- apply(draws_m0, 1, indexmax)
            
            ### Case II-c: Other
            } else {
                if(isGam.m){
                    sigma <- summary(new.fit.M)$scale
                } else {
                    sigma <- summary(new.fit.M)$sigma
                }
                error <- rnorm(n, mean=0, sd=sigma)
                PredictM1 <- predict(new.fit.M, type="response",
                                      newdata=pred.data.t) + error
                PredictM0 <- predict(new.fit.M, type="response",
                                      newdata=pred.data.c) + error
                rm(error)
            }
            
            #####################################
            #  Outcome Predictions
            #####################################
            # ACME(1) Predictions Data
            pred.data.t <- y.data
            pred.data.c <- y.data
            
            if(isFactorT.y){
                pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                pred.data.c[,treat] <- factor(cat.1, levels = t.levels)
                if(!is.null(control)){
                    pred.data.t[,control] <- factor(cat.0, levels = t.levels)
                    pred.data.c[,control] <- factor(cat.0, levels = t.levels)
                }
            } else {
                pred.data.t[,treat] <- cat.1
                pred.data.c[,treat] <- cat.1
                if(!is.null(control)){
                    pred.data.t[,control] <- cat.0
                    pred.data.c[,control] <- cat.0
                }
            }
            
            if(isFactorM) {
                pred.data.t[,mediator] <- factor(PredictM1, labels = m.levels)
                pred.data.c[,mediator] <- factor(PredictM0, labels = m.levels)
            } else {
                pred.data.t[,mediator] <- PredictM1
                pred.data.c[,mediator] <- PredictM0
            }
            
            # ACME(1) Predictions
            probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
            probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
            
            delta.1.tmp <-  probs_p1 - probs_p0
            rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
            
            # ACME(0) Predictions Data
            pred.data.t <- y.data
            pred.data.c <- y.data
            
            if(isFactorT.y){
                pred.data.t[,treat] <- factor(cat.0, levels = t.levels)
                pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                if(!is.null(control)){
                    pred.data.t[,control] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,control] <- factor(cat.1, levels = t.levels)
                }
            } else {
                pred.data.t[,treat] <- cat.0    
                pred.data.c[,treat] <- cat.0
                if(!is.null(control)){
                    pred.data.t[,control] <- cat.1
                    pred.data.c[,control] <- cat.1
                }
            }
            
            if(isFactorM) {
                pred.data.t[,mediator] <- factor(PredictM1, labels = m.levels)
                pred.data.c[,mediator] <- factor(PredictM0, labels = m.levels)
            } else {
                pred.data.t[,mediator] <- PredictM1
                pred.data.c[,mediator] <- PredictM0
            }
            
            # ACME(0) Predictions
            probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
            probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
            
            delta.0.tmp <- probs_p1 - probs_p0
            
            rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
            
            # DE(1) Predictions Data
            pred.data.t <- y.data
            pred.data.c <- y.data
            
            if(isFactorT.y){
                pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                if(!is.null(control)){
                    pred.data.t[,control] <- factor(cat.0, levels = t.levels)
                    pred.data.c[,control] <- factor(cat.1, levels = t.levels)
                }
            } else {
                pred.data.t[,treat] <- cat.1
                pred.data.c[,treat] <- cat.0
                if(!is.null(control)){
                    pred.data.t[,control] <- cat.0
                    pred.data.c[,control] <- cat.1
                }
            }
            
            if(isFactorM) {
                pred.data.t[,mediator] <- factor(PredictM1, labels = m.levels)
                pred.data.c[,mediator] <- factor(PredictM1, labels = m.levels)
            } else {
                pred.data.t[,mediator] <- PredictM1 
                pred.data.c[,mediator] <- PredictM1
            }
            
            # DE(1) Predictions
            probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
            probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
            zeta.1.tmp <- probs_p1 - probs_p0
            
            rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
            
            # DE(0) Predictions Data
            pred.data.t <- y.data
            pred.data.c <- y.data
            
            if(isFactorT.y){
                pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                if(!is.null(control)){
                    pred.data.t[,control] <- factor(cat.0, levels = t.levels)
                    pred.data.c[,control] <- factor(cat.1, levels = t.levels)
                }
            } else {
                pred.data.t[,treat] <- cat.1
                pred.data.c[,treat] <- cat.0
                if(!is.null(control)){
                   pred.data.t[,control] <- cat.0
                   pred.data.c[,control] <- cat.1
                }
            }
            
            if(isFactorM) {
                pred.data.t[,mediator] <- factor(PredictM0, labels = m.levels)
                pred.data.c[,mediator] <- factor(PredictM0, labels = m.levels)
            } else {
                pred.data.t[,mediator] <- PredictM0
                pred.data.c[,mediator] <- PredictM0
            }

            # DE(0) Predictions
            probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
            probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
            
            zeta.0.tmp <- probs_p1 - probs_p0
            rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
            
            # Compute all QoIs
            zeta.1[b,] <- apply(zeta.1.tmp, 2, mean)
            zeta.0[b,] <- apply(zeta.0.tmp, 2, mean)
            delta.1[b,] <- apply(delta.1.tmp, 2, mean)
            delta.0[b,] <- apply(delta.0.tmp, 2, mean)
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
        
        if(long) {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                        d0.sims=delta.0, d1.sims=delta.1,
                        tau.coef=tau.coef, tau.ci=tau.ci, 
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci, 
                        z1.sims=zeta.1, z0.sims=zeta.0, tau.sims=tau,
                        boot=boot, treat=treat, mediator=mediator, 
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m, 
                        control.value=control.value, treat.value=treat.value)
        } else {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                        tau.coef=tau.coef, tau.ci=tau.ci, 
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,  
                        boot=boot, treat=treat, mediator=mediator, 
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m, 
                        control.value=control.value, treat.value=treat.value)
        }
        class(out) <- "mediate.order"
        out
    }
}



print.mediate <- function(x, ...){
    print(unlist(x[1:11]))
    invisible(x)
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
    if (x$INT == FALSE && class(x$model.y)[1] == "lm"){
        # Print only one set of values if lm without interaction
        cat("Mediation Effect: ", format(x$d1, digits=4), clp, "% CI ", 
                format(x$d1.ci, digits=4), "\n")
        cat("Direct Effect: ", format(x$z0, digits=4), clp, "% CI ", 
                format(x$z0.ci, digits=4), "\n")
        cat("Total Effect: ", format(x$tau.coef, digits=4), clp, "% CI ", 
                format(x$tau.ci, digits=4), "\n")
        cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4), 
                clp, "% CI ", format(x$pct.ci, digits=4),"\n")
        #cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"\n")
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
        cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4), 
                clp, "% CI ", format(x$pct.ci, digits=4),"\n")
        #cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"\n")
    } 
    invisible(x)
}



print.mediate.order <- function(x, ...){
    print(unlist(x[1:50]))
    invisible(x)
}



summary.mediate.order <- function(object, ...){
    structure(object, class = c("summary.mediate.order", class(object)))
}



print.summary.mediate.order <- function(x, ...){
    tab.d0 <- rbind(x$d0, x$d0.ci)
    tab.d1 <- rbind(x$d1, x$d1.ci)
    tab.z0 <- rbind(x$d0, x$d0.ci)
    tab.z1 <- rbind(x$d1, x$d1.ci)
    tab.tau <- rbind(x$tau.coef, x$tau.ci)
    
    # Outcome Table Labels
    m.lab <- sort(unique(levels(model.frame(x$model.y)[,1])))
    out.names <- c()
    for(i in 1:length(x$m.lab)){
        out.names.tmp <- paste("Pr(Y=",x$m.lab[i],")",sep="")
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
    
    # TODO: Given the nonlinear models its possible for d0 not = d1 even if
    # no interaction. we have to fix this in the regular code too.
    
    if(x$INT==TRUE){
        cat("\n Causal Mediation Analysis \n\n")
        cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
        print(tab.d0, digits=4)
        print(tab.d1, digits=4)
        print(tab.z0, digits=4)
        print(tab.z1, digits=4)
        print(tab.tau, digits=4)
    } else {
        cat("\n Causal Mediation Analysis \n\n")
        cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
        print(tab.d0, digits=4) 
        cat("\n")
        print(tab.z0, digits=4)
        cat("\n")
        print(tab.tau, digits=4)
    }
    invisible(x)
}



plot.process <- function(model) {
    coef.vec.1 <- c(model$d1, model$z1, model$tau.coef)
    lower.vec.1 <- c(model$d1.ci[1], model$z1.ci[1],model$tau.ci[1])
    upper.vec.1 <- c(model$d1.ci[2], model$z1.ci[2],model$tau.ci[2])
    range.1 <- range(model$d1.ci[1], model$z1.ci[1],model$tau.ci[1],
                      model$d1.ci[2], model$z1.ci[2],model$tau.ci[2])
    
    coef.vec.0 <- c(model$d0, model$z0, model$tau.coef)
    lower.vec.0 <- c(model$d0.ci[1], model$z0.ci[1],model$tau.ci[1])
    upper.vec.0 <- c(model$d0.ci[2], model$z0.ci[2],model$tau.ci[2])
    range.0 <- range(model$d0.ci[1], model$z0.ci[1],model$tau.ci[1],
                      model$d0.ci[2], model$z0.ci[2],model$tau.ci[2])
    
    return(list(coef.vec.1=coef.vec.1, lower.vec.1=lower.vec.1, 
                upper.vec.1=upper.vec.1, coef.vec.0=coef.vec.0,
                lower.vec.0=lower.vec.0, upper.vec.0=upper.vec.0,
                range.1=range.1, range.0=range.0))
}



plot.mediate <- function(x, treatment = NULL,
                        xlim = NULL, xlab = "", ylim = NULL, ylab = "",
                        labels = c("ACME","Direct\nEffect","Total\nEffect"), 
                        main = NULL, lwd = 1.5, cex = .85,
                        col = "black", ...){
    # Determine which graph to plot
    if(is.null(treatment)){
        if(x$INT){
            treatment <- c(0,1)
        } else {
            treatment <- 1
        }
    }
    
    param <- plot.process(x)
    y.axis <- c(length(param$coef.vec.1):.5)
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
        ylim <- c(min(y.axis) - 0.5, max(y.axis) + 0.5)
    }
    
    # Plot
    plot(param$coef.vec.1, y.axis, type = "n", xlab = xlab, ylab = ylab,
            yaxt = "n", xlim = xlim, ylim = ylim, main = main, ...) 
            # create empty plot with labels
    
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
             # coef +/-1.96*se = 95% interval, lwd adjusts line thickness
    }
    if(0 %in% treatment) {
        points(param$coef.vec.0, y.axis - adj, type = "p", pch = 1, cex = cex, col = col)
         # plot coefficients as points, turning off axes and labels. 
        segments(param$lower.vec.0, y.axis - adj, param$upper.vec.0, y.axis - adj, 
                lwd = lwd, lty = 3, col = col)
         # coef +/-1.96*se = 95% interval, lwd adjusts line thickness
    }
    axis(2, at = y.axis, labels = labels, las = 1, tick = TRUE, ...)
         # draw y-axis with tick marks, make labels perpendicular to axis and closer to axis
    abline(v = 0, lty = 2)
}



mediations <- function(datasets, treatment, mediators, outcome, 
                    covariates=NULL, family=c("gaussian", "gaussian"),
                    tau_m=.5,tau_y=.5, LowerY=NULL, UpperY=NULL, interaction=FALSE,
                    conf.level=.95, sims=500) {
    data <- names(datasets)
    labels <- c()
    out <- list()
    count <- 1
        for (i in 1:length(treatment)) {
            d1 <- sprintf("datasets$%s", data[i])
            dataarg <- eval(parse(text=d1))
       for (o in 1:length(outcome)) {
         
            for (j in 1:length(mediators)) {
                #create model formulas
                if(is.null(covariates)) {
                    f1 <- sprintf("%s ~ %s + %s", mediators[j], treatment[i])
                    f2 <- sprintf("%s ~ %s + %s ", outcome[o], treatment[i], mediators[j])
                    if (interaction) {
                        f2 <- sprintf("%s ~ %s*%s ",outcome[o], treatment[i], mediators[j])
                    }
                } else {
                    f1 <- sprintf("%s ~ %s + %s", mediators[j], treatment[i], covariates)
                    f2 <- sprintf("%s ~ %s + %s + %s", outcome[o], treatment[i], 
                                                                mediators[j], covariates)
                    if (interaction) {
                        f2 <- sprintf("%s ~ %s*%s + %s", outcome[o], treatment[i], 
                                                                mediators[j], covariates)
                    }
                }
                fmla <- as.formula(paste(f1))
                fmla2 <- as.formula(paste(f2))
                
                # Depending on the family chosen for the mediator or outcome model, 
                # this generates the correct formulas.
                # NOTE: mediations is a WRAPPER function for mediate--things like 
                # "update" etc. have no use here!
                # Mediations does not currently support the use of GAM's. Logits 
                # are not supported bc. users should use probits so they can do 
                # sensitivity analyses.  
                
                if(family[1] == "binomial") {  # run Mediator model using new data/specification
                    result1 <- glm(fmla, family=binomial("probit"), data=dataarg)
                } else if(family[1] == "quantile") {
                    result1 <- rq(fmla, data=dataarg, tau=tau_m)
                } else {
                    result1 <- glm(fmla, family=family[1], data=dataarg)
                }
                
                if(family[2] == "binomial") {  # run Outcome model using new data/specification
                    result2 <- glm(fmla2, family=binomial("probit"), data=dataarg)
                } else if(family[1] == "quantile") {
                    result2 <- rq(fmla2, data=dataarg, tau=tau_y)
                } else if(family[2] == "tobit") {
                    result2 <- vglm(fmla2, tobit(Lower=LowerY,Upper=UpperY), data=dataarg)
                } else {
                    result2 <- glm(fmla2, family=family[2], data=dataarg)
                }
              
              out[[(count)]] <- mediate(result1, result2, sims=sims, 
                                        treat=treatment[i], mediator=mediators[j],
                                        conf.level=conf.level)
              rm(result1)
              rm(result2)
              labels[(count)] <- sprintf("%s.%s.%s", outcome[o],data[i], mediators[j])
              count <- count + 1
            }
        }
    }    
    names(out) <- labels
    class(out) <- "mediations"
    out
}



plot.mediations <- function(x, which = names(x),
            ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for(i in 1:length(which)){
        plot.mediate(x[[i]], xlab = which[i], ...)
    }
}



summary.mediations <- function(object, ...){
    structure(object, class = c("summary.mediations", class(object)))
}



print.summary.mediations <- function(x, ...){
    clp <- 100 * x[[1]]$conf.level
    name.list <- names(x)
    for(i in 1:length(name.list)){
        cat("\n Causal Mediation Analysis \n\n")
        cat("Specification",name.list[i], "\n\n")
        if(x[[i]]$boot==TRUE){
            cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
        } else {
            cat("Quasi-Bayesian Confidence Intervals\n\n")
        }
        cat("Mediation Effect_0: ", format(x[[i]]$d0, digits=4), clp, "% CI ", 
                format(x[[i]]$d0.ci, digits=4), "\n")
        cat("Mediation Effect_1: ", format(x[[i]]$d1, digits=4), clp, "% CI ", 
                format(x[[i]]$d1.ci, digits=4), "\n")
        cat("Direct Effect_0: ", format(x[[i]]$z0, digits=4), clp, "% CI ", 
                format(x[[i]]$z0.ci, digits=4), "\n")
        cat("Direct Effect_1: ", format(x[[i]]$z1, digits=4), clp, "% CI ", 
                format(x[[i]]$z1.ci, digits=4), "\n")
        cat("Total Effect: ", format(x[[i]]$tau.coef, digits=4), clp, "% CI ", 
                format(x[[i]]$tau.ci, digits=4), "\n")
        cat("Proportion of Total Effect via Mediation: ", format(x[[i]]$pct.coef, digits=4), 
                clp, "% CI ", format(x[[i]]$pct.ci, digits=4),"\n")
        #cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"\n")
    }
    invisible(x)
}
