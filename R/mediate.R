mediate <- function(model.m, model.y, sims=1000, boot=FALSE, treat="treat.name", mediator="med.name", control=NULL, conf.level=.95, treat.0=0, treat.1=1){
    if(isS4(model.y)!=TRUE) {
    INT <- paste(treat,mediator,sep=":") %in% attr(model.y$terms,"term.labels") | 
         paste(mediator,treat,sep=":") %in% attr(model.y$terms,"term.labels") 
    } else {
    INT <- paste(treat,mediator,sep=":") %in% attr(model.y@terms,"term.labels") |
         paste(mediator,treat,sep=":") %in% attr(model.y@terms,"term.labels") 
    }
    
      
########
if (class(model.y)[1]!="polr") {#sense whether an ordered outcome model or not?
########

     
    
    B <- sims
    #model.m <- z
    model.y.t <- model.y
    if(class(model.y)[1]!="tobit"){ 
    n.m <-m.data <- model.frame(model.m)  #Call.M$data
    n.y <-y.t.data <- model.frame(model.y.t) #Call.Y$data
    }

    k.t <- ncol(y.t.data)
    k.y <- ncol(n.y)
    k.m <- ncol(n.m)
    k <- k.y + k.m
    n.m <- length(m.data[,1])
    n.y <- length(y.t.data[,1])
    n <- length(y.t.data[,1]) # Isn't this = n.y?
    m <- length(sort(unique(model.frame(model.m)[,1])))
    m.min <- as.numeric(sort(unique(model.frame(model.m)[,1]))[1]) # Seems unnecessarily complicated
   
    if(class(model.m[1])=="gam"){
        test <- class(model.m)[2]
    } else {
        test <- class(model.m)[1]    
    }
if(class(model.y.t)[1]!="vglm") {
            test2 <- class(model.y.t)[1]
            if(class(model.y.t[1])=="gam"){#note how the[1] is on the object as opposed to outside class. this messes up vglm which uses S4 objects; is this a mistake?
                test2 <- class(model.y.t)[2] 
                }  
            if(class(model.y.t)[1]!="gam"){
                test2 <- class(model.y.t)[1]
                }               
    } else {
    test2 <- class(model.y.t)[1]
    tobit.ind<-model.y.t@family@vfamily #indicates whether the vglm model used is a tobit.
    }

    test3 <- class(model.y)[1]

    if(is.factor(y.t.data[,paste(treat)])==TRUE){
        cat.c <- levels(y.t.data[,treat])[1] 
        cat.t <- levels(y.t.data[,treat])[2]
        T.cat <- paste(treat,cat.t, sep="") 
    } else {
        cat.c <- NULL
        cat.t <- NULL
        T.cat <- paste(treat,cat.t, sep="")
    }

    if(n.m != n.y) stop("Error: Missing Values Present in Data")
    
    # Define indexmax function
    indexmax <- function(x){
        n <- length(x)
        imax <- (1:n)[(order(x))[n]]
        imax
    }
    
    if(boot == FALSE){
################################################################################
# @@@@@@@@@@@@@@ Quasi-Bayesian Monte Carlo @@@@@@@@@@@@@@@@@@@
################################################################################
    cat.0 <- treat.0#defines the values of the treatment; 0 and 1 by default
    cat.1 <- treat.1
    MModel.coef <- model.m$coef
    if(test=="polr"){
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
    MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov)
    TMmodel <- mvrnorm(sims, mu=TMmodel.coef, Sigma=TMmodel.var.cov)
    
    ############################################################################
    ##  Mediator Predictions
    ############################################################################
    if(is.factor(m.data[,paste(treat)])==TRUE){
        pred.data.t <- m.data
        pred.data.t[,treat] <- list(factor(unique(m.data[,treat])[1], levels = levels(m.data[,treat])))
        pred.data.c <- m.data
        pred.data.c[,treat] <- list(factor(unique(m.data[,treat])[2], levels = levels(m.data[,treat])))
    } else {
        pred.data.t <- m.data
        pred.data.t[,treat] <- cat.1
        pred.data.c <- m.data
        pred.data.c[,treat] <- cat.0    
    }

    mmat.t <- model.matrix(terms(model.m), data=pred.data.t)
    mmat.c <- model.matrix(terms(model.m), data=pred.data.c)

    if(test=="glm"){  ### CASE 1: GLM
        PredictM1 <- model.m$family$linkinv(MModel %*% t(mmat.t))
        PredictM0 <- model.m$family$linkinv(MModel %*% t(mmat.c))
        
    } else if(test=="polr"){ ### CASE 2: ORDERED -- Handrolled Predictions for polr
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

            for (j in 1:(m-1)) {         # loop to get category-specific probabilities
                cprobs_m1[,j] <- linkfn(lambda[j]-ystar_m1[i,])
                cprobs_m0[,j] <- linkfn(lambda[j]-ystar_m0[i,])  # cumulative probabilities
                probs_m1[,m] <- 1-cprobs_m1[,m-1] # top category 
                probs_m0[,m] <- 1-cprobs_m0[,m-1] # top category
                probs_m1[,1] <- cprobs_m1[,1]     # bottom category 
                probs_m0[,1] <- cprobs_m0[,1]     # bottom category 
            }

            for (j in 2:(m-1)){          # middle categories
                probs_m1[,j] <- cprobs_m1[,j]-cprobs_m1[,j-1]
                probs_m0[,j] <- cprobs_m0[,j]-cprobs_m0[,j-1]
            }
                
            ## Random Prediction Draws:
            #max.pr_m1 <- apply(probs_m1, 1, max)
            #max.pr_m0 <- apply(probs_m0, 1, max)
            #cat.m1 <- apply(probs_m1,1,indexmax)
            #cat.m0 <- apply(probs_m0,1,indexmax)
            #draws_m1 <- round(runif(n, 1, m), 0)
            #draws_m0 <- round(runif(n, 1, m), 0)

            draws_m1 <- matrix(NA, n, m)
            draws_m0 <- matrix(NA, n, m)

            for(ii in 1:n){
                draws_m1[ii,] <- t(rmultinom(1, 1, prob = probs_m1[ii,]))
                draws_m0[ii,] <- t(rmultinom(1, 1, prob = probs_m0[ii,]))
            }
            
            PredictM1[i,] <- apply(draws_m1, 1, indexmax)
            PredictM0[i,] <- apply(draws_m0, 1, indexmax)
        }
        
    } else { ## CASE 3: LINEAR
        sigma <- summary(model.m)$sigma
        error <- rnorm(n, mean=0, sd=sigma)
        PredictM1 <- MModel %*% t(mmat.t)
        PredictM0 <- MModel %*% t(mmat.c)
        PredictM1 <- PredictM1 + error
        PredictM0 <- PredictM0 + error
        rm(error)
    }
    rm(mmat.t, mmat.c)

    ########################################################################
    ##   Outcome Predictions
    ########################################################################
    #Treatment Predictions Data
    Pr1 <- matrix(,nrow=n, ncol=sims)
    Pr0 <- matrix(,nrow=n, ncol=sims)
    
    for(j in 1:sims){
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
        } else {
            pred.data.t[,treat] <- cat.1
            pred.data.c[,treat] <- cat.1
        }
    
        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM1[j,], labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM0[j,], labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM1[j,]
            pred.data.c[,mediator] <- PredictM0[j,]
        }
        if(class(model.y)[1]!="vglm"){
        ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
        ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
        } else {
        ymat.t <- model.matrix(model.y@terms, data=pred.data.t) #there seems to be a problem here in how things work
        ymat.c <- model.matrix(model.y@terms, data=pred.data.c)        
        }
       
        #Treatment Predictions
        #print(as.matrix(TMmodel[j,])) #the second row of this is a Log(scale) parameter from tobit.
        #print(t(ymat.t))#this only has 3 rows
        #print(dim(ymat.t))        
        #print(as.matrix(TMmodel[j,]))        
        #print(dim(as.matrix(TMmodel[j,])))

      
        if(class(model.y)[1]=="vglm"){
                if(tobit.ind=="tobit") {
                TMmodel.tmp<-TMmodel[,-2]
                Pr1[,j] <- t(as.matrix(TMmodel.tmp[j,])) %*% t(ymat.t)
                Pr0[,j] <- t(as.matrix(TMmodel.tmp[j,])) %*% t(ymat.c)
                
                model.y.upper<-model.y@misc$Upper#THIS SHOULD BE MADE INTO A LITTLE FUNCTION
                model.y.lower<-model.y@misc$Lower
                logical<-(Pr1[,j]>model.y.upper)
                Pr1[logical,j]<-model.y.upper
                logical<-(Pr1[,j]<model.y.lower)
                Pr1[logical,j]<-model.y.lower
                logical<-(Pr0[,j]>model.y.upper)
                Pr0[logical,j]<-model.y.upper
                logical<-(Pr0[,j]<model.y.lower)
                Pr0[logical,j]<-model.y.lower
                }        

        } else {

        Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
        Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
        
        }
        
        
        rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
    }    
    
    if(test2=="glm"){
        Pr1 <- apply(Pr1, 2, model.y$family$linkinv);        Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
        delta.1.tmp <- Pr1 - Pr0
    } else {
        delta.1.tmp <- Pr1 - Pr0 #Binary Mediator
    }
    rm(Pr1, Pr0)
            
    #Control Predictions Data
    Pr1 <- matrix(,nrow=n, ncol=sims)
    Pr0 <- matrix(,nrow=n, ncol=sims)
    
    for(j in 1:sims){
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
        } else{
            pred.data.t[,treat] <- cat.0    
            pred.data.c[,treat] <- cat.0    
        }
            
        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM1[j,], labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM0[j,], labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM1[j,]
            pred.data.c[,mediator] <- PredictM0[j,]
        }

        ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
        ymat.c <- model.matrix(terms(model.y), data=pred.data.c)

        #Control Predictions
        if(class(model.y)[1]=="vglm"){
                        if(tobit.ind=="tobit") {
                        TMmodel.tmp<-TMmodel[,-2]
                        Pr1[,j] <- t(as.matrix(TMmodel.tmp[j,])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(TMmodel.tmp[j,])) %*% t(ymat.c)

                        model.y.upper<-model.y@misc$Upper#THIS SHOULD BE MADE INTO A LITTLE FUNCTION
                        model.y.lower<-model.y@misc$Lower
                        logical<-(Pr1[,j]>model.y.upper)
                        Pr1[logical,j]<-model.y.upper
                        logical<-(Pr1[,j]<model.y.lower)
                        Pr1[logical,j]<-model.y.lower
                        logical<-(Pr0[,j]>model.y.upper)
                        Pr0[logical,j]<-model.y.upper
                        logical<-(Pr0[,j]<model.y.lower)
                        Pr0[logical,j]<-model.y.lower
                    }
        } else {
        Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
        Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
        }
        rm(ymat.t, ymat.c)
    }
    
    if(test2=="glm"){
        Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
        Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
        delta.0.tmp <- Pr1 - Pr0
    } else {
        delta.0.tmp <- Pr1 - Pr0 #Binary Mediator
    }
    rm(Pr1, Pr0)
    
    #####################################################################
    ## Direct Effects
    #####################################################################
    #Zeta_1
    Pr1 <- matrix(,nrow=n, ncol=sims)
    Pr0 <- matrix(,nrow=n, ncol=sims)
    
    for(j in 1:sims){
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
        } else {
            pred.data.t[,treat] <- cat.1 #Treatment
            pred.data.c[,treat] <- cat.0  #Control
        }
        
        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM1[j,], labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM1[j,], labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM1[j,]
            pred.data.c[,mediator] <- PredictM1[j,]
        }
            
        ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
        ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
        
        #Direct Predictions
        if(class(model.y)[1]=="vglm"){
                    if(tobit.ind=="tobit") {
                    TMmodel.tmp<-TMmodel[,-2]
                    Pr1[,j] <- t(as.matrix(TMmodel.tmp[j,])) %*% t(ymat.t)
                    Pr0[,j] <- t(as.matrix(TMmodel.tmp[j,])) %*% t(ymat.c)
            
                    model.y.upper<-model.y@misc$Upper#THIS SHOULD BE MADE INTO A LITTLE FUNCTION
                    model.y.lower<-model.y@misc$Lower
                    logical<-(Pr1[,j]>model.y.upper)
                    Pr1[logical,j]<-model.y.upper
                    logical<-(Pr1[,j]<model.y.lower)
                    Pr1[logical,j]<-model.y.lower
                    logical<-(Pr0[,j]>model.y.upper)
                    Pr0[logical,j]<-model.y.upper
                    logical<-(Pr0[,j]<model.y.lower)
                    Pr0[logical,j]<-model.y.lower
                    }
        } else {
        Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
        Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
        }
        #rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
    }    
    
    if(test2=="glm"){
        Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
        Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
        zeta.1.tmp <- Pr1 - Pr0
    } else {
        zeta.1.tmp <- Pr1 - Pr0 #Binary Mediator
    }
    rm(Pr1, Pr0)
        
    #Zeta_0
    Pr1 <- matrix(,nrow=n, ncol=sims)
    Pr0 <- matrix(,nrow=n, ncol=sims)
    
    for(j in 1:sims){
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
        } else {
            pred.data.t[,treat] <- cat.1
            pred.data.c[,treat] <- cat.0
        }
        
        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM0[j,], labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM0[j,], labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM0[j,]
            pred.data.c[,mediator] <- PredictM0[j,]
        }
            
        ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
        ymat.c <- model.matrix(terms(model.y), data=pred.data.c)

        #Direct Predictions
        if(class(model.y)[1]=="vglm"){
                    if(tobit.ind=="tobit") {
                    TMmodel.tmp<-TMmodel[,-2]
                    Pr1[,j] <- t(as.matrix(TMmodel.tmp[j,])) %*% t(ymat.t)
                    Pr0[,j] <- t(as.matrix(TMmodel.tmp[j,])) %*% t(ymat.c)
                    
                    model.y.upper<-model.y@misc$Upper#THIS SHOULD BE MADE INTO A LITTLE FUNCTION
                    model.y.lower<-model.y@misc$Lower
                    logical<-(Pr1[,j]>model.y.upper)
                    Pr1[logical,j]<-model.y.upper
                    logical<-(Pr1[,j]<model.y.lower)
                    Pr1[logical,j]<-model.y.lower
                    logical<-(Pr0[,j]>model.y.upper)
                    Pr0[logical,j]<-model.y.upper
                    logical<-(Pr0[,j]<model.y.lower)
                    Pr0[logical,j]<-model.y.lower
                    }
        } else {        
        Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
        Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
        }
        rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
    }
    #print(summary(Pr1[,j])) #why only one summary gets printed?
    if(test2=="glm"){
        Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
        Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
        zeta.0.tmp <- Pr1 - Pr0
    } else {
        zeta.0.tmp <- Pr1 - Pr0 #Binary Mediator
    }
    
    rm(Pr1, Pr0, PredictM1, PredictM0, TMmodel, MModel) 
    zeta.1 <- t(as.matrix(apply(zeta.1.tmp, 2, mean)))
    rm(zeta.1.tmp)
    zeta.0 <- t(as.matrix(apply(zeta.0.tmp, 2, mean)))
    rm(zeta.0.tmp)
    delta.1 <- t(as.matrix(apply(delta.1.tmp, 2, mean)))
    #rm(delta.1.tmp)
    delta.0 <- t(as.matrix(apply(delta.0.tmp, 2, mean)))
    #rm(delta.0.tmp)
    tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
    
    
    } else {
################################################################################
# @@@@@@@@@@@@@@ Nonparametric Bootstrap @@@@@@@@@@@@@@@@@@@
################################################################################

    n <- n.m
    Call.M <- model.m$call
    if(isS4(model.y.t)){
    Call.Y.t <- model.y.t@call
    } else {
    Call.Y.t <- model.y.t$call
    }
    
    if(is.factor(y.t.data[,treat])==TRUE){
        cat.0 <- levels(y.t.data[,treat])[1]
        cat.1 <- levels(y.t.data[,treat])[2]
    } else {
        cat.0 <- treat.0
        cat.1 <- treat.1
    }

    #Storage
    delta.1 <- matrix(NA, B, 1)
    delta.0 <- matrix(NA, B, 1)
    zeta.1 <- matrix(NA, B, 1)
    zeta.0 <- matrix(NA, B, 1)
    tau <- matrix(NA, B, 1)
    for(b in 1:B){
        index <- sample(1:n,n, repl=TRUE)
        Call.M$data <- m.data[index,]
        Call.Y.t$data <- y.t.data[index,]
      
        if(test=="polr"){
            if(length(unique(y.t.data[index,mediator]))!=m){
                cat("Insufficient Variation on Mediator")
                break
            }
        }
   #print("here tobit error") #after X iterations, where X can vary, we get the following: Error in lm.fit(X_vlm, z_vlm, ...) :   NA/NaN/Inf in foreign function call (arg 4)
        #printing the call objects doesn't indicate an issue.
        #Refit Models with Resampled Data
        new.fit.M <- eval.parent(Call.M)
        new.fit.t <- eval.parent(Call.Y.t)
        #print(eval.parent(Call.M))
        #print(eval.parent(Call.Y.t))
        #Generate Mediation Model Predictions
        pred.data.t <- m.data
        pred.data.t[,treat] <- cat.1
        #pred.data.t[,control] <- cat.0
        pred.data.c <- m.data
        pred.data.c[,treat] <- cat.0
        #pred.data.c[,control] <- cat.1
      
        if(is.factor(m.data[,treat])==TRUE){
            pred.data.t[,treat] <- as.factor(pred.data.t[,treat])
            pred.data.c[,treat] <- as.factor(pred.data.c[,treat])
        } else { 
            pred.data.t[,treat] <- as.numeric(pred.data.t[,treat])
            pred.data.c[,treat] <- as.numeric(pred.data.c[,treat])
        } 
#print(b)
        if(test=="glm"){
            PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
            PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)
        } else if(test=="polr") {
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
            
            #max.pr_m1 <- apply(probs_m1, 1, max)
            #max.pr_m0 <- apply(probs_m0, 1, max)
            #cat.m1 <- predict(new.fit.M, newdata=pred.data.t)
            #cat.m0 <- predict(new.fit.M, newdata=pred.data.c)
            #draws_m1 <- round(runif(n, m.min, m), 0)
            #draws_m0 <- round(runif(n, m.min, m), 0)
            #PredictM1 <- ifelse(max.pr_m1 > runif(n, 0, .75), draws_m1, cat.m1)
            #PredictM0 <- ifelse(max.pr_m0 > runif(n, 0, .75), draws_m0, cat.m0)
        
        } else {
            if(class(model.m)[1]=="gam"){
                sigma <- summary(new.fit.M)$scale
            } else {
                #print(summary(new.fit.M))#when of class rq (quantile regression) there is not such sigma variable
                #print(names(summary(new.fit.M))), this produces  "call"         "terms"        "coefficients" "rdf"          "tau" 
                sigma <- summary(new.fit.M)$sigma
            }
            if(class(model.m)[1]!="rq"){
            error <- rnorm(n, mean=0, sd=sigma)
            PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t) + error
            PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c) + error
            rm(error)
            } else {
            #for quantile regression we are unable to introduce error into the prediction. Hence sample quantiles from 0 1 and then predict.         
            tau.list<-runif(n)
                    PredictM1<-PredictM0<-vector(length=n)
                    for(zz in 1:n){
                    PredictM1[zz] <- predict(update(new.fit.M, tau=tau.list[zz]), type="response", newdata=pred.data.t)
                    PredictM0[zz] <- predict(update(new.fit.M, tau=tau.list[zz]), type="response", newdata=pred.data.c)
                    
                    }
            
            }
        }

       
        ############################################################################
        #Treatment Predictions Data
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
            
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
                pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
                pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            if(is.null(control)!=TRUE){
                    pred.data.t[,control] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
                    pred.data.c[,control] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            }
        } else{
            pred.data.t[,treat] <- cat.1    
            pred.data.c[,treat] <- cat.1
            if(is.null(control)!=TRUE){
                pred.data.t[,control] <- cat.0
                pred.data.c[,control] <- cat.0
            }
        } 
            
        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM1, labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM0, labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM1
            pred.data.c[,mediator] <- PredictM0
        }
        
               
        #Treatment Predictions
        if(test3=="rq"){#Why not test2?
            pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none")
            pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none")
        } else {
            pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
            pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
        }
        pr.mat <- as.matrix(cbind(pr.1, pr.0))
        delta.1.tmp <- pr.mat[,1] - pr.mat[,2]
        
        rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)

        ############################################################################
        #Control Predictions Data
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
            
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            if(is.null(control)!=TRUE){
                pred.data.t[,control] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
                pred.data.c[,control] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            }
        } else{
            pred.data.t[,treat] <- cat.0
            pred.data.c[,treat] <- cat.0
            if(is.null(control)!=TRUE){
                pred.data.t[,control] <- cat.1
                pred.data.c[,control] <- cat.1
            }
        } 

        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM1, labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM0, labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM1
            pred.data.c[,mediator] <- PredictM0
        }
        
        #Controls Predictions
        if(test3=="rq"){
            pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none")
            pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
        } else {
            pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
            pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
        }
        pr.mat <- as.matrix(cbind(pr.1, pr.0))
        delta.0.tmp <-pr.mat[,1] - pr.mat[,2]
        
        rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
       
        ########################################################################
        #Direct Effects 
        #Zeta.1
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            if(is.null(control)!=TRUE){
                pred.data.t[,control] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
                pred.data.c[,control] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            }
        } else{
            pred.data.t[,treat] <- cat.1
            pred.data.c[,treat] <- cat.0
            if(is.null(control)!=TRUE){
                pred.data.t[,control] <- cat.0
                pred.data.c[,control] <- cat.1
            }
        }
        
        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM1, labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM1, labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM1
            pred.data.c[,mediator] <- PredictM1
        }
            
        if(test3=="rq"){
            pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none")
            pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none")
        } else {
            pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
            pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
        }
        pr.mat <- as.matrix(cbind(pr.1, pr.0))
        zeta.1.tmp <- pr.mat[,1] - pr.mat[,2]
          
        rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)
        
        #Zeta.0
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            if(is.null(control)!=TRUE){
                pred.data.t[,control] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
                pred.data.c[,control] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            }
        } else{
            pred.data.t[,treat] <- cat.1
            pred.data.c[,treat] <- cat.0
            if(is.null(control)!=TRUE){
                pred.data.t[,control] <- cat.0
                pred.data.c[,control] <- cat.1
            }
        }
        
        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM0, labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM0, labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM0
            pred.data.c[,mediator] <- PredictM0
        }
        
        if(test3=="rq"){
            pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none")
            pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none")
        } else {
            pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
            pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
        }
        pr.mat <- as.matrix(cbind(pr.1, pr.0))
        zeta.0.tmp <-pr.mat[,1] - pr.mat[,2]
            
        zeta.1[b] <- mean(zeta.1.tmp)
        zeta.0[b] <- mean(zeta.0.tmp)
        delta.1[b] <- mean(delta.1.tmp)
        delta.0[b] <- mean(delta.0.tmp)
        tau[b] <- (zeta.1[b] + zeta.0[b] + delta.0[b] + delta.1[b])/2
    } #bootstrap loop ends
    } #nonpara boot branch ends
    
    d0 <- mean(delta.0)
    d1 <- mean(delta.1)
    low <- (1 - conf.level)/2; high <- 1 - low
    d0.ci <- quantile(delta.0,c(low,high), na.rm=TRUE)
    d1.ci <- quantile(delta.1,c(low,high), na.rm=TRUE)
    tau.coef <- mean(tau)
    tau.ci <- quantile(tau,c(low,high), na.rm=TRUE)
    avg.delta <- (d0 + d1)/2
    pct.dist <- avg.delta/tau
    pct.coef <- median(pct.dist)
    pct.ci <- quantile(pct.dist,c(low,high), na.rm=TRUE)
    z1 <- mean(zeta.1)    
    z1.ci <- quantile(zeta.1,c(low,high), na.rm=TRUE)
    z0 <- mean(zeta.0)    
    z0.ci <- quantile(zeta.0,c(low,high), na.rm=TRUE)
    
    out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci, d0.sims=delta.0, d1.sims=delta.1,
        pct.coef=pct.coef, pct.ci=pct.ci,
        tau.coef=tau.coef, tau.ci=tau.ci, z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
        boot=boot, treat=treat, mediator=mediator, INT=INT, conf.level=conf.level,
        model.y=model.y, model.m=model.m)
    class(out) <- "mediate"
    out
    
    
    
    
    
    
    
    
    
    #ordered outcome 
    ##########
    } else {
    ##########
    boot<-TRUE
    
      
    
        B <- sims
    model.y.t <- model.y
    m.data <- model.frame(model.m)  #Call.M$data
    y.t.data <- model.frame(model.y.t) #Call.Y$data
    n <- length(y.t.data[,1])
    m <- length(sort(unique(model.frame(model.y)[,1])))
    #m.min <- as.numeric(sort(unique(model.frame(model.m)[,1]))[1]) #What Do I Need This For?
    form.y <- formula(model.y)
    form.m <- formula(model.m)
    cat.0 <- treat.0
    cat.1 <- treat.1
    
    if(class(model.m[1])=="gam"){
        test <- class(model.m)[2]
        } else {
        test <- class(model.m)[1]   
            }
            
    if(class(model.y.t[1])=="gam"){
        test2 <- class(model.y.t)[2]
        } else {
        test2 <- class(model.y.t)[1]    
            }
    test3 <- class(model.y)[1]
    if(is.factor(y.t.data[,paste(treat)])==TRUE){
                cat.c <- levels(y.t.data[,treat])[1] 
                cat.t <- levels(y.t.data[,treat])[2]
                T.cat <- paste(treat,cat.t, sep="") 
            } else {
                cat.c <- NULL
                cat.t <- NULL
                T.cat <- paste(treat,cat.t, sep="")
                }

    #Storage - Now Dynamic
    delta.1 <- matrix(NA, B, m)
    delta.0 <- matrix(NA, B, m)
    zeta.1 <- matrix(NA, B, m)
    zeta.0 <- matrix(NA, B, m)
    tau <- matrix(NA, B, m)

    for(b in 1:B){#Bootstrap Loop
        if(test=="polr"){
            if(length(unique(y.t.data[index,mediator]))!=m){
            cat("Insufficient Variation on Mediator")
            break
            }
            }
        
        #Resampling Step
        data.star.m <- m.data[sample(1:nrow(m.data),n,replace=TRUE),] #Pull off data.star.
        data.star.y <- y.t.data[sample(1:nrow(y.t.data),n,replace=TRUE),] #Pull off data.star.
        new.fit.M <- update(model.m, data=data.star.m)
        new.fit.t <- update(model.y, data=data.star.y)
           

        #Generate Mediation Model Predictions
        pred.data.t <- m.data
        pred.data.t[,treat] <- cat.1
        pred.data.c <- m.data
        pred.data.c[,treat] <- cat.0
        
        
        if(is.factor(m.data[,treat])==TRUE){
        pred.data.t[,treat] <- as.factor(pred.data.t[,treat])
        pred.data.c[,treat] <- as.factor(pred.data.c[,treat])
        } else { 
        pred.data.t[,treat] <- as.numeric(pred.data.t[,treat])
        pred.data.c[,treat] <- as.numeric(pred.data.c[,treat])
        } 
                
        if(test=="glm"){
        PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
        PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)
       
       #above new.fit.M was from a linear model...
        } else if(test=="polr") {
        probs_m1 <- predict(new.fit.M, newdata=pred.data.t, type="probs")
        probs_m0 <- predict(new.fit.M, newdata=pred.data.c, type="probs")
        draws_m1 <- matrix(NA, n, m)
        draws_m0 <- matrix(NA, n, m)

        for(ii in 1:n){ 
            draws_m1[ii,] <- t(rmultinom(1, 1, prob = probs_m1[ii,]))
            draws_m0[ii,] <- t(rmultinom(1, 1, prob = probs_m0[ii,]))
            }
        max.pr_m1 <- apply(probs_m1, 1, max)
        max.pr_m0 <- apply(probs_m0, 1, max)
        cat.m1 <- predict(new.fit.M, newdata=pred.data.t)
        cat.m0 <- predict(new.fit.M, newdata=pred.data.c)
        draws_m1 <- round(runif(n, m.min, m), 0)
        draws_m0 <- round(runif(n, m.min, m), 0)
        #where is .75 from?
        PredictM1 <- ifelse(max.pr_m1 > runif(n, 0, .75), draws_m1, cat.m1)
        PredictM0 <- ifelse(max.pr_m0 > runif(n, 0, .75), draws_m0, cat.m0)
        } else {
        if(class(model.m)[1]=="gam"){
            sigma <- summary(new.fit.M)$scale
            } else {
            sigma <- summary(new.fit.M)$sigma
                }
        error <- rnorm(n, mean=0, sd=sigma)
        PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t) + error
        PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c) + error
        rm(error)
        }
        
    
#####################################################################################       
        #Treatment Predictions Data
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data

    if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
    if(is.null(control)!=TRUE){
            pred.data.t[,control] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            pred.data.c[,control] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            }
        } else{
        #please annotate what these loops are doing...cat.1 etc were just 1 or 0, so what were the above doing? allowing for multivalued treatments?
            pred.data.t[,treat] <- cat.1    
            pred.data.c[,treat] <- cat.1
            if(is.null(control)!=TRUE){
            pred.data.t[,control] <- cat.0
            pred.data.c[,control] <- cat.0
            }
        } 
        
    if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
            pred.data.t[,mediator] <- factor(PredictM1, labels = levels(y.t.data[,mediator]))
            pred.data.c[,mediator] <- factor(PredictM0, labels = levels(y.t.data[,mediator]))
        } else {
            pred.data.t[,mediator] <- PredictM1
            pred.data.c[,mediator] <- PredictM0
    }           
        #Treatment Predictions
        probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
        probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
                        
        delta.1.tmp <-  probs_p1 - probs_p0     
        rm(pred.data.t, pred.data.c, probs_p1, probs_p0)

        #Control Predictions Data
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data

            if(is.factor(y.t.data[,paste(treat)])==TRUE){
            pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
    if(is.null(control)!=TRUE){
            pred.data.t[,control] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            pred.data.c[,control] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            }
    } else{
            pred.data.t[,treat] <- cat.0    
            pred.data.c[,treat] <- cat.0
    if(is.null(control)!=TRUE){
            pred.data.t[,control] <- cat.1
            pred.data.c[,control] <- cat.1
            }
        } 

if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
    pred.data.t[,mediator] <- factor(PredictM1, labels = levels(y.t.data[,mediator]))
    pred.data.c[,mediator] <- factor(PredictM0, labels = levels(y.t.data[,mediator]))
            } else {
    pred.data.t[,mediator] <- PredictM1
    pred.data.c[,mediator] <- PredictM0
    }
            
        #Controls Predictions   
        probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
        probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")

        delta.0.tmp <- probs_p1 - probs_p0

        rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
########################################################################
        #Direct Effects 
########################################################################    
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
        
        #Zeta.1
if(is.factor(y.t.data[,paste(treat)])==TRUE){
    pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
    pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
    if(is.null(control)!=TRUE){
            pred.data.t[,control] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            pred.data.c[,control] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            }
    } else{
    pred.data.t[,treat] <- cat.1    
    pred.data.c[,treat] <- cat.0
    if(is.null(control)!=TRUE){
            pred.data.t[,control] <- cat.0
            pred.data.c[,control] <- cat.1
            }
        } 

        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
    pred.data.t[,mediator] <- factor(PredictM1, labels = levels(y.t.data[,mediator]))
    pred.data.c[,mediator] <- factor(PredictM1, labels = levels(y.t.data[,mediator]))
            } else {
    pred.data.t[,mediator] <- PredictM1
    pred.data.c[,mediator] <- PredictM1
    }

            
        probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
        probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
        zeta.1.tmp <- probs_p1 - probs_p0

        rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
        
        #Zeta.0
        pred.data.t <- y.t.data
        pred.data.c <- y.t.data
        
        if(is.factor(y.t.data[,paste(treat)])==TRUE){
    pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
    pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
    if(is.null(control)!=TRUE){
            pred.data.t[,control] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
            pred.data.c[,control] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
            }
    } else{
    pred.data.t[,treat] <- cat.1    
    pred.data.c[,treat] <- cat.0
    if(is.null(control)!=TRUE){
            pred.data.t[,control] <- cat.0
            pred.data.c[,control] <- cat.1
            }
        } 

        if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
    pred.data.t[,mediator] <- factor(PredictM0, labels = levels(y.t.data[,mediator]))
    pred.data.c[,mediator] <- factor(PredictM0, labels = levels(y.t.data[,mediator]))
            } else {
    pred.data.t[,mediator] <- PredictM0
    pred.data.c[,mediator] <- PredictM0
    }

        #Predictions
        probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
        probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")

        zeta.0.tmp <- probs_p1 - probs_p0
        rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
        
        zeta.1[b,] <- apply(zeta.1.tmp, 2, median)
        zeta.0[b,] <- apply(zeta.0.tmp, 2, median)
        delta.1[b,] <- apply(delta.1.tmp, 2, median)
        delta.0[b,] <- apply(delta.0.tmp, 2, median)
        tau[b,] <- delta.0[b,] + zeta.1[b,]
        } 
        
    #Results 
    low <- (1 - conf.level)/2; high <- 1 - low
    d0 <- apply(delta.0, 2, mean)
    d1 <- apply(delta.1, 2, mean)
    d0.ci <- apply(delta.0, 2, quantile, c(low,high))
    d1.ci <- apply(delta.1, 2, quantile, c(low,high))
    
    z1 <- apply(zeta.1, 2, mean)
    z0 <- apply(zeta.0,2, mean)     
    z1.ci <- apply(zeta.1,2, quantile, c(low,high))
    z0.ci <- apply(zeta.0,2, quantile, c(low,high))
    
    tau.coef <- apply(tau, 2, mean)
    tau.ci <- apply(tau, 2, quantile, c(low,high))
    
    #avg.delta <- (d0 + d1)/2
    #pct.dist <- avg.delta/tau
    #pct.coef <- median(pct.dist)
    #pct.ci <- quantile(pct.dist,c(.025,.975), na.rm=TRUE)
    m.lab <- sort(unique(levels(model.frame(model.y)[,1])))
    
out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci, d0.sims=delta.0, d1.sims=delta.1, z1=z1, z0=z0, z1.ci=z1.ci, z0.ci=z0.ci,tau.coef=tau.coef, tau.ci=tau.ci, treat=treat, mediator=mediator, INT=INT, m.lab=m.lab)
class(out) <- "mediate.order"
out

    ###
    #Else closing loop for ordered outcome
    }
    ###
}

print.mediate <- function(x, ...){
    print(unlist(x[1:11]))
    invisible(x)
}

summary.mediate <- function(object, ...)
    structure(object, class = c("summary.mediate", class(object)))
    
    

print.summary.mediate <- function(x, ...){
    clp <- 100 * x$conf.level
    if(x$INT==TRUE){
        cat("\n Causal Mediation Analysis \n\n")
    if(x$boot==TRUE){
        cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
        } else {
        cat("Quasi-Bayesian Confidence Intervals\n\n")
        }
    cat("Mediation Effect_0: ", format(x$d0, digits=4), clp, "% CI ", format(x$d0.ci, digits=4), "\n")
    cat("Mediation Effect_1: ", format(x$d1, digits=4), clp, "% CI ", format(x$d1.ci, digits=4), "\n")
    cat("Direct Effect_0: ", format(x$z0, digits=4), clp, "% CI ", format(x$z0.ci, digits=4), "\n")
    cat("Direct Effect_1: ", format(x$z1, digits=4), clp, "% CI ", format(x$z1.ci, digits=4), "\n")
    cat("Total Effect: ", format(x$tau.coef, digits=4), clp, "% CI ", format(x$tau.ci, digits=4), "\n")
    cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),clp, "% CI ", format(x$pct.ci, digits=4),"\n")
    #cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"\n")
        } else {
            cat("\n Causal Mediation Analysis \n\n")
            if(x$boot==TRUE){
        cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
        } else {
        cat("Quasi-Bayesian Confidence Intervals\n\n")
        }
    cat("Mediation Effect: ", format(x$d1, digits=4), clp, "% CI ", format(x$d1.ci, digits=4), "\n")
    cat("Direct Effect: ", format(x$z0, digits=4), clp, "% CI ", format(x$z0.ci, digits=4), "\n")
    cat("Total Effect: ", format(x$tau.coef, digits=4), clp, "% CI ", format(x$tau.ci, digits=4), "\n")
    cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),clp, "% CI ", format(x$pct.ci, digits=4),"\n")
    #cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"\n")
                }
    invisible(x)
}





################
#SUMMARY FUNCTIONS FOR ORDERED OUTCOME CASE
################

print.mediate.order <- function(x, ...){
    print(unlist(x[1:50]))
    invisible(x)
    }

summary.mediate.order <- function(object, ...)
    structure(object, class = c("summary.mediate.order", class(object)))

print.summary.mediate.order <- function(x, ...){
    
    tab.d0 <- rbind(x$d0, x$d0.ci)
    tab.d1 <- rbind(x$d1, x$d1.ci)
    tab.z0 <- rbind(x$d0, x$d0.ci)
    tab.z1 <- rbind(x$d1, x$d1.ci)
    tab.tau <- rbind(x$tau.coef, x$tau.ci)
    #Outcome Table Labels
    out.names <- c()

        for(i in 1:length(x$m.lab)){
            out.names.tmp <- paste("Pr(Y=",x$m.lab[i],")",sep="")
            out.names <- c(out.names, out.names.tmp)
            }
    #Label Tables
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


#actually, given the nonlinear models its possible for d0 not = d1 even if no interaction. we have to fix this in the regular code too.
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
    


#function the processes output of mediate model.
plot.process<-function(model) {
coef.vec.1 <- c(model$d1, model$z1, model$tau.coef)
lower.vec.1<- c(model$d1.ci[1], model$z1.ci[1],model$tau.ci[1])
upper.vec.1<- c(model$d1.ci[2], model$z1.ci[2],model$tau.ci[2])
range.1<-range(model$d1.ci[1], model$z1.ci[1],model$tau.ci[1],model$d1.ci[2], model$z1.ci[2],model$tau.ci[2] )

coef.vec.0 <- c(model$d0, model$z0, model$tau.coef)
lower.vec.0<- c(model$d0.ci[1], model$z0.ci[1],model$tau.ci[1])
upper.vec.0<- c(model$d0.ci[2], model$z0.ci[2],model$tau.ci[2])
range.0<-range(model$d0.ci[1], model$z0.ci[1],model$tau.ci[1],model$d0.ci[2], model$z0.ci[2],model$tau.ci[2] )
return(list(coef.vec.1=coef.vec.1, lower.vec.1=lower.vec.1, upper.vec.1=upper.vec.1, coef.vec.0=coef.vec.0, lower.vec.0=lower.vec.0, upper.vec.0=upper.vec.0, range.1=range.1, range.0=range.0) )
}

#plotting function that takes output of plot.process
plot.mediate<-function(param,title=NULL,both=FALSE,reset=TRUE,xlim=NULL, lwd =  1.5, xlab="",cex=1,cex.main=1.5,cex.axis=1.5,cex.lab=1.5) {
var.names<-c("ACME","Direct Effect","Total Effect")
y.axis <- c(length(param$coef.vec.1):.5)#create indicator for y.axis, descending so that R orders vars from top to bottom on y-axis
#pdf("MediationPlots.pdf", width=8.5, height=11, onefile=FALSE, paper="special")
if(reset){
par(mfrow=c(1,1))#reset graphical window
}
par(mar=c(3, 5.5, 2, 0))#set margins for plot, leaving lots of room on left-margin (2nd number in margin command) for variable names

if(is.null(xlim)){
if(both) {
xlim<-range(param$range.1,param$range.0)
} else {
xlim<-param$range.1
}
} else {
xlim<-xlim
}
plot(param$coef.vec.1, y.axis, type = "p", axes = F, xlab = xlab, ylab = "", pch = 19, cex = cex,#plot coefficients as points, turning off axes and labels. 
    xlim = xlim, xaxs = "r", main = "") #set limits of x-axis so that they include mins and maxs of 
segments(param$lower.vec.1, y.axis, param$upper.vec.1, y.axis, lwd=lwd)#coef +/-1.96*se = 95% interval, lwd adjusts line thickness      
if(both) {
points(param$coef.vec.0, y.axis-.1, type = "p",pch = 1, cex = cex#,#plot coefficients as points, turning off axes and labels. 
    )
segments(param$lower.vec.0, y.axis-.1, param$upper.vec.0, y.axis-.1, lwd =  lwd, lty=3)#coef +/-1.96*se = 95% interval, lwd adjusts line thickness      
    
}
axis(1, at = seq(-2,2,by=.5), labels =  seq(-2,2,by=.5), tick = T,#draw x-axis and labels with tick marks
    cex.axis = cex.axis, cex.lab=cex.lab, mgp = c(2,.5,0))#reduce label size, moves labels closer to tick marks
axis(2, at = y.axis, label = var.names, las = 1, tick = T, 
    cex.axis = cex.axis) #draw y-axis with tick marks, make labels perpendicular to axis and closer to axis
segments(0,0,0,length(var.names),lty=2) # draw dotted line through 0
mtext(xlab,side=1,line=1.3,cex=cex.lab)
ti<-title
title(ti)
}



##Assume list of mediators, treatments as character vector
##Assume datasets as named list in same order as treatment
##Covariate specification as a single character vector, e.g., covariates <- c("gender_p+Trust_GSS_pre+ideo7")
mediations <- function(datasets, treatment, mediators, covariates, family=c("gaussian", "gaussian"),LowerY=NULL,UpperY=NULL, interaction=FALSE, conf.level=.95, sims=500) {
  covariates<-covariates
  data <- names(datasets)
  labels <- c()
  out <- list()
  count<-1
  for (i in 1:length(treatment)) {
    d1 <- sprintf("datasets$%s", data[i])
    dataarg<- eval(parse(text=d1))
    for (j in 1:length(mediators)) {
      f1 <- sprintf("%s ~ %s + %s", mediators[j], treatment[i], covariates)
      f2 <- sprintf("Y ~ %s + %s + %s", treatment[i], mediators[j], covariates)
      if (interaction==TRUE) {
        f2 <- sprintf("Y ~ %s*%s + %s", treatment[i], mediators[j], covariates)
      }
      fmla <- as.formula(paste(f1))
      fmla2 <- as.formula(paste(f2))
      if(family[1] == "binomial") {
        result1 <- glm(fmla, family=binomial("probit"), data=dataarg)
      } else {
        result1 <- glm(fmla, family=family[1], data=dataarg)
      }
      if(family[2] == "binomial") {
        result2 <- glm(fmla2, family=binomial("probit"), data=dataarg)
      } else if(family[2] == "tobit") {
        result2 <- vglm(fmla2, tobit(Lower=LowerY,Upper=UpperY), data=dataarg)    
      } else {  
        result2 <- glm(fmla2, family=family[2], data=dataarg)
      }
      #print(result1)
      #print(result2)
      #if (interaction==TRUE) out[[(count)]] <- mediate(result1, result2, sims=sims, treat=treatment[i], mediator=mediators[j], INT=TRUE, conf.level=conf.level)
      #if (interaction==FALSE) out[[(count)]] <- mediate(result1, result2, sims=sims, treat=treatment[i], mediator=mediators[j], conf.level=conf.level)
      out[[(count)]] <- mediate(result1, result2, sims=sims, treat=treatment[i], mediator=mediators[j], conf.level=conf.level)
      rm(result1)
      rm(result2)
      labels[(count)] <- sprintf("%s.%s", data[i], mediators[j])   
      count<-count+1 
    }
  }
  names(out) <- labels

                    #generates what is necessary for plotting and an additional output.  
                    names<-names(out)
                    plot<-list()
                    for(i in 1:length(names)){
                    plot[[i]]<-plot.process(out[[i]])
                    }
                    names(plot)<-names
  
  return(list(out=out,plot=plot))
}
      
    
  
