mediate <- function(model.m, model.y, sims=1000, boot=FALSE, INT=FALSE, T="treat.name", M="med.name", C=NULL){
	B <- sims
	#model.m <- z
	model.y.t <- model.y
	
	m.data <- model.frame(model.m)  #Call.M$data
    y.t.data <- model.frame(model.y.t) #Call.Y$data
    k.t <- ncol(y.t.data)
	k.m <- ncol(m.data)
	n.m <- model.frame(model.m)
	n.y <- model.frame(model.y)	
	k.y <- ncol(n.y)
	k.m <- ncol(n.m)
	k <- k.y + k.m
	n.m <- length(m.data[,1])
	n.y <- length(y.t.data[,1])
	n <- length(y.t.data[,1])
	m <- length(sort(unique(model.frame(model.m)[,1])))
	m.min <- as.numeric(sort(unique(model.frame(model.m)[,1]))[1])
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
	if(is.factor(y.t.data[,paste(T)])==TRUE){
				cat.c <- levels(y.t.data[,T])[1] 
				cat.t <- levels(y.t.data[,T])[2]
				T.cat <- paste(T,cat.t, sep="") 
			} else {
				cat.c <- NULL
				cat.t <- NULL
				T.cat <- paste(T,cat.t, sep="")
				}

	if(n.m != n.y){
	stop("Error: Missing Values Present in Data")
	} else {
	if(boot == FALSE){ 
	cat.0 <- 0
	cat.1 <- 1
	MModel.coef <- model.m$coef
	if(test=="polr"){
	k <- length(model.m$coef)
	MModel.var.cov <- vcov(model.m)[(1:k),(1:k)]
	} else {
	MModel.var.cov <- vcov(model.m)
		}
	TMmodel.coef <- model.y$coef
	TMmodel.var.cov <- vcov(model.y)
	MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov)
	TMmodel <- mvrnorm(sims, mu=TMmodel.coef, Sigma=TMmodel.var.cov)
#####################################################################################
##  Mediator Predictions
#####################################################################################
	if(is.factor(m.data[,paste(T)])==TRUE){
	pred.data.t <- m.data
	pred.data.t[,T] <- list(factor(unique(m.data[,T])[1], levels = levels(m.data[,T])))
	pred.data.c <- m.data
	pred.data.c[,T] <- list(factor(unique(m.data[,T])[2], levels = levels(m.data[,T])))
	} else {
	pred.data.t <- m.data
	pred.data.t[,T] <- cat.1
	pred.data.c <- m.data
	pred.data.c[,T] <- cat.0	
		}
	mmat.t <- model.matrix(terms(model.m), data=pred.data.t)
	mmat.c <- model.matrix(terms(model.m), data=pred.data.c)
	
	if(test=="glm"){
	PredictM1 <- model.m$family$linkinv(MModel %*% t(mmat.t))
	PredictM0 <- model.m$family$linkinv(MModel %*% t(mmat.c))
		} else if(test=="polr"){ #Handrolled Predictions for polr
  if(model.m$method=="logit"){
    linkfn <- plogis
  	} else {
  	linkfn <- pnorm
  	}

indexmax <- function(x){
    n <- length(x)
    imax <- (1:n)[(order(x))[n]]
    imax
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
#Random Prediction Draws
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
########################################################################
##   Outcome Predictions
########################################################################
	#Treatment Predictions Data
	Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
	
	for(j in 1:sims){
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data
	if(is.factor(y.t.data[,paste(T)])==TRUE){
	pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	} else {
	pred.data.t[,T] <- cat.1
	pred.data.c[,T] <- cat.1
		}
	
	if(is.factor(y.t.data[,paste(M)])==TRUE) {
		pred.data.t[,M] <- factor(PredictM1[j,], labels = levels(y.t.data[,M]))
		pred.data.c[,M] <- factor(PredictM0[j,], labels = levels(y.t.data[,M]))
		} else {
	pred.data.t[,M] <- PredictM1[j,]
	pred.data.c[,M] <- PredictM0[j,]
	}
	
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
	
	#Treatment Predictions
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
	if(test2=="glm"){
		Pr1 <- apply(Pr1, 2, model.y$family$linkinv);		Pr0 <- apply(Pr0, 2, model.y$family$linkinv);
		delta.1.tmp <- Pr1 - Pr0;
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
	if(is.factor(y.t.data[,paste(T)])==TRUE){
	pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	} else{
	pred.data.t[,T] <- cat.0	
	pred.data.c[,T] <- cat.0	
		} 
		
	if(is.factor(y.t.data[,paste(M)])==TRUE) {
	pred.data.t[,M] <- factor(PredictM1[j,], labels = levels(y.t.data[,M]))
	pred.data.c[,M] <- factor(PredictM0[j,], labels = levels(y.t.data[,M]))
			} else {
	pred.data.t[,M] <- PredictM1[j,]
	pred.data.c[,M] <- PredictM0[j,]
	}

	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)

	#Control Predictions
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c)
	}
	
	if(test2=="glm"){
		Pr1 <- apply(Pr1, 2, model.y$family$linkinv);		Pr0 <- apply(Pr0, 2, model.y$family$linkinv);
		delta.0.tmp <- Pr1 - Pr0;
		} else {
		delta.0.tmp <- Pr1 - Pr0 #Binary Mediator
			}
	rm(Pr1, Pr0)
#####################################################################
#Direct Effects
#####################################################################
	#Zeta_1
    Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
	
	for(j in 1:sims){
	pred.data.t <- y.t.data
	pred.data.c <- y.t.data
	if(is.factor(y.t.data[,paste(T)])==TRUE){
	pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	} else {
	pred.data.t[,T] <- cat.1	#Treatment
	pred.data.c[,T] <- cat.0    #Control
	}
	
	if(is.factor(y.t.data[,paste(M)])==TRUE) {
	pred.data.t[,M] <- factor(PredictM1[j,], labels = levels(y.t.data[,M]))
	pred.data.c[,M] <- factor(PredictM1[j,], labels = levels(y.t.data[,M]))
	} else {
	pred.data.t[,M] <- PredictM1[j,]
	pred.data.c[,M] <- PredictM1[j,]
		}
		
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
	
	#Direct Predictions
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	#rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
		if(test2=="glm"){
		Pr1 <- apply(Pr1, 2, model.y$family$linkinv);		Pr0 <- apply(Pr0, 2, model.y$family$linkinv);
		zeta.1.tmp <- Pr1 - Pr0;
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
	if(is.factor(y.t.data[,paste(T)])==TRUE){
	pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	} else {
	pred.data.t[,T] <- cat.1
	pred.data.c[,T] <- cat.0	
		}
	
	if(is.factor(y.t.data[,paste(M)])==TRUE) {
	pred.data.t[,M] <- factor(PredictM0[j,], labels = levels(y.t.data[,M]))
	pred.data.c[,M] <- factor(PredictM0[j,], labels = levels(y.t.data[,M]))
	} else {
	pred.data.t[,M] <- PredictM0[j,]
	pred.data.c[,M] <- PredictM0[j,]
		}
		
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)

	#Direct Predictions
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
	if(test2=="glm"){
		Pr1 <- apply(Pr1, 2, model.y$family$linkinv);		Pr0 <- apply(Pr0, 2, model.y$family$linkinv);
		zeta.0.tmp <- Pr1 - Pr0;
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
	delta.0 <- t(as.matrix(apply(delta.0.tmp, 2, mean)));
	#rm(delta.0.tmp)
	tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2

#######################################################################
	} else { #@@@@@@@@@@@@@@Nonparametric Bootstrap@@@@@@@@@@@@@@@@@@@
	n <- n.m
	Call.M <- model.m$call
	Call.Y.t <- model.y.t$call
	
	if(is.factor(y.t.data[,T])==TRUE){

		cat.0 <- levels(y.t.data[,T])[1]
		cat.1 <- levels(y.t.data[,T])[2]
		} else {
		
		cat.0 <- 0
		cat.1 <- 1
		}
		
	indexmax <- function(x){
    n <- length(x)
    imax <- (1:n)[(order(x))[n]]
    imax
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
		if(length(unique(y.t.data[index,M]))!=m){
			cat("Insufficient Variation on Mediator")
			break
			}
			}

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M)
		new.fit.t <- eval.parent(Call.Y.t)

		#Generate Mediation Model Predictions
		pred.data.t <- m.data
		pred.data.t[,T] <- cat.1
		#pred.data.t[,C] <- cat.0
		pred.data.c <- m.data
		pred.data.c[,T] <- cat.0
		#pred.data.c[,C] <- cat.1
		
		if(is.factor(m.data[,T])==TRUE){
		pred.data.t[,T] <- as.factor(pred.data.t[,T])
		pred.data.c[,T] <- as.factor(pred.data.c[,T])
		} else { 
		pred.data.t[,T] <- as.numeric(pred.data.t[,T])
		pred.data.c[,T] <- as.numeric(pred.data.c[,T])
		} 
		
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
		PredictM1 <- apply(draws_m1, 1, indexmax);		PredictM0 <- apply(draws_m0, 1, indexmax);
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
			sigma <- summary(new.fit.M)$sigma
				}
		error <- rnorm(n, mean=0, sd=sigma)
		PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t) + error;		PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c) + error
		rm(error)
		}
		
#####################################################################################		
		#Treatment Predictions Data
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data
		
	if(is.factor(y.t.data[,paste(T)])==TRUE){
			pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
			pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	if(is.null(C)!=TRUE){
			pred.data.t[,C] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
			pred.data.c[,C] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
			}
		} else{
			pred.data.t[,T] <- cat.1	
			pred.data.c[,T] <- cat.1
			if(is.null(C)!=TRUE){
			pred.data.t[,C] <- cat.0
			pred.data.c[,C] <- cat.0
			}
		} 
		
	if(is.factor(y.t.data[,paste(M)])==TRUE) {
			pred.data.t[,M] <- factor(PredictM1, labels = levels(y.t.data[,M]))
			pred.data.c[,M] <- factor(PredictM0, labels = levels(y.t.data[,M]))
		} else {
			pred.data.t[,M] <- PredictM1
			pred.data.c[,M] <- PredictM0
	}

			
		#Treatment Predictions
		if(test3=="rq"){#Why not test2?
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none");		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
			} else {
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
				}
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)

		#Control Predictions Data
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data
		
	if(is.factor(y.t.data[,paste(T)])==TRUE){
			pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
			pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	if(is.null(C)!=TRUE){
			pred.data.t[,C] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
			pred.data.c[,C] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
			}
	} else{
			pred.data.t[,T] <- cat.0	
			pred.data.c[,T] <- cat.0
	if(is.null(C)!=TRUE){
			pred.data.t[,C] <- cat.1
			pred.data.c[,C] <- cat.1
			}
		} 

if(is.factor(y.t.data[,paste(M)])==TRUE) {
	pred.data.t[,M] <- factor(PredictM1, labels = levels(y.t.data[,M]))
	pred.data.c[,M] <- factor(PredictM0, labels = levels(y.t.data[,M]))
			} else {
	pred.data.t[,M] <- PredictM1
	pred.data.c[,M] <- PredictM0
	}
			
		#Controls Predictions	
		if(test3=="rq"){
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none");		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
			} else {
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
				}	
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.0.tmp <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
########################################################################
		#Direct Effects 
########################################################################		
		#Zeta.1
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data
	if(is.factor(y.t.data[,paste(T)])==TRUE){
	pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	if(is.null(C)!=TRUE){
			pred.data.t[,C] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
			pred.data.c[,C] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
			}
	} else{
	pred.data.t[,T] <- cat.1	
	pred.data.c[,T] <- cat.0
	if(is.null(C)!=TRUE){
			pred.data.t[,C] <- cat.0
			pred.data.c[,C] <- cat.1
			}
		} 

		if(is.factor(y.t.data[,paste(M)])==TRUE) {
	pred.data.t[,M] <- factor(PredictM1, labels = levels(y.t.data[,M]))
	pred.data.c[,M] <- factor(PredictM1, labels = levels(y.t.data[,M]))
			} else {
	pred.data.t[,M] <- PredictM1
	pred.data.c[,M] <- PredictM1
	}
			
		if(test3=="rq"){
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none");		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
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
	if(is.factor(y.t.data[,paste(T)])==TRUE){
	pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	if(is.null(C)!=TRUE){
			pred.data.t[,C] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
			pred.data.c[,C] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
			}
	} else{
	pred.data.t[,T] <- cat.1	
	pred.data.c[,T] <- cat.0
	if(is.null(C)!=TRUE){
			pred.data.t[,C] <- cat.0
			pred.data.c[,C] <- cat.1
			}
		} 

		if(is.factor(y.t.data[,paste(M)])==TRUE) {
	pred.data.t[,M] <- factor(PredictM0, labels = levels(y.t.data[,M]))
	pred.data.c[,M] <- factor(PredictM0, labels = levels(y.t.data[,M]))
			} else {
	pred.data.t[,M] <- PredictM0
	pred.data.c[,M] <- PredictM0
	}
		if(test3=="rq"){
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none");		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
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
		} #bootstrap loop
	} #nonpara boot branch
	d0 <- mean(delta.0)
	d1 <- mean(delta.1)
	d0.ci <- quantile(delta.0,c(.025,.975), na.rm=TRUE)
	d1.ci <- quantile(delta.1,c(.025,.975), na.rm=TRUE)
	tau.coef <- mean(tau)
	tau.ci <- quantile(tau,c(.025,.975), na.rm=TRUE)
	avg.delta <- (d0 + d1)/2
	pct.dist <- avg.delta/tau
	pct.coef <- median(pct.dist)
	pct.ci <- quantile(pct.dist,c(.025,.975), na.rm=TRUE)
	z1 <- mean(zeta.1)	
	z1.ci <- quantile(zeta.1,c(.025,.975), na.rm=TRUE)
	z0 <- mean(zeta.0)	
	z0.ci <- quantile(zeta.0,c(.025,.975), na.rm=TRUE)
		
 }

out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci, pct.coef=pct.coef, pct.ci=pct.ci, tau.coef=tau.coef, tau.ci=tau.ci, z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,boot=boot, INT=INT)
class(out) <- "mediate"
out

}

print.mediate <- function(x, ...){
	print(unlist(x[1:11]))
	invisible(x)
	}

summary.mediate <- function(object, ...)
	structure(object, class = c("summary.mediate", class(object)))
 
print.summary.mediate <- function(x, ...){
	if(x$INT==TRUE){
		cat("\n Causal Mediation Analysis \n\n")
	if(x$boot==TRUE){
		cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
		} else {
		cat("Quasi-Bayesian Confidence Intervals\n\n")
		}
	cat("Mediation Effect_0: ", format(x$d0, digits=4), "95% CI ", format(x$d0.ci, digits=4), "\n")
	cat("Mediation Effect_1: ", format(x$d1, digits=4), "95% CI ", format(x$d1.ci, digits=4), "\n")
	cat("Direct Effect_0: ", format(x$z0, digits=4), "95% CI ", format(x$z0.ci, digits=4), "\n")
	cat("Direct Effect_1: ", format(x$z1, digits=4), "95% CI ", format(x$z1.ci, digits=4), "\n")
	cat("Total Effect: ", format(x$tau.coef, digits=4), "95% CI ", format(x$tau.ci, digits=4), "\n")
	cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"95% CI ", format(x$pct.ci, digits=4),"\n")
    #cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"\n")
		} else {
			cat("\n Test For Mediation Effect \n\n")
			if(x$boot==TRUE){
		cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
		} else {
		cat("Quasi-Bayesian Confidence Intervals\n\n")
		}
	cat("Mediation Effect: ", format(x$d1, digits=4), "95% CI ", format(x$d1.ci, digits=4), "\n")
	cat("Direct Effect: ", format(x$z0, digits=4), "95% CI ", format(x$z0.ci, digits=4), "\n")
	cat("Total Effect: ", format(x$tau.coef, digits=4), "95% CI ", format(x$tau.ci, digits=4), "\n")
    cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"95% CI ", format(x$pct.ci, digits=4),"\n")
    #cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"\n")
				}
	invisible(x)
	}
	


