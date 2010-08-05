mediate <- function(model.m, model.y, sims=1000, boot=FALSE, INT=FALSE, treat="treat.name", mediator="med.name", control=NULL){
	
	
#######################################################################
# Code here extracts relevant parts of models and dimensions them.	
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
	
######################################################################
# Extract Model Class  GAM require special extraction
	
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
	
######################################################################
#Set Treatment levels for factor	
	
	if(is.factor(y.t.data[,paste(treat)])==TRUE){
				cat.c <- levels(y.t.data[,treat])[1] 
				cat.t <- levels(y.t.data[,treat])[2]
				T.cat <- paste(treat,cat.t, sep="") 
			} else {
				cat.c <- NULL
				cat.t <- NULL
				T.cat <- paste(treat,cat.t, sep="")
				}
######################################################################
	if(n.m != n.y){ #Error Message for Missing Values
	stop("Error: Missing Values Present in Data")
	} else {
	if(boot == FALSE){  #Nonbootstrap code starts here.
	cat.0 <- 0 #Treatment and control category settings for predictions
	cat.1 <- 1
######################################################################
# Draws from Multivariate Normal - Extract Coefs and Var-Cov Matrices
	MModel.coef <- model.m$coef
	if(test=="polr"){
	k <- length(model.m$coef)
	MModel.var.cov <- vcov(model.m)[(1:k),(1:k)]
	} else {
	MModel.var.cov <- vcov(model.m)
		}
	TMmodel.coef <- model.y$coef
	TMmodel.var.cov <- vcov(model.y)
	MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov) #Actual MvtNormal Draws
	TMmodel <- mvrnorm(sims, mu=TMmodel.coef, Sigma=TMmodel.var.cov)
	
	
#####################################################################################
##  Mediator Predictions
#####################################################################################
#First we create two new model frames. One is the original model frame with observed treatment
#status is replaced with treatment set to treated status and control status.
#Basic complication is doing this when treatment is a factor variable.


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

#####################################################################################
#Predictions For Supported Model Types

#Logit and Probit Predictions	
	if(test=="glm"){
	PredictM1 <- model.m$family$linkinv(MModel %*% t(mmat.t))
	PredictM0 <- model.m$family$linkinv(MModel %*% t(mmat.c))
		} else if(test=="polr"){ #Handrolled Predictions for polr
  if(model.m$method=="logit"){
    linkfn <- plogis
  	} else {
  	linkfn <- pnorm
  	}

#Ordered Probit and Ordered Logit Predictions
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
	} else { #Predictions for Regresion Model
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
#Again we generate two data sets.  One where mediation status
#done as M(t=0) and a second where M(t=1).  In both instances T itself is not set to 1.
#This provides and estimate of delta(t=1).
	#Treatment Predictions Data
	Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
#Again must account for treatment variables that are factors.	
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
	
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
########################################################################
	
	#Treatment Predictions under linear scale
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
	if(test2=="glm"){
		Pr1 <- apply(Pr1, 2, model.y$family$linkinv); #Convert to probability scale for binary outcome.		Pr0 <- apply(Pr0, 2, model.y$family$linkinv);
		delta.1.tmp <- Pr1 - Pr0;
		} else {
		delta.1.tmp <- Pr1 - Pr0 #Take difference of these two predictions for distribution of point estimtes
			}
	rm(Pr1, Pr0)

#######################################################################
#This process now repeats.  Again we gernate two data sets:ne where mediation status
#done as M(t=0) and a second where M(t=1).  In both instances T itself is now set to 0.
#This gives us our estimate of delta(t=0).

	#Control Predictions Data
	Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
	
		for(j in 1:sims){
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data
	if(is.factor(y.t.data[,paste(treat)])==TRUE){ #Take care of factors
	pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
	pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
	} else{
	pred.data.t[,treat] <- cat.0	
	pred.data.c[,treat] <- cat.0	
		} 
		
	if(is.factor(y.t.data[,paste(mediator)])==TRUE) { #Take care of factors
	pred.data.t[,mediator] <- factor(PredictM1[j,], labels = levels(y.t.data[,mediator]))
	pred.data.c[,mediator] <- factor(PredictM0[j,], labels = levels(y.t.data[,mediator]))
			} else {
	pred.data.t[,mediator] <- PredictM1[j,]
	pred.data.c[,mediator] <- PredictM0[j,]
	}

	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)

	#Control Predictions Under Linear Scale
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c)
	}
	
	if(test2=="glm"){
		Pr1 <- apply(Pr1, 2, model.y$family$linkinv); #Apply link function to get predictions in probability scale.		Pr0 <- apply(Pr0, 2, model.y$family$linkinv);
		delta.0.tmp <- Pr1 - Pr0;
		} else {
		delta.0.tmp <- Pr1 - Pr0 #Difference of two predictions is distribution of point estimates.
			}
	rm(Pr1, Pr0)
#####################################################################
#Direct Effects
#####################################################################
#Process is the same, we just change settings of mediatior status and treatment status.
#For direct effect under t=1.  For first direct effect prediction is done under M(t=1)
#while T is set alternately to 1 and 0.  Point estimate distribution is difference in these
#two predictions.

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
	pred.data.t[,treat] <- cat.1	#Key is Here not treatment status varies instead of being held constant.
	pred.data.c[,treat] <- cat.0    #Control
	}
	#Notice here the mediatior is fixed at M(t=1) for both predictions.
	if(is.factor(y.t.data[,paste(mediator)])==TRUE) {  #But not mediator prediction stays constant, in this case at M(t=1).
	pred.data.t[,mediator] <- factor(PredictM1[j,], labels = levels(y.t.data[,mediator]))
	pred.data.c[,mediator] <- factor(PredictM1[j,], labels = levels(y.t.data[,mediator]))
	} else {
	pred.data.t[,mediator] <- PredictM1[j,]  
	pred.data.c[,mediator] <- PredictM1[j,]
		}
	#Datasets for prediction.	
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
	
	#Actual Predictions under each dataset. 
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	#rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
		if(test2=="glm"){ #For binary case transform to probability scale
		Pr1 <- apply(Pr1, 2, model.y$family$linkinv);		Pr0 <- apply(Pr0, 2, model.y$family$linkinv);
		zeta.1.tmp <- Pr1 - Pr0;
		} else {
		zeta.1.tmp <- Pr1 - Pr0 #Take difference in predictions.
			}
	rm(Pr1, Pr0)	
		
	#Zeta_0 For second direct effect prediction is done under M(t=0)
#while T is set alternately to 1 and 0.  Point estimate distribution is difference in these
#two predictions.

	Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
	
	#Take care of factor treatment variables.
	for(j in 1:sims){
	pred.data.t <- y.t.data
	pred.data.c <- y.t.data
	if(is.factor(y.t.data[,paste(treat)])==TRUE){
	pred.data.t[,treat] <- list(factor(unique(y.t.data[,treat])[1], levels = levels(y.t.data[,treat])))
	pred.data.c[,treat] <- list(factor(unique(y.t.data[,treat])[2], levels = levels(y.t.data[,treat])))
	} else {
	pred.data.t[,treat] <- cat.1 #Alternate treatment status.
	pred.data.c[,treat] <- cat.0	
		}
	
	if(is.factor(y.t.data[,paste(mediator)])==TRUE) {
	pred.data.t[,mediator] <- factor(PredictM0[j,], labels = levels(y.t.data[,mediator]))
	pred.data.c[,mediator] <- factor(PredictM0[j,], labels = levels(y.t.data[,mediator]))
	} else {
	pred.data.t[,mediator] <- PredictM0[j,]
	pred.data.c[,mediator] <- PredictM0[j,]
		}
	#Final Data for preiction.
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)

	#Actual predictions.
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
	if(test2=="glm"){#Convert to probability scale for binary outcomes.
		Pr1 <- apply(Pr1, 2, model.y$family$linkinv);		Pr0 <- apply(Pr0, 2, model.y$family$linkinv);
		zeta.0.tmp <- Pr1 - Pr0; #Take difference in outcome predictions.
		} else {
		zeta.0.tmp <- Pr1 - Pr0 #Take difference in outcome predictions.			}
	rm(Pr1, Pr0, PredictM1, PredictM0, TMmodel, MModel)	
	#Now average over individual level predictions.
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
		
	#Principles are exactly the same as before in terms of the predictions.  
	#Key difference is now we do this one step at a time instead of all at once.
	#Each time we resample from the data with replace and do the predictions under
	#this model.
	n <- n.m
	#Model calls save all the information about a specific model object.
	#Using the model call allows the software
	Call.M <- model.m$call
	Call.Y.t <- model.y.t$call
	
#######################################################################
#Again taking care of factors.  Here saving the factor level names.
	if(is.factor(y.t.data[,treat])==TRUE){

		cat.0 <- levels(y.t.data[,treat])[1]
		cat.1 <- levels(y.t.data[,treat])[2]
		} else {
		
		cat.0 <- 0
		cat.1 <- 1
		}

	#Storage
	delta.1 <- matrix(NA, B, 1)
	delta.0 <- matrix(NA, B, 1)
	zeta.1 <- matrix(NA, B, 1)
	zeta.0 <- matrix(NA, B, 1)
	tau <- matrix(NA, B, 1)
	
	#The bootstrap loop starts here.
	for(b in 1:B){
		index <- sample(1:n,n, repl=TRUE) #Sample from data row indicators
		#Replace data in model call with resampled data.
		Call.M$data <- m.data[index,]
		Call.Y.t$data <- y.t.data[index,]
		
		if(test=="polr"){
		if(length(unique(y.t.data[index,mediator]))!=m){
			cat("Insufficient Variation on Mediator")
			break
			}
			}

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M) #This resimates the model with the new resample data set.
		new.fit.t <- eval.parent(Call.Y.t)

		#Everything preceds exactly as before except we can now rely on the predict command to generate the predictions.
		
		#Set up data for mediator model predictions
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
		
		#Data set up complete #Now generate model predictions.
		if(test=="glm"){
		PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
		PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)
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
		
####################################################################################
#Now we generate outcome predictions


		#Set up Treatment Predictions Data
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data
		#Set factors to appropriate levels.
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
 #Data setup Complete
 
 #Generate Outcome model Predictons
		#Treatment Predictions
		if(test3=="rq"){
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none");		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
			} else {
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
				}
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1.tmp <- pr.mat[,1] - pr.mat[,2] #Take difference in predictions.

		rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)

		#Process Repeats for Control Predictions Data
		#Set Up Data
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
			
		#Generate Predictions
		if(test3=="rq"){
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none");		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
			} else {
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
				}	
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.0.tmp <-pr.mat[,1] - pr.mat[,2] #Take Difference in Predictions.

		rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
########################################################################
		#Direct Effects 
########################################################################		#Same Process for Direct Effects.

		#Zeta.1
		#Set Up Data.
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
	pred.data.t[,treat] <- cat.1 #Alternate Treatment status.
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
	pred.data.t[,mediator] <- PredictM1 #Predictions under M(t=1)
	pred.data.c[,mediator] <- PredictM1
	}
	#Data Set up Complete
	
	#Generate Predictions		
		if(test3=="rq"){
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none");		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
			} else {
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
				}
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		zeta.1.tmp <- pr.mat[,1] - pr.mat[,2] #Take Difference in predictions.

		rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)
		
		#Now predictions under M(t=0)
		#Zeta.0
		#Set up Data for Predictions
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
	#Data Setup Complete
	
	#Generate Predictions
		if(test3=="rq"){
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t, interval="none");		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c, interval="none");
			} else {
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
				}
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		zeta.0.tmp <-pr.mat[,1] - pr.mat[,2] #Take differnce in predictions
		
		#Average over individual level predictions
		zeta.1[b] <- mean(zeta.1.tmp)
		zeta.0[b] <- mean(zeta.0.tmp)		
		delta.1[b] <- mean(delta.1.tmp)
		delta.0[b] <- mean(delta.0.tmp)
		tau[b] <- (zeta.1[b] + zeta.0[b] + delta.0[b] + delta.1[b])/2
		} #bootstrap loop
	} #nonpara boot branch
	
	#Now average over simulation distribution of point esimates
	d0 <- mean(delta.0)
	d1 <- mean(delta.1)
	#Quantiles for Confidence Intervals.
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

# Saving model attributes
call.y <- model.y$call
call.m <- model.m$call
env.y <- environment(model.y$terms)
env.m <- environment(model.m$terms)

out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci, pct.coef=pct.coef, pct.ci=pct.ci, 
tau.coef=tau.coef, tau.ci=tau.ci, z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
boot=boot, treat=treat, mediator=mediator, INT=INT, call.m=call.m, call.y=call.y, env.m=env.m, env.y=env.y)
class(out) <- "mediate"
out

}

#Everything from here on out are other functions that print a tables of results
#The key difference is that one gets two tables if an interaction is included.

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
			cat("\n Causal Mediation Analysis \n\n")
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
	


