mediate.binary <- function(z, ...){
	UseMethod("mediate.binary", z)
	}
	
mediate.binary.default <- function(z, model.y, sims=1000, boot=FALSE){
	B <- sims
	model.m <- z
	model.y.t <- model.y
	n.m <- model.frame(model.m)
	n.y <- model.frame(model.y)	
	k.y <- ncol(n.y)
	k.m <- ncol(n.m)
	k <- k.y + k.m
	n.m <- length(n.m[,1])
	n.y <- length(n.y[,1])
	n <- n.m
	sigma <- summary(model.m)$sigma
	test <- class(model.m)[1]
	if(n.m != n.y){
	cat("Error: Missing Values Present in Data")
	} else {
	
	if (test == "glm"){ #Binary Mediator
		if (boot == FALSE){ #Parametric Bootstrap
				MModel.coef <- model.m$coef
	MModel.var.cov <- vcov(model.m)
	TMmodel.coef <- model.y$coef
	TMmodel.var.cov <- vcov(model.y)
	MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov)
	TMmodel <- mvrnorm(sims, mu=TMmodel.coef, Sigma=TMmodel.var.cov)

 if(k.m == 2){ #No Covariates in Mediation Model
			#Generate Predictions From Mediation Model with Random  Error
 	PredictM1temp <- as.matrix(apply(MModel, 1, sum))
 	PredictM1_t  <- matrix(PredictM1temp, n, sims, byrow=TRUE)
 	PredictM1_t <- pnorm(PredictM1_t) 
 	PredictM0temp <- MModel[,1] 
 	PredictM0_t <- matrix(PredictM0temp, n, sims, byrow=TRUE)
 	PredictM0_t <- pnorm(PredictM0_t)	
 	PredictM1 <- matrix(,nrow=n, ncol=sims)
 	PredictM0 <- matrix(,nrow=n, ncol=sims)
			
			for (j in 1:sims) {
        		#error <- runif(n, min=0, max=1)
        		PredictM1[,j] <- as.numeric(PredictM1_t[,j] > runif(n,min=0, max=1))
        		PredictM0[,j] <- as.numeric(PredictM0_t[,j] > runif(n,min=0, max=1))
 				}

 		}
 	if(k.y==3){ #No Covariates
	#Storage Matrices For Outcome Predictions
 	Prob1 <- matrix(,nrow=n, ncol=sims)
	Prob0 <- matrix(,nrow=n, ncol=sims)
	#Prediction For Outcome Based on Mediation Model Preductions
	#T=1
	for (j in 1:sims) {
       				 Prob1[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j] 
       			 	 Prob0[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM0[,j]
					}

	delta.1.tmp <- Prob1 - Prob0 
	rm(Prob1, Prob0)
	
	#Storage Matrices
	Prob1 <- matrix(,nrow=n, ncol=sims)
	Prob0 <- matrix(,nrow=n, ncol=sims)
	#T=0
	for (j in 1:sims) {
       			 	Prob1[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM1[,j]
        		 	Prob0[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM0[,j]
					}

	delta.0.tmp <- Prob1 - Prob0 
	rm(Prob1, Prob0)
			
	#Storage Matrices
	Prob1_temp <- matrix(,nrow=n, ncol=sims)
	Prob0_temp <- matrix(,nrow=n, ncol=sims)

	#Calculate Total Effect
	for (j in 1:sims) {
       			 Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j]  
        		Prob0_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM0[,j]  
        		}

	tau.tmp <- Prob1_temp - Prob0_temp
	} else {#Predictions For Covariates
		PredictM1temp <- as.matrix(apply(MModel, 1, sum))
		PredictM1_t  <- matrix(PredictM1temp, n, sims, byrow=TRUE)
		PredictM1_t <- pnorm(PredictM1_t) 
		PredictM0temp <-  cbind(MModel[,1], MModel[,3:k.m]) 
		PredictM0temp <-  as.matrix(apply(PredictM0temp, 1, sum)) 
		PredictM0_t <- matrix(PredictM0temp, n, sims, byrow=TRUE)
		PredictM0_t <- pnorm(PredictM0_t)	
		PredictM1 <- matrix(,nrow=n, ncol=sims)
		PredictM0 <- matrix(,nrow=n, ncol=sims)
			
			for (j in 1:sims) {
        		#error <- runif(n, min=0, max=1)
        		PredictM1[,j] <- as.numeric(PredictM1_t[,j] > runif(n,min=0, max=1))
        		PredictM0[,j] <- as.numeric(PredictM0_t[,j] > runif(n,min=0, max=1))
 				}
 			
	 	TMmodel.1 <- as.matrix(TMmodel[,1:3])
		TMmodel.2 <- as.matrix(TMmodel[,4:k.y])
		ii <- 4
   		X <- model.frame(model.y)
    	X <- as.matrix(X[,ii:k.y])
    	X.pred <- matrix( , n, sims)
    	k.x <- ncol(X)
   		for (i in 1:k.x){
    			for (j in 1:sims) {
        			X.pred[,j] <- TMmodel.2[j,i]*X[,i]
    				}
    			}
    	
    	#Storage Matrices For Outcome Predictions
		Prob1 <- matrix(,nrow=n, ncol=sims)
		Prob0 <- matrix(,nrow=n, ncol=sims)
		#Prediction For Outcome Based on Mediation Model Preductions
		#T=1
		for (j in 1:sims) {
        		Prob1[,j] <- TMmodel.1[j,1] + TMmodel.1[j,2]*1 + TMmodel.1[j,3]*PredictM1[,j] 
        		Prob0[,j] <- TMmodel.1[j,1] + TMmodel.1[j,2]*1 + TMmodel.1[j,3]*PredictM0[,j]
				}
			
		Prob1 <- Prob1 + X.pred
		Prob0 <- Prob0 + X.pred
		delta.1.tmp <- Prob1 - Prob0 
		rm(Prob1, Prob0)
	
		#Storage Matrices
		Prob1 <- matrix(,nrow=n, ncol=sims)
		Prob0 <- matrix(,nrow=n, ncol=sims)
		#T=0
		for (j in 1:sims) {
       			 Prob1[,j] <- TMmodel.1[j,1] + TMmodel.1[j,2]*0 + TMmodel.1[j,3]*PredictM1[,j]
        		 Prob0[,j] <- TMmodel.1[j,1] + TMmodel.1[j,2]*0 + TMmodel.1[j,3]*PredictM0[,j]
				}
		Prob1 <- Prob1 + X.pred
		Prob0 <- Prob0 + X.pred
		delta.0.tmp <- Prob1 - Prob0 
		rm(Prob1, Prob0)
			
		#Storage Matrices
		Prob1_temp <- matrix(,nrow=n, ncol=sims)
		Prob0_temp <- matrix(,nrow=n, ncol=sims)

		#Calculate Total Effect
			for (j in 1:sims) {
        		Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j]  
        		Prob0_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM0[,j]  
        		}
		Prob1_temp <- Prob1_temp + X.pred
		Prob0_temp <- Prob0_temp + X.pred
		tau.tmp <- Prob1_temp - Prob0_temp
		}
		delta.1 <- t(as.matrix(apply(delta.1.tmp, 2 , mean)))
		delta.0 <- t(as.matrix(apply(delta.0.tmp, 2 , mean)))
		tau <- t(as.matrix(apply(tau.tmp, 2, mean)))
			} else { #@@@@@@@@@@@@@@@@@Nonparametric Bootstrap@@@@@@@@@@@@@@@@@@@@
				n <- n.m
	k.t <- length(names(model.y.t$coef))
	k.m <- length(names(model.m$coef))
	Call.M <- model.m$call
	Call.Y.t <- model.y.t$call
	m.data <- model.frame(model.m)
    y.t.data <- model.frame(model.y.t)
	#Storage
	delta.1 <- matrix(NA, B, 1)
	delta.0 <- matrix(NA, B, 1)
	tau <- matrix(NA, B, 1)

    if(k.m == 2){ #No Covariates in Mediation Model
		for (b in 1:B) {
		#Resample Data
		index <- sample(1:n,n, repl=TRUE)
		Call.M$data <- m.data[index,]
		Call.Y.t$data <- y.t.data[index,]

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M)
		new.fit.t <- eval.parent(Call.Y.t)

		#Generate Mediation Model Predictions Without Covariates
		sigma <- summary(new.fit.M)$sigma
		PredictM1temp <- pnorm(sum(new.fit.M$coef[1:k.m])) 
		PredictM1 <- as.numeric(PredictM1temp > runif(n, min=0, max=1))
		PredictM0temp <- pnorm(new.fit.M$coef[1]) 
		PredictM0 <- as.numeric(PredictM0temp > runif(n, min=0, max=1))
		
		#Treatment Predictions
		pred.data.t <- data.frame(1,PredictM1)
		pred.data.c <- data.frame(1,PredictM0)
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0)

		pred.data.t <- data.frame(0,PredictM1)
		pred.data.c <- data.frame(0,PredictM0)
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		delta.0.tmp <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0)
		#Calculate Total Effect
		pred.data.t <- data.frame(1,PredictM1)
		pred.data.c <- data.frame(0,PredictM0)
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])
		
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		tau.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, PredictM1, PredictM0)
		
		delta.1[b] <- mean(delta.1.tmp)
		delta.0[b] <- mean(delta.0.tmp)
		tau[b] <- mean(tau.tmp)
		}
	} else {
	#With Covariates
	for (b in 1:B) {
		
		index <- sample(1:n,n, repl=TRUE)
		Call.M$data <- m.data[index,]
		Call.Y.t$data <- y.t.data[index,]

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M)
		new.fit.t <- eval.parent(Call.Y.t)

		#Generate Mediation Model Predictions
		sigma <- summary(new.fit.M)$sigma
		PredictM1temp <- pnorm(sum(new.fit.M$coef[1:k.m])) 
		PredictM1 <- as.numeric(PredictM1temp > runif(n, min=0, max=1))
		PredictM0temp <- pnorm(new.fit.M$coef[1] + sum(new.fit.M$coef[3:k.m])) 
		PredictM0 <- as.numeric(PredictM0temp > runif(n, min=0, max=1))
		
		#Treatment Predictions
		pred.data.t <- data.frame(1,PredictM1, y.t.data[index,3:k.t])
		pred.data.c <- data.frame(1,PredictM0, y.t.data[index,3:k.t])

		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)

		pred.data.t <- data.frame(0,PredictM1, y.t.data[index,3:k.t])
		pred.data.c <- data.frame(0,PredictM0, y.t.data[index,3:k.t])
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.0.tmp <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)

		#Calculate Total Effect
		pred.data.t <- data.frame(1,PredictM1, y.t.data[index,3:k.t])
		pred.data.c <- data.frame(0,PredictM0, y.t.data[index,3:k.t])
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])
		
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		tau.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, PredictM1, PredictM0)

		delta.1[b] <- mean(delta.1.tmp)
		delta.0[b] <- mean(delta.0.tmp)
		tau[b] <- mean(tau.tmp)
			}
		  }
		} 
		} else { #@@@@@@@@@@@@@@@@Binary Outcomes
			if(boot == FALSE){ #Parametric Bootstrap
			#Draw From Posterior of Each Model
			MModel.coef <- model.m$coef
			MModel.var.cov <- vcov(model.m)
			TMmodel.coef <- model.y$coef
			TMmodel.var.cov <- vcov(model.y)
			MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov)
			TMmodel <- mvrnorm(sims, mu=TMmodel.coef, Sigma=TMmodel.var.cov)
	
		if(k.m == 2){ #No Covariates in Mediation Model
			#Generate Predictions From Mediation Model with Random  Error
			error <- matrix(,nrow=n, ncol=sims)
			mean <- apply(MModel, 1, sum)
			for (j in 1:sims) {
        		error[,j] <- rnorm(n, mean[j], sd=sigma)
 				}		
 			PredictM1 <- mean + error
 			PredictM0 <- MModel[,1] + error		
			} else { #With Covariates
			#Generate Predictions From Mediation Model with Random  Error
			error <- matrix(,nrow=n, ncol=sims)
			mean <- apply(MModel, 1, sum)
			for (j in 1:sims) {
        		error[,j] <- rnorm(n, mean[j], sd=sigma)
 				}
 				
 			PredictM1 <- mean + error
 			meanM0 <- cbind(MModel[,1], MModel[,3:k.m]) 
 			meanM0 <- apply(meanM0, 1, sum)
 			PredictM0 <- meanM0 + error
			}
			if(k.y==3){ #@@@@@@@@@@@@@@No Covariates@@@@@@@@@@@@@
				#Storage Matrices For Outcome Predictions
				Prob1 <- matrix(,nrow=n, ncol=sims)
				Prob0 <- matrix(,nrow=n, ncol=sims)
				#Prediction For Outcome Based on Mediation Model Preductions
				#T=1
				for (j in 1:sims) {
       				 Prob1[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j] 
       			 	 Prob0[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM0[,j]
					}
			
				Prob1 <- apply(Prob1, 2, pnorm)
				Prob0 <- apply(Prob0, 2, pnorm)
				delta.1.tmp <- Prob1 - Prob0 
				rm(Prob1, Prob0)
	
				#Storage Matrices
				Prob1 <- matrix(,nrow=n, ncol=sims)
				Prob0 <- matrix(,nrow=n, ncol=sims)
				#T=0
				for (j in 1:sims) {
       			 	Prob1[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM1[,j]
        		 	Prob0[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM0[,j]
					}
			Prob1 <- apply(Prob1, 2, pnorm)
			Prob0 <- apply(Prob0, 2, pnorm)
			delta.0.tmp <- Prob1 - Prob0 
			rm(Prob1, Prob0)
			
			#Storage Matrices
			Prob1_temp <- matrix(,nrow=n, ncol=sims)
			Prob0_temp <- matrix(,nrow=n, ncol=sims)

			#Calculate Total Effect
			for (j in 1:sims) {
       			 Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j]  
        		Prob0_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM0[,j]  
        		}
			Prob1_temp <- apply(Prob1_temp, 2, pnorm)
			Prob0_temp <- apply(Prob0_temp, 2, pnorm)
			tau.tmp <- Prob1_temp - Prob0_temp
				} else {#Predictions For Covariates
			TMmodel.1 <- as.matrix(TMmodel[,1:3])
			TMmodel.2 <- as.matrix(TMmodel[,4:k.y])
			ii <- 4
   			X <- model.frame(model.y)
    		X <- as.matrix(X[,ii:k.y])
    		X.pred <- matrix( , n, sims)
    		k.x <- ncol(X)
   			 for (i in 1:k.x){
    			for (j in 1:sims) {
        			X.pred[,j] <- TMmodel.2[j,i]*X[,i]
    				}
    			}
    	
    		#Storage Matrices For Outcome Predictions
			Prob1 <- matrix(,nrow=n, ncol=sims)
			Prob0 <- matrix(,nrow=n, ncol=sims)
			#Prediction For Outcome Based on Mediation Model Preductions
			#T=1
			for (j in 1:sims) {
        		Prob1[,j] <- TMmodel.1[j,1] + TMmodel.1[j,2]*1 + TMmodel.1[j,3]*PredictM1[,j] 
        		Prob0[,j] <- TMmodel.1[j,1] + TMmodel.1[j,2]*1 + TMmodel.1[j,3]*PredictM0[,j]
				}
			
			Prob1 <- Prob1 + X.pred
			Prob0 <- Prob0 + X.pred
			Prob1 <- apply(Prob1, 2, pnorm)
			Prob0 <- apply(Prob0, 2, pnorm)
			delta.1.tmp <- Prob1 - Prob0 
			rm(Prob1, Prob0)
	
			#Storage Matrices
			Prob1 <- matrix(,nrow=n, ncol=sims)
			Prob0 <- matrix(,nrow=n, ncol=sims)
			#T=0
			for (j in 1:sims) {
       			 Prob1[,j] <- TMmodel.1[j,1] + TMmodel.1[j,2]*0 + TMmodel.1[j,3]*PredictM1[,j]
        		 Prob0[,j] <- TMmodel.1[j,1] + TMmodel.1[j,2]*0 + TMmodel.1[j,3]*PredictM0[,j]
				}
			Prob1 <- Prob1 + X.pred
			Prob0 <- Prob0 + X.pred
			Prob1 <- apply(Prob1, 2, pnorm)
			Prob0 <- apply(Prob0, 2, pnorm)	
			delta.0.tmp <- Prob1 - Prob0 
			rm(Prob1, Prob0)
			
			#Storage Matrices
			Prob1_temp <- matrix(,nrow=n, ncol=sims)
			Prob0_temp <- matrix(,nrow=n, ncol=sims)

			#Calculate Total Effect
			for (j in 1:sims) {
        		Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j]  
        		Prob0_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM0[,j]  
        		}
			Prob1_temp <- Prob1_temp + X.pred
			Prob0_temp <- Prob0_temp + X.pred
			Prob1_temp <- apply(Prob1_temp, 2, pnorm)
			Prob0_temp <- apply(Prob0_temp, 2, pnorm)
			tau.tmp <- Prob1_temp - Prob0_temp
			} 
			delta.1 <- t(as.matrix(apply(delta.1.tmp, 2 , mean)))
			delta.0 <- t(as.matrix(apply(delta.0.tmp, 2 , mean)))
			tau <- t(as.matrix(apply(tau.tmp, 2, mean)))
			} else { #@@@@@@@@@@@@@@Nonparametric Bootstrap@@@@@@@@@@@@@@@@@@@
	n <- n.m
	k.t <- length(names(model.y.t$coef))
	k.m <- length(names(model.m$coef))
	Call.M <- model.m$call
	Call.Y.t <- model.y.t$call
	m.data <- model.frame(model.m)
    y.t.data <- model.frame(model.y.t)
	#Storage
	delta.1 <- matrix(NA, B, 1)
	delta.0 <- matrix(NA, B, 1)
	tau <- matrix(NA, B, 1)
	
	if(k.m == 2){ #No Covariates in Mediation Model
		for (b in 1:B) {
		#Resample Data
		index <- sample(1:n,n, repl=TRUE)
		Call.M$data <- m.data[index,]
		Call.Y.t$data <- y.t.data[index,]

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M)
		new.fit.t <- eval.parent(Call.Y.t)

		#Generate Mediation Model Predictions
		sigma <- summary(new.fit.M)$sigma
		mean.sim <- sum(new.fit.M$coef[1:k.m])
		error <- rnorm(n, mean.sim, sd=sigma) 
		PredictM1 <- sum(new.fit.M$coef[1:k.m]) + error
		PredictM0 <- new.fit.M$coef[1] + error 

		#Treatment Predictions
		pred.data.t <- data.frame(1,PredictM1)
		pred.data.c <- data.frame(1,PredictM0)
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0)

		pred.data.t <- data.frame(0,PredictM1)
		pred.data.c <- data.frame(0,PredictM0)
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		delta.0.tmp <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0)

		#Calculate Total Effect
		pred.data.t <- data.frame(1,PredictM1)
		pred.data.c <- data.frame(0,PredictM0)
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])
		
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		tau.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, PredictM1, PredictM0)
		
		delta.1[b] <- mean(delta.1.tmp)
		delta.0[b] <- mean(delta.0.tmp)
		tau[b] <- mean(tau.tmp)
		
		}
	} else {
	#With Covariates
	for (b in 1:B) {
		
		index <- sample(1:n,n, repl=TRUE)
		Call.M$data <- m.data[index,]
		Call.Y.t$data <- y.t.data[index,]

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M)
		new.fit.t <- eval.parent(Call.Y.t)

		#Generate Mediation Model Predictions
		sigma <- summary(new.fit.M)$sigma
		mean.sim <- sum(new.fit.M$coef[1:k.m])
		error <- rnorm(n, mean.sim, sd=sigma) 
		PredictM1 <- sum(new.fit.M$coef[1:k.m]) + error
		PredictM0 <- new.fit.M$coef[1] + sum(new.fit.M$coef[3:k.m]) + error 

		#Treatment Predictions
		pred.data.t <- data.frame(1,PredictM1, y.t.data[index,3:k.t])
		pred.data.c <- data.frame(1,PredictM0, y.t.data[index,3:k.t])

		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)

		pred.data.t <- data.frame(0,PredictM1, y.t.data[index,3:k.t])
		pred.data.c <- data.frame(0,PredictM0, y.t.data[index,3:k.t])
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.0.tmp <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)

		#Calculate Total Effect
		pred.data.t <- data.frame(1,PredictM1, y.t.data[index,3:k.t])
		pred.data.c <- data.frame(0,PredictM0, y.t.data[index,3:k.t])
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])
		
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		tau.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, PredictM1, PredictM0)

		delta.1[b] <- mean(delta.1.tmp)
		delta.0[b] <- mean(delta.0.tmp)
		tau[b] <- mean(tau.tmp)
		}
				
	  }

	}
		
 }
	
	d0 <- mean(delta.0)
	d1 <- mean(delta.1)
	d0.ci <- quantile(delta.0,c(.025,.975), na.rm=TRUE)
	d1.ci <- quantile(delta.1,c(.025,.975), na.rm=TRUE)
	avg.delta <- (delta.1 + delta.0)/2
	pct.dist <- avg.delta/tau
	pct.coef <- median(pct.dist)
	pct.ci <- quantile(pct.dist,c(.025,.975), na.rm=TRUE)
	tau.coef <- mean(tau)
	tau.ci <- quantile(tau,c(.025,.975), na.rm=TRUE)

	
}
out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci, pct.coef=pct.coef, pct.ci=pct.ci, tau.coef=tau.coef, tau.ci=tau.ci, boot=boot)
class(out) <- "my.mediate"
out

}

print.my.mediate <- function(x, ...){
	print(unlist(x[1:5]))
	invisible(x)
	}

summary.my.mediate <- function(object)
	structure(object, class = c("sum.my.mediate", class(object)))
 
print.sum.my.mediate <- function(x, ...){
	cat("\n Test For Mediation Effect \n\n")
	cat("Delta_0: ", format(x$d0, digits=4), "95% CI ", format(x$d0.ci, digits=4), "\n")
	cat("Delta_1: ", format(x$d1, digits=4), "95% CI ", format(x$d1.ci, digits=4), "\n\n")
	cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4), "\n")
	cat("95% Confidence Interval: ", format(x$pct.ci, digits=4), "\n \n")
	if(x$boot==TRUE){
		cat("Confidence Intervals Based on Nonparametric Bootstrap")
		} else {
		cat("Confidence Intervals Based on Parametric Bootstrap")
		}
	invisible(x)
	}
	


