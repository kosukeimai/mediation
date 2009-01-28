
library(MASS)

med.binary <- function(model.1, model.2, sims=1000){
	
	model.m <- model.1
	model.y <- model.2
	n.m <- model.frame(model.m)
	n.y <- model.frame(model.y)	
	k.y <- ncol(n.y)
	k.m <- ncol(n.m)
	k <- k.y + k.m
	n.m <- length(n.m[,1])
	n.y <- length(n.y[,1])
	n <- n.m
	sigma <- summary(model.m)$sigma
	
if(n.m != n.y){
	cat("Error: Missing Values Present in Data")
	} else {
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
			 #PredictM1 <- matrix(,nrow=n, ncol=sims)
 			 #PredictM0 <- matrix(,nrow=n, ncol=sims)
			#	for (j in 1:sims) {
			#	mean <- MModel[j,1] + MModel[j,2]
    		#	error <- mvrnorm(n, mu=mean, Sigma=sigma)  
        	#	#PredictM1[,j] <- MModel[j,1] + MModel[j,2] + error
        #		PredictM0[,j] <- MModel[j,1] + error  
 		#		}
			} else { #With Covariates
			#Generate Predictions From Mediation Model with Random  Error
			#mean.sim <- mean(apply(MModel, 1, sum))
			#error <- rnorm(sims, mean.sim, sd=sigma)  
			#PredictM1temp <- as.matrix(apply(cbind(MModel, error), 1, sum))
			#PredictM1  <- matrix(PredictM1temp, n, sims, byrow=TRUE)
			#PredictM1 <- PredictM1
			#PredictM0temp <- cbind(MModel[,1], MModel[,3:k.m], error)  
 			#PredictM0temp <- as.matrix(apply(PredictM0temp, 1, sum))
 			#PredictM0 <- matrix(PredictM0temp, n, sims, byrow=TRUE)
 			
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
			
				Prob1 <- apply(Prob1, 2, pnorm)
				Prob0 <- apply(Prob0, 2, pnorm)
				Delta1 <- Prob1 - Prob0 
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
			Delta0 <- Prob1 - Prob0 
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
			Tau <- Prob1_temp - Prob0_temp
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
    		#X.pred <- apply(X.pred, 1, sum)
    	
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
			Delta1 <- Prob1 - Prob0 
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
			Delta0 <- Prob1 - Prob0 
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
			Tau <- Prob1_temp - Prob0_temp
			}
		}
		
#Average Across Treatment and Control Predictions
Avg_Delta1 <- t(as.matrix(apply(Delta1, 2 , mean)))
Avg_Delta0 <- t(as.matrix(apply(Delta0, 2 , mean)))
Delta1_Delta0 <- Avg_Delta1 + Avg_Delta0
Delta1_Delta0 <- Delta1_Delta0/2
#Average over the i's 
avg.tau <- t(as.matrix(apply(Tau, 2, mean)))
#To calcuate proportion mediated we take for each draw j the ratio of 1/2(Delta1+Delta0)/Tau
med.eff.dist <- Delta1_Delta0 / avg.tau
med.eff.pct <- mean(med.eff.dist)
pr.med.ci <- quantile(med.eff.dist,c(.025,.975), na.rm=TRUE)


out <- list(pr.med = med.eff.pct, ci = pr.med.ci, dist=med.eff.dist)
}
		






