mediate.cont <- function(z, ...){
	UseMethod("mediate.cont", z)
	}
	
mediate.cont.default <- function(z, model.y, sims=1000, boot=FALSE, INT=FALSE, T="treat.name", M="med.name"){
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
	if(n.m != n.y){
	cat("Error: Missing Values Present in Data")
	} else {

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
				if(INT == TRUE){
					Prob1 <- matrix(,nrow=n, ncol=sims)
				Prob0 <- matrix(,nrow=n, ncol=sims)
				#Prediction For Outcome Based on Mediation Model Preductions
				#T=1
				for (j in 1:sims) {
       				 Prob1[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j] + TMmodel[j,4]*PredictM1[,j]*1
       			 	 Prob0[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM0[,j] + TMmodel[j,4]*PredictM0[,j]*1
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
       			 Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j] + TMmodel[j,4]*PredictM1[,j]*1  
        		 Prob0_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM0[,j] + TMmodel[j,4]*PredictM0[,j]*0 
        		}
			tau.tmp <- Prob1_temp - Prob0_temp
					} else {
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
			}
				} else {#Predictions For Covariates
					if(INT==TRUE){
			#Predictions For Covariates Without No Interaction Assumption
			TMmodel.1 <- as.matrix(TMmodel[,1:4])
			TMmodel.2 <- as.matrix(TMmodel[,5:(k.y+1)])
			ii <- 4
   			X <- model.frame(model.y)
    		X <- as.matrix(X[,ii:k.y])
    		X.pred <- matrix( , n, sims)
    		k.x <- ncol(X)
   			 for (i in 1:k.x){
    			for (j in 1:sims) {
        			X.pred[,j] <- TMmodel.2[j,i]*as.numeric(X[,i])
    				}
    			}
    			#Storage Matrices For Outcome Predictions
			Prob1 <- matrix(,nrow=n, ncol=sims)
			Prob0 <- matrix(,nrow=n, ncol=sims)
			#Prediction For Outcome Based on Mediation Model Preductions
			#T=1
			for (j in 1:sims) {
       				 Prob1[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j] + TMmodel[j,4]*PredictM1[,j]*1
       			 	 Prob0[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM0[,j] + TMmodel[j,4]*PredictM0[,j]*1
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
       			 Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1 + TMmodel[j,3]*PredictM1[,j] + TMmodel[j,4]*PredictM1[,j]*1  
        		 Prob0_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*0 + TMmodel[j,3]*PredictM0[,j] + TMmodel[j,4]*PredictM0[,j]*0 
        		}			
        	Prob1_temp <- Prob1_temp + X.pred
			Prob0_temp <- Prob0_temp + X.pred
			tau.tmp <- Prob1_temp - Prob0_temp
			} else {#End Int Branch
			TMmodel.1 <- as.matrix(TMmodel[,1:3])
			TMmodel.2 <- as.matrix(TMmodel[,4:k.y])
			ii <- 4
   			X <- model.frame(model.y)
    		X <- as.matrix(X[,ii:k.y])
    		X.pred <- matrix( , n, sims)
    		k.x <- ncol(X)
   			 for (i in 1:k.x){
    			for (j in 1:sims) {
        			X.pred[,j] <- TMmodel.2[j,i]*as.numeric(X[,i])
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
			} 
			delta.1 <- t(as.matrix(apply(delta.1.tmp, 2 , mean)))
			delta.0 <- t(as.matrix(apply(delta.0.tmp, 2 , mean)))
			tau <- t(as.matrix(apply(tau.tmp, 2, mean)))
			} else { #@@@@@@@@@@@@@@Nonparametric Bootstrap@@@@@@@@@@@@@@@@@@@
	n <- n.m
	Call.M <- model.m$call
	Call.Y.t <- model.y.t$call
	m.data <- model.frame(model.m)  #Call.M$data
    y.t.data <- model.frame(model.y.t) #Call.Y$data
    k.t <- ncol(y.t.data)
	k.m <- ncol(m.data)
	
	if(is.factor(y.t.data[,T])==TRUE){

		cat.0 <- levels(y.t.data[,T])[1]
		cat.1 <- levels(y.t.data[,T])[2]
		} else {
		
		cat.0 <- 0
		cat.1 <- 1
		}
	T.cat <- paste(T,cat.1, sep="")
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

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M)
		new.fit.t <- eval.parent(Call.Y.t)

		#Direct Effect
		if (INT==TRUE){
			zeta.0[b] <- new.fit.t$coef[paste(T,cat.1,sep="")] + new.fit.t$coef[paste(T.cat,M,sep=":")]*(new.fit.M$coef[1] + 				new.fit.M$coef[paste(T,cat.1,sep="")]*0)
			zeta.1[b] <- new.fit.t$coef[paste(T,cat.1,sep="")] + new.fit.t$coef[paste(T.cat,M,sep=":")]*(new.fit.M$coef[1] + 			new.fit.M$coef[paste(T,cat.1,sep="")]*1)			} else {
			zeta.1[b] <- new.fit.t$coef[paste(T)]
			zeta.0[b] <- new.fit.t$coef[paste(T)] 
				}
				
		#Generate Mediation Model Predictions
		sigma <- summary(new.fit.M)$sigma
		error <- rnorm(n, mean=0, sd=sigma) 
		pred.data.t <- m.data
		pred.data.t[,T] <- cat.1
		pred.data.c <- m.data
		pred.data.c[,T] <- cat.0
		
		if(is.factor(m.data[,T])==TRUE){
		pred.data.t[,T] <- as.factor(pred.data.t[,T])
		pred.data.c[,T] <- as.factor(pred.data.c[,T])
		} else { 
		pred.data.t[,T] <- as.numeric(pred.data.t[,T])
		pred.data.c[,T] <- as.numeric(pred.data.c[,T])
		} 
		
		PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t) + error
		PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c) + error
		
		#Treatment Predictions Data
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.1
		pred.data.t[,M] <- PredictM1
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.1
		pred.data.c[,M] <- PredictM0
			
		if(is.factor(y.t.data[,T])==TRUE){
		pred.data.t[,T] <- as.factor(pred.data.t[,T])
		pred.data.c[,T] <- as.factor(pred.data.c[,T])
		} else { 
		pred.data.t[,T] <- as.numeric(pred.data.t[,T])
		pred.data.c[,T] <- as.numeric(pred.data.c[,T])
		}
		
		#Treatment Predications	
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)

		#Control Predictions Data
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.0
		pred.data.t[,M] <- PredictM1
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.0
		pred.data.c[,M] <- PredictM0
			
		if(is.factor(y.t.data[,T])==TRUE){
		pred.data.t[,T] <- as.factor(pred.data.t[,T])
		pred.data.c[,T] <- as.factor(pred.data.c[,T])
		} else { 
		pred.data.t[,T] <- as.numeric(pred.data.t[,T])
		pred.data.c[,T] <- as.numeric(pred.data.c[,T])
		} 
		
		#Control Predictions	
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.0.tmp <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
		
				
		delta.1[b] <- mean(delta.1.tmp)
		delta.0[b] <- mean(delta.0.tmp)
		if(INT==TRUE){
		tau[b] <- zeta.1[b] + delta.0[b]
		} else {
		tau[b] <- zeta.1[b] + mean(delta.0.tmp)	
			}
		
		} #bootstrap loop
	} #nonpara boot branch
	d0 <- mean(delta.0)
	d1 <- mean(delta.1)
	d0.ci <- quantile(delta.0,c(.025,.975), na.rm=TRUE)
	d1.ci <- quantile(delta.1,c(.025,.975), na.rm=TRUE)
	tau.coef <- mean(tau)
	tau.ci <- quantile(tau,c(.025,.975), na.rm=TRUE)
	#avg.delta <- (d0 + d1)/2
	#pct.coef <- abs(avg.delta/tau.coef)
	#avg.delta <- (delta.0)/tau
	#pct.coef <- abs(avg.delta/tau.coef)
	pct.dist <- delta.0/tau
	pct.coef <- median(pct.dist)
	pct.ci <- quantile(pct.dist,c(.025,.975), na.rm=TRUE)
	z1 <- mean(zeta.1)	
	z1.ci <- quantile(zeta.1,c(.025,.975), na.rm=TRUE)
	z0 <- mean(zeta.0)	
	z0.ci <- quantile(zeta.0,c(.025,.975), na.rm=TRUE)
		
 }

out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci, pct.coef=pct.coef, pct.ci=pct.ci, tau.coef=tau.coef, tau.ci=tau.ci, z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,boot=boot, INT=INT)
class(out) <- "my.mediate"
out

}

print.my.mediate <- function(x, ...){
	print(unlist(x[1:11]))
	invisible(x)
	}

summary.my.mediate <- function(object)
	structure(object, class = c("sum.my.mediate", class(object)))
 
print.sum.my.mediate <- function(x, ...){
	if(x$INT==TRUE){
		cat("\n Test For Mediation Effect \n\n")
	if(x$boot==TRUE){
		cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
		} else {
		cat("Confidence Intervals Based on Parametric Bootstrap\n\n")
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
		cat("Confidence Intervals Based on Parametric Bootstrap\n\n")
		}
	cat("Mediation Effect: ", format(x$d1, digits=4), "95% CI ", format(x$d1.ci, digits=4), "\n")
	cat("Direct Effect: ", format(x$z0, digits=4), "95% CI ", format(x$z0.ci, digits=4), "\n")
	cat("Total Effect: ", format(x$tau.coef, digits=4), "95% CI ", format(x$tau.ci, digits=4), "\n")
    cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"95% CI ", format(x$pct.ci, digits=4),"\n")
    #cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4),"\n")
				}
	invisible(x)
	}
	


