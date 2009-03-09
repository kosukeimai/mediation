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
	k.m <- ncol(m.data)
	
	if(is.factor(y.t.data[,T])==TRUE){

		cat.0 <- levels(y.t.data[,T])[1]
		cat.1 <- levels(y.t.data[,T])[2]
		} else {
		
		cat.0 <- 0
		cat.1 <- 1
		}
	
	#Storage
	delta.1 <- matrix(NA, B, 1)
	delta.0 <- matrix(NA, B, 1)
	tau <- matrix(NA, B, 1)
	
	for(b in 1:B){
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
		PredictM1 <- mean.sim + error
		
		
		if(is.factor(y.t.data[,T])==TRUE){
			PredictM0 <- sum(as.numeric(names(new.fit.M$coef) != paste(T,cat.1,sep="")) * new.fit.M$coef) + error
		} else {
			PredictM0 <- sum(as.numeric(names(new.fit.M$coef) != paste(T)) * new.fit.M$coef) + error
		}

		#Treatment Predictions Data
		if(INT==TRUE){
		#Treatment Predictions
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.1
		pred.data.t[,M] <- PredictM1
		pred.data.t$int <- PredictM1*1
		names(pred.data.t)[(k.t+1)] <- paste(T,M,sep=":") 
		
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.1
		pred.data.c[,M] <- PredictM0
		pred.data.c$int <- PredictM0*1
		names(pred.data.c)[(k.t+1)] <- paste(T,M,sep=":") 
			} else { 
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.1
		pred.data.t[,M] <- PredictM1
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.1
		pred.data.c[,M] <- PredictM0
		}
			
		if(is.factor(y.t.data[,T])==TRUE){
		pred.data.t[,T] <- as.factor(pred.data.t[,T])
		pred.data.c[,T] <- as.factor(pred.data.t[,T])
		} else { 
		pred.data.t[,T] <- as.numeric(pred.data.t[,T])
		pred.data.c[,T] <- as.numeric(pred.data.t[,T])
		}
		
		#Treatment Predications	
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)

		#Control Predictions Data
		if(INT==TRUE){
		#Treatment Predictions
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.0
		pred.data.t[,M] <- PredictM1
		pred.data.t$int <- PredictM1*0
		names(pred.data.t)[(k.t+1)] <- paste(T,M,sep=":") 
		
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.0
		pred.data.c[,M] <- PredictM0
		pred.data.c$int <- PredictM0*0
		names(pred.data.c)[(k.t+1)] <- paste(T,M,sep=":") 
			} else { 
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.0
		pred.data.t[,M] <- PredictM1
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.0
		pred.data.c[,M] <- PredictM0
		}
			
		if(is.factor(y.t.data[,T])==TRUE){
		pred.data.t[,T] <- as.factor(pred.data.t[,T])
		pred.data.c[,T] <- as.factor(pred.data.t[,T])
		} else { 
		pred.data.t[,T] <- as.numeric(pred.data.t[,T])
		pred.data.c[,T] <- as.numeric(pred.data.t[,T])
		} 
		
		#Control Predictions	
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.0.tmp <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
		
		#Calculate Total Effect
if(INT==TRUE){
		#Treatment Predictions
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.1
		pred.data.t[,M] <- PredictM1
		pred.data.t$int <- PredictM1*1
		names(pred.data.t)[(k.t+1)] <- paste(T,M,sep=":") 
		
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.0
		pred.data.c[,M] <- PredictM0
		pred.data.c$int <- PredictM0*0
		names(pred.data.c)[(k.t+1)] <- paste(T,M,sep=":") 
			} else { 
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.1
		pred.data.t[,M] <- PredictM1
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.0
		pred.data.c[,M] <- PredictM0
		}
			
		if(is.factor(y.t.data[,T])==TRUE){
		pred.data.t[,T] <- as.factor(pred.data.t[,T])
		pred.data.c[,T] <- as.factor(pred.data.t[,T])
		} else { 
		pred.data.t[,T] <- as.numeric(pred.data.t[,T])
		pred.data.c[,T] <- as.numeric(pred.data.t[,T])
		}
				
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		tau.tmp <- pr.mat[,1] - pr.mat[,2]
		
				rm(pred.data.t, pred.data.c, pr.1, pr.0, PredictM1, PredictM0)
		
		delta.1[b] <- mean(delta.1.tmp)
		delta.0[b] <- mean(delta.0.tmp)
		tau[b] <- mean(tau.tmp)
		
		} #bootstrap loop
	} #nonpara boot branch
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

out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci, pct.coef=pct.coef, pct.ci=pct.ci, tau.coef=tau.coef, tau.ci=tau.ci, boot=boot, INT=INT)
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
	if(x$INT==TRUE){
		cat("\n Test For Mediation Effect \n\n")
	if(x$boot==TRUE){
		cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
		} else {
		cat("Confidence Intervals Based on Parametric Bootstrap\n\n")
		}
	cat("Delta_0: ", format(x$d0, digits=4), "95% CI ", format(x$d0.ci, digits=4), "\n")
	cat("Delta_1: ", format(x$d1, digits=4), "95% CI ", format(x$d1.ci, digits=4), "\n\n")
	cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4), "\n")
	cat("95% Confidence Interval: ", format(x$pct.ci, digits=4), "\n \n")
		} else {
			cat("\n Test For Mediation Effect \n\n")
			if(x$boot==TRUE){
		cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
		} else {
		cat("Confidence Intervals Based on Parametric Bootstrap\n\n")
		}
	cat("Delta(t): ", format(x$d0, digits=4), "95% CI ", format(x$d0.ci, digits=4), "\n")
	cat("Proportion of Total Effect via Mediation: ", format(x$pct.coef, digits=4), "\n")
	cat("95% Confidence Interval: ", format(x$pct.ci, digits=4), "\n \n")
				}
	invisible(x)
	}
	


