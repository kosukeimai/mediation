mediate.cont <- function(z, ...){
	UseMethod("mediate.cont", z)
	}
	
mediate.cont.default <- function(z, model.y, sims=1000, boot=FALSE, INT=FALSE, T="treat.name", M="med.name"){
	B <- sims
	model.m <- z
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
	sigma <- summary(model.m)$sigma
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
	cat("Error: Missing Values Present in Data")
	} else {
	if(boot == FALSE){ 
	cat.0 <- 0
	cat.1 <- 1

	MModel.coef <- model.m$coef
	MModel.var.cov <- vcov(model.m)
	TMmodel.coef <- model.y$coef
	TMmodel.var.cov <- vcov(model.y)
	MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov)
	TMmodel <- mvrnorm(sims, mu=TMmodel.coef, Sigma=TMmodel.var.cov)
	
	#Mediator Predictions
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
    
	error <- rnorm(n, mean=0, sd=sigma)
	PredictM1 <- MModel %*% t(mmat.t)
	PredictM0 <- MModel %*% t(mmat.c)
	PredictM1 <- PredictM1 + error
	PredictM0 <- PredictM0 + error
	
################################################
# Outcome Predictions
###############################################	
	#Treatment Predictions Data
	Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
	
	for(j in 1:sims){
	if(is.factor(y.t.data[,paste(T)])==TRUE){
	pred.data.t <- y.t.data
	pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	pred.data.t[,M] <- PredictM1[j,]
	pred.data.c <- y.t.data
	pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[2], levels = levels(y.t.data[,T])))
	pred.data.c[,M] <- PredictM0[j,]
			} else {
	pred.data.t <- y.t.data
	pred.data.t[,T] <- cat.1
	pred.data.t[,M] <- PredictM1[j,]
	pred.data.c <- y.t.data
	pred.data.c[,T] <- cat.1
	pred.data.c[,M] <- PredictM0[j,]

	}
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
	
	#Treatment Predictions
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
	delta.1.tmp <- Pr1 - Pr0
	
	#Control Predictions Data
	Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
	
	for(j in 1:sims){
	if(is.factor(y.t.data[,paste(T)])==TRUE){
	pred.data.t <- y.t.data
	pred.data.t[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	pred.data.t[,M] <- PredictM1[j,]
	pred.data.c <- y.t.data
	pred.data.c[,T] <- list(factor(unique(y.t.data[,T])[1], levels = levels(y.t.data[,T])))
	pred.data.c[,M] <- PredictM0[j,]
			} else {
	pred.data.t <- y.t.data
	pred.data.t[,T] <- cat.0
	pred.data.t[,M] <- PredictM1[j,]
	pred.data.c <- y.t.data
	pred.data.c[,T] <- cat.0
	pred.data.c[,M] <- PredictM0[j,]
	}
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)

	#Control Predictions
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c)
	}
	
	delta.0.tmp <- Pr1 - Pr0
	
	#Direct Effects
	Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
	
	for(j in 1:sims){
	pred.data.t <- y.t.data
	pred.data.t[,T] <- cat.1
	pred.data.t[,M] <- PredictM1[j,]
	pred.data.c <- y.t.data
	pred.data.c[,T] <- cat.0
	pred.data.c[,M] <- PredictM1[j,]
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
	
	#Direct Predictions
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
	zeta.1.tmp <- Pr1 - Pr0
	
	#Zeta-0
	Pr1 <- matrix(,nrow=n, ncol=sims)
	Pr0 <- matrix(,nrow=n, ncol=sims)
	
	for(j in 1:sims){
	pred.data.t <- y.t.data
	pred.data.t[,T] <- cat.1
	pred.data.t[,M] <- PredictM0[j,]
	pred.data.c <- y.t.data
	pred.data.c[,T] <- cat.0
	pred.data.c[,M] <- PredictM0[j,]
	ymat.t <- model.matrix(terms(model.y), data=pred.data.t) 
	ymat.c <- model.matrix(terms(model.y), data=pred.data.c)
	
	#Direct Predictions
	Pr1[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.t)
	Pr0[,j] <- t(as.matrix(TMmodel[j,])) %*% t(ymat.c)
	
	rm(ymat.t, ymat.c, pred.data.t,pred.data.c)
	}	
	
	zeta.0.tmp <- Pr1 - Pr0

	
	delta.1 <- t(as.matrix(apply(delta.1.tmp, 2, mean)))
	delta.0 <- t(as.matrix(apply(delta.0.tmp, 2, mean)));
	zeta.1 <- t(as.matrix(apply(zeta.1.tmp, 2, mean)))
	zeta.0 <- t(as.matrix(apply(zeta.0.tmp, 2, mean)))
	
	if(INT==TRUE){
		tau <- zeta.0 + delta.1
		} else {
		tau <- zeta.1 + delta.0	
		}
		
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
		temp.coef <- new.fit.M$coef	
		#X-Predictions
		temp <- model.matrix(new.fit.M)
		temp[,1] <- 0
	   	meanmat <- apply(temp, 2, mean)
	    #new.fit.M$coef[1] <- 0
	    if(is.factor(y.t.data[,T])==TRUE){
        temp.coef[paste(T.cat)] <- 0
        } else {
        temp.coef[paste(T.cat)] <- 0
        	}
        Bm.X <- new.fit.M$coef %*% meanmat
		zeta.0[b] <- new.fit.t$coef[paste(T.cat)] + new.fit.t$coef[paste(T.cat,M,sep=":")]*(new.fit.M$coef[1] + new.fit.M$coef[paste(T.cat)]*0 + Bm.X)
		zeta.1[b] <- new.fit.t$coef[paste(T.cat)] + new.fit.t$coef[paste(T.cat,M,sep=":")]*(new.fit.M$coef[1] + new.fit.M$coef[paste(T.cat)]*1 + Bm.X)			} else {
		if(is.factor(y.t.data[,T])==TRUE){
			zeta.1[b] <- new.fit.t$coef[paste(T.cat)]
			zeta.0[b] <- new.fit.t$coef[paste(T.cat)]
				} else {
			zeta.1[b] <- new.fit.t$coef[paste(T)]
			zeta.0[b] <- new.fit.t$coef[paste(T)] 
			}
				}
				
		#Generate Mediation Model Predictions
		if(class(model.m)[1]=="gam"){
			sigma <- summary(new.fit.M)$scale
			} else {
			sigma <- summary(new.fit.M)$sigma
				}
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
		
		#Direct Effects
		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.1
		pred.data.t[,M] <- PredictM1
		pred.data.c <- y.t.data
		pred.data.c[,T] <- cat.0
		pred.data.c[,M] <- PredictM1
			
		if(is.factor(y.t.data[,T])==TRUE){
		pred.data.t[,T] <- as.factor(pred.data.t[,T])
		pred.data.c[,T] <- as.factor(pred.data.c[,T])
		} else { 
		pred.data.t[,T] <- as.numeric(pred.data.t[,T])
		pred.data.c[,T] <- as.numeric(pred.data.c[,T])
		}
		
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		zeta.1.tmp <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0,pr.mat)

		pred.data.t <- y.t.data
		pred.data.t[,T] <- cat.1
		pred.data.t[,M] <- PredictM0
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
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		zeta.0.tmp <-pr.mat[,1] - pr.mat[,2]
		
		zeta.1[b] <- mean(zeta.1.tmp)
		zeta.0[b] <- mean(zeta.0.tmp)		
		delta.1[b] <- mean(delta.1.tmp)
		delta.0[b] <- mean(delta.0.tmp)
		if(INT==TRUE){
		tau[b] <- zeta.1[b] + delta.0[b]
		} else {
		tau[b] <- zeta.0[b] + delta.1[b]	
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
	pct.dist <- abs(delta.0/tau)
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
	


