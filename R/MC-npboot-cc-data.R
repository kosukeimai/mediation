

#Nonparametric Bootstrap Mediation Function
med.boot <- function(object.1, object.2, object.3, B=500){

	model.m <- object.1
	model.y.t <- object.2
	model.y.c <- object.3
	n.m <- model.frame(model.m)
	n.y <- model.frame(model.y.t)	
	n.m <- length(n.m[,1])
	n.y <- length(n.y[,1])
	
if(n.m != n.y){
	cat("Error: Missing Values Present in Data")
	} else {
		n <- n.m
	k.t <- length(names(model.y.t$coef))
	k.c <- length(names(model.y.c$coef))
	k.m <- length(names(model.m$coef))
	Call.M <- model.m$call
	Call.Y.t <- model.y.t$call
	Call.Y.c <- model.y.c$call
	m.data <- model.frame(model.m)
    y.t.data <- model.frame(model.y.t)
    y.c.data <- model.frame(model.y.c)
	#Storage
	avg.delta.1 <- matrix(NA, B, 1)
	avg.delta.0 <- matrix(NA, B, 1)
	med.eff <- matrix(NA, B, 1)
	avg.tau <- matrix(NA, B, 1)
	pct.med <- matrix(NA, B, 1)
	
	if(k.m == 2){ #No Covariates in Mediation Model
		for (b in 1:B) {
		#Resample Data
		index <- sample(1:n,n, repl=TRUE)
		Call.M$data <- m.data[index,]
		Call.Y.t$data <- y.t.data[index,]
		Call.Y.c$data <- y.c.data[index,]

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M)
		new.fit.t <- eval.parent(Call.Y.t)
		new.fit.c <- eval.parent(Call.Y.c)

		#Generate Mediation Model Predictions
		sigma <- summary(new.fit.M)$sigma
		mean.sim <- sum(new.fit.M$coef[1:k.m])
		error <- rnorm(n, mean.sim, sd=sigma) 
		PredictM1 <- sum(new.fit.M$coef[1:k.m]) + error
		#Without Covariattes
		PredictM0 <- new.fit.M$coef[1] + error 

		#Treatment Predictions
		pred.data.t <- data.frame(1,PredictM1)
		pred.data.c <- data.frame(1,PredictM0)

		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1 <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0)

		pred.data.t <- data.frame(PredictM1)
		pred.data.c <- data.frame(PredictM0)
		names(pred.data.t) <- names(model.y.c$coef[2:k.c])
		names(pred.data.c) <- names(model.y.c$coef[2:k.c])

		pr.1 <- predict(new.fit.c, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.c, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		delta.0 <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0)

		#Calculate Total Effect
		pred.data.t <- data.frame(1,PredictM1)
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)

		pred.data.c <- data.frame(PredictM0)
		names(pred.data.c) <- names(model.y.c$coef[2:k.c])
		pr.0 <- predict(new.fit.c, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
		tau <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, PredictM1, PredictM0)

		avg.delta.1[b] <- mean(delta.1)
		avg.delta.0[b] <- mean(delta.0)
		med.eff[b] <- (avg.delta.1[b] + avg.delta.0[b])/2
		avg.tau[b] <- mean(tau)
		pct.med[b] <- med.eff[b]/avg.tau[b]
		}
	} else {
	#With Covariates
	for (b in 1:B) {
		
		index <- sample(1:n,n, repl=TRUE)
		Call.M$data <- m.data[index,]
		Call.Y.t$data <- y.t.data[index,]
		Call.Y.c$data <- y.c.data[index,]

		#Refit Models with Resampled Data
		new.fit.M <- eval.parent(Call.M)
		new.fit.t <- eval.parent(Call.Y.t)
		new.fit.c <- eval.parent(Call.Y.c)

		#Generate Mediation Model Predictions
		sigma <- summary(new.fit.M)$sigma
		mean.sim <- sum(new.fit.M$coef[1:k.m])
		error <- rnorm(n, mean.sim, sd=sigma) 
		PredictM1 <- sum(new.fit.M$coef[1:k.m]) + error
		#Without Covariattes
		PredictM0 <- new.fit.M$coef[1] + sum(new.fit.M$coef[3:k.m]) + error 

		#Treatment Predictions
		pred.data.t <- data.frame(1,PredictM1, y.t.data[index,3:k.t])
		pred.data.c <- data.frame(1,PredictM0, y.t.data[index,3:k.t])

		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		names(pred.data.c) <- names(model.y.t$coef[2:k.t])

		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.1 <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0)

		pred.data.t <- data.frame(PredictM1, y.c.data[index,3:k.c])
		pred.data.c <- data.frame(PredictM0, y.c.data[index,3:k.c])
		names(pred.data.t) <- names(model.y.c$coef[2:k.c])
		names(pred.data.c) <- names(model.y.c$coef[2:k.c])

		pr.1 <- predict(new.fit.c, type="response", newdata=pred.data.t)
		pr.0 <- predict(new.fit.c, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		delta.0 <-pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0)

		#Calculate Total Effect
		pred.data.t <- data.frame(1,PredictM1, y.t.data[index,3:k.t])
		names(pred.data.t) <- names(model.y.t$coef[2:k.t])
		pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)

		pred.data.c <- data.frame(PredictM0, y.c.data[index,3:k.c])
		names(pred.data.c) <- names(model.y.c$coef[2:k.c])
		pr.0 <- predict(new.fit.c, type="response", newdata=pred.data.c)
		pr.mat <- as.matrix(cbind(pr.1, pr.0))
		tau <- pr.mat[,1] - pr.mat[,2]

		rm(pred.data.t, pred.data.c, pr.1, pr.0, PredictM1, PredictM0)

		avg.delta.1[b] <- mean(delta.1)
		avg.delta.0[b] <- mean(delta.0)
		med.eff[b] <- (avg.delta.1[b] + avg.delta.0[b])/2
		avg.tau[b] <- mean(tau)
		pct.med[b] <- med.eff[b]/avg.tau[b]

			}
		}
	}
			
	#Calculate Quantities of Interest
	ci.d1 <- quantile(avg.delta.1, c(.025,.975))
	ci.d0 <- quantile(avg.delta.0, c(.025,.975))
	Delta1 <- mean(avg.delta.1)
	Delta0 <- mean(avg.delta.0)
	med.eff.dist <- med.eff / avg.tau
	med.eff.pct <- abs(mean(t(med.eff.dist)))*100

	out <- list(med.t = Delta1, med.c = Delta0, conf.1 = ci.d1, conf.0 = ci.d0, pct = med.eff.pct, dist=avg.delta.1)

	
}