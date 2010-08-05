library(MASS)
library(foreign)

rm(list=ls())

set.seed(31453)

setwd("~/Documents/Mediation/example/Anxiety")

im_emo <- read.dta("brader_trim_nomiss.dta")
attach(im_emo)

#Mediation Model
model.m <- lm(emo ~  tone_eth + ppage + ppeducat + ppgender + ppincimp, data=im_emo)

table(english)
table(immigr)

#Outcome Model - English Only Laws
eng <- 0

if(eng==1){
model.y <- polr(as.factor(english) ~ emo + tone_eth + ppage, data=im_emo, method="probit", Hess=TRUE)
	} else {
model.y <- polr(as.factor(immigr) ~ emo + tone_eth + ppage, data=im_emo, method="probit", Hess=TRUE)
		}
		
#Adapting the basic algorithm to the ordered outcome case is in some ways very straightforward.
#There is one problem that I note below. The biggest complication is that one now gets an ACME for each
#category of Y.  This means the summary function needs to be rewritten quite a bit to print out a dynamic
#set of ACMEs that will be a function of the number of categories in Y.

mediator <- "emo"
treat <- "tone_eth"
sims <- 1000
B <- sims
model.y.t <- model.y
	
	m.data <- model.frame(model.m)  #Call.M$data
    y.t.data <- model.frame(model.y.t) #Call.Y$data
	n <- length(y.t.data[,1])
	m <- length(sort(unique(model.frame(model.m)[,1])))
	m.min <- as.numeric(sort(unique(model.frame(model.m)[,1]))[1])

	cat.0 <- 0
	cat.1 <- 1

	#Storage
	#if(eng==1){
	delta.1 <- matrix(NA, B, 4)
	delta.0 <- matrix(NA, B, 4)
	zeta.1 <- matrix(NA, B, 4)
	zeta.0 <- matrix(NA, B, 4)
	tau <- matrix(NA, B, 4)
	#} else {
	#delta.1 <- matrix(NA, B, 5)
	#delta.0 <- matrix(NA, B, 5)
	#zeta.1 <- matrix(NA, B, 5)
	#zeta.0 <- matrix(NA, B, 5)
	#tau <- matrix(NA, B, 5)
	#	}
	
	for(b in 1:B){
		
		#The other problem with the ordered outcome is that the Call function used before will not
		#work with the polr() function. The problem is that this means one needs to make sure this
		#code is general enough to be the same specification as the original model.  
		#This is probably best done via using the formula() and paste() commands.
		data.star <- im_emo[sample(1:nrow(im_emo),n,replace=TRUE),]
		new.fit.M <- lm(emo ~  tone_eth + ppage + ppeducat + ppgender + ppincimp, data=data.star)
		if(eng==1){
		new.fit.t <- polr(as.factor(english) ~ emo + tone_eth + ppage, data=data.star, method="probit", Hess=TRUE)
		} else {
		new.fit.t <-  polr(as.factor(immigr) ~ emo + tone_eth + ppage, data=data.star, method="probit", Hess=TRUE)
			}

#The rest of the algorithm is the same as the other finished code.
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
		
		sigma <- summary(new.fit.M)$sigma
		error <- rnorm(n, mean=0, sd=sigma)
		PredictM1 <- predict(new.fit.M, type="response", newdata=pred.data.t) + error
		PredictM0 <- predict(new.fit.M, type="response", newdata=pred.data.c) + error
		rm(error)
		
	
#####################################################################################		
		#Treatment Predictions Data
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data

		pred.data.t[,treat] <- cat.1	
		pred.data.c[,treat] <- cat.1

		pred.data.t[,mediator] <- PredictM1
		pred.data.c[,mediator] <- PredictM0
			
		#Treatment Predictions
		probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
		probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
						
		delta.1.tmp <-  probs_p1 - probs_p0		
		rm(pred.data.t, pred.data.c, probs_p1, probs_p0)

		#Control Predictions Data
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data

		pred.data.t[,treat] <- cat.0	
		pred.data.c[,treat] <- cat.0

		pred.data.t[,mediator] <- PredictM1
		pred.data.c[,mediator] <- PredictM0
			
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
		pred.data.t[,treat] <- cat.1	
		pred.data.c[,treat] <- cat.0

		pred.data.t[,mediator] <- PredictM1
		pred.data.c[,mediator] <- PredictM1
			
		probs_p1 <- predict(new.fit.t, newdata=pred.data.t, type="probs")
		probs_p0 <- predict(new.fit.t, newdata=pred.data.c, type="probs")
		zeta.1.tmp <- probs_p1 - probs_p0

		rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
		
		#Zeta.0
		pred.data.t <- y.t.data
		pred.data.c <- y.t.data
		pred.data.t[,treat] <- cat.1	
		pred.data.c[,treat] <- cat.0

		pred.data.t[,mediator] <- PredictM0
		pred.data.c[,mediator] <- PredictM0

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
	apply(delta.0, 2, mean)
	apply(delta.0, 2, quantile, c(.025,.975))
	
	apply(zeta.1, 2, mean)
	apply(zeta.1, 2, quantile, c(.025,.975))
	
	apply(tau, 2, mean)
	apply(tau, 2, quantile, c(.025,.975))
	
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
	
	