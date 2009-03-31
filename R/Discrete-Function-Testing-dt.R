#Testing of Discrete Mediator
#Author: Luke Keele/Dustin
#Date: 032609

####################################################
# Discrete Mediator With Covariates                  #
####################################################
library(MASS)
library(Zelig)
rm(list=ls())
set.seed(3)

setwd("~/Documents/Mediation Analysis/Local/mediation/mediation/R")
#setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/R")
source("binary.med.R")

n <-  1000
#Population Values
alpha.2 <- .25
alpha.3 <- .25
beta.2 <- .25
beta.3 <- .25
gamma <- .25
kappa <- .25

#Treatment Indicator
T <- round(runif(n), 0)
X.1 <- rnorm(n)
X.2 <- rnorm(n)

latentm <- alpha.2 + beta.2*T + beta.2*X.1 + beta.2*X.2 + rnorm(n)
M <- latentm
M[latentm < -.3] <- 1
M[latentm >= -.3 & latentm < .45] <- 2
M[latentm >= .45 & latentm < 1.1] <- 3
M[latentm >= 1.1] <- 4
Y <- alpha.3 + beta.3*T + gamma*M + beta.2*X.1 + beta.2*X.2 + rnorm(n)

M <- as.factor(M)
data<- data.frame(M,T,X.1,X.2)

#calculated different in probability at mediator level 1-4 under treatment and under control
mmodel <- polr(M ~ T + X.1 + X.2, method="probit", Hess=TRUE, data=data)
treat <- data
treat[,2] <- 1
m.treat<-predict(mmodel, new=treat, type="probs")
control <- data
control[,2] <- 0
m.control <- predict(mmodel, new=control, type="probs")
m.def <- m.treat-m.control

#calculate predicted values of the outcome variable under the observed treatment and under the 4 categories of the mediator
ymodel <- lm(Y ~ T + M + X.1 + X.2)

y.predict.m1<-ymodel$coef["T"]*T+ymodel$coef["X.1"]*X.1+ymodel$coef["X.2"]*X.2
y.predict.m2<-ymodel$coef["T"]*T+ymodel$coef["M2"]+ymodel$coef["X.1"]*X.1+ymodel$coef["X.2"]*X.2
y.predict.m3<-ymodel$coef["T"]*T+ymodel$coef["M3"]+ymodel$coef["X.1"]*X.1+ymodel$coef["X.2"]*X.2
y.predict.m4<-ymodel$coef["T"]*T+ymodel$coef["M4"]+ymodel$coef["X.1"]*X.1+ymodel$coef["X.2"]*X.2

#Ypred*ProbDiff for each level of the mediator
cme<-matrix(, nrow=n, ncol=4)
for(i in 1:n){
cme[i,1]<-y.predict.m1[i]*m.def[i,1]
cme[i,2]<-y.predict.m2[i]*m.def[i,2]
cme[i,3]<-y.predict.m3[i]*m.def[i,3]
cme[i,4]<-y.predict.m4[i]*m.def[i,4]
}

#Average over the rows and then sum across the columns
acme<-apply(cme,2,mean, na.rm=TRUE)
delta <- sum(acme)
delta


#Simulation Approach
M <- as.factor(M)
mmodel <- polr(M ~ T + X.1 + X.2, method="probit", Hess=TRUE)
ymodel <- lm(Y ~ T + M + X.1 + X.2)

time.start <- Sys.time()
mod.1 <- mediate.binary(mmodel, ymodel, sims=1000, boot=TRUE, T="T", M="M")
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(mmodel, ymodel, sims=1000, T="T", M="M")
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)


#Relax No-Interaction Assumption
Y <- alpha.3 + beta.3*T + gamma*M + kappa*T*M + beta.2*X.1 + beta.2*X.2 + rnorm(n)

M <- as.factor(M)
mmodel <- polr(M ~ T + X.1 + X.2, method="probit", Hess=TRUE)
ymodel <- lm(Y ~ T + M + T:M + X.1 + X.2)


time.start <- Sys.time()
mod.1 <- mediate.binary(mmodel, ymodel, sims=1000, boot=TRUE, T="T", M="M", INT=TRUE)
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(mmodel, ymodel, sims=1000, T="T", M="M", INT=TRUE)
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)
