#Testing of Binary Mediation Function
#Author: Luke Keele
#Date: 2/13/2009

####################################################
# Binary Mediator With Covariates                  #
####################################################
library(MASS)
rm(list=ls())
set.seed(3)

setwd("~/Documents/Mediation Analysis/Local/mediation/mediation/R")
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
