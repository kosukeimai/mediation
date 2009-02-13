#Testing of Binary Mediation Function
#Author: Luke Keele
#Date: 2/13/2009

#####################################
#          Binary Mediator          #
#####################################

library(MASS)
rm(list=ls())
set.seed(373456)

setwd("~/Documents/Mediation Analysis/Local/mediation/mediation/R")
source("binary.med.R")

n <-  1000
#Population Values
alpha.2 <- .25
alpha.3 <- .25
beta.2 <- .25
beta.3 <- .25
gamma <- .25

#Treatment Indicator
T <- round(runif(n), 0)
latentm <- alpha.2 + beta.2*T + rnorm(n)
M <- latentm
M[latentm <  1] <- 0
M[latentm >= 1] <- 1
Y <- alpha.3 + beta.3*T + gamma*M + rnorm(n)

#Standard Method - Without Covariates
mmodel <- glm(M ~ T, family=binomial(link=probit))
ymodel <- lm(Y ~ T + M)

time.start <- Sys.time()
mod.1 <- mediate.binary(mmodel, ymodel, sims=1000, boot=TRUE)
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(mmodel, ymodel, sims=1000)
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)

####################################################
#                 With Covariates                  #
####################################################

rm(list=ls())
set.seed(3)
source("binary.med.R")
n <-  1000
#Population Values
alpha.2 <- .25
alpha.3 <- .25
beta.2 <- .25
beta.3 <- .25
gamma <- .25

#Treatment Indicator
T <- round(runif(n), 0)
X.1 <- rnorm(n)
X.2 <- rnorm(n)

latentm <- alpha.2 + beta.2*T + beta.2*X.1 + beta.2*X.2 + rnorm(n)
M <- latentm
M[latentm <  1] <- 0
M[latentm >= 1] <- 1
Y <- alpha.3 + beta.3*T + gamma*M + beta.2*X.1 + beta.2*X.2 + rnorm(n)

mmodel <- glm(M ~ T + X.1 + X.2, family=binomial(link=probit))
ymodel <- lm(Y ~ T + M + X.1 + X.2)

#New Methods
time.start <- Sys.time()
mod.1 <- mediate.binary(mmodel, ymodel, sims=1000, boot=TRUE)
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(mmodel, ymodel, sims=1000)
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)

#############################################
#             Binary Outcome                #
#############################################
rm(list=ls())
set.seed(373456)
source("binary.med.R")

#Simulate X Vars
n <-  1000
T <- round(runif(n), 0)
#Population Values
alpha.2 <- .25
alpha.3 <- .25
beta.2 <- .25
beta.3 <- .25
gamma <- .25

M <- alpha.2 + beta.2*T + rnorm(n)
latenty <-  alpha.3 + beta.3*T + gamma*M + rnorm(n)
Y <- latenty
Y[latenty <  1] <- 0
Y[latenty >= 1] <- 1

#Standard Method - Without Covariates
ymodel <- glm(Y ~ T + M, family=binomial(link=probit))
mmodel <- lm(M ~ T)

#New Methods
time.start <- Sys.time()
mod.1 <- mediate.binary(mmodel, ymodel, sims=1000, boot=TRUE)
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(mmodel, ymodel, sims=1000)
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)

#############################################################
#                     With Covariates                       #
#############################################################
rm(list=ls())
set.seed(3)
source("binary.med.R")

n <-  1000
#Simulate X Vars
T <- round(runif(n), 0)
X.1 <- rnorm(n)
X.2 <- rnorm(n)
#Population Values
alpha.2 <- .25
alpha.3 <- .25
beta.2 <- .25
beta.3 <- .25
gamma <- .25

M <- alpha.2 + beta.2*T + beta.2*X.1 + beta.2*X.2 + rnorm(n)
latenty <-  alpha.3 + beta.3*T + gamma*M + beta.2*X.1 + beta.2*X.2 + rnorm(n)
Y <- latenty
Y[latenty <  1] <- 0
Y[latenty >= 1] <- 1

ymodel <- glm(Y ~ T + M + X.1 + X.2, family=binomial(link=probit))
mmodel <- lm(M ~ T + X.1 + X.2)

#New Method
time.start <- Sys.time()
mod.1 <- mediate.binary(mmodel, ymodel, sims=1000, boot=TRUE)
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(mmodel, ymodel, sims=1000)
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)

