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

#Standard Method - With Covariates
mmodel <- glm(M ~ T + X.1 + X.2, family=binomial(link=probit))
alpha2.est<- mmodel$coef[1]
beta2.est <- mmodel$coef[2]
ymodel <- lm(Y ~ T + M + X.1 + X.2)
gamma.est <- ymodel$coef[3]
beta3.est <- ymodel$coef[2]
delta.true <- gamma*(pnorm(alpha.2 + beta.2) - pnorm(alpha.2))
tau.true <- beta.3 + gamma*(pnorm(alpha.2 + beta.2) - pnorm(alpha.2))
pr.true <- delta.true/tau.true
delta.t.est <- gamma.est*(pnorm(alpha2.est + beta2.est) - pnorm(alpha2.est))
scalar_Li <- pnorm(mmodel$coef[1] + mmodel$coef[2])
delta.est <- gamma.est*beta2.est*scalar_Li
tau.est <- beta3.est + gamma.est*(pnorm(alpha2.est + beta2.est) - pnorm(alpha2.est))
pr.est <-  delta.est/tau.est

#New Methods
time.start <- Sys.time()
mod.1 <- mediate.binary(mmodel, ymodel, sims=1000, boot=TRUE)
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(mmodel, ymodel, sims=1000)
print (round((Sys.time()-time.start),1))

cat("True Point Estimate: ", delta.true, "\n")
cat("True Proportion Mediated: ", pr.true, "\n" )
cat("Point Estimate Li Method: ", delta.est, "\n")
cat("Proportion Est Li Method: ", pr.est, "\n")
summary(mod.1)
summary(mod.2)


