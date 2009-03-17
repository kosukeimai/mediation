
###########################################
#          With Covariates                #
##########################################
library(MASS)
rm(list=ls())
set.seed(3)

setwd("~/Documents/Mediation Analysis/Local/mediation/mediation/R")
source("cont.med.R")
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
X.2 <- round(runif(n, 1, 4), 0)

M <- alpha.2 + beta.2*T + beta.2*X.1 + beta.2*X.2 + rnorm(n)
Y <- alpha.3 + beta.3*T + gamma*M + beta.2*X.1 + beta.2*X.2 + rnorm(n)

mmodel <- lm(M ~ T + X.1 + X.2)
model.y <- lm(Y ~ T + M + X.1 + X.2)
lm(Y ~ T + X.1 + X.2)

mod.1 <- mediate.cont(mmodel, model.y, sims=1000, boot=TRUE, T="T", M="M")
mod.2 <- mediate.cont(mmodel, model.y, sims=1000, T="T", M="M")

summary(mod.1)
summary(mod.2)

#Mediation Effect
mmodel$coef[2]*model.y$coef[3] 

#Total Effect
mmodel$coef[2]*model.y$coef[3] + model.y$coef[2]

#Proportion
(mmodel$coef[2]*model.y$coef[3])/(mmodel$coef[2]*model.y$coef[3] + model.y$coef[2])

#Population Values
#Indirect
beta.2*gamma
#Direct
beta.3
#Total
beta.3 + beta.2*gamma


#Test Relaxing No Interaction Assumption
T <- round(runif(n), 0)
M <- alpha.2 + beta.2*T + beta.2*X.1 + beta.2*X.2 + rnorm(n)
Y <-  alpha.3 + beta.3*T + gamma*M + kappa*T*M + beta.2*X.1 + beta.2*X.2 + rnorm(n)

mmodel <- lm(M ~ T + X.1 + X.2)
model.y <- lm(Y ~ T + M + T:M + X.1 + X.2)


mod.1 <- mediate.cont(mmodel, model.y, sims=1000, boot=TRUE, INT=TRUE, T="T", M="M")
mod.2 <- mediate.cont(mmodel, model.y, sims=1000, INT=TRUE, T="T", M="M")

summary(mod.1)
summary(mod.2)

#Mediation Effect Under Control
d0 <- mmodel$coef[2]*model.y$coef[3] 
#Mediation Effect Under Treatment
d1 <- mmodel$coef[2]*(model.y$coef[3] + model.y$coef[4])
#Direct Effects
Bm.X <- model.y$coef[4]*mean(X.1) + model.y$coef[5]*mean(X.2)

z0 <- model.y$coef[2] + model.y$coef[6]*(mmodel$coef[1] + Bm.X)
z1 <- model.y$coef[2] + model.y$coef[6]*(mmodel$coef[1] + mmodel$coef[2] + Bm.X)


d0
d1
z0
z1
z1 + d0
z0 + d1


#Total Effect
mmodel$coef[2]*model.y$coef[3] + model.y$coef[2] + model.y$coef[6]*(mmodel$coef[1] + mmodel$coef[2])

mean(c(d0,d1))/(mmodel$coef[2]*model.y$coef[3] + model.y$coef[2] + model.y$coef[6]*(mmodel$coef[1] + mmodel$coef[2]))

#Pop Values
#Indirect Effects
#Delta_0
beta.2*gamma
#Delta_1
beta.2*(gamma + kappa)
#Zeta_0
beta.3 + kappa*(alpha.2)
#Zeta_1
beta.3 + kappa*(alpha.2 + beta.2)
#Tau
beta.2*gamma + beta.3 + kappa*(alpha.2 + beta.2)


############################################################
# Without Covariates
#############################################################
library(MASS)
rm(list=ls())
set.seed(3)

setwd("~/Documents/Mediation Analysis/Local/mediation/mediation/R")
source("cont.med.R")

n <- 1000
#Population Values
alpha.2 <- .25
alpha.3 <- .25
beta.2 <- .25
beta.3 <- .25
gamma <- .25
kappa <- .25

#Treatment Indicator
T <- round(runif(n), 0)
M <- alpha.2 + beta.2*T + rnorm(n)
Y <-  alpha.3 + beta.3*T + gamma*M + rnorm(n)

#Step 1: Fit the model Y=T,M using OLS
ymodel <- lm(Y ~ T + M)
gamma.est <- ymodel$coef[3]
beta3.est <- ymodel$coef[2]
#Step 2: Fit the model M=T
mmodel <- lm(M ~ T)
beta2.est <- mmodel$coef[2]


mod.1 <- mediate.cont(mmodel, ymodel, sims=1000, boot=TRUE, T="T", M="M")
mod.2 <- mediate.cont(mmodel, ymodel, sims=1000, T="T", M="M")

summary(mod.1)
summary(mod.2)

#Now show that this is equal to the product of coefficients
gamma.est*beta2.est
beta3.est
(gamma.est*beta2.est) + beta3.est
#The DGP estimate is
gamma*beta.2
beta.3
(gamma*beta.2) + beta.3

n <- 1000
#Test Relaxing No Interaction Assumption
#Population Values
alpha.2 <- .25
alpha.3 <- .25
beta.2 <- -.25
beta.3 <- .25
gamma <- .25
kappa <- .25

T <- round(runif(n), 0)
M <- alpha.2 + beta.2*T + rnorm(n)
Y <-  alpha.3 + beta.3*T + gamma*M + kappa*T*M + rnorm(n)

mmodel <- lm(M ~ T)
ymodel <- lm(Y ~ T + M + T:M)

z <- lm(M ~ T)
model.y <- lm(Y ~ T + M + T:M)

mod.1 <- mediate.cont(mmodel, ymodel, sims=100, boot=TRUE, INT=TRUE, T="T", M="M")
mod.2 <- mediate.cont(mmodel, ymodel, sims=1000, INT=TRUE, T="T", M="M")

summary(mod.1)
summary(mod.2)

#Mediation Effect Under Control
mmodel$coef[2]*ymodel$coef[3] 
#Mediation Effect Under Treatment
mmodel$coef[2]*(ymodel$coef[3] + ymodel$coef[4])
#Direct Effects
ymodel$coef[2] + ymodel$coef[4]*(mmodel$coef[1] + mmodel$coef[2])
ymodel$coef[2] + ymodel$coef[4]*(mmodel$coef[1])
#Total Effect
(mmodel$coef[2]*ymodel$coef[3]) + (ymodel$coef[2] + ymodel$coef[4]*(mmodel$coef[1] + mmodel$coef[2]))
(mmodel$coef[2]*(ymodel$coef[3] + ymodel$coef[4])) + (ymodel$coef[2] + ymodel$coef[4]*(mmodel$coef[1]))
#Pop Values
beta.2*gamma
beta.2*(beta.3 + kappa)

B <- 100

