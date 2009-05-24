library(foreign)
library(MASS)

rm(list=ls())
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/analysis/example/JOBS")
job <- read.dta("Jobs-Combined-Data.dta", convert.underscore=FALSE)
job$treat<-as.numeric(job$treat=="exp")
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/data")
save(job,file="job.RData")



#Generate simulation data for mediation function
rm(list=ls())

n<-5000

ALPHA.2 <- .25
ALPHA.3 <- .25
BETA.2 <- .25
BETA.3 <- -.25
GAMMA <- .25
GAMMA.2 <- .25
GAMMA.3 <- .5
GAMMA.4 <- .75

KAPPA<-.25
ETA.3<-.25
ETA.2<-.25
SIGMA.3 <- 1
SIGMA.3.SQ<-SIGMA.3^2
SIGMA.2 <- 1
SIGMA.2.SQ<-1

cov <- 0
Sigma <- matrix(c(SIGMA.2^2,cov,cov,SIGMA.3^2), 2,2)
e <- mvrnorm(n, rep(0,2), Sigma)


#Generate Data
X.1 <- rnorm(n)
T <- round(runif(n), 0)
M.cont <- ALPHA.2 + BETA.2*T + ETA.2*X.1 + e[,1]
latentm <-  ALPHA.2 + BETA.2*T + ETA.2*X.1 + e[,2]
M.dich <- latentm
M.dich[latentm <  0] <- 0
M.dich[latentm >= 0] <- 1
M.ord <- latentm
M.ord[latentm < -.3] <- 1
M.ord[latentm >= -.3 & latentm < .45] <- 2
M.ord[latentm >= .45 & latentm < 1.1] <- 3
M.ord[latentm >= 1.1] <- 4

M.ord.1<-as.numeric(M.ord==1)
M.ord.2<-as.numeric(M.ord==2)
M.ord.3<-as.numeric(M.ord==3)
M.ord.4<-as.numeric(M.ord==4)


Y.cont <-  ALPHA.3 + BETA.3*T + GAMMA*M.cont + ETA.3*X.1+ e[,2]
Y.ord <-  ALPHA.3 + BETA.3*T + GAMMA.2*M.ord.2+GAMMA.3*M.ord.3+GAMMA.4*M.ord.4 + ETA.3*X.1+ e[,2]

latenty.dich <-  ALPHA.3 + BETA.3*T + GAMMA*M.dich + ETA.3*X.1+ e[,2]
Y.dich <- latenty.dich
Y.dich[latenty.dich <  0] <- 0
Y.dich[latenty.dich >= 0] <- 1


#Analytic answers
#Dichotomous outcomes
IKY_DELTA<-function(alpha.2,alpha.3,beta.2,beta.3,gamma, eta.2, eta.3, X, sigma.2.sq) {
std.scalar<-sqrt(sigma.2.sq*gamma^2+1)
t<-0
temp1<-pnorm( (alpha.3 + beta.3*t + eta.3*X+ gamma*(alpha.2+beta.2+eta.2*X))/(std.scalar) )
temp2<-pnorm( (alpha.3 + beta.3*t + eta.3*X+ gamma*(alpha.2+eta.2*X))/(std.scalar) )
delta.0 <-mean(temp1-temp2)
t<-1
temp1<-pnorm( (alpha.3 + beta.3*t + eta.3*X+ gamma*(alpha.2+beta.2+eta.2*X))/(std.scalar) )
temp2<-pnorm( (alpha.3 + beta.3*t + eta.3*X+ gamma*(alpha.2+eta.2*X))/(std.scalar) )
delta.1 <-mean(temp1-temp2)
DELTA<-(delta.0+delta.1)/2
return(DELTA)
}

#USES EQUATION 85 and 86 of IKT
IKY_TAU<-function(alpha2,alpha3,beta2,beta3,gamma, eta2, eta3, X, sigma.2.sq) {
std.scalar <- sqrt(gamma^2*sigma.2.sq + 1)
alpha1<-(alpha3+gamma*alpha2)/std.scalar
beta1<-(beta3+gamma*beta2)/std.scalar
eta1<-(eta3+eta2*gamma)/(std.scalar)
TAU<-mean(pnorm(alpha1+beta1+eta1*X)-pnorm(alpha1+eta1*X))
return(TAU)
}

ACME<-IKY_DELTA(ALPHA.2,ALPHA.3,BETA.2,BETA.3,GAMMA, ETA.2, ETA.3, X.1, SIGMA.2.SQ)
TAU<-IKY_TAU(ALPHA.2,ALPHA.3,BETA.2,BETA.3,GAMMA, ETA.2, ETA.3, X.1, SIGMA.2.SQ)
PROPMED<-ACME/TAU
DichOutcomesAnswers<-cbind(ACME,TAU,PROPMED)

ACME<-GAMMA*BETA.2
TAU<-BETA.3
ContContAnswers<-cbind(ACME,TAU)

sim<-as.data.frame(cbind(T,X.1,M.cont,M.dich,M.ord,Y.cont,Y.dich,Y.ord))
sim.answers<-cbind(DichOutcomesAnswers,ContContAnswers)
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/data")
save(sim,file="sim.RData")
save(sim.answers,file="sim.answers.RData")
