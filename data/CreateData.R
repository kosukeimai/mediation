library(foreign)
library(MASS)

rm(list=ls())
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/analysis/example/JOBS")
job <- read.dta("Jobs-Combined-Data.dta", convert.underscore=FALSE)
job$treat<-as.numeric(job$treat=="exp")
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/data")
save(job,file="job.RData")


#Create global parameters and storage
rm(list=ls())
set.seed(1)
n<-4000

ALPHA.2 <- .25
ALPHA.3 <- .25
BETA.2 <- .25
BETA.3 <- .25
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

#Non-Interaction
#Generate Data
X.1 <- rnorm(n)
T <- round(runif(n), 0)
M.cont <- ALPHA.2 + BETA.2*T + ETA.2*X.1 + e[,1]
latentm <-  ALPHA.2 + BETA.2*T + ETA.2*X.1 + e[,1]
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

Y.dich.temp <-  ALPHA.3 + BETA.3*T + GAMMA*M.cont + ETA.3*X.1+ e[,2]
Y.dich<-0
Y.dich[Y.dich.temp  <  0] <- 0
Y.dich[Y.dich.temp  >= 0] <- 1
Y.cont.dich<-ALPHA.3 + BETA.3*T + GAMMA*M.dich + ETA.3*X.1+ e[,2]



#Generate interaction data
M.cont.T.int<-M.cont*T 
M.dich.T.int<-M.dich*T
M.ord.T.int.1<-M.ord.1*T
M.ord.T.int.2<-M.ord.2*T
M.ord.T.int.3<-M.ord.3*T
M.ord.T.int.4<-M.ord.4*T


Y.cont.int <-  ALPHA.3 + BETA.3*T + GAMMA*M.cont+ KAPPA*M.cont.T.int + ETA.3*X.1+ e[,2]
Y.ord.int <-  ALPHA.3 + BETA.3*T + GAMMA.2*M.ord.2+GAMMA.3*M.ord.3+GAMMA.4*M.ord.4 + ETA.3*X.1+ e[,2]
Y.dich.temp <-  ALPHA.3 + BETA.3*T + GAMMA*M.cont +KAPPA*M.cont.T.int+ ETA.3*X.1+ e[,2]
Y.dich.int<-0
Y.dich.int[Y.dich.temp  <  0] <- 0
Y.dich.int[Y.dich.temp  >= 0] <- 1
Y.cont.dich.int<-ALPHA.3 + BETA.3*T + GAMMA*M.dich + KAPPA*M.dich.T.int+ ETA.3*X.1+ e[,2]


setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/data")
sim<-as.data.frame(cbind(T,X.1,M.cont,M.dich,M.ord,Y.cont,Y.cont.dich,Y.dich,Y.ord,Y.cont.int,Y.dich.int,Y.ord.int,Y.cont.dich.int))
save(sim,file="sim.RData")




#Analytic answers
####################
#Cont Mediator - Cont outcome 
####################


ACME<-GAMMA*BETA.2
TAU<-BETA.3+GAMMA*BETA.2
PROPMED<-ACME/TAU
ContContAnswers<-cbind(ACME,TAU,PROPMED)


####################
#Dichotomous outcomes cont mediator
####################
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


####################
#Dichotomous Mediator - Cont outcome
####################

ACME<-mean(GAMMA*(pnorm(ALPHA.2+BETA.2+ETA.2*X.1)-pnorm(ALPHA.2+ETA.2*X.1)))
#ACME<-mean(GAMMA*(pnorm(-ALPHA.2-ETA.2*X.1)-pnorm(-ALPHA.2-BETA.2-ETA.2*X.1))) #Li's expression
TAU<-BETA.3+ACME
PROPMED<-ACME/TAU
DichMediatorAnswers<-cbind(ACME,TAU,PROPMED)


#INTERACTIONS

####################
#Cont Mediator - Cont outcome - Interaction
####################

delta.0<-BETA.2*(GAMMA+KAPPA)
delta.1<-BETA.2*(GAMMA+KAPPA)

TAU<-BETA.2*GAMMA+BETA.3+KAPPA*(ALPHA.2+BETA.2+ETA.2*mean(X.1))
ContContAnswers.int<-cbind(delta.0,delta.1,TAU)



####################
#Dichotomous Mediator - Cont outcome - Interaction
####################

t<-0
delta.0<- mean((GAMMA+KAPPA*t)*(pnorm(ALPHA.2+BETA.2+ETA.2*X.1)-pnorm(ALPHA.2+ETA.2*X.1)))
t<-1
delta.1<- mean((GAMMA+KAPPA*t)*(pnorm(ALPHA.2+BETA.2+ETA.2*X.1)-pnorm(ALPHA.2+ETA.2*X.1)))
ACME<-(delta.0+delta.1)/2
TAU<-BETA.3+ACME
PROPMED<-ACME/TAU
DichMediatorAnswers.int<-cbind(delta.0,delta.1,TAU,PROPMED)


####################
#Cont Mediator - Dict outcome - Interaction
####################
#must be done

ACME<-1
TAU<-1
delta.0<-0
delta.1<-1
PROPMED<-ACME/TAU
DichOutcomesAnswers.int<-cbind(delta.0,delta.1,TAU,PROPMED)


answertag<-"Order: cont-cont, dichO-contM, ContO-DichM, repeat but w interactions."


sim.answers<-t(cbind(answertag, ContContAnswers,DichOutcomesAnswers,DichMediatorAnswers,ContContAnswers.int,DichOutcomesAnswers.int,DichMediatorAnswers.int))
save(sim.answers,file="sim.answers.RData")
