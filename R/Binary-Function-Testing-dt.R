#This file checks the Luke mediation function for binary outcomes 

library(MASS)
library(foreign)
rm(list=ls())

#setwd("~/Documents/Mediation Analysis/Local/mediation/mediation/R")
#setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/R")

source("binary.med.R")

MacKinPropMed<-function(a,b,c) {
alpha2 <- b$coefficients[1]
beta2 <-b$coefficients[2]
gamma <- c$coefficients[3]
beta3<- c$coefficients[2]
alpha3<- c$coefficients[1]
sigma <- summary(b)$sigma
beta1 <- a$coefficients[2]
var.res_mmodel <- sigma^2
std.scalar <- sqrt(gamma^2*var.res_mmodel + 1)
pr.med_scale_1 <- beta2*gamma/(beta1*std.scalar)
pr.med_scale_2 <- beta2*gamma/(beta2*gamma+beta3)
return(pr.med_scale_1,pr.med_scale_2 )
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





ALPHA.2 <- .25
ALPHA.3 <- .25
BETA.2 <- .25
BETA.3 <- 2
GAMMA <- .5
ETA.3<-.5
ETA.2<-.5
SIGMA.2<-1
SIGMA.2.SQ<-SIGMA.2^2
RHO<-0
sims<-3000
n<-100000
X.1 <- rnorm(n)


sigma22 <- SIGMA.2
sigma33 <- 1
cov <- RHO*sigma22*sigma33
Sigma <- matrix(c(sigma22^2,cov,cov,sigma33^2), 2,2)
e <- mvrnorm(n, rep(0,2), Sigma)
var(e)
cor(e[,1], e[,2])

#Generate Data
T <- round(runif(n), 0)
M <- ALPHA.2 + BETA.2*T + ETA.2*X.1 + e[,1]
latenty <-  ALPHA.3 + BETA.3*T + GAMMA*M + ETA.3*X.1+ e[,2]
#M <- ALPHA.2 + BETA.2*T  + e[,1]
#latenty <-  ALPHA.3 + BETA.3*T + GAMMA*M + e[,2]
Y <- latenty
Y[latenty <  0] <- 0
Y[latenty >= 0] <- 1


#Estimating regressions
a <- glm(Y ~ T +X.1, family=binomial(link=probit))
b <- lm(M ~ T+X.1)
c <- glm(Y ~ T + M+X.1, family=binomial(link=probit))

#a <- glm(Y ~ T , family=binomial(link=probit))
#b <- lm(M ~ T)
#c <- glm(Y ~ T + M, family=binomial(link=probit))


alpha2 <- b$coefficients[1]
beta2 <-b$coefficients[2]
gamma <- c$coefficients[3]
beta3<- c$coefficients[2]
alpha3<- c$coefficients[1]
sigma <- summary(b)$sigma
beta1 <- a$coefficients[2]
var.res_mmodel <- sigma^2
std.scalar <- sqrt(gamma^2*var.res_mmodel + 1)
Beta2Gamma<-b$coef["T"]*c$coef["M"]
Beta1Scalar<-a$coefficients[2]*std.scalar 
Beta2GammaBeta3<-Beta2Gamma+beta3
eta2<-b$coef[3]
eta3<-c$coef[4]

DELTA.01<-IKY_DELTA(ALPHA.2,ALPHA.3,BETA.2,BETA.3,GAMMA, ETA.2, ETA.3, X.1, SIGMA.2.SQ)
     TAU<-IKY_TAU(ALPHA.2,ALPHA.3,BETA.2,BETA.3,GAMMA, ETA.2, ETA.3, X.1, SIGMA.2.SQ)
PROPMED<-DELTA.01/TAU
WangProp<-BETA.2*GAMMA/(BETA.2*GAMMA+BETA.3)

DELTA.01_sample<-IKY_DELTA(alpha2,alpha3,beta2,beta3,gamma, eta2, eta3, X.1, var.res_mmodel)
     TAU_sample<-IKY_TAU(alpha2,alpha3,beta2,beta3,gamma, eta2, eta3, X.1, var.res_mmodel)
PropMed_sample<-DELTA.01_sample/TAU_sample
PropMed_sample

MacKinPropMed(a,b,c)

LukeCode<- mediate.binary(b, c, sims=sims, T="T", M="M")
LukeACME<-(LukeCode$d0+LukeCode$d1)/2
LukePropMed<-LukeACME/LukeCode$tau.coef


#(LukeCode$d0+LukeCode$d1+LukeCode$z0+LukeCode$z1)/4

#LukeCode$z1+LukeCode$d1

t<-0
temp1<-pnorm( (alpha3 + beta3*t + eta3*X.1+ gamma*(alpha2+beta2+eta2*X.1))/(std.scalar) )
temp2<-pnorm( (alpha3 + beta3*t + eta3*X.1+ gamma*(alpha2+eta2*X.1))/(std.scalar) )
delta.0 <-mean(temp1-temp2)
t<-1
temp1<-pnorm( (alpha3 + beta3*t + eta3*X.1+ gamma*(alpha2+beta2+eta2*X.1))/(std.scalar) )
temp2<-pnorm( (alpha3 + beta3*t + eta3*X.1+ gamma*(alpha2+eta2*X.1))/(std.scalar) )
delta.1 <-mean(temp1-temp2)
delta.0 
delta.1


print("IKT ACME")
print(DELTA.01)
print("LUKE ACME")
print(LukeACME)
print("Beta2Gamma")
print(Beta2Gamma)
print("IKT TAU")
print(TAU)


print("IKY Prop Mediated")
print(PROPMED)
print("Luke Point Estimate")
print(LukeCode)
print("Luke d0+d1/tau point estimate")
print(LukePropMed)
print("MacKinon Prop Mediated")
print(MacKinPropMed(a,b,c))
print("Wang Proportion Mediated")
print(WangProp)




#################################################################################################
#Generate coefficients to use

#We don't run this code right now.
if (F) {

setwd("H:/imai_methods/Mediation/")
CigData<-read.csv("MacKinnonCigUsage.csv", header=TRUE)

a <- glm(usage ~ treatment, data=CigData, family=binomial(link=probit))
b <- lm(intent ~ treatment, data=CigData)
c <- glm(usage ~ treatment + intent, data=CigData, family=binomial(link=probit))



MacKinPropMed(a,b,c)
mod.2 <- mediate.binary(b, c, sims=1000, T="treatment", M="intent")




#Jobs data
setwd("~/Documents/Mediation Analysis/Local/mediation/analysis/example/JOBS")
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/analysis/example/JOBS")

job <- read.dta("Jobs-NoMiss-Binary.dta")




#Binary Outcome
#Models



a <- glm(work1 ~ treat  + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)
b <- lm(job_seek ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
c <- glm(work1 ~ treat + job_seek + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)

NI_M <- "treat + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp"
NI_O <- "treat +job_seek+ depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp"
I_O <- "treat+job_seek+ treat:job_seek + depress1 +age+educ+sex+nonwhite+econ_hard+marital+occp"

NI_M <- "treat + depress1 "
NI_O <- "treat +job_seek+ depress1 "
I_O <- "treat+job_seek+ treat:job_seek + depress1 "



#Models
a<-glm(as.formula(paste("work1 ~", NI_M)),family=binomial(link="probit"), data=job)
b <- lm(as.formula(paste("job_seek ~", NI_M)), data=job)
c <- glm(as.formula(paste("work1 ~", NI_O)),family=binomial(link="probit"), data=job)
d <- glm(as.formula(paste("work1 ~", I_O)),family=binomial(link="probit"), data=job)




MacKinPropMed(a,b,c)


time.start <- Sys.time()
mod.1 <- mediate.binary(b, c, sims=1000, boot=TRUE, T="treat", M="job_seek")
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(b, c, sims=500, T="treat", M="job_seek")
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)





#Double Binary - No Direct Effect
#Models
b <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)
c <- glm(work1 ~ treat + job_dich + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, family=binomial(link="probit"), data=job)

time.start <- Sys.time()
mod.1 <- mediate.binary(b, c, sims=1000, boot=TRUE, T="treat", M="job_dich")
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(b, c, sims=1000, T="treat", M="job_dich")
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)

#Binary Mediator - Seems Fine
rm(list=ls())

setwd("~/Documents/Mediation Analysis/Local/mediation/analysis/example/JOBS")

job <- read.dta("Jobs-NoMiss-Cont.dta")

setwd("~/Documents/Mediation Analysis/Local/mediation/mediation/R")
source("binary.med.R")

b <- glm(job_dich ~ treat + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job, family=binomial(link="probit"))
c <- lm(depress2 ~ treat + job_dich + depress1 + econ_hard + sex + age + occp + marital + nonwhite + educ + income, data=job)

time.start <- Sys.time()
mod.1 <- mediate.binary(b, c, sims=1000, boot=TRUE, T="treat", M="job_dich")
print (round((Sys.time()-time.start),1))

time.start <- Sys.time()
mod.2 <- mediate.binary(b, c, sims=1000, T="treat", M="job_dich")
print (round((Sys.time()-time.start),1))

summary(mod.1)
summary(mod.2)

c$coef[3]*(pnorm(b$coef[1] + b$coef[2]) - pnorm(b$coef[1]))

c$coef[2] + c$coef[3]*(pnorm(b$coef[1] + b$coef[2]) - pnorm(b$coef[1]))


}
