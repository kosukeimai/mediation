library(foreign)

rm(list=ls())
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/analysis/example/JOBS")
job <- read.dta("Jobs-Combined-Data.dta", convert.underscore=FALSE)
job$treat<-as.numeric(job$treat=="exp")
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/data")
save(job,file="job.RData")



#Generate simulation data
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
SIGMA.SQ<-SIGMA.3^2

sigma22 <- 1
sigma33 <- 1
cov <- 0
Sigma <- matrix(c(sigma22^2,cov,cov,sigma33^2), 2,2)
e <- mvrnorm(n, rep(0,2), Sigma)


#Generate Data
X.1 <- rnorm(n)
T <- round(runif(n), 0)
M.cont <- ALPHA.2 + BETA.2*T + ETA.2*X.1 + e[,1]
latentm <-  ALPHA.3 + BETA.3*T + GAMMA*M.cont + ETA.3*X.1+ e[,2]
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



sim<-cbind(T,X.1,M.cont,M.dich,M.ord,Y.cont,Y.dich,Y.ord)
setwd("H:/imai_methods/Mediation/CheckoutCVS/mediation/mediation/data")
save(sim,file="sim.RData")
