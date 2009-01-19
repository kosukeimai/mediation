library(foreign)
library(MASS)

setwd("/Users/Luke/Documents/Mediation Analysis/Local/mediation/analysis/example/JOBS")

setwd("/Users/keele/Documents/Mediation Analysis/Local/mediation/analysis/example/JOBS")

rm(list=ls())
set.seed(314567)
job <- read.dta("JOBS-cont.dta")
attach(job)
n <- length(job$work)
time.start <- Sys.time()

logit <- function(xb){
    1/(1+exp(-xb))
    }
    
####################################################################################
# Draw From Multivariate Normal Posterior
#Step 1
#Fit the model Y=T,M using logit
TMmodel <- glm(work ~ treat + job_seek + depress1  , data=job, family=binomial(link=logit))

#Step 2
#Fit the model M=T
MModel <- lm(job_seek ~ treat + depress1 , data=job)
#Need to save the model's error variance term below. Extract the "residual standard error" from the model
sigma <- summary(MModel)$sigma

#Step 3
#generate vectors that contain simulated draws of the coefficients
#extract coefficients and variance expressions in order to do simulation draws
MModel.coef <- MModel$coef
MModel.var.cov <- vcov(MModel)

TMmodel.coef <- TMmodel$coef
TMmodel.var.cov <- vcov(TMmodel)

sims <- 1000
MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov)
TMmodel <- mvrnorm(sims, mu=TMmodel.coef, Sigma=TMmodel.var.cov)

#####################################################################################
#Step 4
 mean.sim <- mean(apply(MModel, 1, sum))
 error <- rnorm(sims, mean.sim, sd=sigma)  
 PredictM1temp <- as.matrix(apply(MModel, 1, sum))
 PredictM1  <- matrix(PredictM1temp, n, sims, byrow=TRUE)
 PredictM1 <- PredictM1 + error
 PredictM0temp <- cbind(MModel[,1], MModel[,2] + error)
 PredictM0temp <- as.matrix(apply(PredictM0temp, 1, sum))
 PredictM0 <- matrix(PredictM0temp, n, sims, byrow=TRUE)
 
#gen matrices required to hold necessary values for making predictions of the outcome variable
Prob1_temp <- matrix(,nrow=n, ncol=sims)
Prob0_temp <- matrix(,nrow=n, ncol=sims)

#generate the linear combinations of the coefficients and data values and store in _temp variable

#T=1
for (j in 1:sims) {
        Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1+ TMmodel[j,3]*PredictM1[,j] + TMmodel[j,4]*job$depress1
        Prob0_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1+ TMmodel[j,3]*PredictM0[,j] + TMmodel[j,4]*job$depress1
		}

    
Prob1 <- apply(Prob1_temp, 2, logit)
Prob0 <- apply(Prob0_temp, 2, logit)
Delta1 <- Prob1 -Prob0 
rm(Prob1, Prob0, Prob1_temp, Prob0_temp)
Prob1_temp <- matrix(,nrow=n, ncol=sims)
Prob0_temp <- matrix(,nrow=n, ncol=sims)

#T=0
for (j in 1:sims) {
        Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*0+ TMmodel[j,3]*PredictM1[,j] + TMmodel[j,4]*job$depress1
        Prob0_temp[,j] <- TMmodel[j,1]+ TMmodel[j,2]*0+TMmodel[j,3]*PredictM0[,j]  + TMmodel[j,4]*job$depress1
		}

Prob1 <- apply(Prob1_temp, 2, logit)
Prob0 <- apply(Prob0_temp, 2, logit)
Delta0 <- Prob1 - Prob0 
rm(Prob1, Prob0, Prob1_temp, Prob0_temp)


#Step 4.3, average over the observations of i
Avg_Delta1 <- t(as.matrix(apply(Delta1, 2 , mean)))
Avg_Delta0 <- t(as.matrix(apply(Delta0, 2 , mean)))

#Now we add delta 1 to delta 0 and divide by 2
Delta1_Delta0 <- Avg_Delta1 + Avg_Delta0
Delta1_Delta0 <- Delta1_Delta0/2

Prob1_temp <- matrix(,nrow=n, ncol=sims)
Prob0_temp <- matrix(,nrow=n, ncol=sims)

#Calculate the vector for the denominator for proportion mediated equation
#Must be careful with missing values in covariates.
for (j in 1:sims) {
        Prob1_temp[,j] <- TMmodel[j,1] + TMmodel[j,2]*1+ TMmodel[j,3]*PredictM1[,j] + TMmodel[j,4]*job$depress1
        Prob0_temp[,j] <- TMmodel[j,1]+ TMmodel[j,2]*0+TMmodel[j,3]*PredictM0[,j]   + TMmodel[j,4]*job$depress1
}

Prob1 <- apply(Prob1_temp, 2, logit)
Prob0 <- apply(Prob0_temp, 2, logit)
Tau <- Prob1 - Prob0 

#Now average over the i's 
Avg_Tau <- t(as.matrix(apply(Tau, 2, mean)))

#To calcuate proportion mediated we take for each draw j the ratio of 1/2(Delta1+Delta0)/Tau
IKY_ProportionMediated_temp <- Delta1_Delta0 / Avg_Tau

IKY_ProportionMediated <- mean(t(IKY_ProportionMediated_temp))
IKY_ProportionMediated_median <- median(t(IKY_ProportionMediated_temp))

CausalEffectVectors <- cbind(t(IKY_ProportionMediated_temp), t(Delta1_Delta0), t(Avg_Tau))

#The vector IKY_ProportionMediated_temp can contain extreme values depending on the draws from the distributions
#The median of the above measure was regularly almost exactly that of Mac's estimate, but the mean could be pulled away
#the below excludes some of these outliers
#IKY_ProportionMediated_temp_out <- matrix(,nrow=1, ncol=sims)
#for (j in 1:sims) {
# if (abs(IKY_ProportionMediated_temp[,j])>1) {
# IKY_ProportionMediated_temp_out[,j]<-NA
# }
# else {
# IKY_ProportionMediated_temp_out[,j]<-IKY_ProportionMediated_temp[,j]
# }
# }

print (round((Sys.time()-time.start),1))

#################################################################################################


IKY_ProportionMediated
IKY_ProportionMediated_median
#Confidence Interval
quantile(IKY_ProportionMediated_temp_out,c(.025,.975), na.rm=TRUE)
quantile(IKY_ProportionMediated_temp,c(.025,.975), na.rm=TRUE)

pdf("pdens.pdf", width=5, height=5, onefile=FALSE, paper="special")
plot(density(IKY_ProportionMediated_temp), main="Posterior Density of Proportion Mediated")
dev.off()