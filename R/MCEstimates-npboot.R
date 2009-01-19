library(foreign)

setwd("/Users/Luke/Documents/Mediation Analysis/Local/mediation/analysis/example/JOBS")

setwd("/Users/keele/Documents/Mediation Analysis/Local/mediation/analysis/example/JOBS")

rm(list=ls())
set.seed(314567)
job <- read.dta("JOBS-cont.dta")
attach(job)

#Fit Models - work with generic model objects
MModel <- lm(job_seek ~ treat + depress1 , data=job)
TMmodel.t <- glm(work ~ treat + job_seek + depress1  , data=job, family=binomial(link=logit))
TMmodel.c <- glm(work ~ job_seek + depress1  , data=job, family=binomial(link=logit))

#Start Function Here

model.m <- MModel
model.y.t <- TMmodel.t
model.y.c <- TMmodel.c
k.t <- length(names(model.y.t$coef))
k.c <- length(names(model.y.c$coef))

n <- length(job$work)
B <- 25
pr.med <- matrix(NA, B, 1)
time.start <- Sys.time()

Call.M <- model.m$call
Call.Y.t <- model.y.t$call
Call.Y.c <- model.y.c$call
	
for (b in 1:B) {
	
index <- sample(1:n,n,repl=TRUE)
Call.M$data <- job[index,]	
Call.Y.t$data <- job[index,]	
Call.Y.c$data <- job[index,]

new.fit.M <- eval.parent(Call.M)

#Need to save the model's error variance term below. Extract the "residual standard error" from the model
sigma <- summary(new.fit.M)$sigma
mean.sim <- new.fit.M$coef[1] + new.fit.M$coef[2]
error <- rnorm(n, 0, sd=sigma) 
PredictM1 <- new.fit.M$coef[1] + new.fit.M$coef[2] + error
PredictM0 <- new.fit.M$coef[1] + error 

#PredictM1 <- rep(PredictM1, n)
#PredictM0 <- rep(PredictM0, n)

new.fit.t <- eval.parent(Call.Y.t)
new.fit.c <- eval.parent(Call.Y.c)

#Treatment Predictions
pred.data.t <- data.frame(1,PredictM1, job[index,]$depress1)
pred.data.c <- data.frame(1,PredictM0, job[index,]$depress1)
names(pred.data.t) <- names(model.y.t$coef[2:k.t])
names(pred.data.c) <- names(model.y.c$coef[2:k.t])

pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)
pr.0 <- predict(new.fit.t, type="response", newdata=pred.data.c)
pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
delta.1 <- pr.mat[,1] - pr.mat[,2]

rm(pred.data.t, pred.data.c, pr.1, pr.0)

#Control Predictions
pred.data.t <- data.frame(PredictM1, job[index,]$depress1)
pred.data.c <- data.frame(PredictM0, job[index,]$depress1)
names(pred.data.t) <- names(model.y.t$coef[2:k.c])
names(pred.data.c) <- names(model.y.c$coef[2:k.c])

pr.1 <- predict(new.fit.c, type="response", newdata=pred.data.t)
pr.0 <- predict(new.fit.c, type="response", newdata=pred.data.c)
pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
delta.0 <-pr.mat[,1] - pr.mat[,2]

avg_delta.1 <- mean(delta.1)
avg_delta.0 <- mean(delta.0)
med.eff <- (avg_delta.1 + avg_delta.0)/2

rm(pred.data.t, pred.data.c, pr.1, pr.0)

#Calculate Total Effect
pred.data.t <- data.frame(1,PredictM1, job[index,]$job_seek, job[index,]$depress1)
names(pred.data.t) <- names(model.y.t$coef[2:k.t])
pr.1 <- predict(new.fit.t, type="response", newdata=pred.data.t)

pred.data.c <- data.frame(PredictM0, job[index,]$job_seek, job[index,]$depress1)
names(pred.data.c) <- names(model.y.c$coef[2:k.c])
pr.0 <- predict(new.fit.c, type="response", newdata=pred.data.c)
pr.mat <- as.matrix(na.omit(cbind(pr.1, pr.0)))
tau <- pr.mat[,1] - pr.mat[,2]

rm(pred.data.t, pred.data.c, pr.1, pr.0)

avg.tau <- mean(tau)
pr.med[b] <- med.eff/avg.tau


}

#Confidence Interval
ci <- quantile(pr.med,c(.025,.975), na.rm=TRUE)

cat("Proportion Mediated: ", mean(pr.med), "95% Confidence Interval", ci, "\n")



pdf("pdens.pdf", width=5, height=5, onefile=FALSE, paper="special")
plot(density(pr.med), main="Posterior Density of Proportion Mediated")
dev.off()