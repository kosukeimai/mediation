
mediate.npdiscrete <- function(M, Y, T, INT=TRUE,boot = FALSE, sims = 1000, conf.level = .95){
	
	samp <- data.frame(na.omit(cbind(M,Y,T)))
	m.cat <- sort(unique(samp$M))
    n <- length(samp$Y)
    n1 <- sum(samp$T)
    n0 <- sum(1-samp$T)
    t <- (1 - conf.level)/2
##########################################################
# Point Estimates
##########################################################

if(boot==FALSE){
#Indirect Effects
	delta_0m <- matrix(NA, length(m.cat), 1)
    for(i in 1:length(m.cat)){
    	delta_0m[i] <- (sum(samp$T[samp$M==m.cat[i]]) * sum(samp$Y[samp$M==m.cat[i] & samp$T==0]))/(n1*(sum((1 - samp$T[samp$M==m.cat[i]]))))
        }
    
    d0 <- sum(delta_0m) - mean(samp$Y[samp$T==0])
    delta_1m <- matrix(NA, length(m.cat), 1)
    
    for(i in 1:length(m.cat)){
    	delta_1m[i] <-  (sum((1 - samp$T[samp$M==m.cat[i]])) * sum(samp$Y[samp$M==m.cat[i] & samp$T==1]))/(n0*(sum(samp$T[samp$M==m.cat[i]])))    
        }
    
   d1 <- mean(samp$Y[samp$T==1]) - sum(delta_1m)
  } else {
  	idx <- seq(1,n,1)
    #A Storage Matrix
	d0.bs <- matrix(NA,sims,1)
	d1.bs <- matrix(NA,sims,1)
	#Now Let's do the Resampling
 	for(b in 1:sims){ 
 	#First take a resample
 	resample <- sample(idx,n,replace=TRUE)
  	samp.star <- samp[resample,]
  	
  		delta_0m <- matrix(NA, length(m.cat), 1)
    for(i in 1:length(m.cat)){
    	delta_0m[i] <- (sum(samp.star$T[samp.star$M==m.cat[i]]) * sum(samp.star$Y[samp.star$M==m.cat[i] & samp.star$T==0]))/(n1*(sum((1 - samp.star$T[samp.star$M==m.cat[i]]))))
        }
    
		d0.bs[b,] <- sum(delta_0m) - mean(samp.star$Y[samp.star$T==0])

    	delta_1m <- matrix(NA, length(m.cat), 1)
    
    for(i in 1:length(m.cat)){
    	delta_1m[i] <-  (sum((1 - samp.star$T[samp.star$M==m.cat[i]])) * sum(samp.star$Y[samp.star$M==m.cat[i] & samp.star$T==1]))/(n0*(sum(samp.star$T[samp.star$M==m.cat[i]])))    
        }
		d1.bs[b,] <- mean(samp.star$Y[samp.star$T==1]) - sum(delta_1m)  	
  }
}
    
##########################################################
# Variance Estimates
##########################################################

pr.t.1 <- sum(samp$T)/n
pr.t.0 <- sum(1-samp$T)/n

lambda.1m <- matrix(NA, length(m.cat), 1)
lambda.0m <- matrix(NA, length(m.cat), 1)
for(i in 1:length(m.cat)){
lambda.1m[i] <- (sum(length(samp$M[samp$M==m.cat[i] & samp$T==1])))/n /pr.t.1
lambda.0m[i] <- (sum(length(samp$M[samp$M==m.cat[i] & samp$T==0])))/n /pr.t.0
}

mu.1m <- matrix(NA, length(m.cat), 1)
mu.0m <- matrix(NA, length(m.cat), 1)

for(i in 1:length(m.cat)){
	mu.1m[i] <- mean(samp$Y[samp$M==m.cat[i] & samp$T==1])
	mu.0m[i] <- mean(samp$Y[samp$M==m.cat[i] & samp$T==0])
	}

# Variance of Delta_0 
    var.delta.0m <- matrix(NA, length(m.cat), 1)
    for(i in 1:length(m.cat)){
        var.delta.0m[i] <-  lambda.1m[i]*(((lambda.1m[i]/lambda.0m[i]) - 2) * var(samp$Y[samp$M==m.cat[i] & samp$T==0]) + ((n0*(1-lambda.1m[i])*mu.0m[i]^2)/n1))
        }

    m.leng <- length(m.cat)
    mterm.d0 <- c()
    for(i in 1:m.leng-1){
        mterm.d0 <- c(mterm.d0, lambda.1m[i] * lambda.1m[(i+1):m.leng] * mu.0m[i] * mu.0m[(i+1):m.leng])
        }
var.delta_0 <- (1/n0) * sum(var.delta.0m) - (2/n1) * sum(mterm.d0) + var(samp$Y[samp$T == 0])/n0 

# Variance of Delta_1
    var.delta.1m <- matrix(NA, length(m.cat), 1)
    
    for(i in 1:length(m.cat)){
        var.delta.1m[i] <-  lambda.0m[i]*(((lambda.0m[i]/lambda.1m[i]) - 2) * var(samp$Y[samp$M==m.cat[i] & samp$T==1]) + ((n1*(1-lambda.0m[i])*mu.1m[i]^2)/n0)) 
        }

    mterm.d1 <- c()
    for(i in 1:m.leng-1){
        mterm.d1 <- c(mterm.d1, lambda.0m[i] * lambda.0m[(i+1):m.leng] * mu.1m[i] * mu.1m[(i+1):m.leng])
        }
var.delta_1 <- (1/n1) * sum(var.delta.1m) - (2/n0)* sum(mterm.d1) + var(samp$Y[samp$T == 1])/n1

        if(boot==FALSE){
        d0.ci <- c(d0 - (-qnorm(t)*sqrt(var.delta_0)), d0 + (-qnorm(t)*sqrt(var.delta_0)))
        d1.ci <- c(d1 - (-qnorm(t)*sqrt(var.delta_1)), d1 + (-qnorm(t)*sqrt(var.delta_1)))
        } else {
        d0 <- mean(d0.bs, na.rm=TRUE)
        d1 <- mean(d1.bs, na.rm=TRUE)
        low <- (1 - conf.level)/2
        high <- 1 - low
        d0.ci <- quantile(d0.bs,c(low,high), na.rm=TRUE)
        d1.ci <- quantile(d1.bs,c(low,high), na.rm=TRUE)
        }

        out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci, INT=INT, conf.level=conf.level,nobs=n,boot=boot)
	    class(out) <- "mediate.npdiscrete"
        out
 }
        
print.mediate.npdiscrete <- function(x, ...){
	print(unlist(x[1:11]))
	invisible(x)
	}

summary.mediate.npdiscrete <- function(object, ...)
    structure(object, class = c("summary.mediate.npdiscrete", class(object)))

print.summary.mediate.npdiscrete <- function(x, ...){
    clp <- 100 * x$conf.level
    cat("\n NonParametric Causal Mediation Analysis for Discrete Mediators \n\n")
    
    if(x$boot==TRUE){
        cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
    } else {
        cat("Large Sample Confidence Intervals\n\n")
    }
    
    printone <- x$INT == FALSE 
    
    if (printone){
        cat("Mediation Effect: ", format(x$d1, digits=4), clp, "% CI ", 
                format(x$d1.ci, digits=4), "\n")
        cat("Sample Size Used:", x$nobs,"\n\n")        
    } else {
        cat("Mediation Effect_0: ", format(x$d0, digits=4), clp, "% CI ", 
                format(x$d0.ci, digits=4), "\n")
        cat("Mediation Effect_1: ", format(x$d1, digits=4), clp, "% CI ", 
                format(x$d1.ci, digits=4), "\n")
        cat("Sample Size Used:", x$nobs,"\n\n") 
    } 
    invisible(x)
}


