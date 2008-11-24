med.nonpara <- function(Y,M,T){
	
samp <- data.frame(Y,M,T)
samp <- na.omit(samp)
ftab <- ftable(samp$M, samp$T)
test.1 <- unique(c(ftab[,1]==0, ftab[,2]==0))
test.2 <- unique(c(ftab[,1]==1, ftab[,2]==1))

if(length(test.1)==1){test.1 <- c(test.1, test.1)}
if(length(test.2)==1){test.2 <- c(test.2, test.2)}


if(test.1[1]==FALSE & test.1[2]==FALSE & test.2[1]==FALSE & test.2[2]==FALSE){
	
m.cat <- sort(unique(samp$M))
n <- length(samp$Y)
n1 <- sum(samp$T)
n0 <- sum(1-samp$T)
################################################################
# Point Estimates
################################################################

    delta_0m <- matrix(NA, length(m.cat), 1)
    for(i in 1:length(m.cat)){
    	delta_0m[i] <- (sum(samp$T[samp$M==m.cat[i]]) * sum(samp$Y[samp$M==m.cat[i] & samp$T==0]))/(n1*(sum((1 - samp$T[samp$M==m.cat[i]]))))
    }
    
delta_0.est <- sum(delta_0m) - mean(samp$Y[samp$T==0])

    delta_1m <- matrix(NA, length(m.cat), 1)
    
    for(i in 1:length(m.cat)){
    	delta_1m[i] <-  (sum((1 - samp$T[samp$M==m.cat[i]])) * sum(samp$Y[samp$M==m.cat[i] & samp$T==1]))/(n0*(sum(samp$T[samp$M==m.cat[i]])))    
    }
    
delta_1.est <- mean(samp$Y[samp$T==1]) - sum(delta_1m)

    
##########################################################
# Variance Estimates
#########################################################

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


#Print Results
est <- c(delta_0.est, delta_1.est)
se <- c(sqrt(var.delta_0), sqrt(var.delta_1))
zscore <- est/se
pz <-  2*pnorm(-abs(est/se))

coef.table <- cbind(est, se, zscore, pz)
dimnames(coef.table) <- list(c("Delta_0", "Delta_1"), c("Estimate", "Std. Error", "z-value", "Pr(>|z|)"))
print(coef.table)
cat("Delta_0 Represents the Average Mediation Effect Under Control \n")
cat("Delta_1 Represents the Average Mediation Effect Under Treatment \n")		
		}

else if(test.1[1]==TRUE){print("Error: Empty Cells in Data Under Control")
} else if(test.1[2]==TRUE){print("Error: Empty Cells in Data Under Treatment")
} else if(test.2[1]==TRUE){print("Error Some Cells Have No Variation Under Control")
} else {print("Error Some Cells Have No Variation Under Treatment")} 
	
}
