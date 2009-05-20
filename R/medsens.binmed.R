

medsens.binmed <- function(z, ...){
	UseMethod("medsens.binmed", z)
	}

medsens.binmed.default <- function(model.m, model.y, T="treat.name", M="med.name", INT=FALSE, DETAIL=TRUE)
{

#Variable labels (raw uppercase letters)
y.t.data <- model.frame(model.y)
Y <- colnames(y.t.data)[1]
if(is.factor(y.t.data[,T])==TRUE){
	cat.c <- levels(y.t.data[,T])[1] 
	cat.t <- levels(y.t.data[,T])[2]
	T.out <- paste(T,cat.t, sep="") 
	} else {
	cat.c <- NULL
	cat.t <- NULL
	T.out <- paste(T,cat.t, sep="")
	}

if(is.factor(y.t.data[,M])==TRUE){
	cat.m0 <- levels(y.t.data[,M])[1] 
	cat.m1 <- levels(y.t.data[,M])[2]
	M.out <- paste(M,cat.m1, sep="") 
	} else {
	cat.m0 <- NULL
	cat.m1 <- NULL
	M.out <- paste(M,cat.m1, sep="")
	}

if(INT==TRUE){
TM <- paste(T,M, sep=":")
TM.out <- paste(T.out,M.out, sep=":")
	}

#Variable values (LABEL.value)
Y.value <- y.t.data[,1]
T.value <- y.t.data[,T]
M.value <- y.t.data[,M]

		
#Function for computing lambdas
lambda <- function(mmodel) {
  mu <- model.matrix(mmodel) %*% coef(mmodel)
  m <- mmodel$y
  return((m*dnorm(-mu)-(1-m)*dnorm(-mu))/(m*pnorm(mu)+(1-m)*pnorm(-mu)))
}

#Sensitivity Parameter
if(DETAIL==TRUE){
	rho <- round(seq(-.99, .99, by = .01),2)
	} else {
	rho <- round(seq(-.90, .90, by = .1),2)
		}

d0 <- matrix(NA, length(rho), 1)
d1 <- matrix(NA, length(rho), 1)
d0.var <- matrix(NA, length(rho), 1)
d1.var <- matrix(NA, length(rho), 1)

#Loop begins here
for(i in 1:length(rho)){
# Refit outcome model
adj <- rho[i] * lambda(model.m)
model.y.adj <- update(model.y, as.formula(paste(". ~ . + adj")))
y.coefs <- model.y.adj$coef

#Save Estimates
m.mat.1 <- model.matrix(model.m)
m.mat.1[,T.out] <- 1
m.mat.0 <- model.matrix(model.m)
m.mat.0[,T.out] <- 0
mu.1 <- m.mat.1 %*% coef(model.m)
mu.0 <- m.mat.0 %*% coef(model.m)
lambda11 <- dnorm(-mu.1) / pnorm(mu.1) #m=1,t=1
lambda10 <- dnorm(-mu.0) / pnorm(mu.0) #m=1,t=0
lambda01 <- -dnorm(-mu.1) / pnorm(-mu.1) #m=0,t=1
lambda00 <- -dnorm(-mu.0) / pnorm(-mu.0) #m=0,t=0

d0[i,] <- mean((y.coefs[M.out] + rho[i] * y.coefs["adj"] * (lambda10 - lambda00)) * (dnorm(mu.1) - dnorm(mu.0)))
if(INT==TRUE){
    d1[i,] <- mean((y.coefs[M.out] + y.coefs[TM.out] + 
        rho[i] * y.coefs["adj"] * (lambda11 - lambda01)) * (dnorm(mu.1) - dnorm(mu.0)))
    } else {
    d1[i,] <- mean((y.coefs[M.out] + 
        rho[i] * y.coefs["adj"] * (lambda11 - lambda01)) * (dnorm(mu.1) - dnorm(mu.0)))
    }

}

#Name Elements - Extract Quantities
m.names <- names(model.m$coef)
y.names <- names(model.y$coef)
b.names <- c(m.names, y.names)
row.names(b.sur) <- b.names
m.coefs <- as.matrix(b.sur[1:m.k])
y.coefs <- as.matrix(b.sur[(m.k+1):(m.k+y.k)])
row.names(m.coefs) <- m.names
row.names(y.coefs) <- y.names
rownames(v.cov) <- b.names
colnames(v.cov) <- b.names
v.m <- v.cov[1:m.k,1:m.k]
v.y <- v.cov[(m.k+1):(m.k+y.k),(m.k+1):(m.k+y.k)]

#Save Estimates
if(INT==TRUE){
	d0[i,] <- m.coefs[paste(T.cat),]*y.coefs[paste(M),] 
	d1[i,] <- m.coefs[paste(T.cat),]*(y.coefs[paste(M),] + y.coefs[paste(int.lab),])
	} else {
		d0[i,] <- m.coefs[paste(T.cat),]*y.coefs[paste(M),]
		d1[i,] <- m.coefs[paste(T.cat),]*y.coefs[paste(M),] 
		}

#Save Variance Estimates
if(INT==TRUE){
	d0.var[i,] <- (y.coefs[paste(M),] + 0*y.coefs[paste(int.lab),])^2*v.m[T.cat,T.cat] + 	m.coefs[paste(T.cat),]^2*(v.y[M,M] + 0*v.y[int.lab, int.lab] + 0*2*v.y[M, int.lab])
	d1.var[i,] <- (y.coefs[paste(M),] + y.coefs[paste(int.lab),])^2*v.m[T.cat,T.cat] + 	m.coefs[paste(T.cat),]^2*(v.y[M,M] + v.y[int.lab, int.lab] + 2*v.y[M, int.lab])
	} else {
	d0.var[i,] <- (m.coefs[paste(T.cat),]^2*v.y[M,M]) + (y.coefs[paste(M),]^2*v.m[T.cat,T.cat])
	d1.var[i,] <- (m.coefs[paste(T.cat),]^2*v.y[M,M]) + (y.coefs[paste(M),]^2*v.m[T.cat,T.cat])
		}
		
rm(b.sur, m.coefs, y.coefs, v.cov, v.m, v.y)

}

if(INT==TRUE){
upper.d0 <- d0 + qnorm(0.975) * sqrt(d0.var)
lower.d0 <- d0 - qnorm(0.975) * sqrt(d0.var)
upper.d1 <- NULL
lower.d1 <- NULL
ind.d0 <- as.numeric(lower.d0 < 0 & upper.d0 > 0)
ind.d1 <- NULL	
	} else {
upper.d0 <- d0 + qnorm(0.975) * sqrt(d0.var)
lower.d0 <- d0 - qnorm(0.975) * sqrt(d0.var)
upper.d1 <- d1 + qnorm(0.975) * sqrt(d1.var)
lower.d1 <- d1 - qnorm(0.975) * sqrt(d1.var)	
ind.d0 <- as.numeric(lower.d0 < 0 & upper.d0 > 0)
ind.d1 <- as.numeric(lower.d1 < 0 & upper.d1 > 0)
		}

out <- list(rho = rho, err.cr=err.cr, d0=d0, d1=d1, upper.d0=upper.d0, lower.d0=lower.d0, upper.d1=upper.d1, lower.d1=lower.d1, ind.d0=ind.d0, ind.d1=ind.d1, INT=INT, DETAIL=DETAIL)
class(out) <- "sens.c"
out
	}

print.sens.c <- function(x, ...){
	print(unlist(x[1:16]))
	invisible(x)
	}

summary.sens.c <- function(object)
	structure(object, class = c("sum.sens.c", class(object)))
 
print.sum.sens.c <- function(x, ...){
	if(x$INT==FALSE){
	tab <- cbind(x$rho, round(x$err.cr,4), round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
	tab <- tab[x$ind.d0==1, -6]
    colnames(tab) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
    rownames(tab) <- NULL
    cat("\nMediation Sensitivity Analysis\n")
    cat("\nSensitivity Region\n\n")
    print(tab)
	invisible(x)	
		} else {
	tab.d0 <- cbind(x$rho, round(x$err.cr,4), round(x$d0,4), round(x$lower.d0,4), round(x$upper.d0, 4), x$ind.d0)
	tab.d0 <- tab.d0[x$ind==1, -6]
    colnames(tab.d0) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
    rownames(tab.d0) <- NULL
    tab.d1 <- cbind(x$rho, round(x$err.cr,4), round(x$d1,4), round(x$lower.d1,4), round(x$upper.d1, 4), x$ind.d1)
	tab.d1 <- tab[x$ind.d1==1, -6]
    colnames(tab.d1) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
    rownames(tab.d1) <- NULL
    cat("\nMediation Sensitivity Analysis\n")
    cat("\nSensitivity Region: d0\n\n")
    print(tab.d0)
    cat("\nSensitivity Region: d1\n\n")
    print(tab.d1)
	invisible(x)		
			}
	}

plot.sens.c <- function(x, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL,...){
	if(x$INT==FALSE){
		plot(x$rho, x$d0, type="n", xlab="", ylab = "", main=main, xlim=xlim, ylim=ylim)
		polygon(x=c(x$rho, rev(x$rho)), y=c(x$lower.d0, rev(x$upper.d0)), border=FALSE, col=8, lty=2)
		lines(x$rho, x$d0, lty=1)
		abline(h=0)
		abline(v=0)
		title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
		title(ylab = expression(paste("Average Mediation Effect: ", bar(delta))), cex.lab=.9)
		if(x$DETAIL==TRUE){
		abline(h=x$d0[96], lty=2)
			} else {
		abline(h=x$d0[10], lty=2)		
				}
		} else {
		par(mfrow=c(1,2))
		plot(x$rho, x$d0, type="n", xlab="", ylab = "", ylim = c(-.2,.2))
		polygon(x=c(x$rho, rev(x$rho)), y=c(x$lower.d0, rev(x$upper.d0)), border=FALSE, col=8, lty=2)
		lines(x$rho, x$d0, lty=1)
		abline(h=0)
		if(x$DETAIL==TRUE){
		abline(h=x$d0[96], lty=2)
			} else {
		abline(h=x$d0[10], lty=2)		
				}
		abline(v=0)
		title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
		title(ylab = expression(paste("Average Mediation Effect: ", bar(delta[0]))), cex.lab=.9)

		#Delta_1
		plot(x$rho, x$d1, type="n", xlab="", ylab = "", ylim = c(-.2,.2))
		polygon(x=c(x$rho, rev(x$rho)), y=c(x$lower.d1, rev(x$upper.d1)), border=FALSE, col=8, lty=2)
		lines(rho, d1, lty=1)
		abline(h=0)
				if(x$DETAIL==TRUE){
		abline(h=x$d1[96], lty=2)
			} else {
		abline(h=x$d1[10], lty=2)		
				}
		abline(v=0)
		title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
		title(ylab = expression(paste("Average Mediation Effect: ", bar(delta[1]))), cex.lab=.9)
			}

	}
