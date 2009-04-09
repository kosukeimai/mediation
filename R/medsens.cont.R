

medsens.cont <- function(z, ...){
	UseMethod("medsens.cont", z)
	}

medsens.cont.default <- function(z, model.y, T="treat.name", M="med.name")
{

y.t.data <- model.frame(model.y)	
if(is.factor(y.t.data[,paste(T)])==TRUE){
	cat.c <- levels(y.t.data[,T])[1] 
	cat.t <- levels(y.t.data[,T])[2]
	T.cat <- paste(T,cat.t, sep="") 
	} else {
	cat.c <- NULL
	cat.t <- NULL
	T.cat <- paste(T,cat.t, sep="")
	}
	
#Estimate Error Correlation
mod.y <- update(model.y,as.formula(paste(". ~ . -", M)))
err.cr <- cor(model.m$resid, mod.y$resid)

#Sensitivity Parameter
rho <- round(seq(-.90, .90, by = .1),2)
med.eff <- matrix(NA, length(rho), 1)
med.var <- matrix(NA, length(rho), 1)

for(i in 1:length(rho)){

e.cor <- rho[i]

b.dif <- 1
eps <- .001 #.Machine$double.eps

#Stacked Equations
m.mat <- model.matrix(model.m)
y.mat <- model.matrix(model.y)
m.k <- ncol(m.mat)
m.n <- nrow(m.mat)
y.k <- ncol(y.mat)
y.n <- nrow(y.mat)
n <- y.n
m.zero <- matrix(0, m.n, y.k)
y.zero <- matrix(0, y.n, m.k)
X.1 <- cbind(m.mat, m.zero)
X.2 <- cbind(y.zero, y.mat)
X <- rbind(X.1, X.2)

m.frame <- model.frame(model.m)
y.frame <- model.frame(model.y)
Y.c <- rbind(as.matrix(m.frame[,1]), as.matrix(y.frame[,1]))

#Estimates of OLS Start Values
inxx <- solve(crossprod(X))
b.ols <- inxx %*% crossprod(X,Y.c)
b.tmp <- b.ols

while(abs(b.dif) > eps){

e.hat <- as.matrix(Y.c - (X %*% b.tmp))

e.1 <- e.hat[1:n]
e.2 <- e.hat[(n+1): (2*n)]

sd.1 <- sd(e.1)
sd.2 <- sd(e.2)

omega <- matrix(NA, 2,2)

omega[1,1] <- crossprod(e.1)/(n-1)
omega[2,2] <- crossprod(e.2)/(n-1) 
omega[2,1] <- e.cor*sd.1*sd.2
omega[1,2] <- e.cor*sd.1*sd.2

I <- diag(1,n)
omega.i <- solve(omega)
v.i <-  kronecker(omega.i, I)

X.sur <- crossprod(X, v.i) %*% X 
b.sur <- solve(X.sur) %*% crossprod(X, v.i) %*% Y.c

#Variance-Covariance Matrix
v.cov <- crossprod(X, v.i) %*% X
v.cov <- solve(v.cov)

b.old <- b.tmp
b.dif <- sum((b.sur - b.old)^2)
b.tmp <- b.sur

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

med.eff[i,] <- m.coefs[paste(T.cat),]*y.coefs[paste(M),] 
med.var[i,] <- (m.coefs[paste(T.cat),]^2*v.y[M,M]) + (y.coefs[paste(M),]^2*v.m[T.cat,T.cat])
rm(b.sur, m.coefs, y.coefs, v.cov, v.m, v.y)

}

upper <- med.eff + qnorm(0.975) * sqrt(med.var)
lower <- med.eff - qnorm(0.975) * sqrt(med.var)
ind <- as.numeric(lower < 0 & upper > 0)

out <- list(rho = rho, err.cr=err.cr, delta=med.eff, d.var=med.var, lo=lower, up=upper, ind=ind)
class(out) <- "sens.c"
out
	}

print.sens.c <- function(x, ...){
	print(unlist(x[1:6]))
	invisible(x)
	}

summary.sens.c <- function(object)
	structure(object, class = c("sum.sens.c", class(object)))
 
print.sum.sens.c <- function(x, ...){
	tab <- cbind(x$rho, round(x$err.cr,4), round(x$delta,4), round(x$lo,4), round(x$up, 4), x$ind)
	tab <- tab[x$ind==1, -6]
    colnames(tab) <-  c("Rho","Error Cor.", "Med. Eff.", "95% CI Lower", "95% CI Upper")
    rownames(tab) <- NULL
    cat("\nMediation Sensitivity Analysis\n")
    cat("\nSensitivity Region\n\n")
    print(tab)
	invisible(x)
	}

plot.sens.c <- function(x, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, main=NULL, ...){
	plot(x$rho, x$delta, type="n", xlab="", ylab = "", main=main, xlim=xlim, ylim=ylim)
	polygon(x=c(x$rho, rev(x$rho)), y=c(x$lo, rev(x$up)), border=FALSE, col=8, lty=2)
	lines(x$rho, x$delta, lty=1)
	abline(h=0)
	abline(h=x$delta[96], lty=2)
	abline(v=0)
	title(xlab=expression(paste("Sensitivity Parameter: ", rho)), line=2.5, cex.lab=.9)
	title(ylab = expression(paste("Average Mediation Effect: ", bar(delta))), cex.lab=.9)
	}