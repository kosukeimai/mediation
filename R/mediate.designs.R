## User interface functions

#single experiment design
mediate.sed <- function(outcome, mediator, treat, data,
                        SI = FALSE, sims = 1000, conf.level = 0.95, boot = FALSE) {

    varnames <- c(outcome, mediator, treat)
    data <- data[,varnames]
    if(sum(is.na(data))){
        warning("NA's in data. Observations with missing data removed.")
        data <- na.omit(data)
    }
    data <- sapply(data, function(x) if(is.factor(x)) as.numeric(x)-1 else as.numeric(as.character(x)))

    if(SI) {
        out <- mediate.np(data[,outcome], data[,mediator], data[,treat],
                           sims, conf.level, boot)
        out$design <- "SED.NP.SI"
    } else {
		check <- apply(data, 2, function(x) identical(sort(unique(x)), c(0,1)))
        if(sum(check) != ncol(data)) {
	        stop("Invalid values in variables.")
        }

        out <- mechanism.bounds(data[,outcome], data[,mediator], data[,treat],
                                NULL, design = "SED")
        out$design <- "SED.NP.NOSI"
    }
    out
}


#parallel design
mediate.pd <- function(outcome, mediator, treat, manipulated, data,
                        NINT = TRUE, sims = 1000, conf.level = 0.95) {

    varnames <- c(outcome, mediator, treat, manipulated)
    data <- data[,varnames]
    if(sum(is.na(data))){
        warning("NA's in data. Observations with missing data removed.")
        data <- na.omit(data)
    }
    data <- sapply(data, function(x) as.numeric(as.character(x)))
    
    check <- apply(data, 2, function(x) identical(sort(unique(x)), c(0,1)))
    if(sum(check) != ncol(data)) {
        stop("Invalid values in variables.")
    }

    if(NINT) {
       out <- boot.pd(data[,outcome], data[,mediator], data[,treat],
                       data[,manipulated], sims, conf.level)
    } else {
       out <- mechanism.bounds(data[,outcome], data[,mediator], data[,treat],
                              data[,manipulated], design = "PD")
    }
    out
}


#parallel encouragement design
mediate.ped <- function(outcome, mediator, treat, encourage, data) {

    varnames <- c(outcome, mediator, treat, encourage)
    data <- data[,varnames]
    if(sum(is.na(data))){
        warning("NA's in data. Observations with missing data removed.")
        data <- na.omit(data)
    }
    data <- sapply(data, function(x) as.numeric(as.character(x)))
    		
	check.1 <- apply(data[,1:3], 2, function(x) identical(sort(unique(x)), c(0,1)))
	check.2 <- identical(sort(unique(data[,4])), c(-1,0,1))
    if(sum(check.1, check.2) != ncol(data)) {
        stop("Invalid values in variables.")
    }

    out <- mechanism.bounds(data[,outcome], data[,mediator], data[,treat],
                             data[,encourage], design = "PED")
    out
}


#crossover encouragement design
mediate.ced <- function(outcome, med.1, med.2, treat, encourage, data,
                        sims = 1000, conf.level = .95){

    varnames <- c(outcome, treat, med.1, med.2, encourage)
    data <- data[,varnames]
    if(sum(is.na(data))){
        warning("NA's in data. Observations with missing data removed.")
        data <- na.omit(data)
    }
    data <- sapply(data, function(x) as.numeric(as.character(x)))
    data <- as.data.frame(data)
    		
	check <- apply(data, 2, function(x) identical(sort(unique(x)), c(0,1)))
    if(sum(check) != ncol(data)) {
        stop("Invalid values in variables.")
    }

    names(data) <- c("Y2","T","M1","M2","V")
    n <- nrow(data)

    # Storage
        d.p.c <- d.p.t <- matrix(NA, sims, 1)

    # Bootstrap function.
    for(b in 1:(sims+1)){
        index <- sample(1:n, n, replace = TRUE)
        if(b == sims + 1){
        	index <- 1:n
        }
        d <- data[index,]

        #Atmv
        A111 <- sum(d$M2==1 & d$T==1 & d$M1==1 & d$V==1)/sum(d$T==1 & d$M1==1 & d$V==1)
        A011 <- sum(d$M2==1 & d$T==0 & d$M1==1 & d$V==1)/sum(d$T==0 & d$M1==1 & d$V==1)
        A101 <- sum(d$M2==1 & d$T==1 & d$M1==0 & d$V==1)/sum(d$T==1 & d$M1==0 & d$V==1)
        A001 <- sum(d$M2==1 & d$T==0 & d$M1==0 & d$V==1)/sum(d$T==0 & d$M1==0 & d$V==1)

        A110 <- sum(d$M2==1 & d$T==1 & d$M1==1 & d$V==0)/sum(d$T==1 & d$M1==1 & d$V==0)
        A010 <- sum(d$M2==1 & d$T==0 & d$M1==1 & d$V==0)/sum(d$T==0 & d$M1==1 & d$V==0)
        A100 <- sum(d$M2==1 & d$T==1 & d$M1==0 & d$V==0)/sum(d$T==1 & d$M1==0 & d$V==0)
        A000 <- sum(d$M2==1 & d$T==0 & d$M1==0 & d$V==0)/sum(d$T==0 & d$M1==0 & d$V==0)

        #Gtm1m2
        G1111 <- mean(d$Y2[d$T==1 & d$M1==1 & d$M2==1 & d$V==1])
        G0111 <- mean(d$Y2[d$T==0 & d$M1==1 & d$M2==1 & d$V==1])
        G1011 <- mean(d$Y2[d$T==1 & d$M1==0 & d$M2==1 & d$V==1])
        G0011 <- mean(d$Y2[d$T==0 & d$M1==0 & d$M2==1 & d$V==1])

        G1101 <- mean(d$Y2[d$T==1 & d$M1==1 & d$M2==0 & d$V==1])
        G0101 <- mean(d$Y2[d$T==0 & d$M1==1 & d$M2==0 & d$V==1])
        G1001 <- mean(d$Y2[d$T==1 & d$M1==0 & d$M2==0 & d$V==1])
        G0001 <- mean(d$Y2[d$T==0 & d$M1==0 & d$M2==0 & d$V==1])

        G1110 <- mean(d$Y2[d$T==1 & d$M1==1 & d$M2==1 & d$V==0])
        G0110 <- mean(d$Y2[d$T==0 & d$M1==1 & d$M2==1 & d$V==0])
        G1010 <- mean(d$Y2[d$T==1 & d$M1==0 & d$M2==1 & d$V==0])
        G0010 <- mean(d$Y2[d$T==0 & d$M1==0 & d$M2==1 & d$V==0])

        G1100 <- mean(d$Y2[d$T==1 & d$M1==1 & d$M2==0 & d$V==0])
        G0100 <- mean(d$Y2[d$T==0 & d$M1==1 & d$M2==0 & d$V==0])
        G1000 <- mean(d$Y2[d$T==1 & d$M1==0 & d$M2==0 & d$V==0])
        G0000 <- mean(d$Y2[d$T==0 & d$M1==0 & d$M2==0 & d$V==0])

        Q11 <- sum(d$T==1 & d$M1==1)/sum(d$T==1)
        Q10 <- sum(d$T==1 & d$M1==0)/sum(d$T==1)
        Q01 <- sum(d$T==0 & d$M1==1)/sum(d$T==0)
        Q00 <- sum(d$T==0 & d$M1==0)/sum(d$T==0)
        B11 <- Q11/((A100 - A101)*Q10 + (A111 - A110)*Q11)
        B10 <- Q10/((A100 - A101)*Q10 + (A111 - A110)*Q11)
        B01 <- Q01/((A000 - A001)*Q00 + (A011 - A010)*Q01)
        B00 <- Q00/((A000 - A001)*Q00 + (A011 - A010)*Q01)

		if(b == sims + 1){
            d.p.c.mu <- (-A110*G1110 - (1-A110)*G1100 + A111*G1111 + (1-A111)*G1101) * B11 +
                            (-A100*G1010 - (1-A100)*G1000 + A101*G1011 + (1-A101)*G1001) * B10
            d.p.t.mu <- (A000*G0010 + (1-A000)*G0000 - A001*G0011 - (1-A001)*G0001) * B00 +
                            (A010*G0110 + (1-A010)*G0100 - A011*G0111 - (1-A011)*G0101) * B01
		} else {
            d.p.c[b] <- (-A110*G1110 - (1-A110)*G1100 + A111*G1111 + (1-A111)*G1101) * B11 +
                            (-A100*G1010 - (1-A100)*G1000 + A101*G1011 + (1-A101)*G1001) * B10
            d.p.t[b] <- (A000*G0010 + (1-A000)*G0000 - A001*G0011 - (1-A001)*G0001) * B00 +
                            (A010*G0110 + (1-A010)*G0100 - A011*G0111 - (1-A011)*G0101) * B01
		}
    }#bootstraploop
    
    if(is.nan(d.p.c.mu) | is.nan(d.p.t.mu)){
    	warning("NaN produced; distribution of observed variables may be too sparse")
    }

    d.p.c[d.p.c==-Inf] <- NA
    d.p.t[d.p.t==-Inf] <- NA
    d.p.c[d.p.c==Inf] <- NA
    d.p.t[d.p.t==Inf] <- NA

    low <- (1 - conf.level)/2
    high <- 1 - low

    d.p.c.ci <- quantile(d.p.c, c(low,high), na.rm=TRUE)
    d.p.t.ci <- quantile(d.p.t, c(low,high), na.rm=TRUE)

    out <- list(d0 = d.p.c.mu, d1 = d.p.t.mu, d0.ci = d.p.c.ci, d1.ci = d.p.t.ci,
                nobs = n, sims = sims, conf.level = conf.level, design = "CED")

    class(out) <- "mediate.design"
    out
}


## Internal functions

#nonparametric under SI assumption
mediate.np <- function(Y, M, T, sims, conf.level, boot){
#    samp <- data.frame(na.omit(cbind(M,Y,T)))
    samp <- data.frame(Y=Y, M=M, T=T)
    m.cat <- sort(unique(samp$M))
    n <- length(samp$Y)
    n1 <- sum(samp$T)
    n0 <- sum(1-samp$T)
    t <- (1 - conf.level)/2

    if(boot){
        idx <- seq(1,n,1)
        d0.bs <- matrix(NA,sims,1)
        d1.bs <- matrix(NA,sims,1)
        for(b in 1:(sims+1)){
            resample <- sample(idx,n,replace=TRUE)
            if(b == sims + 1){
            	resample <- 1:n
            }
            samp.star <- samp[resample,]

            delta_0m <- delta_1m <- matrix(NA, length(m.cat), 1)
            for(i in 1:length(m.cat)){
                delta_0m[i] <- (sum(samp.star$T[samp.star$M==m.cat[i]]) *

                                sum(samp.star$Y[samp.star$M==m.cat[i] & samp.star$T==0])) /
                                (n1*(sum((1 - samp.star$T[samp.star$M==m.cat[i]]))))
                delta_1m[i] <- (sum((1 - samp.star$T[samp.star$M==m.cat[i]])) *
                            sum(samp.star$Y[samp.star$M==m.cat[i] & samp.star$T==1])) /
                            (n0*(sum(samp.star$T[samp.star$M==m.cat[i]])))
            }
            if(b == sims + 1){
	            d0 <- sum(delta_0m) - mean(samp.star$Y[samp.star$T==0])
    	        d1 <- mean(samp.star$Y[samp.star$T==1]) - sum(delta_1m)
            } else {
	            d0.bs[b,] <- sum(delta_0m) - mean(samp.star$Y[samp.star$T==0])
    	        d1.bs[b,] <- mean(samp.star$Y[samp.star$T==1]) - sum(delta_1m)
            }
        }

        low <- (1 - conf.level)/2
        high <- 1 - low
        d0.ci <- quantile(d0.bs,c(low,high), na.rm=TRUE)
        d1.ci <- quantile(d1.bs,c(low,high), na.rm=TRUE)

    } else {

        delta_0m <- delta_1m <- matrix(NA, length(m.cat), 1)

        for(i in 1:length(m.cat)){
            delta_0m[i] <- (sum(samp$T[samp$M==m.cat[i]]) *
                                sum(samp$Y[samp$M==m.cat[i] & samp$T==0])) /
                                (n1*(sum((1 - samp$T[samp$M==m.cat[i]]))))
            delta_1m[i] <- (sum((1 - samp$T[samp$M==m.cat[i]])) *
                                sum(samp$Y[samp$M==m.cat[i] & samp$T==1])) /
                                (n0*(sum(samp$T[samp$M==m.cat[i]])))
        }

        d0 <- sum(delta_0m) - mean(samp$Y[samp$T==0])
        d1 <- mean(samp$Y[samp$T==1]) - sum(delta_1m)

        pr.t.1 <- sum(samp$T)/n
        pr.t.0 <- sum(1-samp$T)/n

        lambda.1m <- lambda.0m <- matrix(NA, length(m.cat), 1)
        for(i in 1:length(m.cat)){
            lambda.1m[i] <- (sum(length(samp$M[samp$M==m.cat[i] & samp$T==1])))/n /pr.t.1
            lambda.0m[i] <- (sum(length(samp$M[samp$M==m.cat[i] & samp$T==0])))/n /pr.t.0
        }

        mu.1m <- mu.0m <- matrix(NA, length(m.cat), 1)
        for(i in 1:length(m.cat)){
            mu.1m[i] <- mean(samp$Y[samp$M==m.cat[i] & samp$T==1])
            mu.0m[i] <- mean(samp$Y[samp$M==m.cat[i] & samp$T==0])
        }

        # Variance of Delta_0
        var.delta.0m <- matrix(NA, length(m.cat), 1)
        for(i in 1:length(m.cat)){
            var.delta.0m[i] <-  lambda.1m[i]*(((lambda.1m[i]/lambda.0m[i]) - 2) *
                                var(samp$Y[samp$M==m.cat[i] & samp$T==0]) +
                                ((n0*(1-lambda.1m[i])*mu.0m[i]^2)/n1))
        }

        m.leng <- length(m.cat)
        mterm.d0 <- c()
        for(i in 1:m.leng-1){
            mterm.d0 <- c(mterm.d0, lambda.1m[i] * lambda.1m[(i+1):m.leng] *
                            mu.0m[i] * mu.0m[(i+1):m.leng])
        }
        var.delta_0 <- (1/n0) * sum(var.delta.0m) - (2/n1) * sum(mterm.d0) +
                        var(samp$Y[samp$T == 0])/n0

        # Variance of Delta_1
        var.delta.1m <- matrix(NA, length(m.cat), 1)
        for(i in 1:length(m.cat)){
            var.delta.1m[i] <-  lambda.0m[i]*(((lambda.0m[i]/lambda.1m[i]) - 2) *
                                var(samp$Y[samp$M==m.cat[i] & samp$T==1]) +
                                ((n1*(1-lambda.0m[i])*mu.1m[i]^2)/n0))
        }

        mterm.d1 <- c()
        for(i in 1:m.leng-1){
            mterm.d1 <- c(mterm.d1, lambda.0m[i] * lambda.0m[(i+1):m.leng] *
                            mu.1m[i] * mu.1m[(i+1):m.leng])
        }
        var.delta_1 <- (1/n1) * sum(var.delta.1m) - (2/n0)* sum(mterm.d1) +
                        var(samp$Y[samp$T == 1])/n1

        d0.ci <- c(d0 - (-qnorm(t)*sqrt(var.delta_0)), d0 + (-qnorm(t)*sqrt(var.delta_0)))
        d1.ci <- c(d1 - (-qnorm(t)*sqrt(var.delta_1)), d1 + (-qnorm(t)*sqrt(var.delta_1)))
    }

    out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
                conf.level=conf.level, nobs=n, boot=boot)
    class(out) <- "mediate.design"
    out

}


#parallel design under no interaction assumption
boot.pd <- function(outcome, mediator, treatment, manipulated,
                    sims, conf.level) {
    n <- length(outcome)
    data <- matrix(, nrow = n, ncol = 4)
    data[,1] <- outcome
    data[,2] <- treatment
    data[,3] <- mediator
    data[,4] <- manipulated
    data <- as.data.frame(data)
    names(data) <- c("Y","T","M","D")

    acme <- matrix(NA, sims, 1)
    for(b in 1:(sims+1)){  # bootstrap
        index <- sample(1:n, n, replace = TRUE)
        if(b == sims + 1){
        	index <- 1:n
        }
        d <- data[index,]
        d0.temp <- mean(d$Y[d$T==1 & d$D==0]) - mean(d$Y[d$T==0 & d$D==0])
        weight.m1 <- sum(d$M==1 & d$D==1 )/sum(d$D==1)
        weight.m0 <- sum(d$M==0 & d$D==1 )/sum(d$D==1)
        m1 <- mean(d$Y[d$T==1 & d$M==1 & d$D==1]) - mean(d$Y[d$T==0 & d$M==1 & d$D==1])
        m0 <- mean(d$Y[d$T==1 & d$M==0 & d$D==1]) - mean(d$Y[d$T==0 & d$M==0 & d$D==1])
        if(b == sims + 1){
	        acme.mu <- weight.m1*m1 + weight.m0*m0
        } else {
	        acme[b] <- weight.m1*m1 + weight.m0*m0
    	}
    }

    acme[acme==-Inf] <- NA
    acme[acme==Inf] <- NA

    low <- (1 - conf.level)/2
    high <- 1 - low
    acme.ci <- quantile(acme, c(low,high), na.rm=TRUE)
    out <- list(d0 = acme.mu, d1 = acme.mu, d0.ci = acme.ci, d1.ci = acme.ci,
                nobs = n, conf.level = conf.level, sims = sims, design = "PD.NINT")
    class(out) <- "mediate.design"
    out

}



#Bounds Function for parallel, parallel encouragement, and single experiment designs
mechanism.bounds <- function(outcome, mediator, treatment, DorZ, design) {
    d <- matrix(, nrow=length(outcome), ncol=3)

    d[,1] <- outcome
    d[,2] <- treatment
    d[,3] <- mediator
    d <- as.data.frame(d)
    names(d) <- c("Y","T","M")
    nobs <- nrow(d)

    ## "DorZ" indicates the variable that assigns manipulation direction.
    ## Z in the JRSSA paper. This is left to null for the single experiment design.

    # Single Experiment
    if(design=="SED") {
        d$Z <- 0
        P111 <- sum(d$Y==1 & d$M==1 & d$Z==0 & d$T==1)/sum(d$T==1 & d$Z==0)
        P101 <- sum(d$Y==1 & d$M==0 & d$Z==0 & d$T==1)/sum(d$T==1 & d$Z==0)
        P001 <- sum(d$Y==0 & d$M==0 & d$Z==0 & d$T==1)/sum(d$T==1 & d$Z==0)
        P011 <- sum(d$Y==0 & d$M==1 & d$Z==0 & d$T==1)/sum(d$T==1 & d$Z==0)

        P110 <- sum(d$Y==1 & d$M==1 & d$Z==0 & d$T==0)/sum(d$T==0 & d$Z==0)
        P100 <- sum(d$Y==1 & d$M==0 & d$Z==0 & d$T==0)/sum(d$T==0 & d$Z==0)
        P000 <- sum(d$Y==0 & d$M==0 & d$Z==0 & d$T==0)/sum(d$T==0 & d$Z==0)
        P010 <- sum(d$Y==0 & d$M==1 & d$Z==0 & d$T==0)/sum(d$T==0 & d$Z==0)

        l.delta.t <- max(-P001-P011, -P000-P001-P100, -P011-P010-P110)
        u.delta.t <- min(P101+P111, P000+P100+P101, P010+P110+P111)
        l.delta.c <- max(-P100-P110, -P001-P101-P100, -P011-P111-P110)
        u.delta.c <- min(P000+P010, P011+P111+P010, P000+P001+P101)
        S.t <- cbind(l.delta.t,u.delta.t)
        S.c <- cbind(l.delta.c,u.delta.c)

    # Parallel
    } else if(design=="PD") {
        ## D indicates if there was manipulation
        d$D <- DorZ
        Z111 <- sum(d$Y==1 & d$T==1 & d$M==1 & d$D==1)/sum(d$T==1 & d$M==1 & d$D==1)
        Z011 <- sum(d$Y==0 & d$T==1 & d$M==1 & d$D==1)/sum(d$T==1 & d$M==1 & d$D==1)
        Z001 <- sum(d$Y==0 & d$T==1 & d$M==0 & d$D==1)/sum(d$T==1 & d$M==0 & d$D==1)
        Z110 <- sum(d$Y==1 & d$T==0 & d$M==1 & d$D==1)/sum(d$T==0 & d$M==1 & d$D==1)
        Z010 <- sum(d$Y==0 & d$T==0 & d$M==1 & d$D==1)/sum(d$T==0 & d$M==1 & d$D==1)
        Z000 <- sum(d$Y==0 & d$T==0 & d$M==0 & d$D==1)/sum(d$T==0 & d$M==0 & d$D==1)
        Z101 <- sum(d$Y==1 & d$T==1 & d$M==0 & d$D==1)/sum(d$T==1 & d$M==0 & d$D==1)
        Z100 <- sum(d$Y==1 & d$T==0 & d$M==0 & d$D==1)/sum(d$T==0 & d$M==0 & d$D==1)

        P111 <- sum(d$Y==1 & d$M==1 & d$T==1 & d$D==0)/sum(d$T==1 & d$D==0)
        P101 <- sum(d$Y==1 & d$M==0 & d$T==1 & d$D==0)/sum(d$T==1 & d$D==0)
        P001 <- sum(d$Y==0 & d$M==0 & d$T==1 & d$D==0)/sum(d$T==1 & d$D==0)
        P011 <- sum(d$Y==0 & d$M==1 & d$T==1 & d$D==0)/sum(d$T==1 & d$D==0)
        P110 <- sum(d$Y==1 & d$M==1 & d$T==0 & d$D==0)/sum(d$T==0 & d$D==0)
        P100 <- sum(d$Y==1 & d$M==0 & d$T==0 & d$D==0)/sum(d$T==0 & d$D==0)
        P000 <- sum(d$Y==0 & d$M==0 & d$T==0 & d$D==0)/sum(d$T==0 & d$D==0)
        P010 <- sum(d$Y==0 & d$M==1 & d$T==0 & d$D==0)/sum(d$T==0 & d$D==0)

        P.l.delta.t <- max(-P001-P011, -P011-P010-P110-P001+Z001, -P000-P001-P100-P011+Z011, -P001-P011+Z001-Z111, -P001+P101-Z101, -P011+P111-Z111)
        P.u.delta.t <- min(P101+P111, P010+P110+P101+P111-Z101, P000+P100+P101+P111-Z111, P101+P111+Z001-Z111, P111-P011+Z011, P101-P001+Z001)
        P.l.delta.c <- max(-P100-P110, -P011-P111-P110-P100+Z000, -P001-P101-P100-P110+Z110,-P100-P110+Z100-Z010, -P110+P010-Z010, -P100+P000-Z100)
        P.u.delta.c <- min(P000+P010, P011+P111+P010+P000-Z100, P000+P001+P101+P010-Z010, P000+P010+Z100-Z010, P010-P110+Z110, P000-P100+Z100)
        P.t <- cbind(P.l.delta.t,P.u.delta.t)
        P.c <- cbind(P.l.delta.c,P.u.delta.c)

    #Parallel Encouragement
    } else if(design=="PED") {
        d$Z <- DorZ
        dP000 <- sum(d$Y==0 & d$M==0 & d$T==0 & d$Z==-1)/sum(d$T==0 & d$Z==-1)
        dP001 <- sum(d$Y==0 & d$M==0 & d$T==1 & d$Z==-1)/sum(d$T==1 & d$Z==-1)
        dP010 <- sum(d$Y==0 & d$M==1 & d$T==0 & d$Z==-1)/sum(d$T==0 & d$Z==-1)
        dP011 <- sum(d$Y==0 & d$M==1 & d$T==1 & d$Z==-1)/sum(d$T==1 & d$Z==-1)
        dP100 <- sum(d$Y==1 & d$M==0 & d$T==0 & d$Z==-1)/sum(d$T==0 & d$Z==-1)
        dP101 <- sum(d$Y==1 & d$M==0 & d$T==1 & d$Z==-1)/sum(d$T==1 & d$Z==-1)
        dP110 <- sum(d$Y==1 & d$M==1 & d$T==0 & d$Z==-1)/sum(d$T==0 & d$Z==-1)
        dP111 <- sum(d$Y==1 & d$M==1 & d$T==1 & d$Z==-1)/sum(d$T==1 & d$Z==-1)

        tP000 <- sum(d$Y==0 & d$M==0 & d$T==0 & d$Z==0)/sum(d$T==0 & d$Z==0)
        tP001 <- sum(d$Y==0 & d$M==0 & d$T==1 & d$Z==0)/sum(d$T==1 & d$Z==0)
        tP010 <- sum(d$Y==0 & d$M==1 & d$T==0 & d$Z==0)/sum(d$T==0 & d$Z==0)
        tP011 <- sum(d$Y==0 & d$M==1 & d$T==1 & d$Z==0)/sum(d$T==1 & d$Z==0)
        tP100 <- sum(d$Y==1 & d$M==0 & d$T==0 & d$Z==0)/sum(d$T==0 & d$Z==0)
        tP101 <- sum(d$Y==1 & d$M==0 & d$T==1 & d$Z==0)/sum(d$T==1 & d$Z==0)
        tP110 <- sum(d$Y==1 & d$M==1 & d$T==0 & d$Z==0)/sum(d$T==0 & d$Z==0)
        tP111 <- sum(d$Y==1 & d$M==1 & d$T==1 & d$Z==0)/sum(d$T==1 & d$Z==0)

        sP000 <- sum(d$Y==0 & d$M==0 & d$T==0 & d$Z==1)/sum(d$T==0 & d$Z==1)
        sP001 <- sum(d$Y==0 & d$M==0 & d$T==1 & d$Z==1)/sum(d$T==1 & d$Z==1)
        sP010 <- sum(d$Y==0 & d$M==1 & d$T==0 & d$Z==1)/sum(d$T==0 & d$Z==1)
        sP011 <- sum(d$Y==0 & d$M==1 & d$T==1 & d$Z==1)/sum(d$T==1 & d$Z==1)
        sP100 <- sum(d$Y==1 & d$M==0 & d$T==0 & d$Z==1)/sum(d$T==0 & d$Z==1)
        sP101 <- sum(d$Y==1 & d$M==0 & d$T==1 & d$Z==1)/sum(d$T==1 & d$Z==1)
        sP110 <- sum(d$Y==1 & d$M==1 & d$T==0 & d$Z==1)/sum(d$T==0 & d$Z==1)
        sP111 <- sum(d$Y==1 & d$M==1 & d$T==1 & d$Z==1)/sum(d$T==1 & d$Z==1)


        P0 <- c(dP000,dP010,dP100,dP110,tP000,tP010,tP100,tP110,sP000,sP010,sP100,sP110)
        P1 <- c(dP001,dP011,dP101,dP111,tP001,tP011,tP101,tP111,sP001,sP011,sP101,sP111)
        Q0 <- c(dP010+dP110,tP010+tP110,sP010+sP110)
        Q1 <- c(dP011+dP111,tP011+tP111,sP011+sP111)

        #Linear programming functions
        bounds.bidir.cmpl <- function(P,Q,dir){
            f.obj <- c(rep(0,16), 0,0,1,0, 0,0,1,0, 0,-1,0,0, 0,-1,0,0,
                    0,0,-1,0, 0,0,-1,0, 0,1,0,0, 0,1,0,0, rep(0,16))
            f.con <- matrix(c(
                    rep(c(1,1,1,0),8), rep(0,32), #dP00t
                    rep(c(rep(c(0,0,0,1),4), rep(0,16)),2), #dP01t
                    rep(0,32), rep(c(1,1,1,0),8), #dP10t
                    rep(c(rep(0,16), rep(c(0,0,0,1),4)),2), #dP11t
                    rep(c(1,1,0,0),8), rep(0,32), #tP00t
                    rep(c(rep(c(0,0,1,1),4), rep(0,16)),2), #tP01t
                    rep(0,32), rep(c(1,1,0,0),8), #tP10t
                    rep(c(rep(0,16), rep(c(0,0,1,1),4)),2), #tP11t
                    rep(c(1,0,0,0),8), rep(0,32), #sP00t
                    rep(c(rep(c(0,1,1,1),4), rep(0,16)),2), #sP01t
                    rep(0,32), rep(c(1,0,0,0),8), #sP10t
                    rep(c(rep(0,16), rep(c(0,1,1,1),4)),2), #sP11t
                    rep(c(rep(0,12), rep(1,4)),4), #dQ
                    rep(c(rep(0,8), rep(1,8)),4), #tQ
                    rep(c(rep(0,4), rep(1,12)),4), #sQ
                    rep(1,64) #sum=1
                    ), nrow=16, byrow=T)
            f.dir <- rep("=", 16)
            f.rhs <- c(P,Q,1)
            r <- lp(dir, f.obj, f.con, f.dir, f.rhs)
            r$objval
        }

        bounds.bidir.popl <- function(P,Q,dir){
            f.obj <- c(rep(0,16), 0,0,1,1, 0,0,1,1, -1,-1,0,0, -1,-1,0,0,
                    0,0,-1,-1, 0,0,-1,-1, 1,1,0,0, 1,1,0,0, rep(0,16))
            f.con <- matrix(c(
                    rep(c(1,1,1,0),8), rep(0,32), #dP00t
                    rep(c(rep(c(0,0,0,1),4), rep(0,16)),2), #dP01t
                    rep(0,32), rep(c(1,1,1,0),8), #dP10t
                    rep(c(rep(0,16), rep(c(0,0,0,1),4)),2), #dP11t
                    rep(c(1,1,0,0),8), rep(0,32), #tP00t
                    rep(c(rep(c(0,0,1,1),4), rep(0,16)),2), #tP01t
                    rep(0,32), rep(c(1,1,0,0),8), #tP10t
                    rep(c(rep(0,16), rep(c(0,0,1,1),4)),2), #tP11t
                    rep(c(1,0,0,0),8), rep(0,32), #sP00t
                    rep(c(rep(c(0,1,1,1),4), rep(0,16)),2), #sP01t
                    rep(0,32), rep(c(1,0,0,0),8), #sP10t
                    rep(c(rep(0,16), rep(c(0,1,1,1),4)),2), #sP11t
                    rep(c(rep(0,12), rep(1,4)),4), #dQ
                    rep(c(rep(0,8), rep(1,8)),4), #tQ
                    rep(c(rep(0,4), rep(1,12)),4), #sQ
                    rep(1,64) #sum=1
                    ), nrow=16, byrow=T)
            f.dir <- rep("=", 16)
            f.rhs <- c(P,Q,1)
            r <- lp(dir, f.obj, f.con, f.dir, f.rhs)
            r$objval
        }

        BE.t <- c(bounds.bidir.popl(P1, Q0, "min"), bounds.bidir.popl(P1, Q0, "max"))
        BE.c <- c(-bounds.bidir.popl(P0 ,Q1, "max"), -bounds.bidir.popl(P0, Q1, "min"))
        num.BET.t.lo <- bounds.bidir.cmpl(P1, Q0, "min")
        num.BET.t.up <- bounds.bidir.cmpl(P1, Q0, "max")
        num.BET.c.lo <- -bounds.bidir.cmpl(P0, Q1, "max")
        num.BET.c.up <- -bounds.bidir.cmpl(P0, Q1, "min")
        denom.BET.t <- P1[1] + P1[3] - P1[9] - P1[11]
        denom.BET.c <- P0[1] + P0[3] - P0[9] - P0[11]
        BET.t <- c(num.BET.t.lo/denom.BET.t, num.BET.t.up/denom.BET.t)
        BET.c <- c(num.BET.c.lo/denom.BET.c, num.BET.c.up/denom.BET.c)

    }#end computation section

    # Output
    if(design=="SED"){
        out <- list(d1 = S.t, d0 = S.c, design=design, nobs=nobs)
    } else if(design=="PD"){
        out <- list(d1 = P.t, d0 = P.c, design=design, nobs=nobs)
    } else if(design=="PED") {
        out <- list(d1=BE.t, d0=BE.c,  d1.p=BET.t, d0.p=BET.c, design=design, nobs=nobs)
    }

    rm(d)
    class(out) <- "mediate.design"
    out
}



## Summary methods

summary.mediate.design <- function(object, ...){
    structure(object, class = c("summary.mediate.design", class(object)))
}



print.summary.mediate.design <- function(x, ...){
    cat("\n")
    cat("Design-Based Causal Mediation Analysis\n\n")

    # Print values
    if(x$design=="SED.NP.NOSI" | x$design=="PD"){
        if(x$design=="SED.NP.NOSI"){
            cat("Single Experiment Design (without Sequential Ignorability) \n\n")
        }
        if(x$design=="PD"){
            cat("Parallel Design (Interaction Allowed) \n\n")
        }
        smat <- rbind(x$d0, x$d1)
        colnames(smat) <- c("Lower Bound", "Upper Bound")
        rownames(smat) <- c("ACME (control)", "ACME (treated)")

    } else if(x$design=="PD.NINT"){
            clp <- 100 * x$conf.level
            cat("Parallel Design (with No Interaction Assumption) \n\n")
            smat <- t(as.matrix(c(x$d0, x$d0.ci)))
            colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                          paste(clp, "% CI Upper", sep=""))
            rownames(smat) <- "ACME"

    } else if(x$design=="PED"){
        cat("Parallel Encouragement Design\n\n")
        smat <- rbind(x$d0, x$d0.p, x$d1, x$d1.p)
        colnames(smat) <- c("Lower Bound", "Upper Bound")
        rownames(smat) <- c("Population ACME (control)", "Complier ACME (control)",
                             "Population ACME (treated)", "Complier ACME (treated)")

    } else if(x$design=="CED"){
        clp <- 100 * x$conf.level
        cat("Crossover Encouragement Design \n\n")
        smat <- rbind(c(x$d0, x$d0.ci), c(x$d1, x$d1.ci))
        colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                          paste(clp, "% CI Upper", sep=""))
        rownames(smat) <- c("Pliable ACME (control)", "Pliable ACME (treated)")

    } else if(x$design=="SED.NP.SI"){
        clp <- 100 * x$conf.level
        cat("Single Experiment Design with Sequential Ignorability \n\n")

        if(x$boot==TRUE){
            cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
        } else {
            cat("Confidence Intervals Based on Asymptotic Variance\n\n")
        }
        
        smat <- rbind(c(x$d0, x$d0.ci), c(x$d1, x$d1.ci))
        colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
                          paste(clp, "% CI Upper", sep=""))
        rownames(smat) <- c("ACME (control)", "ACME (treated)")
        
    }
    
    printCoefmat(smat, digits=4)
    cat("\n")
    cat("Sample Size Used: ", x$nobs,"\n\n")
    invisible(x)
}
