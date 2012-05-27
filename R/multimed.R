##################################################################
# Sensitivity Analysis for Multiple Mediators
##################################################################

multimed <- function(outcome, med.main, med.alt = NULL, treat, 
					 covariates = NULL, experiment = NULL, 
					 data, design = c("single", "parallel"),
					 sims = 1000, R2.by = 0.01, conf.level = 0.95){

    varnames <- c(outcome, treat, med.main, med.alt, experiment, covariates)
    data <- na.omit(data[,varnames])
    
    design <- match.arg(design)
    
    if(design == "single"){
	    if (is.null(covariates)){
	        f.y <- as.formula(paste(outcome, "~", treat, "+", med.main, "+",
	                                paste(treat, med.main, sep=":"), "+", med.alt, "+",
	                                paste(treat, med.alt, sep=":")))
	        f.ytot <- as.formula(paste(outcome, "~", treat))
	        f.m <- as.formula(paste(med.main, "~", treat))
	        f.w <- as.formula(paste(med.alt, "~", treat))
	    } else {
	        f.y <- as.formula(paste(outcome, "~", treat, "+", med.main, "+",
	                                paste(treat, med.main, sep=":"), "+", med.alt, "+",
	                                paste(treat, med.alt, sep=":"), "+",
	                                paste(covariates, collapse = " + ")))
	        f.ytot <- as.formula(paste(outcome, "~", treat, "+", 
	        						paste(covariates, collapse = " + ")))
	        f.m <- as.formula(paste(med.main, "~", treat, "+", paste(covariates, collapse = " + ")))
	        f.w <- as.formula(paste(med.alt, "~", treat, "+", paste(covariates, collapse = " + ")))
	    }
    } else if (design == "parallel"){
    	if (!is.null(covariates)){
    		warning("covariates currently cannot be used for the parallel design; option ignored")
    	}
        f.y <- as.formula(paste(outcome, "~", treat, "+", med.main, "+",
                                paste(treat, med.main, sep=":")))
        f.ytot <- as.formula(paste(outcome, "~", treat))
        f.m <- as.formula(paste(med.main, "~", treat))
    } else {
    	stop("design unsupported by the multimed procedure")
    }
    
    # Compute sensitivity parameters
    R2.s <- seq(0, 1, by = R2.by)

	data.1 <- switch(design,
					 single = data,
					 parallel = subset(data, data[[experiment]] == 1))
					   
    ETM2 <- mean(data.1[,treat] * data.1[,med.main]^2)
    model.y <- lm(f.y, data=data.1)
    sigma <- summary(model.y)$sigma * sqrt(R2.s/ETM2)

    VY <- var(data.1[,outcome])
    R2.t <- ETM2 * sigma^2/VY
    
    # Bootstrap ACME values
    
    ACME.1.lo <- ACME.0.lo <- ACME.1.up <- ACME.0.up <- 
    ACME.ave.lo <- ACME.ave.up <- matrix(NA, nrow=length(sigma), ncol=sims)
    tau <- rep(NA, length = sims)
    for(b in 1:(sims+1)){
        # Resample
        data.b <- data[sample(1:nrow(data), nrow(data), replace=TRUE),]
        if(b == sims + 1){
        	data.b <- data
        }

		# Subset data for parallel design
		data.b.0 <- switch(design,
						   single = data.b,
						   parallel = subset(data.b, data.b[[experiment]] == 0))
		data.b.1 <- switch(design,
						   single = data.b,
						   parallel = subset(data.b, data.b[[experiment]] == 1))
		
        # Fit models
        model.y <- lm(f.y, data=data.b.1)
        model.ytot <- lm(f.ytot, data=data.b.0)
        model.m <- lm(f.m, data=data.b.0)
        if(design == "single"){
	        model.w <- lm(f.w, data=data.b)
        }
        
        beta3 <- coef(model.y)[treat]
        kappa <- coef(model.y)[paste(treat, med.main, sep=":")]
        if(design == "single"){
	        xi3 <- coef(model.y)[med.alt]
    	    mu3 <- coef(model.y)[paste(treat, med.alt, sep=":")]
		}
		
        # E(M|T=t)
        mf.m1 <- mf.m0 <- mf.m <- model.frame(model.m)
        mf.m1[,treat] <- 1
        mf.m0[,treat] <- 0
        EM.1 <- mean(predict(model.m, mf.m1))
        EM.0 <- mean(predict(model.m, mf.m0))

        # V(M|T=t)
        VM.1 <- sum(model.m$residuals[mf.m[,treat]==1]^2)/(sum(mf.m[,treat]) - length(coef(model.m)))
        VM.0 <- sum(model.m$residuals[mf.m[,treat]==0]^2)/(sum(1-mf.m[,treat]) - length(coef(model.m)))

        # E(W|T=t)
        if(design == "single"){
	        mf.w1 <- mf.w0 <- model.frame(model.w)
	        mf.w1[,treat] <- 1
	        mf.w0[,treat] <- 0
	        EW.1 <- mean(predict(model.w, mf.w1))
	        EW.0 <- mean(predict(model.w, mf.w0))
	    }

        ## Bounds
        # ACME
        if(b == sims + 1){
        	tau.o <- coef(model.ytot)[treat]
        	WXterms <- switch(design,
        					  single = (xi3 + mu3)*EW.1 - xi3*EW.0,
        					  parallel = 0)
	        ACME.1.lo.o <- tau.o - beta3 - kappa*EM.0 - sigma*sqrt(VM.0) - WXterms
	        ACME.0.lo.o <- tau.o - beta3 - kappa*EM.1 - sigma*sqrt(VM.1) - WXterms
	        ACME.1.up.o <- tau.o - beta3 - kappa*EM.0 + sigma*sqrt(VM.0) - WXterms
        	ACME.0.up.o <- tau.o - beta3 - kappa*EM.1 + sigma*sqrt(VM.1) - WXterms
   	        P <- mean(data.b[,treat])
    	    ACME.ave.lo.o <- P * ACME.1.lo.o + (1-P) * ACME.0.lo.o
        	ACME.ave.up.o <- P * ACME.1.up.o + (1-P) * ACME.0.up.o
        } else {
    	    tau[b] <- coef(model.ytot)[treat]
        	WXterms <- switch(design,
        					  single = (xi3 + mu3)*EW.1 - xi3*EW.0,
        					  parallel = 0)
	        ACME.1.lo[,b] <- tau[b] - beta3 - kappa*EM.0 - sigma*sqrt(VM.0) - WXterms
    	    ACME.0.lo[,b] <- tau[b] - beta3 - kappa*EM.1 - sigma*sqrt(VM.1) - WXterms
        	ACME.1.up[,b] <- tau[b] - beta3 - kappa*EM.0 + sigma*sqrt(VM.0) - WXterms
        	ACME.0.up[,b] <- tau[b] - beta3 - kappa*EM.1 + sigma*sqrt(VM.1) - WXterms
	        P <- mean(data.b[,treat])
    	    ACME.ave.lo[,b] <- P * ACME.1.lo[,b] + (1-P) * ACME.0.lo[,b]
        	ACME.ave.up[,b] <- P * ACME.1.up[,b] + (1-P) * ACME.0.up[,b]
        }
    }

    ACME.ave.lo.var <- apply(ACME.ave.lo, 1, var)
    ACME.1.lo.var <- apply(ACME.1.lo, 1, var)
    ACME.0.lo.var <- apply(ACME.0.lo, 1, var)
    ACME.ave.up.var <- apply(ACME.ave.up, 1, var)
    ACME.1.up.var <- apply(ACME.1.up, 1, var)
    ACME.0.up.var <- apply(ACME.0.up, 1, var)

    ACME.ave.CI <- ACME.1.CI <- ACME.0.CI <- matrix(NA, nrow=2, ncol=length(sigma))
    for(i in 1:length(sigma)){
        ACME.ave.CI[,i] <- IMCI(ACME.ave.up.o[i], ACME.ave.lo.o[i], 
                                ACME.ave.up.var[i], ACME.ave.lo.var[i], conf.level = conf.level)$ci
        ACME.1.CI[,i] <- IMCI(ACME.1.up.o[i], ACME.1.lo.o[i], 
                              ACME.1.up.var[i], ACME.1.lo.var[i], conf.level = conf.level)$ci
        ACME.0.CI[,i] <- IMCI(ACME.0.up.o[i], ACME.0.lo.o[i], 
                              ACME.0.up.var[i], ACME.0.lo.var[i], conf.level = conf.level)$ci
    }
    tau.CI <- quantile(tau, probs = c((1-conf.level)/2, (1+conf.level)/2), na.rm = TRUE)
    out <- list(sigma = sigma, R2tilde = R2.t, R2star = R2.s, tau = tau.o, tau.ci = tau.CI,
         d1.ci = ACME.1.CI, d0.ci = ACME.0.CI, d.ave.ci = ACME.ave.CI,
         d1.lb = ACME.1.lo.o, d0.lb = ACME.0.lo.o, d.ave.lb = ACME.ave.lo.o,
         d1.ub = ACME.1.up.o, d0.ub = ACME.0.up.o, d.ave.ub = ACME.ave.up.o)
    class(out) <- "multimed"
    out
}


## Calculates Imbens-Manski confidence set for nonparametric bounds
IMCI <- function(upper, lower, var.upper, var.lower, 
        		 conf.level){
    A <- (upper-lower)/sqrt(max(var.upper,var.lower))
    C <- seq(0,10,by=.001)
    const <- abs(pnorm(C + A) - pnorm(-C) - conf.level)
    Cn <- C[const==min(const)]
    ci <- c(0,0)
    names(ci) <- c("lower","upper")
    ci <- c(lower - Cn*sqrt(var.lower), upper + Cn*sqrt(var.upper))
    if (is.na(ci[1])) ci <- c(NA,NA)
    list(ci=ci, conf.level=conf.level)
}

## Summary
summary.multimed <- function(object, ...){
    structure(object, class = c("summary.multimed", class(object)))
}

print.summary.multimed <- function(x, ...){
    cat("\n")
    cat("Causal Mediation Analysis with Confounding by an Alternative Mechanism\n\n")
    cat("Estimates under the Homogeneous Interaction Assumption:\n")
    
    cmat <- c(x$d1.lb[1], x$d0.lb[1], x$d.ave.lb[1], x$tau)
    cmat <- cbind(cmat, rbind(x$d1.ci[,1], x$d0.ci[,1], x$d.ave.ci[,1], x$tau.ci))
    colnames(cmat) <- c("Estimate", "CI lower", "CI upper")
    rownames(cmat) <- c("ACME(treated)", "ACME(control)", "ACME(average)", "Total")
    printCoefmat(cmat[,1:3], digits=3)
    cat("\n")
    
    cat("Sensitivity Analysis: \n")
    cat("Values of the sensitivity parameters at which ACME first crosses zero:\n")
    ind.d1.b <- sum(sign(x$d1.lb) * sign(x$d1.ub) > 0) + 1
    ind.d1.c <- sum(sign(x$d1.ci[1,]) * sign(x$d1.ci[2,]) > 0) + 1
    ind.d0.b <- sum(sign(x$d0.lb) * sign(x$d0.ub) > 0) + 1
    ind.d0.c <- sum(sign(x$d0.ci[1,]) * sign(x$d0.ci[2,]) > 0) + 1
    ind.d.ave.b <- sum(sign(x$d.ave.lb) * sign(x$d.ave.ub) > 0) + 1
    ind.d.ave.c <- sum(sign(x$d.ave.ci[1,]) * sign(x$d.ave.ci[2,]) > 0) + 1
    smat <- c(x$sigma[ind.d1.b], x$sigma[ind.d1.c], 
               x$R2star[ind.d1.b], x$R2star[ind.d1.c], x$R2tilde[ind.d1.b], x$R2tilde[ind.d1.c])
    smat <- rbind(smat, c(x$sigma[ind.d0.b], x$sigma[ind.d0.c], 
               x$R2star[ind.d0.b], x$R2star[ind.d0.c], x$R2tilde[ind.d0.b], x$R2tilde[ind.d0.c]))
    smat <- rbind(smat, c(x$sigma[ind.d.ave.b], x$sigma[ind.d.ave.c], 
               x$R2star[ind.d.ave.b], x$R2star[ind.d.ave.c], x$R2tilde[ind.d.ave.b], x$R2tilde[ind.d.ave.c]))
    colnames(smat) <- c("sigma(bounds)", "sigma(CI)", "R2s(bounds)", "R2s(CI)", "R2t(bounds)", "R2t(CI)")
    rownames(smat) <- c("ACME(treated)", "ACME(control)", "ACME(average)")
    printCoefmat(smat[,1:6], digits=3)
    cat("\n")
    
    invisible(x)
}


## Plot
plot.multimed <- function(x, type = c("point", "sigma", "R2-residual", "R2-total"),
						  tgroup = c("average", "treated", "control"),
						  ask = prod(par("mfcol")) < nplots,
						  xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, main = NULL,
						  lwd = par("lwd"), pch = par("pch"), cex = par("cex"), las = par("las"), 
						  col.eff = "black", col.cbar = "black", col.creg = "gray", ...){
	
	type <- match.arg(type, several.ok = TRUE)
	tgroup <- match.arg(tgroup, several.ok = TRUE)
	
	show.point <- "point" %in% type
	nplots <- show.point + (length(type) - show.point) * length(tgroup)
	if(ask){
		oask <- devAskNewPage(TRUE)
		on.exit(devAskNewPage(oask))
	}
	
	eff.up <- eff.lo <- ci.up <- ci.lo <- c()
    if("control" %in% tgroup){
    	eff.lo <- cbind(eff.lo, x$d0.lb)
    	eff.up <- cbind(eff.up, x$d0.ub)
    	ci.lo <- cbind(ci.lo, x$d0.ci[1,])
    	ci.up <- cbind(ci.up, x$d0.ci[2,])
    }
    if("treated" %in% tgroup){
    	eff.lo <- cbind(eff.lo, x$d1.lb)
    	eff.up <- cbind(eff.up, x$d1.ub)
    	ci.lo <- cbind(ci.lo, x$d1.ci[1,])
    	ci.up <- cbind(ci.up, x$d1.ci[2,])
    }
    if("average" %in% tgroup){
    	eff.lo <- cbind(eff.lo, x$d.ave.lb)
    	eff.up <- cbind(eff.up, x$d.ave.ub)
    	ci.lo <- cbind(ci.lo, x$d.ave.ci[1,])
    	ci.up <- cbind(ci.up, x$d.ave.ci[2,])
	}
			
	## 1. Point Estimate under Homogeneous Interaction Assumption	  	
	if(show.point){
		if(is.null(main)){
			ma <- "Point Estimate"
		} else ma <- main
		if(is.null(xlab)){
			xla <- "Average Causal Mediation Effects"
		} else xla <- xlab
		if(is.null(ylab)){
            yla <- expression(paste("Total (", bar(tau), ")"))
            if("control" %in% tgroup){
            	 yla <- c(yla, expression(paste("Control (", bar(delta)[0], ")")))
            }
		    if("treated" %in% tgroup){
                 yla <- c(yla, expression(paste("Treated (", bar(delta)[1], ")")))
            }
			if("average" %in% tgroup){
				 yla <- c(yla, expression(paste("Average (", bar(bar(delta)), ")")))
		    }
		} else yla <- ylab
		
        eff <- c(x$tau, eff.lo[1,])
        ci <- cbind(x$tau.ci, rbind(ci.lo[1,], ci.up[1,]))
		
		if(is.null(xlim)){
			xli <- c(min(ci), max(ci))
		} else xli <- xlim
		if(is.null(ylim)){
			yli <- c(0, length(eff)) + 0.5
		} else yli <- ylim
		
		plot(0, 0, type = "n", main = ma, xlab = xla, ylab = "",
			 xlim = xli, ylim = yli, yaxt = "n", ...)
	    for(i in 1:length(eff)){
	    	segments(ci[1,i], i, ci[2,i], i, lwd = lwd, col = col.cbar)
	    	points(eff[i], i, pch = pch, cex = cex, col = col.eff)
	    }
	    abline(v = 0)
	    axis(side = 2, labels = yla, at = 1:length(eff), las = las)
	}
	
	## 2. Sensitivity analysis
	
	if(is.null(ylab)){
        yla <- as.list(rep(NA, 3))
        if("control" %in% tgroup){
           	 yla[[1]] <- c(expression(paste(bar(delta)[0], "(", sigma, ")")),
           	 			   expression(paste(bar(delta)[0], "(", R^{2}, "*)")),
           	 			   expression(paste(bar(delta)[0], "(", tilde(R)^{2}, ")")))
        }
	    if("treated" %in% tgroup){
           	 yla[[2]] <- c(expression(paste(bar(delta)[1], "(", sigma, ")")),
           	 			   expression(paste(bar(delta)[1], "(", R^{2}, "*)")),
           	 			   expression(paste(bar(delta)[1], "(", tilde(R)^{2}, ")")))
        }
		if("average" %in% tgroup){
           	 yla[[3]] <- c(expression(paste(bar(bar(delta)), "(", sigma, ")")),
           	 			   expression(paste(bar(bar(delta)), "(", R^{2}, "*)")),
           	 			   expression(paste(bar(bar(delta)), "(", tilde(R)^{2}, ")")))
	    }
	    yla <- yla[!is.na(yla)]
	} else yla <- ylab
	
	wh <- c("sigma", "R2-residual", "R2-total") %in% type
	
	for(j in 1:length(wh)){
		if(!wh[j]) next
		if(is.null(main)){
			ma <- ifelse(j == 1, "Sensitivity with Respect to \n Interaction Heterogeneity",
								 "Sensitivity with Respect to \n Importance of Interaction")
		} else ma <- main
		if(is.null(xlab)){
			xla <- switch(j, expression(sigma),
							 expression(paste(R^{2}, "*")),
							 expression(tilde(R)^{2}))
		} else xla <- xlab
		if(is.null(xlim)){
			if(j == 1) xli <- range(x$sigma) else xli <- c(0,1)
		} else xli <- xlim
		if(is.null(ylim)){
			yli <- c(min(ci.lo), max(ci.up))
		} else yli <- ylim
		
		spar <- switch(j, x$sigma, x$R2star, x$R2tilde)
		
		for(i in 1:ncol(eff.lo)){
			plot(0, 0, type = "n", main = ma, xlab = xla, ylab = yla[[i]][j],
				 xlim = xli, ylim = yli, ...)
			polygon(c(spar, rev(spar)), c(ci.lo[,i], rev(ci.up[,i])),
				    border = FALSE, col = col.creg)
			lines(spar, eff.lo[,i], lwd = lwd, col = col.eff)
			lines(spar, eff.up[,i], lwd = lwd, col = col.eff)
			abline(h = 0)
			abline(h = eff.lo[1,i], lty = "dashed")
		}
	}
}
