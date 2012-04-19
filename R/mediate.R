mediate <- function(model.m, model.y, sims = 1000, boot = FALSE,
                    treat = "treat.name", mediator = "med.name",
                    covariates = NULL, outcome = NULL,
                    control = NULL, conf.level = .95,
                    control.value = 0, treat.value = 1,
                    long = TRUE, dropobs = FALSE,
                    robustSE = FALSE, cluster = NULL, ...){

    # Warn users who still use INT option
    if(match("INT", names(match.call()), 0L)){
        warning("'INT' is deprecated - existence of interaction terms is now automatically detected from model formulas")
    }

    # Warning for robustSE and cluster used with boot
    if(robustSE && boot){
        warning("'robustSE' is ignored for nonparametric bootstrap")
    }

    if(!is.null(cluster) && boot){
        warning("'cluster' is ignored for nonparametric bootstrap")
    }

    if(robustSE & !is.null(cluster)){
        stop("choose either `robustSE' or `cluster' option, not both")
    }

    # Drop observations not common to both mediator and outcome models
    if(dropobs){
        odata.m <- model.frame(model.m)
        odata.y <- model.frame(model.y)
        newdata <- merge(odata.m, odata.y, sort=FALSE,
                    by=c("row.names", intersect(names(odata.m), names(odata.y))))
        rownames(newdata) <- newdata$Row.names
        newdata <- newdata[,-1L]
        rm(odata.m, odata.y)

        call.m <- getCall(model.m)
        call.y <- getCall(model.y)

        call.m$data <- call.y$data <- newdata
        if(c("(weights)") %in% names(newdata)){
            call.m$weights <- call.y$weights <- model.weights(newdata)
        }
        model.m <- eval.parent(call.m)
        model.y <- eval.parent(call.y)
    }

    # Model type indicators
    isGam.y <- inherits(model.y, "gam")
    isGam.m <- inherits(model.m, "gam")
    isGlm.y <- inherits(model.y, "glm")  # Note gam and bayesglm also inherits "glm"
    isGlm.m <- inherits(model.m, "glm")  # Note gam and bayesglm also inherits "glm"
    isLm.y <- inherits(model.y, "lm")    # Note gam, glm and bayesglm also inherit "lm"
    isLm.m <- inherits(model.m, "lm")    # Note gam, glm and bayesglm also inherit "lm"
    isVglm.y <- inherits(model.y, "vglm")
    isRq.y <- inherits(model.y, "rq")
    isRq.m <- inherits(model.m, "rq")
    isOrdered.y <- inherits(model.y, "polr")  # Note bayespolr also inherits "polr"
    isOrdered.m <- inherits(model.m, "polr")  # Note bayespolr also inherits "polr"
    isSurvreg.y <- inherits(model.y, "survreg")
    isSurvreg.m <- inherits(model.m, "survreg")

    # Record family of model.m if glm
    if(isGlm.m){
        FamilyM <- model.m$family$family
    }

    # Record vfamily of model.y if vglm (currently only tobit)
    if(isVglm.y){
        VfamilyY <- model.y@family@vfamily
    }

    # Warning for unused options
    if(!is.null(control) && !isGam.y){
        warning("'control' is only used for GAM outcome models - ignored")
        control <- NULL
    }
    if(!is.null(outcome) && !(isSurvreg.y && boot)){
        warning("'outcome' is only relevant for survival outcome models with bootstrap - ignored")
    }

    # Model frames for M and Y models
    m.data <- model.frame(model.m)  # Call.M$data
    y.data <- model.frame(model.y)  # Call.Y$data

    # Numbers of observations and categories
    n.m <- nrow(m.data)
    n.y <- nrow(y.data)
    if(n.m != n.y){
        stop("number of observations do not match between mediator and outcome models")
    } else{
        n <- n.m
    }
    m <- length(sort(unique(model.frame(model.m)[,1])))

    # Extracting weights from models
    weights.m <- model.weights(m.data)
    weights.y <- model.weights(y.data)

    if(!is.null(weights.m) && isGlm.m && FamilyM == "binomial"){
        message("weights taken as sampling weights, not total number of trials")
    }

    if(is.null(weights.m)){
        weights.m <- rep(1,nrow(m.data))
    }
    if(is.null(weights.y)){
        weights.y <- rep(1,nrow(y.data))
    }
    if(!all(weights.m == weights.y)) {
        stop("weights on outcome and mediator models not identical")
    } else {
        weights <- weights.m
    }

    # Convert character treatment to factor
    if(is.character(m.data[,treat])){
        m.data[,treat] <- factor(m.data[,treat])
    }
    if(is.character(y.data[,treat])){
        y.data[,treat] <- factor(y.data[,treat])
    }

    # Convert character mediator to factor
    if(is.character(y.data[,mediator])){
        y.data[,mediator] <- factor(y.data[,mediator])
    }

    # Factor treatment indicator
    isFactorT.m <- is.factor(m.data[,treat])
    isFactorT.y <- is.factor(y.data[,treat])
    if(isFactorT.m != isFactorT.y){
        stop("treatment variable types differ in mediator and outcome models")
    } else {
        isFactorT <- isFactorT.y
    }

    if(isFactorT){
        t.levels <- levels(y.data[,treat])
        if(treat.value %in% t.levels & control.value %in% t.levels){
            cat.0 <- control.value
            cat.1 <- treat.value
        } else {
            cat.0 <- t.levels[1]
            cat.1 <- t.levels[2]
            warning("treatment and control values do not match factor levels; using ", cat.0, " and ", cat.1, " as control and treatment, respectively")
        }
    } else {
        cat.0 <- control.value
        cat.1 <- treat.value
    }

    # Factor mediator indicator
    isFactorM <- is.factor(y.data[,mediator])

    if(isFactorM){
        m.levels <- levels(y.data[,mediator])
    }

    #####################################
    ## Define functions
    #####################################

    indexmax <- function(x){
        ## Return position of largest element in vector x
        order(x)[length(x)]
    }

    getvcov <- function(dat, fm, cluster){
        ## Compute cluster robust standard errors
        ## fm is the model object
        attach(dat, warn.conflicts = F)
        M <- length(unique(cluster))
        N <- length(cluster)
        K <- fm$rank
        dfc <- (M/(M-1))*((N-1)/(N-K))
        uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
        dfc*sandwich(fm, meat. = crossprod(uj)/N)
    }

    ############################################################################
    ############################################################################
    ### CASE I: EVERYTHING EXCEPT ORDERED OUTCOME
    ############################################################################
    ############################################################################
    if (!isOrdered.y) {

        ########################################################################
        ## Case I-1: Quasi-Bayesian Monte Carlo
        ########################################################################
        if(!boot){
            # Error if gam outcome or quantile mediator
            if(isGam.m | isGam.y | isRq.m){
                stop("'boot' must be 'TRUE' for models used")
            }

            # Get mean and variance parameters for mediator simulations
            if(isSurvreg.m && is.null(survreg.distributions[[model.m$dist]]$scale)){
            	MModel.coef <- c(coef(model.m), log(model.m$scale))
            	scalesim.m <- TRUE
            } else {
	            MModel.coef <- coef(model.m)
	            scalesim.m <- FALSE
            }

            if(isOrdered.m){
                if(is.null(model.m$Hess)){
                    cat("Mediator model object does not contain 'Hessian';")
                }
                k <- length(MModel.coef)
                MModel.var.cov <- vcov(model.m)[(1:k),(1:k)]
            } else if(isSurvreg.m){
            	MModel.var.cov <- vcov(model.m)
            } else {
                if(robustSE){
                    MModel.var.cov <- vcovHC(model.m, ...)
                } else if(!is.null(cluster)){
                    MModel.var.cov <- getvcov(m.data, model.m, cluster)
                } else {
                    MModel.var.cov <- vcov(model.m)
                }
            }

            # Get mean and variance parameters for outcome simulations
            if(isSurvreg.y && is.null(survreg.distributions[[model.y$dist]]$scale)){
            	YModel.coef <- c(coef(model.y), log(model.y$scale))
            	scalesim.y <- TRUE  # indicates if survreg scale parameter is simulated
            } else {
	            YModel.coef <- coef(model.y)
	            scalesim.y <- FALSE
            }

            if(isRq.y){
                YModel.var.cov <- summary(model.y, covariance=TRUE)$cov
            } else if(isSurvreg.y){
            	YModel.var.cov <- vcov(model.y)
            } else {
                if(robustSE){
                    YModel.var.cov <- vcovHC(model.y, ...)
                } else if(!is.null(cluster)){
                    YModel.var.cov <- getvcov(y.data, model.y, cluster)
                } else {
                    YModel.var.cov <- vcov(model.y)
                }
            }

            # Draw model coefficients from normal
            if(sum(is.na(MModel.coef), is.na(YModel.coef)) > 0){
            	stop("NA in model coefficients; rerun models with nonsingular design matrix")
            }
            MModel <- mvrnorm(sims, mu=MModel.coef, Sigma=MModel.var.cov)
            YModel <- mvrnorm(sims, mu=YModel.coef, Sigma=YModel.var.cov)

            if(robustSE && (isSurvreg.m | isSurvreg.y)){
           		warning("`robustSE' ignored for survival models; fit the model with `robust' option instead\n")
           	}
           	if(!is.null(cluster) && (isSurvreg.m | isSurvreg.y)){
           		warning("`cluster' ignored for survival models; fit the model with 'cluster()' term in the formula\n")
           	}
           	
            #####################################
            #  Mediator Predictions
            #####################################
            pred.data.t <- pred.data.c <- m.data

            if(isFactorT){
                pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
            } else {
                pred.data.t[,treat] <- cat.1
                pred.data.c[,treat] <- cat.0
            }

            if(!is.null(covariates)){
            	for(p in 1:length(covariates)){
            		vl <- names(covariates[p])
            		x <- ifelse(is.factor(pred.data.t[,vl]),
            						factor(covariates[[p]], levels = levels(m.data[,vl])),
            						covariates[[p]])
            		pred.data.t[,vl] <- pred.data.c[,vl] <- x
            	}
            }

            mmat.t <- model.matrix(terms(model.m), data=pred.data.t)
            mmat.c <- model.matrix(terms(model.m), data=pred.data.c)

            ### Case I-1-a: GLM Mediator
            if(isGlm.m){

                muM1 <- model.m$family$linkinv(tcrossprod(MModel, mmat.t))
                muM0 <- model.m$family$linkinv(tcrossprod(MModel, mmat.c))

                if(FamilyM == "poisson"){
                    PredictM1 <- matrix(rpois(sims*n, lambda = muM1), nrow = sims)
                    PredictM0 <- matrix(rpois(sims*n, lambda = muM0), nrow = sims)
                } else if (FamilyM == "Gamma") {
                    shape <- gamma.shape(model.m)$alpha
                    PredictM1 <- matrix(rgamma(n*sims, shape = shape,
                                        scale = muM1/shape), nrow = sims)
                    PredictM0 <- matrix(rgamma(n*sims, shape = shape,
                                        scale = muM0/shape), nrow = sims)
                } else if (FamilyM == "binomial"){
                    PredictM1 <- matrix(rbinom(n*sims, size = 1,
                                        prob = muM1), nrow = sims)
                    PredictM0 <- matrix(rbinom(n*sims, size = 1,
                                        prob = muM0), nrow = sims)
                } else if (FamilyM == "gaussian"){
                    sigma <- sqrt(summary(model.m)$dispersion)
                    error <- rnorm(sims*n, mean=0, sd=sigma)
                    PredictM1 <- muM1 + matrix(error, nrow=sims)
                    PredictM0 <- muM0 + matrix(error, nrow=sims)
                } else if (FamilyM == "inverse.gaussian"){
                    disp <- summary(model.m)$dispersion
                    PredictM1 <- matrix(SuppDists::rinvGauss(n*sims, nu = muM1,
                                        lambda = 1/disp), nrow = sims)
                    PredictM0 <- matrix(SuppDists::rinvGauss(n*sims, nu = muM0,
                                        lambda = 1/disp), nrow = sims)
                } else {
                    stop("unsupported glm family")
                }

            ### Case I-1-b: Ordered mediator
            } else if(isOrdered.m){
                if(model.m$method=="logistic"){
                    linkfn <- plogis
                } else if(model.m$method=="probit") {
                    linkfn <- pnorm
                } else {
                    stop("unsupported polr method; use 'logistic' or 'probit'")
                }

                m.cat <- sort(unique(model.frame(model.m)[,1]))
                lambda <- model.m$zeta

                mmat.t <- mmat.t[,-1]
                mmat.c <- mmat.c[,-1]

                ystar_m1 <- tcrossprod(MModel, mmat.t)
                ystar_m0 <- tcrossprod(MModel, mmat.c)

                PredictM1 <- matrix(,nrow=sims, ncol=n)
                PredictM0 <- matrix(,nrow=sims, ncol=n)

                for(i in 1:sims){

                    cprobs_m1 <- matrix(NA,n,m)
                    cprobs_m0 <- matrix(NA,n,m)
                    probs_m1 <- matrix(NA,n,m)
                    probs_m0 <- matrix(NA,n,m)

                    for (j in 1:(m-1)) {  # loop to get category-specific probabilities
                        cprobs_m1[,j] <- linkfn(lambda[j]-ystar_m1[i,])
                        cprobs_m0[,j] <- linkfn(lambda[j]-ystar_m0[i,])
                                                           # cumulative probabilities
                        probs_m1[,m] <- 1-cprobs_m1[,m-1] # top category
                        probs_m0[,m] <- 1-cprobs_m0[,m-1] # top category
                        probs_m1[,1] <- cprobs_m1[,1]     # bottom category
                        probs_m0[,1] <- cprobs_m0[,1]     # bottom category
                    }

                    for (j in 2:(m-1)){  # middle categories
                        probs_m1[,j] <- cprobs_m1[,j]-cprobs_m1[,j-1]
                        probs_m0[,j] <- cprobs_m0[,j]-cprobs_m0[,j-1]
                    }

                    draws_m1 <- matrix(NA, n, m)
                    draws_m0 <- matrix(NA, n, m)

                    for(ii in 1:n){
                        draws_m1[ii,] <- t(rmultinom(1, 1, prob = probs_m1[ii,]))
                        draws_m0[ii,] <- t(rmultinom(1, 1, prob = probs_m0[ii,]))
                    }

                    PredictM1[i,] <- apply(draws_m1, 1, indexmax)
                    PredictM0[i,] <- apply(draws_m0, 1, indexmax)
                }

            ### Case I-1-c: Linear
            } else if(isLm.m){
                sigma <- summary(model.m)$sigma
                error <- rnorm(sims*n, mean=0, sd=sigma)
                muM1 <- tcrossprod(MModel, mmat.t)
                muM0 <- tcrossprod(MModel, mmat.c)
                PredictM1 <- muM1 + matrix(error, nrow=sims)
                PredictM0 <- muM0 + matrix(error, nrow=sims)
                rm(error)

            ### Case I-1-d: Survreg
            } else if(isSurvreg.m){
       			dd <- survreg.distributions[[model.m$dist]]
				if (is.null(dd$itrans)){
					itrans <- function(x) x
				} else {
					itrans <- dd$itrans
				}
				dname <- dd$dist
				if(is.null(dname)){
					dname <- model.m$dist
				}
				if(scalesim.m){
    				scale <- exp(MModel[,ncol(MModel)])
					lpM1 <- tcrossprod(MModel[,1:(ncol(MModel)-1)], mmat.t)
					lpM0 <- tcrossprod(MModel[,1:(ncol(MModel)-1)], mmat.c)
				} else {
					scale <- dd$scale
					lpM1 <- tcrossprod(MModel, mmat.t)
					lpM0 <- tcrossprod(MModel, mmat.c)
				}
				error <- switch(dname,
								extreme = log(rweibull(sims*n, shape=1, scale=1)),
								gaussian = rnorm(sims*n),
								logistic = rlogis(sims*n),
								t = rt(sims*n, df=dd$parms))
			    PredictM1 <- itrans(lpM1 + scale * matrix(error, nrow=sims))
			    PredictM0 <- itrans(lpM0 + scale * matrix(error, nrow=sims))
			    rm(error)
			
            } else {
                stop("mediator model is not yet implemented")
            }
            rm(mmat.t, mmat.c)

            #####################################
            ##  Outcome Predictions
            #####################################
            effects.tmp <- array(NA, dim = c(n, sims, 4))
            for(e in 1:4){
                tt <- switch(e, c(1,1,1,0), c(0,0,1,0), c(1,0,1,1), c(1,0,0,0))
                Pr1 <- matrix(, nrow=n, ncol=sims)
                Pr0 <- matrix(, nrow=n, ncol=sims)

                for(j in 1:sims){
                    pred.data.t <- pred.data.c <- y.data

		            if(!is.null(covariates)){
        		    	for(p in 1:length(covariates)){
            				vl <- names(covariates[p])
            				x <- ifelse(is.factor(pred.data.t[,vl]),
            								factor(covariates[[p]], levels = levels(y.data[,vl])),
            								covariates[[p]])
     			       		pred.data.t[,vl] <- pred.data.c[,vl] <- x
            			}
		            }
		
                    # Set treatment values
                    cat.t <- ifelse(tt[1], cat.1, cat.0)
                    cat.c <- ifelse(tt[2], cat.1, cat.0)
                    cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
                    cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
                    if(isFactorT){
                        pred.data.t[,treat] <- factor(cat.t, levels = t.levels)
                        pred.data.c[,treat] <- factor(cat.c, levels = t.levels)
                        if(!is.null(control)){
                            pred.data.t[,control] <- factor(cat.t.ctrl, levels = t.levels)
                            pred.data.c[,control] <- factor(cat.c.ctrl, levels = t.levels)
                        }
                    } else {
                        pred.data.t[,treat] <- cat.t
                        pred.data.c[,treat] <- cat.c
                        if(!is.null(control)){
                            pred.data.t[,control] <- cat.t.ctrl
                            pred.data.c[,control] <- cat.c.ctrl
                        }
                    }

                    # Set mediator values
                    PredictMt <- PredictM1[j,] * tt[3] + PredictM0[j,] * (1 - tt[3])
                    PredictMc <- PredictM1[j,] * tt[4] + PredictM0[j,] * (1 - tt[4])
                    if(isFactorM) {
                        pred.data.t[,mediator] <- factor(PredictMt, levels=1:m, labels=m.levels)
                        pred.data.c[,mediator] <- factor(PredictMc, levels=1:m, labels=m.levels)
                    } else {
                        pred.data.t[,mediator] <- PredictMt
                        pred.data.c[,mediator] <- PredictMc
                    }

                    ymat.t <- model.matrix(terms(model.y), data=pred.data.t)
                    ymat.c <- model.matrix(terms(model.y), data=pred.data.c)

                    if(isVglm.y){
                        if(VfamilyY=="tobit") {
                            Pr1.tmp <- ymat.t %*% YModel[j,-2]
                            Pr0.tmp <- ymat.c %*% YModel[j,-2]
                            Pr1[,j] <- pmin(pmax(Pr1.tmp, model.y@misc$Lower), model.y@misc$Upper)
                            Pr0[,j] <- pmin(pmax(Pr0.tmp, model.y@misc$Lower), model.y@misc$Upper)
                        } else {
                            stop("outcome model is in unsupported vglm family")
                        }
                    } else if(scalesim.y){
                        Pr1[,j] <- t(as.matrix(YModel[j,1:(ncol(YModel)-1)])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(YModel[j,1:(ncol(YModel)-1)])) %*% t(ymat.c)
                    } else {
                        Pr1[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.t)
                        Pr0[,j] <- t(as.matrix(YModel[j,])) %*% t(ymat.c)
                    }
                    rm(ymat.t, ymat.c, pred.data.t, pred.data.c)
                }

                if(isGlm.y){
                    Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
                    Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
                } else if(isSurvreg.y){
	       			dd <- survreg.distributions[[model.y$dist]]
					if (is.null(dd$itrans)){
						itrans <- function(x) x
					} else {
						itrans <- dd$itrans
					}
                	Pr1 <- apply(Pr1, 2, itrans)
                	Pr0 <- apply(Pr0, 2, itrans)
                }

                effects.tmp[,,e] <- Pr1 - Pr0
                rm(Pr1, Pr0)
            }
            rm(PredictM1, PredictM0, YModel, MModel)

            delta.1 <- t(as.matrix(apply(effects.tmp[,,1], 2, weighted.mean, w=weights)))
            delta.0 <- t(as.matrix(apply(effects.tmp[,,2], 2, weighted.mean, w=weights)))
            zeta.1 <- t(as.matrix(apply(effects.tmp[,,3], 2, weighted.mean, w=weights)))
            zeta.0 <- t(as.matrix(apply(effects.tmp[,,4], 2, weighted.mean, w=weights)))
            rm(effects.tmp)
            tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
	        nu.0 <- delta.0/tau
    	    nu.1 <- delta.1/tau
    	    delta.avg <- (delta.1 + delta.0)/2
    	    zeta.avg <- (zeta.1 + zeta.0)/2
    	    nu.avg <- (nu.1 + nu.0)/2

	        d0 <- mean(delta.0)
    	    d1 <- mean(delta.1)
	        z1 <- mean(zeta.1)
	        z0 <- mean(zeta.0)
	        tau.coef <- mean(tau)
	        n0 <- median(nu.0)
    	    n1 <- median(nu.1)
    	    d.avg <- (d0 + d1)/2
    	    z.avg <- (z0 + z1)/2
    	    n.avg <- (n0 + n1)/2

        ########################################################################
        ## Case I-2: Nonparametric Bootstrap
        ########################################################################
        } else {

            Call.M <- getCall(model.m)
            Call.Y <- getCall(model.y)

            # Storage
            delta.1 <- matrix(NA, sims, 1)
            delta.0 <- matrix(NA, sims, 1)
            zeta.1 <- matrix(NA, sims, 1)
            zeta.0 <- matrix(NA, sims, 1)
            tau <- matrix(NA, sims, 1)

            # Bootstrap loop begins
            for(b in 1:(sims+1)){
                index <- sample(1:n, n, replace = TRUE)

                if(b == sims+1){  # in the final run, use the original data
                	index <- 1:n
                }

                if(isSurvreg.m){
                	if(ncol(model.m$y) > 2){
                		stop("unsupported censoring type")
                	}
                	mname <- names(m.data)[1]
                	if(substr(mname, 1, 4) != "Surv"){
                		stop("refit the survival model with `Surv' used directly in model formula")
                	}
                	nc <- nchar(mediator)
                	eventname <- substr(mname, 5 + nc + 3, nchar(mname) - 1)
                	if(nchar(eventname) == 0){
	                	m.data.tmp <- data.frame(m.data,
    	            							as.numeric(m.data[,1L][,1L]))
            	    	names(m.data.tmp)[c(1L, ncol(m.data)+1)] <- c(mname, mediator)
                	} else {
	                	m.data.tmp <- data.frame(m.data,
    	            							as.numeric(m.data[,1L][,1L]),
        	        							as.numeric(model.m$y[,2]))
            	    	names(m.data.tmp)[c(1L, ncol(m.data)+(1:2))] <- c(mname, mediator, eventname)
                	}
                	Call.M$data <- m.data.tmp[index,]
                } else {
                	Call.M$data <- m.data[index,]
                }

                if(isSurvreg.y){
                	if(ncol(model.y$y) > 2){
                		stop("unsupported censoring type")
                	}
                	yname <- names(y.data)[1]
                	if(substr(yname, 1, 4) != "Surv"){
                		stop("refit the survival model with `Surv' used directly in model formula")
                	}
                	if(is.null(outcome)){
                		stop("`outcome' must be supplied for survreg outcome with boot")
                	}
                	nc <- nchar(outcome)
                	eventname <- substr(yname, 5 + nc + 3, nchar(yname) - 1)
                	if(nchar(eventname) == 0){
	                	y.data.tmp <- data.frame(y.data,
    	            							as.numeric(y.data[,1L][,1L]))
            	    	names(y.data.tmp)[c(1L, ncol(y.data)+1)] <- c(yname, outcome)
                	} else {
	                	y.data.tmp <- data.frame(y.data,
    	            							as.numeric(y.data[,1L][,1L]),
        	        							as.numeric(model.y$y[,2]))
            	    	names(y.data.tmp)[c(1L, ncol(y.data)+(1:2))] <- c(yname, outcome, eventname)
                	}
                	Call.Y$data <- y.data.tmp[index,]
                } else {
                	Call.Y$data <- y.data[index,]
                }

                Call.M$weights <- m.data[index,"(weights)"]
                Call.Y$weights  <- y.data[index,"(weights)"]

                if(isOrdered.m && length(unique(y.data[index,mediator]))!=m){
                        stop("insufficient variation on mediator")
                }

                # Refit Models with Resampled Data
                new.fit.M <- eval.parent(Call.M)
                new.fit.Y <- eval.parent(Call.Y)

                #####################################
                #  Mediator Predictions
                #####################################
                pred.data.t <- pred.data.c <- m.data

                if(isFactorT){
                    pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
                } else {
                    pred.data.t[,treat] <- cat.1
                    pred.data.c[,treat] <- cat.0
                }

	            if(!is.null(covariates)){
    	        	for(p in 1:length(covariates)){
        	    		vl <- names(covariates[p])
            			x <- ifelse(is.factor(pred.data.t[,vl]),
            							factor(covariates[[p]], levels = levels(m.data[,vl])),
            							covariates[[p]])
  		          		pred.data.t[,vl] <- pred.data.c[,vl] <- x
        	    	}
           		}

                ### Case I-2-a: GLM Mediator (including GAMs)
                if(isGlm.m){

                    muM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
                    muM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)

                    if(FamilyM == "poisson"){
                        PredictM1 <- rpois(n, lambda = muM1)
                        PredictM0 <- rpois(n, lambda = muM0)
                    } else if (FamilyM == "Gamma") {
                        shape <- gamma.shape(new.fit.M)$alpha
                        PredictM1 <- rgamma(n, shape = shape, scale = muM1/shape)
                        PredictM0 <- rgamma(n, shape = shape, scale = muM0/shape)
                    } else if (FamilyM == "binomial"){
                        PredictM1 <- rbinom(n, size = 1, prob = muM1)
                        PredictM0 <- rbinom(n, size = 1, prob = muM0)
                    } else if (FamilyM == "gaussian"){
                        sigma <- sqrt(summary(new.fit.M)$dispersion)
                        error <- rnorm(n, mean=0, sd=sigma)
                        PredictM1 <- muM1 + error
                        PredictM0 <- muM0 + error
                    } else if (FamilyM == "inverse.gaussian"){
                        disp <- summary(new.fit.M)$dispersion
                        PredictM1 <- SuppDists::rinvGauss(n, nu = muM1, lambda = 1/disp)
                        PredictM0 <- SuppDists::rinvGauss(n, nu = muM0, lambda = 1/disp)
                    } else {
                        stop("unsupported glm family")
                    }

                ### Case I-2-b: Ordered Mediator
                } else if(isOrdered.m) {
                    probs_m1 <- predict(new.fit.M, newdata=pred.data.t, type="probs")
                    probs_m0 <- predict(new.fit.M, newdata=pred.data.c, type="probs")
                    draws_m1 <- matrix(NA, n, m)
                    draws_m0 <- matrix(NA, n, m)
                    for(ii in 1:n){
                        draws_m1[ii,] <- t(rmultinom(1, 1, prob = probs_m1[ii,]))
                        draws_m0[ii,] <- t(rmultinom(1, 1, prob = probs_m0[ii,]))
                    }
                    PredictM1 <- apply(draws_m1, 1, indexmax)
                    PredictM0 <- apply(draws_m0, 1, indexmax)

                ### Case I-2-c: Quantile Regression for Mediator
                } else if(isRq.m){
                    # Use inverse transform sampling to predict M
                    call.new <- new.fit.M$call
                    call.new$tau <- runif(n)
                    newfits <- eval.parent(call.new)
                    tt <- delete.response(terms(new.fit.M))
                    m.t <- model.frame(tt, pred.data.t, xlev = new.fit.M$xlevels)
                    m.c <- model.frame(tt, pred.data.c, xlev = new.fit.M$xlevels)
                    X.t <- model.matrix(tt, m.t, contrasts = new.fit.M$contrasts)
                    X.c <- model.matrix(tt, m.c, contrasts = new.fit.M$contrasts)
                    rm(tt, m.t, m.c)
                    PredictM1 <- rowSums(X.t * t(newfits$coefficients))
                    PredictM0 <- rowSums(X.c * t(newfits$coefficients))
                    rm(newfits, X.t, X.c)

                ### Case I-2-d: Linear
                } else if(isLm.m){
                    sigma <- summary(new.fit.M)$sigma
                    error <- rnorm(n, mean=0, sd=sigma)
                    PredictM1 <- predict(new.fit.M, type="response",
                                          newdata=pred.data.t) + error
                    PredictM0 <- predict(new.fit.M, type="response",
                                          newdata=pred.data.c) + error
                    rm(error)

	            ### Case I-2-e: Survreg
   		        } else if(isSurvreg.m){
       				dd <- survreg.distributions[[new.fit.M$dist]]
					if (is.null(dd$itrans)){
						itrans <- function(x) x
					} else {
						itrans <- dd$itrans
					}
					dname <- dd$dist
					if(is.null(dname)){
						dname <- new.fit.M$dist
					}
					scale <- new.fit.M$scale
					lpM1 <- predict(new.fit.M, newdata=pred.data.t, type="linear")
					lpM0 <- predict(new.fit.M, newdata=pred.data.c, type="linear")
					error <- switch(dname,
									extreme = log(rweibull(n, shape=1, scale=1)),
									gaussian = rnorm(n),
									logistic = rlogis(n),
									t = rt(n, df=dd$parms))
				    PredictM1 <- as.numeric(itrans(lpM1 + scale * error))
				    PredictM0 <- as.numeric(itrans(lpM0 + scale * error))
				    rm(error)

                } else {
                    stop("mediator model is not yet implemented")
                }

                #####################################
                #  Outcome Predictions
                #####################################
                effects.tmp <- matrix(NA, nrow = n, ncol = 4)
                for(e in 1:4){
                    tt <- switch(e, c(1,1,1,0), c(0,0,1,0), c(1,0,1,1), c(1,0,0,0))
                    pred.data.t <- pred.data.c <- y.data

		            if(!is.null(covariates)){
   	 	        		for(p in 1:length(covariates)){
        	    			vl <- names(covariates[p])
            				x <- ifelse(is.factor(pred.data.t[,vl]),
            								factor(covariates[[p]], levels = levels(y.data[,vl])),
            								covariates[[p]])
        	    			pred.data.t[,vl] <- pred.data.c[,vl] <- x
            			}
            		}

                    # Set treatment values
                    cat.t <- ifelse(tt[1], cat.1, cat.0)
                    cat.c <- ifelse(tt[2], cat.1, cat.0)
                    cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
                    cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
                    if(isFactorT){
                        pred.data.t[,treat] <- factor(cat.t, levels = t.levels)
                        pred.data.c[,treat] <- factor(cat.c, levels = t.levels)
                        if(!is.null(control)){
                            pred.data.t[,control] <- factor(cat.t.ctrl, levels = t.levels)
                            pred.data.c[,control] <- factor(cat.c.ctrl, levels = t.levels)
                        }
                    } else {
                        pred.data.t[,treat] <- cat.t
                        pred.data.c[,treat] <- cat.c
                        if(!is.null(control)){
                            pred.data.t[,control] <- cat.t.ctrl
                            pred.data.c[,control] <- cat.c.ctrl
                        }
                    }

                    # Set mediator values
                    PredictM1.tmp <- PredictM1
                    PredictM0.tmp <- PredictM0
                    PredictMt <- PredictM1 * tt[3] + PredictM0 * (1 - tt[3])
                    PredictMc <- PredictM1 * tt[4] + PredictM0 * (1 - tt[4])
                    if(isFactorM) {
                        pred.data.t[,mediator] <- factor(PredictMt, levels=1:m, labels=m.levels)
                        pred.data.c[,mediator] <- factor(PredictMc, levels=1:m, labels=m.levels)
                    } else {
                        pred.data.t[,mediator] <- PredictMt
                        pred.data.c[,mediator] <- PredictMc
                    }

                    if(isRq.y){
                        pr.1 <- predict(new.fit.Y, type="response",
                                        newdata=pred.data.t, interval="none")
                        pr.0 <- predict(new.fit.Y, type="response",
                                        newdata=pred.data.c, interval="none")
                    } else {
                        pr.1 <- predict(new.fit.Y, type="response",
                                        newdata=pred.data.t)
                        pr.0 <- predict(new.fit.Y, type="response",
                                        newdata=pred.data.c)
                    }
                    pr.mat <- as.matrix(cbind(pr.1, pr.0))
                    effects.tmp[,e] <- pr.mat[,1] - pr.mat[,2]

                    rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
                }

                # Compute all QoIs
                if(b == sims+1){
	                d1 <- weighted.mean(effects.tmp[,1], weights)
    	            d0 <- weighted.mean(effects.tmp[,2], weights)
        	        z1 <- weighted.mean(effects.tmp[,3], weights)
            	    z0 <- weighted.mean(effects.tmp[,4], weights)
                } else {
	                delta.1[b] <- weighted.mean(effects.tmp[,1], weights)
    	            delta.0[b] <- weighted.mean(effects.tmp[,2], weights)
        	        zeta.1[b] <- weighted.mean(effects.tmp[,3], weights)
            	    zeta.0[b] <- weighted.mean(effects.tmp[,4], weights)
                }
            }  # bootstrap loop ends

            tau.coef <- (d1 + d0 + z1 + z0)/2
            n0 <- d0/tau.coef
            n1 <- d1/tau.coef
            d.avg <- (d1 + d0)/2
            z.avg <- (z1 + z0)/2
            n.avg <- (n0 + n1)/2

            tau <- (delta.1 + delta.0 + zeta.1 + zeta.0)/2
            nu.0 <- delta.0/tau
            nu.1 <- delta.1/tau
            delta.avg <- (delta.0 + delta.1)/2
            zeta.avg <- (zeta.0 + zeta.1)/2
            nu.avg <- (nu.0 + nu.1)/2

        }  # nonpara boot branch ends

        ########################################################################
        ## Compute Outputs and Put Them Together
        ########################################################################

        low <- (1 - conf.level)/2
        high <- 1 - low
        d0.ci <- quantile(delta.0, c(low,high), na.rm=TRUE)
        d1.ci <- quantile(delta.1, c(low,high), na.rm=TRUE)
        tau.ci <- quantile(tau, c(low,high), na.rm=TRUE)
        z1.ci <- quantile(zeta.1, c(low,high), na.rm=TRUE)
        z0.ci <- quantile(zeta.0, c(low,high), na.rm=TRUE)
        n0.ci <- quantile(nu.0, c(low,high), na.rm=TRUE)
        n1.ci <- quantile(nu.1, c(low,high), na.rm=TRUE)
        d.avg.ci <- quantile(delta.avg, c(low,high), na.rm=TRUE)
        z.avg.ci <- quantile(zeta.avg, c(low,high), na.rm=TRUE)
        n.avg.ci <- quantile(nu.avg, c(low,high), na.rm=TRUE)

        # p-values
        d0.p <- 2 * sum(sign(delta.0) != sign(median(delta.0)))/sims
        d1.p <- 2 * sum(sign(delta.1) != sign(median(delta.1)))/sims
        d.avg.p <- 2 * sum(sign(delta.avg) != sign(median(delta.avg)))/sims
        z0.p <- 2 * sum(sign(zeta.0) != sign(median(zeta.0)))/sims
        z1.p <- 2 * sum(sign(zeta.1) != sign(median(zeta.1)))/sims
        z.avg.p <- 2 * sum(sign(zeta.avg) != sign(median(zeta.avg)))/sims
        n0.p <- 2 * sum(sign(nu.0) != sign(median(nu.0)))/sims
        n1.p <- 2 * sum(sign(nu.1) != sign(median(nu.1)))/sims
        n.avg.p <- 2 * sum(sign(nu.avg) != sign(median(nu.avg)))/sims
        tau.p <- 2 * sum(sign(tau) != sign(median(tau)))/sims

        # Detect whether models include T-M interaction
        INT <- paste(treat,mediator,sep=":") %in% attr(terms(model.y),"term.labels") |
               paste(mediator,treat,sep=":") %in% attr(terms(model.y),"term.labels")
        if(!INT & isGam.y){
            INT <- !isTRUE(all.equal(d0, d1))  # if gam, determine empirically
        }

        if(long) {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
            			d0.p=d0.p, d1.p=d1.p,
                        d0.sims=delta.0, d1.sims=delta.1,
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
            			z0.p=z0.p, z1.p=z1.p,
                        z0.sims=zeta.0, z1.sims=zeta.1,
                        n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
            			n0.p=n0.p, n1.p=n1.p,
                        n0.sims=nu.0, n1.sims=nu.1,
                        tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                        tau.sims=tau,
                        d.avg=d.avg, d.avg.p=d.avg.p, d.avg.ci=d.avg.ci, d.avg.sims=delta.avg,
                        z.avg=z.avg, z.avg.p=z.avg.p, z.avg.ci=z.avg.ci, z.avg.sims=zeta.avg,
                        n.avg=n.avg, n.avg.p=n.avg.p, n.avg.ci=n.avg.ci, n.avg.sims=nu.avg,
                        boot=boot, treat=treat, mediator=mediator,
                        covariates=covariates,
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m,
                        control.value=control.value, treat.value=treat.value,
                        nobs=n, sims=sims)
        } else {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
            			d0.p=d0.p, d1.p=d1.p,
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
            			z0.p=z0.p, z1.p=z1.p,
                        n0=n0, n1=n1, n0.ci=n0.ci, n1.ci=n1.ci,
            			n0.p=n0.p, n1.p=n1.p,
                        tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                        d.avg=d.avg, d.avg.p=d.avg.p, d.avg.ci=d.avg.ci,
                        z.avg=z.avg, z.avg.p=z.avg.p, z.avg.ci=z.avg.ci,
                        n.avg=n.avg, n.avg.p=n.avg.p, n.avg.ci=n.avg.ci,
                        boot=boot, treat=treat, mediator=mediator,
                        covariates=covariates,
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m,
                        control.value=control.value, treat.value=treat.value,
                        nobs=n, sims=sims)
        }
        class(out) <- "mediate"
        out

    ############################################################################
    ############################################################################
    ### CASE II: ORDERED OUTCOME
    ############################################################################
    ############################################################################
    } else {
        if(boot != TRUE){
            warning("ordered outcome model can only be used with nonparametric bootstrap - option forced")
            boot <- TRUE
        }

        n.ycat <- length(unique(model.response(y.data)))

        # Storage
        delta.1 <- matrix(NA, sims, n.ycat)
        delta.0 <- matrix(NA, sims, n.ycat)
        zeta.1 <- matrix(NA, sims, n.ycat)
        zeta.0 <- matrix(NA, sims, n.ycat)
        tau <- matrix(NA, sims, n.ycat)

        # Bootstrap loop begins
        for(b in 1:(sims+1)){

            # Resampling Step
            index <- sample(1:n, n, replace = TRUE)
			if(b == sims + 1){  # use original data for the last iteration
				index <- 1:n
			}
			
            Call.M <- model.m$call
            Call.Y <- model.y$call

            if(isSurvreg.m){
            	if(ncol(model.m$y) > 2){
                	stop("unsupported censoring type")
                }
                mname <- names(m.data)[1]
                if(substr(mname, 1, 4) != "Surv"){
                	stop("refit the survival model with `Surv' used directly in model formula")
                }
                nc <- nchar(mediator)
                eventname <- substr(mname, 5 + nc + 3, nchar(mname) - 1)
                if(nchar(eventname) == 0){
	               	m.data.tmp <- data.frame(m.data,
    	           							as.numeric(m.data[,1L][,1L]))
            	   	names(m.data.tmp)[c(1L, ncol(m.data)+1)] <- c(mname, mediator)
                } else {
	               	m.data.tmp <- data.frame(m.data,
    	           							as.numeric(m.data[,1L][,1L]),
        	       							as.numeric(model.m$y[,2]))
            	   	names(m.data.tmp)[c(1L, ncol(m.data)+(1:2))] <- c(mname, mediator, eventname)
                }
                Call.M$data <- m.data.tmp[index,]
            } else {
                Call.M$data <- m.data[index,]
            }

            Call.Y$data <- y.data[index,]
            Call.M$weights <- m.data[index,"(weights)"]
            Call.Y$weights  <- y.data[index,"(weights)"]
            new.fit.M <- eval.parent(Call.M)
            new.fit.Y <- eval.parent(Call.Y)

            if(isOrdered.m && length(unique(y.data[index,mediator]))!=m){
                # Modify the coefficients when mediator has empty cells
                coefnames.y <- names(model.y$coefficients)
                coefnames.new.y <- names(new.fit.Y$coefficients)
                new.fit.Y.coef <- rep(0, length(coefnames.y))
                names(new.fit.Y.coef) <- coefnames.y
                new.fit.Y.coef[coefnames.new.y] <- new.fit.Y$coefficients
                new.fit.Y$coefficients <- new.fit.Y.coef
            }

            #####################################
            # Mediator Predictions
            #####################################
            pred.data.t <- pred.data.c <- m.data

            if(isFactorT){
                 pred.data.t[,treat] <- factor(cat.1, levels = t.levels)
                 pred.data.c[,treat] <- factor(cat.0, levels = t.levels)
            } else {
                 pred.data.t[,treat] <- cat.1
                 pred.data.c[,treat] <- cat.0
            }

            if(!is.null(covariates)){
   	        	for(p in 1:length(covariates)){
       	    		vl <- names(covariates[p])
           			x <- ifelse(is.factor(pred.data.t[,vl]),
           							factor(covariates[[p]], levels = levels(m.data[,vl])),
           							covariates[[p]])
 	          		pred.data.t[,vl] <- pred.data.c[,vl] <- x
       	    	}
       		}

            ### Case II-a: GLM Mediator (including GAMs)
            if(isGlm.m){

                muM1 <- predict(new.fit.M, type="response", newdata=pred.data.t)
                muM0 <- predict(new.fit.M, type="response", newdata=pred.data.c)

                if(FamilyM == "poisson"){
                    PredictM1 <- rpois(n, lambda = muM1)
                    PredictM0 <- rpois(n, lambda = muM0)
                } else if (FamilyM == "Gamma") {
                    shape <- gamma.shape(model.m)$alpha
                    PredictM1 <- rgamma(n, shape = shape, scale = muM1/shape)
                    PredictM0 <- rgamma(n, shape = shape, scale = muM0/shape)
                } else if (FamilyM == "binomial"){
                    PredictM1 <- rbinom(n, size = 1, prob = muM1)
                    PredictM0 <- rbinom(n, size = 1, prob = muM0)
                } else if (FamilyM == "gaussian"){
                    sigma <- sqrt(summary(model.m)$dispersion)
                    error <- rnorm(n, mean=0, sd=sigma)
                    PredictM1 <- muM1 + error
                    PredictM0 <- muM0 + error
                } else if (FamilyM == "inverse.gaussian"){
                    disp <- summary(model.m)$dispersion
                    PredictM1 <- SuppDists::rinvGauss(n, nu = muM1, lambda = 1/disp)
                    PredictM0 <- SuppDists::rinvGauss(n, nu = muM0, lambda = 1/disp)
                } else {
                    stop("unsupported glm family")
                }

            ### Case II-b: Ordered Mediator
            } else if(isOrdered.m) {
                probs_m1 <- predict(new.fit.M, type="probs", newdata=pred.data.t)
                probs_m0 <- predict(new.fit.M, type="probs", newdata=pred.data.c)
                draws_m1 <- matrix(NA, n, m)
                draws_m0 <- matrix(NA, n, m)

                for(ii in 1:n){
                    draws_m1[ii,] <- t(rmultinom(1, 1, prob = probs_m1[ii,]))
                    draws_m0[ii,] <- t(rmultinom(1, 1, prob = probs_m0[ii,]))
                }
                PredictM1 <- apply(draws_m1, 1, indexmax)
                PredictM0 <- apply(draws_m0, 1, indexmax)

            ### Case II-c: Quantile Regression for Mediator
            } else if(isRq.m){
                # Use inverse transform sampling to predict M
                call.new <- new.fit.M$call
                call.new$tau <- runif(n)
                newfits <- eval.parent(call.new)
                tt <- delete.response(terms(new.fit.M))
                m.t <- model.frame(tt, pred.data.t, xlev = new.fit.M$xlevels)
                m.c <- model.frame(tt, pred.data.c, xlev = new.fit.M$xlevels)
                X.t <- model.matrix(tt, m.t, contrasts = new.fit.M$contrasts)
                X.c <- model.matrix(tt, m.c, contrasts = new.fit.M$contrasts)
                rm(tt, m.t, m.c)
                PredictM1 <- rowSums(X.t * t(newfits$coefficients))
                PredictM0 <- rowSums(X.c * t(newfits$coefficients))
                rm(newfits, X.t, X.c)

            ### Case II-d: Linear
            } else if(isLm.m){
                sigma <- summary(new.fit.M)$sigma
                error <- rnorm(n, mean=0, sd=sigma)
                PredictM1 <- predict(new.fit.M, type="response",
                                      newdata=pred.data.t) + error
                PredictM0 <- predict(new.fit.M, type="response",
                                      newdata=pred.data.c) + error
                rm(error)

	        ### Case I-2-e: Survreg
   		    } else if(isSurvreg.m){
       			dd <- survreg.distributions[[new.fit.M$dist]]
				if (is.null(dd$itrans)){
					itrans <- function(x) x
				} else {
					itrans <- dd$itrans
				}
				dname <- dd$dist
				if(is.null(dname)){
					dname <- new.fit.M$dist
				}
				scale <- new.fit.M$scale
				lpM1 <- predict(new.fit.M, newdata=pred.data.t, type="linear")
				lpM0 <- predict(new.fit.M, newdata=pred.data.c, type="linear")
				error <- switch(dname,
								extreme = log(rweibull(n, shape=1, scale=1)),
								gaussian = rnorm(n),
								logistic = rlogis(n),
								t = rt(n, df=dd$parms))
				PredictM1 <- as.numeric(itrans(lpM1 + scale * error))
				PredictM0 <- as.numeric(itrans(lpM0 + scale * error))
				rm(error)

            } else {
                stop("mediator model is not yet implemented")
            }

            #####################################
            #  Outcome Predictions
            #####################################
            effects.tmp <- array(NA, dim = c(n, n.ycat, 4))
            for(e in 1:4){
                tt <- switch(e, c(1,1,1,0), c(0,0,1,0), c(1,0,1,1), c(1,0,0,0))
                pred.data.t <- pred.data.c <- y.data

	            if(!is.null(covariates)){
	        		for(p in 1:length(covariates)){
       	    			vl <- names(covariates[p])
           				x <- ifelse(is.factor(pred.data.t[,vl]),
          								factor(covariates[[p]], levels = levels(y.data[,vl])),
           								covariates[[p]])
       	    			pred.data.t[,vl] <- pred.data.c[,vl] <- x
           			}
           		}

                # Set treatment values
                cat.t <- ifelse(tt[1], cat.1, cat.0)
                cat.c <- ifelse(tt[2], cat.1, cat.0)
                cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
                cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
                if(isFactorT){
                    pred.data.t[,treat] <- factor(cat.t, levels = t.levels)
                    pred.data.c[,treat] <- factor(cat.c, levels = t.levels)
                    if(!is.null(control)){
                        pred.data.t[,control] <- factor(cat.t.ctrl, levels = t.levels)
                        pred.data.c[,control] <- factor(cat.c.ctrl, levels = t.levels)
                    }
                } else {
                    pred.data.t[,treat] <- cat.t
                    pred.data.c[,treat] <- cat.c
                    if(!is.null(control)){
                        pred.data.t[,control] <- cat.t.ctrl
                        pred.data.c[,control] <- cat.c.ctrl
                    }
                }

                # Set mediator values
                PredictM1.tmp <- PredictM1
                PredictM0.tmp <- PredictM0
                PredictMt <- PredictM1 * tt[3] + PredictM0 * (1 - tt[3])
                PredictMc <- PredictM1 * tt[4] + PredictM0 * (1 - tt[4])
                if(isFactorM) {
                    pred.data.t[,mediator] <- factor(PredictMt, levels=1:m, labels=m.levels)
                    pred.data.c[,mediator] <- factor(PredictMc, levels=1:m, labels=m.levels)
                } else {
                    pred.data.t[,mediator] <- PredictMt
                    pred.data.c[,mediator] <- PredictMc
                }
                probs_p1 <- predict(new.fit.Y, newdata=pred.data.t, type="probs")
                probs_p0 <- predict(new.fit.Y, newdata=pred.data.c, type="probs")
                effects.tmp[,,e] <- probs_p1 - probs_p0
                rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
            }

            # Compute all QoIs
            if(b == sims+1){
	            d1 <- apply(effects.tmp[,,1], 2, weighted.mean, w=weights)
    	        d0 <- apply(effects.tmp[,,2], 2, weighted.mean, w=weights)
        	    z1 <- apply(effects.tmp[,,3], 2, weighted.mean, w=weights)
            	z0 <- apply(effects.tmp[,,4], 2, weighted.mean, w=weights)
            } else {
	            delta.1[b,] <- apply(effects.tmp[,,1], 2, weighted.mean, w=weights)
    	        delta.0[b,] <- apply(effects.tmp[,,2], 2, weighted.mean, w=weights)
        	    zeta.1[b,] <- apply(effects.tmp[,,3], 2, weighted.mean, w=weights)
            	zeta.0[b,] <- apply(effects.tmp[,,4], 2, weighted.mean, w=weights)
            }

        }  # Bootstrap loop ends

        tau.coef <- (d1 + d0 + z1 + z0)/2
  		tau <- (zeta.1 + zeta.0 + delta.0 + delta.1)/2

        ########################################################################
        ## Compute Outputs and Put Them Together
        ########################################################################
        low <- (1 - conf.level)/2
        high <- 1 - low
        d0.ci <- apply(delta.0, 2, quantile, c(low,high))
        d1.ci <- apply(delta.1, 2, quantile, c(low,high))
        tau.ci <- apply(tau, 2, quantile, c(low,high))
        z1.ci <- apply(zeta.1, 2, quantile, c(low,high))
        z0.ci <- apply(zeta.0, 2, quantile, c(low,high))

        # p-values
        d0.p <- 2 * apply(delta.0, 2, function(x) sum(sign(x) != sign(median(x)))/sims)
        d1.p <- 2 * apply(delta.1, 2, function(x) sum(sign(x) != sign(median(x)))/sims)
        z0.p <- 2 * apply(zeta.0, 2, function(x) sum(sign(x) != sign(median(x)))/sims)
        z1.p <- 2 * apply(zeta.1, 2, function(x) sum(sign(x) != sign(median(x)))/sims)
        tau.p <- 2 * apply(tau, 2, function(x) sum(sign(x) != sign(median(x)))/sims)

        # Detect whether models include T-M interaction
        INT <- paste(treat,mediator,sep=":") %in% attr(model.y$terms,"term.labels") |
             paste(mediator,treat,sep=":") %in% attr(model.y$terms,"term.labels")

        if(long) {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
            			d0.p=d0.p, d1.p=d1.p,
                        d0.sims=delta.0, d1.sims=delta.1,
                        tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
            			z0.p=z0.p, z1.p=z1.p,
                        z1.sims=zeta.1, z0.sims=zeta.0, tau.sims=tau,
                        boot=boot, treat=treat, mediator=mediator,
                        covariates=covariates,
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m,
                        control.value=control.value, treat.value=treat.value, nobs=n, sims=sims)
        } else {
            out <- list(d0=d0, d1=d1, d0.ci=d0.ci, d1.ci=d1.ci,
            			d0.p=d0.p, d1.p=d1.p,
                        tau.coef=tau.coef, tau.ci=tau.ci, tau.p=tau.p,
                        z0=z0, z1=z1, z0.ci=z0.ci, z1.ci=z1.ci,
            			z0.p=z0.p, z1.p=z1.p,
                        boot=boot, treat=treat, mediator=mediator,
                        covariates=covariates,
                        INT=INT, conf.level=conf.level,
                        model.y=model.y, model.m=model.m,
                        control.value=control.value, treat.value=treat.value, nobs=n, sims=sims)
        }
        class(out) <- "mediate.order"
        out
    }
}



summary.mediate <- function(object, ...){
    structure(object, class = c("summary.mediate", class(object)))
}



print.summary.mediate <- function(x, ...){
    clp <- 100 * x$conf.level
    cat("\n")
    cat("Causal Mediation Analysis \n\n")
    if(x$boot){
        cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
    } else {
        cat("Quasi-Bayesian Confidence Intervals\n\n")
    }

    if(!is.null(x$covariates)){
    	cat("(Inference Conditional on the Covariate Values Specified in `covariates')\n\n")
    }

    isLinear.y <- (	(class(x$model.y)[1] %in% c("lm", "rq")) ||
        			(inherits(x$model.y, "glm") &&
        				x$model.y$family$family == "gaussian" &&
        				x$model.y$family$link == "identity") ||
				    (inherits(x$model.y, "survreg") &&
				       	x$model.y$dist == "gaussian") )

    printone <- !x$INT && isLinear.y

    if (printone){
        # Print only one set of values if lmY/quanY/linear gamY without interaction
        smat <- c(x$d1, x$d1.ci, x$d1.p)
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        rownames(smat) <- c("Mediation Effect", "Direct Effect",
        					"Total Effect", "Proportion via Mediation")
    } else {
        smat <- c(x$d0, x$d0.ci, x$d0.p)
        smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
        smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
        smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
        smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
        rownames(smat) <- c("Mediation Effect_0", "Mediation Effect_1",
        					"Direct Effect_0", "Direct Effect_1",
        					"Total Effect",
        					"Proportion via Mediation_0",
        					"Proportion via Mediation_1",
        					"Mediation Effect (Ave.)",
        					"Direct Effect (Ave.)",
        					"Proportion via Mediation (Ave.)")
    }
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep=""),
    					paste(clp, "% CI Upper", sep=""), "p-value")
    printCoefmat(smat, digits=3)
    cat("\n")
    cat("Sample Size Used:", x$nobs,"\n\n")
    cat("\n")
    cat("Simulations:", x$sims,"\n\n")
    invisible(x)
}



summary.mediate.order <- function(object, ...){
    structure(object, class = c("summary.mediate.order", class(object)))
}



print.summary.mediate.order <- function(x, ...){
    tab.d0 <- rbind(x$d0, x$d0.ci, x$d0.p)
    tab.d1 <- rbind(x$d1, x$d1.ci, x$d1.p)
    tab.z0 <- rbind(x$z0, x$z0.ci, x$z0.p)
    tab.z1 <- rbind(x$z1, x$z1.ci, x$z1.p)
    tab.tau <- rbind(x$tau.coef, x$tau.ci, x$tau.p)

    # Outcome Table Labels
    y.lab <- sort(unique(levels(model.frame(x$model.y)[,1])))
    out.names <- c()
    for(i in 1:length(y.lab)){
        out.names.tmp <- paste("Pr(Y=",y.lab[i],")",sep="")
        out.names <- c(out.names, out.names.tmp)
    }

    # Label Tables
    rownames(tab.d0)[1] <- "Mediation Effect_0 "
    rownames(tab.d0)[4] <- "p-value "
    colnames(tab.d0) <- out.names
    rownames(tab.d1)[1] <- "Mediation Effect_1 "
    rownames(tab.d1)[4] <- "p-value "
    colnames(tab.d1) <- out.names
    rownames(tab.z0)[1] <- "Direct Effect_0    "
    rownames(tab.z0)[4] <- "p-value "
    colnames(tab.z0) <- out.names
    rownames(tab.z1)[1] <- "Direct Effect_1    "
    rownames(tab.z1)[4] <- "p-value "
    colnames(tab.z1) <- out.names
    rownames(tab.tau)[1] <- "Total Effect      "
    rownames(tab.tau)[4] <- "p-value "
    colnames(tab.tau) <- out.names

    cat("\n")
    cat("Causal Mediation Analysis \n\n")
    cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
    if(!is.null(x$covariates)){
    	cat("(Inference Conditional on the Covariate Values Specified in `covariates')\n\n")
    }
    print(tab.d0, digits=3)
    cat("\n")
    print(tab.d1, digits=3)
    cat("\n")
    print(tab.z0, digits=3)
    cat("\n")
    print(tab.z1, digits=3)
    cat("\n")
    print(tab.tau, digits=3)
    cat("\n\n")
    cat("Sample Size Used:", x$nobs,"\n\n")
    cat("\n\n")
    cat("Simulations:", x$sims,"\n\n")
    invisible(x)
}



plot.process <- function(model) {
    coef.vec.1 <- c(model$d1, model$z1)
    lower.vec.1 <- c(model$d1.ci[1], model$z1.ci[1])
    upper.vec.1 <- c(model$d1.ci[2], model$z1.ci[2])
    tau.vec <- c(model$tau.coef,model$tau.ci[1],model$tau.ci[2])
    range.1 <- range(model$d1.ci[1], model$z1.ci[1],model$tau.ci[1],
                      model$d1.ci[2], model$z1.ci[2],model$tau.ci[2])

    coef.vec.0 <- c(model$d0, model$z0)
    lower.vec.0 <- c(model$d0.ci[1], model$z0.ci[1])
    upper.vec.0 <- c(model$d0.ci[2], model$z0.ci[2])
    range.0 <- range(model$d0.ci[1], model$z0.ci[1],model$tau.ci[1],
                      model$d0.ci[2], model$z0.ci[2],model$tau.ci[2])

    return(list(coef.vec.1=coef.vec.1, lower.vec.1=lower.vec.1,
                upper.vec.1=upper.vec.1, coef.vec.0=coef.vec.0,
                lower.vec.0=lower.vec.0, upper.vec.0=upper.vec.0, tau.vec=tau.vec,
                range.1=range.1, range.0=range.0))
}



plot.mediate <- function(x, treatment = NULL,
                        labels = c("ACME","Direct\nEffect","Total\nEffect"),
                        xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                        main = NULL, lwd = 1.5, cex = .85,
                        col = "black", ...){
    # Determine which graph to plot
    if(is.null(treatment)){
        if(x$INT){
            treatment <- c(0,1)
        } else {
            treatment <- 1
        }
    } else {
        treatment <- switch(treatment,
                                control = 0,
                                treated = 1,
                                both = c(0,1))
    }

    param <- plot.process(x)
    y.axis <- c(length(param$coef.vec.1):.5)
    y.axis <- y.axis + 1
        # create indicator for y.axis, descending so labels go from top to bottom

    # Set xlim
    if(is.null(xlim)){
        if(length(treatment) > 1) {
            xlim <- range(param$range.1, param$range.0) * 1.2
        } else if (treatment == 1){
            xlim <- param$range.1 * 1.2
        } else {
            xlim <- param$range.0 * 1.2
        }
    }

    # Set ylim
    if(is.null(ylim)){
        ylim <- c(min(y.axis) -1- 0.5, max(y.axis) + 0.5)
    }

    # Plot
    plot(param$coef.vec.1, y.axis, type = "n", xlab = xlab, ylab = ylab,
            yaxt = "n", xlim = xlim, ylim = ylim, main = main, ...)

    # Set offset values depending on number of bars to plot
    if(length(treatment) == 1){
        adj <- 0
    } else {
        adj <- 0.05
    }

    if(1 %in% treatment){
        points(param$coef.vec.1, y.axis + adj, type = "p", pch = 19, cex = cex, col = col)
        segments(param$lower.vec.1, y.axis + adj, param$upper.vec.1, y.axis + adj,
                lwd = lwd, col = col)
        points(param$tau.vec[1], 1, type = "p", pch = 19, cex = cex, col = col)
        segments(param$tau.vec[2], 1 , param$tau.vec[3], 1 ,
                lwd = lwd, col = col)
    }
    if(0 %in% treatment) {
        points(param$coef.vec.0, y.axis - adj, type = "p", pch = 1, cex = cex, col = col)
        segments(param$lower.vec.0, y.axis - adj, param$upper.vec.0, y.axis - adj,
                lwd = lwd, lty = 3, col = col)
    }
    if(treatment[1]==0 & length(treatment)==1) {
        points(param$tau.vec[1], 1 , type = "p", pch = 19, cex = cex, col = col)
        segments(param$tau.vec[2], 1 , param$tau.vec[3], 1 ,
                lwd = lwd, col = col)
    }
    y.axis.new <- c(3,2,1)
    axis(2, at = y.axis.new, labels = labels, las = 1, tick = TRUE, ...)
    abline(v = 0, lty = 2)
}



plot.process.order <- function(model){
    length <- length(model$d1)
    coef.vec.1 <- lower.vec.1 <- upper.vec.1 <-
        coef.vec.0 <- lower.vec.0 <- upper.vec.0 <- matrix(NA,ncol=2,nrow=length)
    tau.vec<-matrix(NA,ncol=3,nrow=length)
    for(j in 1:length){
        coef.vec.1[j,] <- c(model$d1[j], model$z1[j])
        lower.vec.1[j,] <- c(model$d1.ci[1,j], model$z1.ci[1,j])
        upper.vec.1[j,] <- c(model$d1.ci[2,j], model$z1.ci[2,j])

        coef.vec.0[j,] <- c(model$d0[j], model$z0[j])
        lower.vec.0[j,] <- c(model$d0.ci[1,j], model$z0.ci[1,j])
        upper.vec.0[j,] <- c(model$d0.ci[2,j], model$z0.ci[2,j])

        tau.vec[j,] <- c(model$tau.coef[j], model$tau.ci[1,j], model$tau.ci[2,j])

    }

    range.1 <- range(model$d1.ci[1,], model$z1.ci[1,],model$tau.ci[1,],
                      model$d1.ci[2,], model$z1.ci[2,],model$tau.ci[2,])
    range.0 <- range(model$d0.ci[1,], model$z0.ci[1,],model$tau.ci[1,],
                      model$d0.ci[2,], model$z0.ci[2,],model$tau.ci[2,])

    return(list(coef.vec.1=coef.vec.1, lower.vec.1=lower.vec.1,
                upper.vec.1=upper.vec.1, coef.vec.0=coef.vec.0,
                lower.vec.0=lower.vec.0, upper.vec.0=upper.vec.0,
                tau.vec=tau.vec,
                range.1=range.1, range.0=range.0, length=length))
}



plot.mediate.order <- function(x, treatment = NULL,
                        labels = c("ACME","Direct\nEffect","Total\nEffect"),
                        xlim = NULL, ylim = NULL, xlab = "", ylab = "",
                        main = NULL, lwd = 1.5, cex = .85,
                        col = "black", ...){
    # Determine which graph to plot
    if(is.null(treatment)){
        if(x$INT){
            treatment <- c(0,1)
        } else {
            treatment <- 1
        }
    } else {
        treatment <- switch(treatment,
                                control = 0,
                                treated = 1,
                                both = c(0,1))
    }

    param <- plot.process.order(x)
    y.axis <- c(ncol(param$coef.vec.1):.5)
    y.axis <- y.axis + 1
    # create indicator for y.axis, descending so labels go from top to bottom

    # Set xlim
    if(is.null(xlim)){
        if(length(treatment) > 1) {
            xlim <- range(param$range.1, param$range.0) * 1.2
        } else if (treatment == 1){
            xlim <- param$range.1 * 1.2
        } else {
            xlim <- param$range.0 * 1.2
        }
    }

    # Set ylim
    if(is.null(ylim)){
        ylim <- c(min(y.axis) - 1 - 0.5, max(y.axis) + 0.5)
    }

    # Plot
    plot(param$coef.vec.1[1,], y.axis, type = "n", xlab = xlab, ylab = ylab,
            yaxt = "n", xlim = xlim, ylim = ylim, main = main, ...)

    # Set offset values depending on number of bars to plot
    if(length(treatment) == 1){
        adj <- 0
    } else {
        adj <- 0.05
    }

    if(1 %in% treatment){
        adj.1 <- adj * nrow(param$coef.vec.1)
        for(z in 1:nrow(param$coef.vec.1)){
            points(param$coef.vec.1[z,], y.axis + adj.1,
                    type = "p", pch = 19, cex = cex, col = col)
            segments(param$lower.vec.1[z,], y.axis + adj.1,
                    param$upper.vec.1[z,], y.axis + adj.1,
                    lwd = lwd, col = col)
            points(param$tau.vec[z,1], 1 + adj.1 ,
                    type = "p", pch = 19, cex = cex, col = col)
            segments(param$tau.vec[z,2], 1 + adj.1 ,
                    param$tau.vec[z,3], 1 + adj.1 ,
                    lwd = lwd, col = col)
            adj.1 <- adj.1 - 0.05
        }

    }
    if(0 %in% treatment) {
        adj.0 <- adj
        for(z in 1:nrow(param$coef.vec.0)){
            points(param$coef.vec.0[z,], y.axis - adj.0,
                    type = "p", pch = 1, cex = cex, col = col)
            segments(param$lower.vec.0[z,], y.axis - adj.0,
                    param$upper.vec.0[z,], y.axis - adj.0,
                    lwd = lwd, lty = 3, col = col)
            adj.0 <- adj.0 + 0.05
        }
    }
        if (treatment[1]==0 & length(treatment)==1){
        print("test")
        adj.1 <- adj * nrow(param$coef.vec.1)
        for(z in 1:nrow(param$tau.vec)){
            points(param$tau.vec[z,1], 1 + adj.1 ,
                    type = "p", pch = 19, cex = cex, col = col)
            segments(param$tau.vec[z,2], 1 + adj.1 ,
                    param$tau.vec[z,3], 1 +adj.1 ,
                    lwd = lwd, col = col)
                    adj.1 <- adj.1 - 0.05
                    }
            }

    y.axis.new <- c(3,2,1)
    axis(2, at = y.axis.new, labels = labels, las = 1, tick = TRUE, ...)
    abline(v = 0, lty = 2)
}
