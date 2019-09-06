#' Sensitivity Analysis for Causal Mediation Effects
#' 
#' 'medsens' is used to perform sensitivity analysis on the average causal 
#' mediation effects and direct effects for violations of the sequential 
#' ignorability assumption.  The function takes output from '\code{mediate}' and
#' calculates the true average causal mediation effects and direct effects for
#' different values of the sensitivity parameter representing the degree of the
#' sequential ignorability violation.
#' 
#' @details This is the workhorse function for sensitivity analyses for average
#'   causal mediation effects. The sensitivity analysis can be used to assess
#'   the robustness of the findings from \code{mediate} to the violation of 
#'   sequential ignorability, the crucial identification assumption necessary
#'   for the estimates to be valid. The analysis proceeds by quantifying the
#'   degree of sequential ignorability violation as the correlation between the
#'   error terms of the mediator and outcome models, and then calculating the
#'   true values of the average causal mediation effect for given values of this
#'   sensitivity parameter, rho. The original findings are deemed sensitive if 
#'   the true effects are found to vary widely as function of rho.
#'   
#'   The sensitivity analysis is only implemented for the following three model 
#'   combinations: linear mediator and outcome models (both of class 'lm'), 
#'   binary probit mediator (fitted via 'glm' with family "binomial" and link 
#'   "probit") and linear outcome models, and linear mediator and binary probit 
#'   outcome models.  In addition, the binary outcome model cannot include a 
#'   treatment-mediator interaction term. An error is returned if the 'mediate' 
#'   object in 'x' is based on other model combinations.  As of version 3.0, the
#'   sensitivity analysis can also be conducted with respect to the average 
#'   direct effect by setting 'effect.type' to "direct" (or "both" if results
#'   for the average causal mediation effect are also desired).
#'   
#'   Users should note that computation can take significant time for 
#'   \code{medsens}. Setting 'rho.by' to a larger number significantly decreases
#'   computational time, as does decreasing 'eps' (for the linear-linear case)
#'   or the number of simulations 'sims' (for the binary-linear and
#'   linear-binary cases).
#'   
#' @param x an object of class 'mediate', typically an output from the 
#'   \code{mediate} function.
#' @param rho.by a numeric value between 0 and 1 indicating the increment for 
#'   the sensitivity parameter, rho.
#' @param sims the number of Monte Carlo draws for the calculation of confidence
#'   intervals. Only used in cases where either the mediator or outcome variable
#'   is binary.
#' @param eps convergence tolerance parameter for the iterative FGLS. Only used 
#'   when both the mediator and outcome models are linear.
#' @param effect.type a character string indicating which effect(s) to be 
#'   analyzed.  Default is "indirect".
#'   
#' @return \code{medsens} returns an object of class "\code{medsens}", a list 
#'   containing the following elements. Some of these elements are not available
#'   depending on the 'effect.type' argument specified by the user. The output 
#'   can then be passed to the \code{\link{summary}} (i.e., 
#'   \code{\link{summary.medsens}}) and \code{\link{plot}} (i.e., 
#'   \code{\link{plot.medsens}}) functions to produce tabular and graphical 
#'   summaries of the results.
#'   
#'   \item{d0, d1}{vectors of point estimates for average causal mediation 
#'   effects under the control and treatment conditions for each value of 
#'   sensitivity parameter rho.}
#'   \item{upper.d0, lower.d0, upper.d1, lower.d1}{vectors of upper and lower 
#'   confidence limits for average causal mediation effect under the control 
#'   and treatment conditions for each value of rho.}
#'   \item{z0, z1}{vectors of point estimates for average direct effect under 
#'   the control and treatment conditions for each value of sensitivity 
#'   parameter rho.}
#'   \item{upper.z0, lower.z0, upper.z1, lower.z1}{vectors of upper and lower 
#'   confidence limits for average direct effect under the control and 
#'   treatment conditions for each value of rho.}
#'   \item{tau}{a vector of point estimates for total effect for each value of 
#'   rho.  Only present when the outcome model is binary.}
#'   \item{upper.tau, lower.tau}{vectors of upper and lower confidence limits 
#'   for total effect.  Only present when the outcome model is binary.}
#'   \item{nu}{a vector of point estimates for the proportion mediated for each 
#'   value of rho.  Only present when the outcome model is binary.}
#'   \item{upper.nu, lower.nu}{vectors of upper and lower confidence limits 
#'   for the proportion mediated.  Only present when the outcome model is binary.}
#'   \item{rho}{a numeric vector containing the values of sensitivity 
#'   parameter rho used.}
#'   \item{rho.by}{a numeric value indicating the increment of rho used.}
#'   \item{sims}{a numeric value indicating the number of Monte Carlo draws 
#'   used.}
#'   \item{err.cr.d, err.cr.z}{the values of rho with which the average causal 
#'   mediation and direct effects are zero.  Vectors of length two if 'INT' is 
#'   'TRUE'; numeric values otherwise.}
#'   \item{ind.d0, ind.d1, ind.z0, ind.z1}{vectors of 0s/1s, indicating whether 
#'   the confidence intervals of d0, d1, z0 and z1 do not cover zero for each 
#'   value of rho.}
#'   \item{R2star.prod}{a numeric vector containing the values of the products 
#'   of the two "R square stars", representing the proportions of residual 
#'   variance in the mediator and outcome explained by the hypothesized 
#'   unobserved confounder. The values correspond to those of rho. See 
#'   \code{\link{plot.medsens}} for details.}
#'   \item{R2tilde.prod}{a numeric vector containing the values of the products 
#'   of the two "R square tildes", representing the proportions of total 
#'   variance in the mediator and outcome explained by the hypothesized 
#'   unobserved confounder. The values correspond to those of rho. See 
#'   \code{\link{plot.medsens}} for details.}
#'   \item{R2star.d.thresh, R2star.z.thresh}{the values of the product of "R 
#'   square stars" for which the average causal mediation and direct effects are 
#'   zero, respectively.}
#'   \item{R2tilde.d.thresh, R2tilde.z.thresh}{the values of the product of "R 
#'   square tildes" for which the average causal mediation and direct effects 
#'   are zero, respectively.}
#'   \item{r.square.y, r.square.m}{the usual R square statistics for the 
#'   outcome and mediator models.}
#'   \item{INT}{a logical value indicating whether interaction between the 
#'   treatment and mediator is allowed in the original mediate object.}
#'   \item{conf.level}{the confidence level used.}
#'   \item{effect.type}{the 'effect.type' argument used.}
#'   \item{type}{a character string indicating the type of the mediator and 
#'   outcome models used. Currently either "ct" (linear mediator and outcome 
#'   models), 'bm' (binary mediator and linear outcome models) or 'bo' (linear 
#'   mediator and binary outcome models).}
#'   \item{robustSE}{`TRUE' or `FALSE'.}
#'   \item{cluster}{the clusters used.}

#' @author Dustin Tingley, Harvard University, 
#'   \email{dtingley@@gov.harvard.edu}; Teppei Yamamoto, Massachusetts Institute
#'   of Technology, \email{teppei@@mit.edu}; Jaquilyn Waddell-Boie, Princeton 
#'   University, \email{jwaddell@@princeton.edu}; Kentaro Hirose, Princeton 
#'   University, \email{hirose@@princeton.edu}; Luke Keele, Penn State 
#'   University, \email{ljk20@@psu.edu}; Kosuke Imai, Princeton University, 
#'   \email{kimai@@princeton.edu}.
#'   
#' @seealso \code{\link{mediate}}, \code{\link{summary.medsens}}, 
#'   \code{\link{plot.medsens}}.
#'   
#' @references Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. 
#'   (2014). "mediation: R package for Causal Mediation Analysis", Journal of 
#'   Statistical Software, Vol. 59, No. 5, pp. 1-38.
#'   
#'   Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the 
#'   Black Box of Causality: Learning about Causal Mechanisms from Experimental 
#'   and Observational Studies, American Political Science Review, Vol. 105, No.
#'   4 (November), pp. 765-789.
#'   
#'   Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal 
#'   Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 
#'   309-334.
#'   
#'   Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and 
#'   Sensitivity Analysis for Causal Mediation Effects, Statistical Science,
#'   Vol. 25, No. 1 (February), pp. 51-71.
#'   
#'   Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation 
#'   Analysis Using R" in Advances in Social Science Research Using R, ed. H. D.
#'   Vinod New York: Springer.
#'   
#' @export
#' @examples
#' # Examples with JOBS II Field Experiment
#' 
#' # **For illustration purposes a small number of simulations are used**
#' 
#' data(jobs)
#' 
#' ####################################################
#' # Example 1: Binary treatment 
#' ####################################################
#' # Fit parametric models
#' b <- lm(job_seek ~ treat + econ_hard + sex + age, data=jobs)
#' c <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data=jobs)
#' 
#' # Pass model objects through mediate function
#' med.cont <- mediate(b, c, treat="treat", mediator="job_seek", sims=50)
#' # med.cont <- mediate(b, c, treat="treat", mediator="job_seek", sims=50, robustSE = T)
#' # jobs$cluster <- rep(1:30, each = 30)[-1]
#' # med.cont <- mediate(b, c, treat="treat", mediator="job_seek", sims=50, cluster = jobs$cluster)
#' 
#' # Pass mediate output through medsens function
#' sens.cont <- medsens(med.cont, rho.by=.1, eps=.01, effect.type="both")
#' 
#' # Use summary function to display results
#' summary(sens.cont)
#' 
#' # Plot true ACMEs and ADEs as functions of rho
#' par.orig <- par(mfrow = c(2,2))
#' plot(sens.cont, main="JOBS", ylim=c(-.2,.2))
#' 
#' # Plot true ACMEs and ADEs as functions of "R square tildes"
#' plot(sens.cont, sens.par="R2", r.type="total", sign.prod="positive")
#' par(par.orig)
#' 
#' ####################################################
#' # Example 2: Categorical treatment 
#' ####################################################
#' \dontrun{
#' 
#' # Purely for illustration, think of educ as a ``treatment''
#' 
#' b <- lm(job_seek ~ educ + sex, data=jobs)
#' c <- lm(depress2 ~ educ + job_seek + sex, data=jobs)
#' 
#' # compare two categories of educ --- gradwk and somcol
#' med.cont <- mediate(b, c, treat="educ", mediator="job_seek", sims=50, 
#'                     control.value = "gradwk", treat.value = "somcol")
#' sens.cont <- medsens(med.cont, rho.by=.1, eps=.01, effect.type="both")
#' summary(sens.cont)
#' }
#' 
medsens <- function(x, rho.by = 0.1, sims = 1000, eps = sqrt(.Machine$double.eps),
        effect.type = c("indirect", "direct", "both"))
{
    effect.type <- match.arg(effect.type)
    if(effect.type == "both"){
        etype.vec <- c("indirect", "direct")
    } else {
        etype.vec <- effect.type
    }
    
    model.y <- x$model.y
    model.m <- x$model.m
    class.y <- class(model.y)[1]
    class.m <- class(model.m)[1]
    treat <- x$treat
    mediator <- x$mediator
    INT <- x$INT
    low <- (1 - x$conf.level)/2
    high <- 1 - low
    robustSE <- x$robustSE
    cluster <- x$cluster
    control.value <- x$control.value
    treat.value <- x$treat.value

    # Setting Variable labels
    ## Uppercase letters (e.g. T) = labels in the input matrix
    ## Uppercase letters + ".out" (e.g. T.out) = labels in the regression output

    y.data <- model.frame(model.y)
    Y <- colnames(y.data)[1]
    if(is.factor(y.data[,treat])){
        t.levels <- levels(y.data[,treat])
        n.t.levels <- length(t.levels)
        ## When treatment = more than two categories 
        ## control.value => 1st factor level (baseline category)
        ## treat.value => 2nd factor level
        ## run regression again with this new data set
        if(n.t.levels != 2){
            number <- 1:n.t.levels
            t.levels.mat <- data.frame(number, t.levels)
            c.num <- t.levels.mat$number[t.levels.mat$t.levels == control.value]
            t.num <- t.levels.mat$number[t.levels.mat$t.levels == treat.value]            
            number.new <- c(c.num, t.num, number[-c(c.num, t.num)])
            y.data[, treat] <- factor(y.data[, treat], levels(y.data[, treat])[number.new])
            model.y <- update(model.y, data = y.data)     
        }
        t.levels <- levels(y.data[,treat])
        cat.c <- t.levels[1]
        cat.t <- t.levels[2]
    } else {
        cat.c <- NULL
        cat.t <- NULL
    }
    T.out <- paste(treat, cat.t, sep="")

    m.data <- model.frame(model.m)
    if(is.factor(m.data[,treat])){
        t.levels <- levels(m.data[,treat])
        n.t.levels <- length(t.levels)
        ## When treatment = more than two categories 
        ## control.value => 1st factor level (baseline category)
        ## treat.value => 2nd factor level
        ## run regression again with this new data set 
        if(n.t.levels != 2){
            number <- 1:n.t.levels
            t.levels.mat <- data.frame(number, t.levels)
            c.num <- t.levels.mat$number[t.levels.mat$t.levels == control.value]
            t.num <- t.levels.mat$number[t.levels.mat$t.levels == treat.value]            
            number.new <- c(c.num, t.num, number[-c(c.num, t.num)])
            m.data[, treat] <- factor(m.data[, treat], levels(m.data[, treat])[number.new])
            model.m <- update(model.m, data = m.data)     
        }
    } 
    
    if(is.factor(y.data[,mediator])){
        m.levels <- levels(y.data[,mediator])
        if(length(m.levels) != 2){
            stop("mediator with more than two categories currently not supported")
        }
        cat.m0 <- m.levels[1]
        cat.m1 <- m.levels[2]
    } else {
        cat.m0 <- NULL
        cat.m1 <- NULL
    }
    M.out <- paste(mediator, cat.m1, sep="")

    if(INT){
        if(paste(treat,mediator,sep=":") %in% attr(model.y$terms,"term.labels")){ # T:M
            TM.out <- paste(T.out, M.out, sep=":")
            TM <- paste(treat, mediator, sep=":")
        } else { # M:T
            TM.out <- paste(M.out, T.out, sep=":")
            TM <- paste(mediator, treat, sep=":")
        }
    }

    ## Setting Up Sensitivity Parameters
    rho <- seq(-1 + rho.by, 1 - rho.by, rho.by)
    R2star.prod <- rho^2
    
    ## Extracting weights
    weights <- model.weights(y.data)
    if(is.null(weights)){
        weights <- rep(1, nrow(y.data))
    }
    
    ## Initializing Containers
    d0 <- d1 <- z0 <- z1 <- nu0 <- nu1 <- tau <-
    upper.d0 <- upper.d1 <- upper.z0 <- upper.z1 <- upper.nu0 <- upper.nu1 <- upper.tau <-
    lower.d0 <- lower.d1 <- lower.z0 <- lower.z1 <- lower.nu0 <- lower.nu1 <- lower.tau <-
    ind.d0 <- ind.d1 <- ind.z0 <- ind.z1 <- err.cr.d <- err.cr.z <- 
    R2star.d.thresh <- R2tilde.d.thresh <- R2star.z.thresh <- R2tilde.z.thresh <- NULL

    getvcov <- function(dat, fm, cluster){
    ## Compute cluster robust standard errors
    ## fm is the model object
        cluster <- factor(cluster)  # remove missing levels and NA
        M <- nlevels(cluster)
        N <- sum(!is.na(cluster))
        K <- fm$rank
        dfc <- (M/(M-1))*((N-1)/(N-K))
        uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
        dfc*sandwich(fm, meat. = crossprod(uj)/N)
    }
    
    #########################################################
    ## CASE 1: Continuous Outcome + Continuous Mediator
    #########################################################
    if(class.y=="lm" & class.m=="lm") {
        type <- "ct"
        d0 <- d1 <- z0 <- z1 <- matrix(NA, length(rho), 1)
        d0.var <- d1.var <- z0.var <- z1.var <- matrix(NA, length(rho), 1)

        for(i in 1:length(rho)){
            e.cor <- rho[i]
            b.dif <- 1

            # Stacked Equations
            m.mat <- model.matrix(model.m) * weights^(1/2)
            y.mat <- model.matrix(model.y) * weights^(1/2)
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
            
            Y.c <- as.matrix(c(model.response(model.frame(model.m)), 
                                model.response(model.frame(model.y)))) * weights^(1/2)
            
            # Estimates of OLS Start Values
            inxx <- solve(crossprod(X))
            b.ols <- inxx %*% crossprod(X,Y.c)
            b.tmp <- b.ols

            while(abs(b.dif) > eps){
                e.hat <- as.matrix(Y.c - (X %*% b.tmp)) 
                e.1 <- e.hat[1:n]
                e.2 <- e.hat[(n+1):(2*n)]

                sd.1 <- sd(e.1)
                sd.2 <- sd(e.2)

                omega <- matrix(NA, 2,2)
                omega[1,1] <- crossprod(e.1)/(n-1)
                omega[2,2] <- crossprod(e.2)/(n-1)
                omega[2,1] <- e.cor*sd.1*sd.2
                omega[1,2] <- e.cor*sd.1*sd.2

                I <- Diagonal(n)
                omega.i <- solve(omega)
                v.i <-  kronecker(omega.i, I)
                Xv.i <- t(X) %*% v.i
                X.sur <- Xv.i %*% X
                b.sur <- solve(X.sur) %*% Xv.i %*% Y.c

                # Variance-Covariance Matrix of SUR 
                if(robustSE){
                    meat <- 0
                    for (ii in 1:n){
                        Xi.med <- X.1[ii, ]
                        Xi.out <- X.2[ii, ]
                        Xi <- rbind(Xi.med, Xi.out)
                        ui <- as.matrix(c(e.1[ii], e.2[ii]))
                        meat <- t(Xi) %*% solve(omega) %*% ui %*% t(ui) %*% solve(omega) %*% Xi + meat
                    }
                    v.cov <- solve(X.sur) %*% meat %*% solve(X.sur)  
                } else if(!is.null(cluster)){
                    meat <- 0
                    for (mm in unique(cluster)) {
                        X.med.m <- matrix(X.1[cluster == mm, ],, ncol(X.1))
                        X.out.m <- matrix(X.2[cluster == mm, ],, ncol(X.2))
                        e.1.m <- e.1[cluster == mm]
                        e.2.m <- e.2[cluster == mm]
                        half.meat <- 0
                        for (ii in 1:nrow(X.med.m)){
                            Xm <- rbind(X.med.m[ii, ], X.out.m[ii, ])
                            half.meat <- t(Xm) %*% solve(omega) %*% as.matrix(c(e.1.m[ii], e.2.m[ii])) + half.meat
                        }
                        meat <- half.meat %*% t(half.meat) + meat
                    }
                    v.cov <- solve(X.sur) %*% meat %*% solve(X.sur)  
                } else {
                    v.cov <- solve(X.sur)
                }
                b.old <- b.tmp
                b.dif <- sum((b.sur - b.old)^2)
                b.tmp <- b.sur
            }

            # Name Elements - Extract Quantities
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
            m.coef.treat <- m.coefs[T.out, ]
            y.coef.treat <- y.coefs[T.out, ]
            y.coef.med <- y.coefs[M.out, ]
            v.m.treat <- v.m[T.out, T.out]
            v.y.treat <- v.y[T.out, T.out]
            v.y.med <- v.y[M.out, M.out]
            if(INT){
                y.coef.TM <- y.coefs[TM.out, ]
                v.y.TM <- v.y[TM.out, TM.out]
            }
            
            if(INT && ("direct" %in% etype.vec)){
                m.bar <- apply(model.matrix(model.m), 2, weighted.mean, w=weights)
                m.covts.ind <- -(match(c("(Intercept)", paste(T.out)), names(model.m$coefficients), NULL))
                m.coef.covts <- m.coefs[m.covts.ind, ]
                m.bar.covts <- m.bar[m.covts.ind]
                m.intercept <- m.coefs["(Intercept)", ]
                v.m.covts <- v.m[m.covts.ind, m.covts.ind]
                v.m.intercept <- v.m["(Intercept)", "(Intercept)"]
                v.m.int.treat <- v.m["(Intercept)", T.out]
                
                if(length(m.coef.covts) != 0){
                    model.m.bar.z0 <- m.intercept + sum(crossprod(m.coef.covts, m.bar.covts))
                    v.model.m.bar.z0 <- v.m.intercept + sum(crossprod(m.bar.covts, v.m.covts) %*% m.bar.covts)
                    model.m.bar.z1 <- m.intercept + m.coef.treat + sum(m.coef.covts * m.bar.covts)
                    v.m.treat.covts <- v.m[T.out, m.covts.ind]
                    v.m.int.covts <- v.m["(Intercept)", m.covts.ind]
                    v.model.m.bar.z1 <- ( v.m.intercept + v.m.treat +
                                sum(crossprod(m.bar.covts, v.m.covts) %*% m.bar.covts) + 
                                2 * ( v.m.int.treat + 
                                    sum(crossprod(m.bar.covts, v.m.treat.covts)) +
                                    sum(crossprod(m.bar.covts, v.m.int.covts)) ) )
                } else {
                    model.m.bar.z0 <- m.intercept
                    v.model.m.bar.z0 <- v.m.intercept
                    model.m.bar.z1 <- m.intercept + m.coef.treat
                    v.model.m.bar.z1 <- v.m.intercept + v.m.treat + 2 * v.m.int.treat
                }
            }

            # Save Estimates
            if(INT){
                if("indirect" %in% etype.vec){
                    d0[i,] <- m.coef.treat * y.coef.med
                    d1[i,] <- m.coef.treat * (y.coef.med + y.coef.TM)
                }
                if("direct" %in% etype.vec){ # direct effects with interaction
                    if(length(m.coef.covts) != 0){ # with covariates
                        z0[i,] <- ( y.coef.treat + 
                                     y.coef.TM * (m.intercept + sum(crossprod(m.coef.covts, m.bar.covts))) )
                        z1[i,] <- ( y.coef.treat + 
                                     y.coef.TM * (m.intercept + m.coef.treat +
                                                  sum(crossprod(m.coef.covts, m.bar.covts))) )
                    } else { # no covariates
                        z0[i,] <- y.coef.treat + y.coef.TM * m.intercept
                        z1[i,] <- y.coef.treat + y.coef.TM * (m.intercept + m.coef.treat)
                    }
                }
            } else {
                if("indirect" %in% etype.vec){
                    d0[i,] <- d1[i,] <- m.coef.treat * y.coef.med
                }
                if("direct" %in% etype.vec) {
                    z0[i,] <- z1[i,] <- y.coef.treat
                }
            }

            # Save Variance Estimates
            if(INT){
                if("indirect" %in% etype.vec){
                    d0.var[i,] <- y.coef.med^2 * v.m.treat + m.coef.treat^2 * v.y.med
                    d1.var[i,] <- ( (y.coef.med + y.coef.TM)^2 * v.m.treat +
                        m.coef.treat^2 * (v.y.med + v.y.TM + 2 * v.y[M.out, TM.out]) )
                }
                if("direct" %in% etype.vec){
                    z0.var[i,] <- ( v.y.treat + model.m.bar.z0^2 * v.y.TM + 
                                     y.coef.TM^2 * v.model.m.bar.z0 +
                                     v.y.TM * v.model.m.bar.z0 + 
                                     2 * model.m.bar.z0 * v.y[T.out, TM.out] )
                    z1.var[i,] <- ( v.y.treat + model.m.bar.z1^2 * v.y.TM +
                                     y.coef.TM^2 * v.model.m.bar.z1 +
                                     v.y.TM * v.model.m.bar.z1 + 
                                     2 * model.m.bar.z1 * v.y[T.out,TM.out] )
                }
            } else {
                if("indirect" %in% etype.vec){
                    d0.var[i,] <- d1.var[i,] <- (m.coef.treat^2 * v.y.med) + (y.coef.med^2 * v.m.treat)
                }
                if("direct" %in% etype.vec){
                    z0.var[i,] <- z1.var[i,] <- v.y.treat
                }
            }

            rm(b.sur, m.coefs, y.coefs, v.cov, v.m, v.y)

        }  # END of Rho Loop

        if("indirect" %in% etype.vec){
            upper.d0 <- d0 + qnorm(high) * sqrt(d0.var)
            lower.d0 <- d0 + qnorm(low) * sqrt(d0.var)
            upper.d1 <- d1 + qnorm(high) * sqrt(d1.var)
            lower.d1 <- d1 + qnorm(low) * sqrt(d1.var)
            ind.d0 <- as.numeric(lower.d0 < 0 & upper.d0 > 0)
            ind.d1 <- as.numeric(lower.d1 < 0 & upper.d1 > 0)
        }
        if("direct" %in% etype.vec){
            upper.z0 <- z0 + qnorm(high) * sqrt(z0.var)
            lower.z0 <- z0 + qnorm(low) * sqrt(z0.var)
            upper.z1 <- z1 + qnorm(high) * sqrt(z1.var)
            lower.z1 <- z1 + qnorm(low) * sqrt(z1.var)
            ind.z0 <- as.numeric(lower.z0 < 0 & upper.z0 > 0)
            ind.z1 <- as.numeric(lower.z1 < 0 & upper.z1 > 0)
        }

        # Save R2 tilde values
        r.sq.m <- summary(model.m)$r.squared
        r.sq.y <- summary(model.y)$r.squared
        R2tilde.prod <- rho^2 * (1 - r.sq.m) * (1 - r.sq.y)

        
    ## END OF CASE 1: Continuous Outcome + Continuous Mediator
    } else

    #########################################################
    ## CASE 2: Continuous Outcome + Binary Mediator
    #########################################################
    if (class.y == "lm" & class.m == "glm") {
        type <- "bm"

        if(model.m$family$link == "logit"){
            stop("medsens is only valid for the probit link")
        }
        
        ## Variable values (LABEL.value)
        Y.value <- y.data[,1]
        y.k <- length(model.y$coef)

        # Step 1: Pre-loop computations
        # Step 1-1: Sample M model parameters
        Mmodel.coef <- model.m$coef
        if(robustSE){
            Mmodel.var.cov <- vcovHC(model.m)
        } else if(!is.null(cluster)) {
            dta <- merge(m.data, as.data.frame(cluster), sort=FALSE, by="row.names")
            fm <- update(model.m, data=dta)
            MModel.var.cov <- getvcov(dta, fm, dta[,ncol(dta)])
        } else {
            Mmodel.var.cov <- vcov(model.m)
        }
        Mmodel.coef.sim <- rmvnorm(sims, mean=Mmodel.coef, sigma=Mmodel.var.cov)
            # simulate M-model parameters

        # Step 1-2: Sample lambda_0 and lambda_1; lambdas are (n x sims) matrix
        m.mat <- m.mat.1 <- m.mat.0 <- model.matrix(model.m)
        m.mat.1[,T.out] <- 1 # M-model matrix with t=1
        m.mat.0[,T.out] <- 0 # M-model matrix with t=0
        mu.sim <- tcrossprod(m.mat, Mmodel.coef.sim) # E(M|T,X)
        mu.1.sim <- tcrossprod(m.mat.1, Mmodel.coef.sim) # E(M|T=1,X)
        mu.0.sim <- tcrossprod(m.mat.0, Mmodel.coef.sim) # E(M|T=0,X)

        # Step 1-3: Define lambda function (inverse Mill's ratio)
        lambda <- function(mmodel, mcoef) {
            mu <- model.matrix(mmodel) %*% mcoef
            m <- mmodel$y # this is M
            return((m*dnorm(-mu)-(1-m)*dnorm(-mu))/(m*pnorm(mu)+(1-m)*pnorm(-mu)))
        }

        # Step 2: Rho loop
        ## Step 2-0: Initialize containers
        if("indirect" %in% etype.vec){
            d0 <- d1 <- upper.d0 <- upper.d1 <- 
            lower.d0 <- lower.d1 <- ind.d0 <- ind.d1 <- matrix(NA, length(rho), 1)
        }
        if("direct" %in% etype.vec){
            z0 <- z1 <- upper.z0 <- upper.z1 <- 
            lower.z0 <- lower.z1 <- ind.z0 <- ind.z1 <- matrix(NA, length(rho), 1)
        }
        Ymodel.coef.sim <- matrix(NA, sims, y.k)
        colnames(Ymodel.coef.sim) <- names(model.y$coef)
        sigma.3.sim <- rep(NA, sims)
        d0.sim <- d1.sim <- z0.sim <- z1.sim <- rep(NA, sims)
        ## START OF RHO LOOP
        for(i in 1:length(rho)){

            ## START OF SIMULATION LOOP
            for(k in 1:sims){
            ## Step 2-1: Obtain the initial Y model with the correction term
            adj <- lambda(model.m, Mmodel.coef.sim[k,]) * rho[i] # the adjustment term
            w <- 1 - rho[i]^2*lambda(model.m, Mmodel.coef.sim[k,]) * 
                (lambda(model.m, Mmodel.coef.sim[k,]) + mu.sim[,k])
            y.data.adj <- data.frame(y.data, w, adj, weights)
            model.y.adj <- update(model.y, as.formula(paste(". ~ . + adj")), 
                                weights = w*weights, data = y.data.adj)
            sigma.3 <- summary(model.y.adj)$sigma

            ## Step 2-2: Update the Y model via Iterative FGLS
            sigma.dif <- 1
            while(abs(sigma.dif) > eps){
                Y.star <- Y.value - sigma.3 * adj
                y.data.star <- data.frame(Y.star, y.data.adj)
                model.y.update <- update(model.y, as.formula(paste("Y.star ~ .")), 
                                weights = w*weights, data = y.data.star)
                sigma.3.temp <- summary(model.y.update)$sigma
                sigma.dif <- sigma.3.temp - sigma.3
                sigma.3 <- sigma.3.temp
            }

            ## Step 2-3: Simulate Y model parameters
            Ymodel.coef <- model.y.update$coef
            if(robustSE){
                Ymodel.var.cov <- vcovHC(model.y.update)
            } else if(!is.null(cluster)) {
                dta <- merge(y.data, as.data.frame(cluster), sort=FALSE, by="row.names")
                fm <- update(model.y.update, data=dta)
                YModel.var.cov <- getvcov(dta, fm, dta[,ncol(dta)])
            } else {
                Ymodel.var.cov <- vcov(model.y.update)
            }
            Ymodel.coef.sim[k,] <- rmvnorm(1, mean=Ymodel.coef, sigma=Ymodel.var.cov) 
                # draw one simulation sample of Y-model parameters for each k

            ## Step 2-4: Simulate ACMEs; means are over observations
            if(INT){
                if("indirect" %in% etype.vec){
                    d1.sim[k] <- weighted.mean( (Ymodel.coef.sim[k,M.out] + Ymodel.coef.sim[k,TM.out]) *
                                            (pnorm(mu.1.sim[,k]) - pnorm(mu.0.sim[,k])), w = weights )
                    d0.sim[k] <- weighted.mean( (Ymodel.coef.sim[k,M.out]) * (pnorm(mu.1.sim[,k]) - 
                                            pnorm(mu.0.sim[,k])), w = weights )
                }
                if("direct" %in% etype.vec){
                    z1.sim[k] <- weighted.mean( sum(Ymodel.coef.sim[k,c(T.out,TM.out)]) * 
                                        pnorm(mu.1.sim[,k]) + Ymodel.coef.sim[k,T.out] * 
                                        pnorm(-mu.1.sim[k]), w = weights )
                    z0.sim[k] <- weighted.mean( sum(Ymodel.coef.sim[k,c(T.out,TM.out)]) * 
                                        pnorm(mu.0.sim[,k]) + Ymodel.coef.sim[k,T.out] * 
                                        pnorm(-mu.0.sim[k]), w = weights )
                }
            } else {
                if("indirect" %in% etype.vec){
                    d0.sim[k] <- d1.sim[k] <- weighted.mean( Ymodel.coef.sim[k,M.out] *
                        (pnorm(mu.1.sim[,k]) - pnorm(mu.0.sim[,k])), w = weights )
                }
                if("direct" %in% etype.vec){
                    z0.sim[k] <- z1.sim[k] <- Ymodel.coef.sim[k,T.out]
                }
            }

            ## END OF SIMULATION LOOP
            }

        ## Step 2-5: Compute Outputs
        if("indirect" %in% etype.vec){
            d0[i] <- mean(d0.sim)
            d1[i] <- mean(d1.sim)
            upper.d0[i] <- quantile(d0.sim, high)
            upper.d1[i] <- quantile(d1.sim, high)
            lower.d0[i] <- quantile(d0.sim, low)
            lower.d1[i] <- quantile(d1.sim, low)
            ind.d0[i] <- as.numeric(lower.d0[i] < 0 & upper.d0[i] > 0)
            ind.d1[i] <- as.numeric(lower.d1[i] < 0 & upper.d1[i] > 0)
        }
        if("direct" %in% etype.vec){
            z0[i] <- mean(z0.sim)
            z1[i] <- mean(z1.sim)
            upper.z0[i] <- quantile(z0.sim, high)
            upper.z1[i] <- quantile(z1.sim, high)
            lower.z0[i] <- quantile(z0.sim, low)
            lower.z1[i] <- quantile(z1.sim, low)
            ind.z0[i] <- as.numeric(lower.z0[i] < 0 & upper.z0[i] > 0)
            ind.z1[i] <- as.numeric(lower.z1[i] < 0 & upper.z1[i] > 0)
        }

        ## END OF RHO LOOP
        }

        # Save R2 tilde values
        fitted <- m.mat %*% Mmodel.coef
        var.mstar <- var(fitted)
        r.sq.m <- var.mstar/(1+var.mstar)
        r.sq.y <- summary(model.y)$r.squared
        R2tilde.prod <- rho^2 * (1 - r.sq.m) * (1 - r.sq.y)

    ## END OF CASE 2: Continuous Outcome + Binary Mediator
    } else

    #########################################################
    ## CASE 3: Binary Outcome + Continuous Mediator
    #########################################################
    if(class.y=="glm" & class.m=="lm") {
        if(INT){
            stop("sensitivity analysis is not available for binary outcome with interactions")
        }

        if(model.y$family$link == "logit"){
            stop("medsens is only valid for the probit link")
        }

        type <- "bo"
        
        # Step 1: Obtain Model Parameters
        ## Step 1-1: Simulate M model parameters
        Mmodel.coef <- model.m$coef
        if(robustSE){
            Mmodel.var.cov <- vcovHC(model.m)
        } else if(!is.null(cluster)) {
            dta <- merge(m.data, as.data.frame(cluster), sort=FALSE, by="row.names")
            fm <- update(model.m, data=dta)
            MModel.var.cov <- getvcov(dta, fm, dta[,ncol(dta)])
        } else {
            Mmodel.var.cov <- vcov(model.m)
        }
        m.k <- length(Mmodel.coef)
        Mmodel.coef.sim <- rmvnorm(sims, mean = Mmodel.coef, sigma = Mmodel.var.cov)
        beta2.sim <- Mmodel.coef.sim[, T.out]
        
        sigma.2 <- summary(model.m)$sigma
        sig2.shape <- model.m$df/2
        sig2.invscale <- (model.m$df/2) * sigma.2^2
        sigma.2.sim <- sqrt(1 / rgamma(sims, shape = sig2.shape, scale = 1/sig2.invscale))

        ## Step 1-2: Simulate Y model parameters
        Ymodel.coef <- model.y$coef
        if(robustSE){
            Ymodel.var.cov <- vcovHC(model.y)
        } else if(!is.null(cluster)) {
            dta <- merge(y.data, as.data.frame(cluster), sort=FALSE, by="row.names")
            fm <- update(model.y, data=dta)
            YModel.var.cov <- getvcov(dta, fm, dta[,ncol(dta)])
        } else {
            Ymodel.var.cov <- vcov(model.y)
        }
        y.k <- length(Ymodel.coef)
        Ymodel.coef.sim <- rmvnorm(sims, mean = Ymodel.coef, sigma = Ymodel.var.cov)
        colnames(Ymodel.coef.sim) <- names(Ymodel.coef)
        gamma.tilde <- Ymodel.coef.sim[, M.out]

        # Step 2: Compute ACME via the procedure in IKT
        ## Step 2-1: Estimate Error Correlation from inconsistent estimate of Y on M
        rho12.sim <- (sigma.2.sim * gamma.tilde) / (1 + sqrt(sigma.2.sim^2 * gamma.tilde^2))

        ## Step 2-2: Calculate alpha_1, beta_1 and xi_1
        YTmodel.coef.sim <- ( Ymodel.coef.sim[, !colnames(Ymodel.coef.sim) %in% M.out] * 
                              sqrt(1-rho12.sim^2) %x% t(rep(1,y.k-1)) + 
                              Mmodel.coef.sim * (rho12.sim/sigma.2.sim) %x% t(rep(1, y.k-1)) )

        ## Step 2-3: Calculate Gamma
        ## Data matrices for the Y model less M
        y.mat.1 <- y.mat.0 <- model.matrix(model.y)[, !colnames(model.matrix(model.y)) %in% M.out]
        y.mat.1[,T.out] <- 1
        y.mat.0[,T.out] <- 0
        

        ## Initialize objects before the rho loop
        d0 <- d1 <- z0 <- z1 <-
        upper.d0 <- upper.d1 <- lower.d0 <- lower.d1 <- 
        upper.z0 <- upper.z1 <- lower.z0 <- lower.z1 <- 
        ind.d0 <- ind.d1 <- ind.z0 <- ind.z1 <- 
        tau <- nu0 <- nu1 <- 
        upper.tau <- upper.nu0 <- upper.nu1 <-
        lower.tau <- lower.nu0 <- lower.nu1 <- rep(NA, length(rho))
        d0.sim <- d1.sim <- z0.sim <- z1.sim <- 
        sigma.1star.sim <- beta3.sim <- 
        tau.sim <- nu0.sim <- nu1.sim <- rep(NA, sims)

        ## START OF RHO LOOP
        for(i in 1:length(rho)){
            gamma.sim <- (-rho[i] + rho12.sim * sqrt((1 - rho[i]^2)/(1 - rho12.sim^2)))/sigma.2.sim
            for(k in 1:sims){
                sigma.1star.sim[k] <- sqrt(gamma.sim[k]^2 * 
                    sigma.2.sim[k]^2 + 2 * gamma.sim[k] * rho[i] * sigma.2.sim[k] + 1)
                YTcoef.k <- YTmodel.coef.sim[k,]
                if("indirect" %in% etype.vec){
                    d0.sim[k] <- weighted.mean( pnorm(y.mat.0 %*% YTcoef.k +
                                                 gamma.sim[k] * beta2.sim[k] / sigma.1star.sim[k]) - 
                                                 pnorm(y.mat.0 %*% YTcoef.k), w = weights )
                    d1.sim[k] <- weighted.mean( pnorm(y.mat.1 %*% YTcoef.k) - 
                                                 pnorm(y.mat.1 %*% YTcoef.k -
                                                 gamma.sim[k] * beta2.sim[k] / sigma.1star.sim[k]), 
                                                 w = weights )
                    tau.sim[k] <- weighted.mean( pnorm(y.mat.1 %*% YTcoef.k) - 
                                                  pnorm(y.mat.0 %*% YTcoef.k), w = weights )
                    nu0.sim[k] <- d0.sim[k]/tau.sim[k]
                    nu1.sim[k] <- d1.sim[k]/tau.sim[k]
                }
                if("direct" %in% etype.vec){
                    beta3.sim[k] <- sigma.1star.sim[k] * YTcoef.k[T.out] -
                                     beta2.sim[k] * gamma.sim[k]
                    z0.sim[k] <- weighted.mean( pnorm(y.mat.0 %*% YTcoef.k +
                                                 beta3.sim[k] / sigma.1star.sim[k]) -
                                                 pnorm(y.mat.0 %*% YTcoef.k), w = weights )
                    z1.sim[k] <- weighted.mean( pnorm(y.mat.1 %*% YTcoef.k) -
                                                 pnorm(y.mat.1 %*% YTcoef.k - 
                                                 beta3.sim[k] / sigma.1star.sim[k]), w = weights )
                }
            }

            ## Step 2-4: Compute Outputs
            if("indirect" %in% etype.vec){
                d0[i] <- mean(d0.sim) # ACME(t=0)
                d1[i] <- mean(d1.sim) # ACME(t=1)
                upper.d0[i] <- quantile(d0.sim, high)
                upper.d1[i] <- quantile(d1.sim, high)
                lower.d0[i] <- quantile(d0.sim, low)
                lower.d1[i] <- quantile(d1.sim, low)
                tau[i] <- mean(tau.sim) # ATE
                nu0[i] <- mean(nu0.sim) # Proportion Mediated (t=0)
                nu1[i] <- mean(nu1.sim) # Proportion Mediated (t=0)
                upper.tau[i] <- quantile(tau.sim, high)
                upper.nu0[i] <- quantile(nu0.sim, high)
                upper.nu1[i] <- quantile(nu1.sim, high)
                lower.tau[i] <- quantile(tau.sim, low)
                lower.nu0[i] <- quantile(nu0.sim, low)
                lower.nu1[i] <- quantile(nu1.sim, low)
                ind.d0[i] <- as.numeric(lower.d0[i] < 0 & upper.d0[i] > 0)
                ind.d1[i] <- as.numeric(lower.d1[i] < 0 & upper.d1[i] > 0)
            }
            if("direct" %in% etype.vec){
                z0[i] <- mean(z0.sim) # ACME(t=0)
                z1[i] <- mean(z1.sim) # ACME(t=1)
                upper.z0[i] <- quantile(z0.sim, high)
                upper.z1[i] <- quantile(z1.sim, high)
                lower.z0[i] <- quantile(z0.sim, low)
                lower.z1[i] <- quantile(z1.sim, low)
                ind.z0[i] <- as.numeric(lower.z0[i] < 0 & upper.z0[i] > 0)
                ind.z1[i] <- as.numeric(lower.z1[i] < 0 & upper.z1[i] > 0)
            }
            
        ## END OF RHO LOOP
        }

        # Save R2 tilde values
        r.sq.m <- summary(model.m)$r.squared
        y.mat <- model.matrix(model.y)
        fitted <- y.mat %*% Ymodel.coef
        var.ystar <- var(fitted)
        r.sq.y <- var.ystar/(1 + var.ystar)
        R2tilde.prod <- rho^2 * (1 - r.sq.m) * (1 - r.sq.y)

    ## END OF CASE 3: Binary Outcome + Continuous Mediator
    } else {
        stop("mediate object fitted with non-supported model combinations")
    }
    
    #########################################################
    ## Compute Output
    #########################################################
    # Calculate rho at which ACME=0
    
    if(INT || type=="bo"){
        if("indirect" %in% etype.vec){
            ii <- which(abs(d0 - 0) == min(abs(d0 - 0)))
            kk <- which(abs(d1 - 0) == min(abs(d1 - 0)))
            err.cr.1.d <- rho[ii]
            err.cr.2.d <- rho[kk]
            err.cr.d <- c(err.cr.1.d, err.cr.2.d)
        }
        if("direct" %in% etype.vec){
            ii.z <- which(abs(z0 - 0) == min(abs(z0 - 0)))
            kk.z <- which(abs(z1 - 0) == min(abs(z1 - 0)))
            err.cr.1.z <- rho[ii.z]
            err.cr.2.z <- rho[kk.z]
            err.cr.z <- c(err.cr.1.z, err.cr.2.z)
        }
    } else {
        if("indirect" %in% etype.vec){
            ii <- which(abs(d0 - 0) == min(abs(d0 - 0)))
            err.cr.d <- rho[ii]
        }
        if("direct" %in% etype.vec){
            ii.z <- which(abs(z0 - 0) == min(abs(z0 - 0)))
            err.cr.z <- rho[ii.z]
        }
    }

    if("indirect" %in% etype.vec){
        R2star.d.thresh <- err.cr.d^2
        R2tilde.d.thresh <- err.cr.d^2 * (1 - r.sq.m) * (1 - r.sq.y)
    }
    if("direct" %in% etype.vec){
        R2star.z.thresh <- err.cr.z^2
        R2tilde.z.thresh <- err.cr.z^2 * (1 - r.sq.m) * (1 - r.sq.y)
    }
        
    out <- list(
        d0=d0, d1=d1,
        upper.d0=upper.d0, lower.d0=lower.d0, 
        upper.d1=upper.d1, lower.d1=lower.d1, 
        z0=z0, z1=z1, 
        upper.z0=upper.z0, lower.z0=lower.z0,
        upper.z1=upper.z1, lower.z1=lower.z1, 
        tau=tau, upper.tau=upper.tau, lower.tau=lower.tau, 
        nu0=nu0, nu1=nu1, upper.nu0=upper.nu0, upper.nu1=upper.nu1,
        lower.nu0=lower.nu0, lower.nu1=lower.nu1,
        rho = rho, rho.by=rho.by, sims=sims,
        err.cr.d=err.cr.d, err.cr.z=err.cr.z, 
        ind.d0=ind.d0, ind.d1=ind.d1, ind.z0=ind.z0, ind.z1=ind.z1,
        R2star.prod=R2star.prod, R2tilde.prod=R2tilde.prod,
        R2star.d.thresh=R2star.d.thresh, R2tilde.d.thresh=R2tilde.d.thresh,
        R2star.z.thresh=R2star.z.thresh, R2tilde.z.thresh=R2tilde.z.thresh,
        r.square.y=r.sq.y, r.square.m=r.sq.m,
        INT=INT, conf.level = x$conf.level, effect.type=effect.type, type=type,
        robustSE = robustSE, cluster = cluster)
    class(out) <- "medsens"
    out
}


#' Summarizing Results from Sensitivity Analysis for Causal Mediation Effects
#' 
#' Functions to report results from the sensitivity analysis for causal 
#' mediation effects via \code{\link{medsens}} in a tabular form.
#' 
#' @aliases summary.medsens print.summary.medsens
#' 
#' @param object output from medsens function.
#' @param x output from summary.medsens function.
#' @param ...  additional arguments affecting the summary produced.
#' 
#' @author Dustin Tingley, Harvard University, 
#'   \email{dtingley@@gov.harvard.edu}; Teppei Yamamoto, Massachusetts Institute
#'   of Technology, \email{teppei@@mit.edu}; Jaquilyn Waddell-Boie, Princeton 
#'   University, \email{jwaddell@@princeton.edu}; Luke Keele, Penn State 
#'   University, \email{ljk20@@psu.edu}; Kosuke Imai, Princeton University, 
#'   \email{kimai@@princeton.edu}.
#'   
#' @seealso \code{\link{medsens}}, \code{\link{summary}}.
#' 
#' @references Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. 
#'   (2014). "mediation: R package for Causal Mediation Analysis", Journal of 
#'   Statistical Software, Vol. 59, No. 5, pp. 1-38.
#'   
#'   Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2011). Unpacking the 
#'   Black Box of Causality: Learning about Causal Mechanisms from Experimental 
#'   and Observational Studies, American Political Science Review, Vol. 105, No.
#'   4 (November), pp. 765-789.
#'   
#'   Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal 
#'   Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 
#'   309-334.
#'   
#'   Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and 
#'   Sensitivity Analysis for Causal Mediation Effects, Statistical Science,
#'   Vol. 25, No. 1 (February), pp. 51-71.
#'   
#'   Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation 
#'   Analysis Using R" in Advances in Social Science Research Using R, ed. H. D.
#'   Vinod New York: Springer.
#'   
#' @export
summary.medsens <- function(object, ...)
    structure(object, class = c("summary.medsens", class(object)))

#' @rdname summary.medsens
#' @export
print.summary.medsens <- function(x, ...){
    etype.vec <- x$effect.type
    if(etype.vec == "both"){
        etype.vec <- c("indirect","direct")
    }
    if(x$INT || x$type=="bo"){
        if("indirect" %in% etype.vec){
            tab.d0 <- cbind(round(x$rho,4), round(x$d0,4), round(x$lower.d0,4), 
                round(x$upper.d0, 4), x$ind.d0, round(x$R2star.prod,4), round(x$R2tilde.prod, 4))
            if(sum(x$ind.d0)==1){
                tab.d0 <- as.matrix(tab.d0[x$ind.d0==1, -5])
                tab.d0 <- t(tab.d0)
            } else {
                tab.d0 <- tab.d0[x$ind.d0==1, -5]
            }
            colnames(tab.d0) <-  c("Rho", "ACME(control)", paste(100*x$conf.level, "% CI Lower", sep=""), 
                paste(100*x$conf.level, "% CI Upper", sep=""), "R^2_M*R^2_Y*", "R^2_M~R^2_Y~")
            rownames(tab.d0) <- NULL
            tab.d1 <- cbind(round(x$rho,4), round(x$d1,4), round(x$lower.d1,4), 
                round(x$upper.d1, 4), x$ind.d1, round(x$R2star.prod,4), round(x$R2tilde.prod, 4))
            if(sum(x$ind.d1)==1){
                tab.d1 <- as.matrix(tab.d1[x$ind.d1==1, -5])
                tab.d1 <- t(tab.d1)
            } else {
                tab.d1 <- tab.d1[x$ind.d1==1, -5]
            }
            colnames(tab.d1) <-  c("Rho", "ACME(treated)", paste(100*x$conf.level, "% CI Lower", sep=""), 
                paste(100*x$conf.level, "% CI Upper", sep=""), "R^2_M*R^2_Y*", "R^2_M~R^2_Y~")
            rownames(tab.d1) <- NULL
            cat("\nMediation Sensitivity Analysis: Average Mediation Effect\n")
            if(sum(x$ind.d0) > 0){
	            cat("\nSensitivity Region: ACME for Control Group\n\n")
    	        print(tab.d0)
    	    }
            cat("\nRho at which ACME for Control Group = 0:", round(x$err.cr.d[1], 4))
            cat("\nR^2_M*R^2_Y* at which ACME for Control Group = 0:", round(x$R2star.d.thresh[1], 4))
            cat("\nR^2_M~R^2_Y~ at which ACME for Control Group = 0:", round(x$R2tilde.d.thresh[1], 4), "\n\n")
            if(sum(x$ind.d1) > 0){
	            cat("\nSensitivity Region: ACME for Treatment Group\n\n")
    	        print(tab.d1)
            }
            cat("\nRho at which ACME for Treatment Group = 0:", round(x$err.cr.d[2], 4))
            cat("\nR^2_M*R^2_Y* at which ACME for Treatment Group = 0:", round(x$R2star.d.thresh[2], 4))
            cat("\nR^2_M~R^2_Y~ at which ACME for Treatment Group = 0:", round(x$R2tilde.d.thresh[2], 4), "\n\n")
            invisible(x)
        }
        if("direct" %in% etype.vec){
            tab.z0 <- cbind(round(x$rho,4), round(x$z0,4), round(x$lower.z0,4), 
                round(x$upper.z0, 4), x$ind.z0, round(x$R2star.prod,4), round(x$R2tilde.prod, 4))
            if(sum(x$ind.z0)==1){
                tab.z0 <- as.matrix(tab.z0[x$ind.z0==1, -5])
                tab.z0 <- t(tab.z0)
            } else {
                tab.z0 <- tab.z0[x$ind.z0==1, -5]
            }
            colnames(tab.z0) <-  c("Rho", "ADE(control)", paste(100*x$conf.level, "% CI Lower", sep=""), 
                paste(100*x$conf.level, "% CI Upper", sep=""), "R^2_M*R^2_Y*", "R^2_M~R^2_Y~")
            rownames(tab.z0) <- NULL
            tab.z1 <- cbind(round(x$rho,4), round(x$z1,4), round(x$lower.z1,4), 
                round(x$upper.z1, 4), x$ind.z1, round(x$R2star.prod,4), round(x$R2tilde.prod, 4))
            if(sum(x$ind.z1)==1){
                tab.z1 <- as.matrix(tab.z1[x$ind.z1==1, -5])
                tab.z1 <- t(tab.z1)
            } else {
                tab.z1 <- tab.z1[x$ind.z1==1, -5]
            }
            colnames(tab.z1) <-  c("Rho", "ADE(treated)", paste(100*x$conf.level, "% CI Lower", sep=""), 
                paste(100*x$conf.level, "% CI Upper", sep=""), "R^2_M*R^2_Y*", "R^2_M~R^2_Y~")
            rownames(tab.z1) <- NULL
            cat("\nMediation Sensitivity Analysis: Average Direct Effect\n")
            if(sum(x$ind.z0) > 0){
	            cat("\nSensitivity Region: ADE for Control Group\n\n")
    	        print(tab.z0)
            }
            cat("\nRho at which ADE for Control Group = 0:", round(x$err.cr.z[1], 4))
            cat("\nR^2_M*R^2_Y* at which ADE for Control Group = 0:", round(x$R2star.z.thresh[1], 4))
            cat("\nR^2_M~R^2_Y~ at which ADE for Control Group = 0:", round(x$R2tilde.z.thresh[1], 4), "\n\n")
            if(sum(x$ind.z1) > 0){
	            cat("\nSensitivity Region: ADE for Treatment Group\n\n")
    	        print(tab.z1)
            }
            cat("\nRho at which ADE for Treatment Group = 0:", round(x$err.cr.z[2], 4))
            cat("\nR^2_M*R^2_Y* at which ADE for Treatment Group = 0:", round(x$R2star.z.thresh[2], 4))
            cat("\nR^2_M~R^2_Y~ at which ADE for Treatment Group = 0:", round(x$R2tilde.z.thresh[2], 4), "\n\n")
            invisible(x)
        }
    } else {
        if("indirect" %in% etype.vec){
            tab <- cbind(round(x$rho,4), round(x$d0, 4), round(x$lower.d0, 4), 
                round(x$upper.d0, 4), x$ind.d0, round(x$R2star.prod,4), round(x$R2tilde.prod,4))
            if(sum(x$ind.d0)==1){
                tab <- as.matrix(tab[x$ind.d0==1, -5])
                tab <- t(tab)
            } else {
                tab <- tab[x$ind.d0==1, -5]
            }
            colnames(tab) <-  c("Rho", "ACME", paste(100*x$conf.level, "% CI Lower", sep=""), 
                paste(100*x$conf.level, "% CI Upper", sep=""), "R^2_M*R^2_Y*", "R^2_M~R^2_Y~")
            rownames(tab) <- NULL
            cat("\nMediation Sensitivity Analysis for Average Causal Mediation Effect\n")
            if(sum(x$ind.d0) > 0){
	            cat("\nSensitivity Region\n\n")
    	        print(tab)
            }
            cat("\nRho at which ACME = 0:", round(x$err.cr.d, 4))
            cat("\nR^2_M*R^2_Y* at which ACME = 0:", round(x$R2star.d.thresh, 4))
            cat("\nR^2_M~R^2_Y~ at which ACME = 0:", round(x$R2tilde.d.thresh, 4), "\n\n")
            invisible(x)
        }
        if("direct" %in% etype.vec){
            tab.z <- cbind(round(x$rho,4), round(x$z0, 4), round(x$lower.z0, 4), 
                round(x$upper.z0, 4), x$ind.z0, round(x$R2star.prod,4), round(x$R2tilde.prod,4))
            if(sum(x$ind.z0)==1){
                tab.z <- as.matrix(tab.z[x$ind.z0==1, -5])
                tab.z <- t(tab.z)
            } else {
                tab.z <- tab.z[x$ind.z0==1, -5]
            }
            colnames(tab.z) <- c("Rho", "ADE", paste(100*x$conf.level, "% CI Lower", sep=""), 
                paste(100*x$conf.level, "% CI Upper", sep=""), "R^2_M*R^2_Y*", "R^2_M~R^2_Y~")
            rownames(tab.z) <- NULL
            cat("\nMediation Sensitivity Analysis for Average Direct Effect\n")
            if(sum(x$ind.z0) > 0){
	            cat("\nSensitivity Region\n\n")
    	        print(tab.z)
            }
            cat("\nRho at which ADE = 0:", round(x$err.cr.z, 4))
            cat("\nR^2_M*R^2_Y* at which ADE = 0:", round(x$R2star.z.thresh, 4))
            cat("\nR^2_M~R^2_Y~ at which ADE = 0:", round(x$R2tilde.z.thresh, 4), "\n\n")
            invisible(x)
        }
    }
}

#' Plotting Results from Sensitivity Analysis for Causal Mediation Effects
#' 
#' This function is used to plot results from the 'medsens' function. Causal 
#' average mediation effects (as well as average direct effects and proportions 
#' mediated for selected models) can be plotted against two alternative 
#' sensitivity parameters.
#' 
#' @details The sensitivity analysis for causal mediation effects can be
#'   conducted in terms of two alternative sensitivity parameters, which both
#'   quantify the degree of violation of the sequential ignorability assumption.
#'   The "rho" parameter represents the correlation between the two error terms
#'   of the (latent) linear models for the mediator and outcome variables. A
#'   large value of rho indicates the existence of important common unobserved
#'   predictors for both the mediator and outcome and therefore a high degree of
#'   sequential ignorability violation, while a value close to zero implies
#'   there is no such confounders.
#'   
#'   The resulting "rho" figures plot the estimated true values of ACME (or ADE,
#'   proportion mediated) against rho, along with the confidence intervals. When
#'   rho is zero, sequantial ignorability holds, so the estimated value at that 
#'   point will be equal to the estimate returned by the \code{\link{mediate}}. 
#'   The confidence level is determined by the 'conf.level' value of the
#'   original \code{\link{mediate}} object.
#'   
#'   The "R2" parameters represent the proportions of the mediator and outcome 
#'   variances that are explained by an unobserved pre-treatment confounder, 
#'   thereby indicating the importance of such a confounder in each model.  When
#'   'r.type' is "residual", the R2 parameters represent the proportions of the 
#'   residual variances of the mediator and outcome models that become explained
#'   by the inclusion of the hypothetical pre-treatment confounder.  These are 
#'   denoted as "R square stars" in Imai, Keele and Yamamoto (2010) and can also
#'   be specified as "star" or using a numeric value 1 in \code{medsens.plot}. 
#'   When 'r.type' is "total", the R2s represent the total mediator and outcome 
#'   variances the unobserved confounder would explain. This option can also be 
#'   specified using "tilde" or a numeric value 2.
#'   
#'   For both types of the "R2" parameters, 'sign.prod' indicates the
#'   hypothesized direction in which the unobserved confounder affects the
#'   mediator and outcome. (The name derives from the fact that this direction
#'   is mathematically represented by the sign of the product of two regression 
#'   coefficients.) If "positive" (or a numeric value 1) is given, the
#'   confounder is assumed to affect the mediator and outcome in the same
#'   direction.  If "negative" (or a numeric value -1), the effect is assumed to
#'   be in opposite directions.
#'   
#'   The resulting contours in the "R2" plots represent the values of the ACME
#'   (or ADE) for different combinations of the mediator R2 and outcome R2
#'   values. When both values are zero (the lower-left corner of the plot), the
#'   unobserved pre-treatment confounder has no effect on either mediator or
#'   outcome and therefore sequantial ignorability is satisfied.
#'   
#' @param x 'medsens' object, typically output from \code{medsens}.
#' @param sens.par a character string indicating the sensitivity parameter to be
#'   used. Default plots effects as functions of "rho". See Details.
#' @param r.type type of the R square parameter to be used in "R2" plots. If 
#'   "residual", effects are plotted against the proportions of the residual 
#'   variances that are explained by the unobserved confounder.  If "total", the
#'   proportions of the total variances are used as sensitivity parameters. Only
#'   relevant if 'sens.par' is "R2".
#' @param sign.prod a value indicating the direction of hypothesized confounding
#'   in the sensitivity analysis. If "positive", the confounder is assumed to 
#'   affect the mediator and outcome variable in the same direction; if 
#'   "negative" the effects are assumed to be in opposite directions.  Only 
#'   relevant if sens.par is set to "R2".
#' @param pr.plot a logical value. If 'TRUE', the "proportions mediated" will be
#'   plotted instead of the average causal mediation effects or direct effects. 
#'   Currently only available if the object 'medsens' is based on the linear 
#'   mediator and binary probit outcome models. Default is 'FALSE'.
#' @param smooth.effect a logical value indicating whether the estimated 
#'   mediation effects are smoothed via \code{\link{lowess}} before being 
#'   plotted. Default is 'FALSE'.
#' @param smooth.ci a logical value indicating whether the confidence bands are 
#'   smoothed via \code{\link{lowess}} before being plotted. Default is 'FALSE'.
#' @param ask a logical value. If 'TRUE', the user is asked for input before a 
#'   new figure is plotted.  Default is to ask only if the number of plots on 
#'   current screen is fewer than necessary.
#' @param levels vector of levels at which to draw contour lines. Only relevant 
#'   if 'sens.par' is set to "R2". If 'NULL', default values in 
#'   \code{\link{contour.default}} are used.
#' @param xlab label for the x axis. Default labels are used if 'NULL'.
#' @param ylab label for the y axis. Default labels are used if 'NULL'.
#' @param xlim limits of the x axis. If 'NULL' default values are used.
#' @param ylim limits of the y axis. If 'NULL' default values are used.
#' @param main main title for the plot. If 'NULL', default titles are used.
#' @param lwd width of the lines used in graphs.
#' @param ...  additional arguments to be passed to plotting functions.
#' 
#' @section Warning: The 'smooth.effect' and 'smooth.ci' options should be used 
#'   with caution since the smoothing could affect substantive implications of 
#'   the graphical analysis in a significant way.
#'   
#' @author Dustin Tingley, Harvard University, 
#'   \email{dtingley@@gov.harvard.edu}; Teppei Yamamoto, Massachusetts Institute
#'   of Technology, \email{teppei@@mit.edu}; Jaquilyn Waddell-Boie, Princeton 
#'   University, \email{jwaddell@@princeton.edu}; Luke Keele, Penn State 
#'   University, \email{ljk20@@psu.edu}; Kosuke Imai, Princeton University, 
#'   \email{kimai@@princeton.edu}.
#'   
#' @seealso \code{\link{medsens}}, \code{\link{plot}}, \code{\link{contour}}.
#' 
#' @references Tingley, D., Yamamoto, T., Hirose, K., Imai, K. and Keele, L. 
#'   (2014). "mediation: R package for Causal Mediation Analysis", Journal of 
#'   Statistical Software, Vol. 59, No. 5, pp. 1-38.
#'   
#'   Imai, K., Keele, L. and Tingley, D. (2010) A General Approach to Causal 
#'   Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 
#'   309-334.
#'   
#'   Imai, K., Keele, L. and Yamamoto, T. (2010) Identification, Inference, and 
#'   Sensitivity Analysis for Causal Mediation Effects, Statistical Science, 
#'   Vol. 25, No. 1 (February), pp. 51-71.
#'   
#'   Imai, K., Keele, L., Tingley, D. and Yamamoto, T. (2009) "Causal Mediation 
#'   Analysis Using R" in Advances in Social Science Research Using R, ed. H. D.
#'   Vinod New York: Springer.
#' @export
plot.medsens <- function(x, sens.par = c("rho", "R2"), 
        r.type = c("residual", "total"), sign.prod = c("positive", "negative"),
        pr.plot = FALSE, smooth.effect = FALSE, smooth.ci = FALSE, 
        ask = prod(par("mfcol")) < nplots, levels = NULL, 
        xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
        main = NULL, lwd = par("lwd"), ...){
    
    # Compatibility with previous version
    if(r.type[1L] == 1 | r.type[1L] == "star") r.type <- "residual"
    if(r.type[1L] == 2 | r.type[1L] == "tilde") r.type <- "total"
    if(sign.prod[1L] == 1) sign.prod <- "positive"
    if(sign.prod[1L] == -1) sign.prod <- "negative"
    
    # Match arguments
    sens.par <- match.arg(sens.par)
    r.type <- match.arg(r.type)
    sign.prod <- match.arg(sign.prod)
    
    # Effect types
    etype.vec <- x$effect.type
    if(x$effect.type == "both"){
        etype.vec <- c("indirect","direct")
    }
    
    # Determine number of plots necessary
    nplots <- ifelse(pr.plot, 2,
                ifelse(x$INT || x$type == "bo", 2, 1) * ifelse(x$effect.type == "both", 2, 1))
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    ##################################################
    ## A. Plots for proportions mediated
    ##################################################
    if(pr.plot){
        if(sens.par != "rho"){
            stop("proportion mediated can only be plotted against rho")
        }
        if(x$type != "bo"){
            stop("proportion mediated is only implemented for binary outcome")
        }
        
        out <- as.list(rep(NA, 2))
        out[[1]] <- list(nu = x$nu0, lo = x$lower.nu0, up = x$upper.nu0)
        out[[2]] <- list(nu = x$nu1, lo = x$lower.nu1, up = x$upper.nu1)
        r <- x$rho
        
        for(i in 1:2){
        	o <- out[[i]]
            if(is.null(ylim)){
                yl <- c(min(o$nu), max(o$nu))
            } else yl <- ylim
            if(is.null(main))
                ma <- bquote(paste("Proportion Mediated"[.(i-1)], (rho)))
                else ma <- main
            if(smooth.effect)
                nu <- lowess(r, o$nu)$y
                else nu <- o$nu
            if(smooth.ci){
                lower <- lowess(r, o$lo)$y
                upper <- lowess(r, o$up)$y
            } else {
                lower <- o$lo
                upper <- o$up
            }
            plot.default(r, nu, type = "n", xlab = "", ylab = "", 
            				main = ma, xlim = xlim, ylim = yl, ...)
            polygon(x = c(r, rev(r)), y = c(lower, rev(upper)), 
            		border = FALSE, col = 8, lty = 2, lwd = lwd)
            lines(r, nu, lty = 1, lwd = lwd)
            abline(h = 0)
            abline(v = 0)
            if(is.null(xlab))
                title(xlab = expression(paste("Sensitivity Parameter: ", rho)), 
                		line = 2.5, cex.lab = .9)
                else title(xlab = xlab)
            if(is.null(ylab))
                title(ylab = bquote(paste("Proportion Mediated: ", bar(nu), (.(i-1)))), 
                		line = 2.5, cex.lab = .9)
                else title(ylab = ylab)
        }
            
    ##################################################
    ## B. Plots for ACME/ADE wrt rho
    ##################################################
    } else if(sens.par == "rho"){
		
		out <- as.list(NA, 4)
		out[[1]] <- list(eff = x$d0, lo = x$lower.d0, up = x$upper.d0,
						 labs = c("ACME", "Average Mediation Effect", 0))
		out[[2]] <- list(eff = x$d1, lo = x$lower.d1, up = x$upper.d1, 
						 labs = c("ACME", "Average Mediation Effect", 1))
		out[[3]] <- list(eff = x$z0, lo = x$lower.z0, up = x$upper.z0,
						 labs = c("ADE", "Average Direct Effect", 0))
		out[[4]] <- list(eff = x$z1, lo = x$lower.z1, up = x$upper.z1,
						 labs = c("ADE", "Average Direct Effect", 1))
		r <- x$rho
		rby <- x$rho.by
		
	  	wh <- x$INT || x$type == "bo"
	   	out[[1]]$show <- "indirect" %in% etype.vec
	   	out[[2]]$show <- wh && "indirect" %in% etype.vec
	   	out[[3]]$show <- "direct" %in% etype.vec
	   	out[[4]]$show <- wh && "direct" %in% etype.vec

        for(i in 1:length(out)){
       		o <- out[[i]]
       		if(!o$show) next
            if(is.null(ylim)){
	            yl <- c(min(o$eff), max(o$eff))
	        } else yl <- ylim
            if(is.null(main)){
            	if(wh){
	                ma <- bquote(paste(.(o$labs[1])[.(o$labs[3])], "(", rho, ")"))
            	} else {
	                ma <- bquote(paste(.(o$labs[1]), "(", rho, ")"))
	            }
	        } else ma <- main
            if(smooth.effect)
                eff <- lowess(r, o$eff)$y
                else eff <- o$eff
            if(smooth.ci){
                lower <- lowess(r, o$lo)$y
                upper <- lowess(r, o$up)$y
            } else {
                lower <- o$lo
                upper <- o$up
            }
            plot.default(r, eff, type = "n", xlab = "", ylab = "", 
            			    main = ma, xlim = xlim, ylim = yl, ...)
            polygon(x = c(r, rev(r)), y = c(lower, rev(upper)), 
                    border = FALSE, col = 8, lty = 2, lwd = lwd)
            lines(r, eff, lty = 1, lwd = lwd)
            abline(h = 0)
            abline(v = 0)
            abline(h = weighted.mean(c(eff[floor(1/rby)], eff[ceiling(1/rby)]),
                				     c(1 - 1/rby + floor(1/rby), 1/rby - floor(1/rby))), 
                   lty = 2, lwd = lwd)
            if(is.null(xlab))
                title(xlab = expression(paste("Sensitivity Parameter: ", rho)), 
                      line = 2.5, cex.lab = .9)
                else title(xlab = xlab)
            if(is.null(ylab)){
            	if(wh){
	                title(ylab = bquote(.(o$labs[2])[.(o$labs[3])]), 
	                	line = 2.5, cex.lab = .9)
            	} else {
	                title(ylab = o$labs[2], line = 2.5, cex.lab = .9)
            	}
            } else title(ylab = ylab)
		}
        
    ##################################################
    ## C. Plots for ACME/ADE wrt R2
    ##################################################
    } else {
        R2Mstar <- seq(0, 1 - x$rho.by, 0.01)
        R2Ystar <- seq(0, 1 - x$rho.by, 0.01)
        R2Mtilde <- (1 - x$r.square.m) * seq(0, 1 - x$rho.by, 0.01)
        R2Ytilde <- (1 - x$r.square.y) * seq(0, 1 - x$rho.by, 0.01)

        if(r.type == "residual"){
            R2M <- R2Mstar; R2Y <- R2Ystar
        } else {
            R2M <- R2Mtilde; R2Y <- R2Ytilde
        }
        dlength <- length(seq(0, (1 - x$rho.by)^2, 0.0001))
        R2prod.mat <- outer(R2Mstar, R2Ystar)
        Rprod.mat <- round(sqrt(R2prod.mat), digits=4)

        if(is.null(xlim))
            xlim <- c(0,1)
        if(is.null(ylim))
            ylim <- c(0,1)

		out <- as.list(NA, 4)
		out[[1]] <- list(eff = x$d0, labs = c("ACME", "Average Mediation Effect", 0))
		out[[2]] <- list(eff = x$d1, labs = c("ACME", "Average Mediation Effect", 1))
		out[[3]] <- list(eff = x$z0, labs = c("ADE", "Average Direct Effect", 0))
		out[[4]] <- list(eff = x$z1, labs = c("ADE", "Average Direct Effect", 1))
		
	  	wh <- x$INT || x$type == "bo"
	   	out[[1]]$show <- "indirect" %in% etype.vec
	   	out[[2]]$show <- wh && "indirect" %in% etype.vec
	   	out[[3]]$show <- "direct" %in% etype.vec
	   	out[[4]]$show <- wh && "direct" %in% etype.vec
	   	
		for(i in 1:length(out)){
       		o <- out[[i]]
       		if(!o$show) next
	        if(is.null(levels)){
    	        levels <- pretty(quantile(o$eff, probs=c(0.1,0.9)), 10)
			}
 	      	if(sign.prod == "positive"){
    	    	eff <- approx(o$eff[((length(o$eff)+1)/2):length(o$eff)], n=dlength)$y
    	    	sgnlab <- "1"
    	    } else {
	    	    eff <- rev(approx(o$eff[1:((length(o$eff)+1)/2)], n=dlength)$y)
	    	    sgnlab <- "-1"
    	    }
   	        effmat <- matrix(eff[Rprod.mat/0.0001+1], nrow=length(R2M))
    	        
           	ma <- switch(wh + 1, bquote(paste(.(o$labs[1]))),				 # wh == F
           						 bquote(paste(.(o$labs[1])[.(o$labs[3])])))  # wh == T
      	    if(is.null(main) & r.type == "residual")
      	    	ma <- bquote(paste(.(ma), "(", R[M]^{2},"*,", R[Y]^2,"*), sgn", 
                                    (lambda[2]*lambda[3])==.(sgnlab)))
            else if(is.null(main) & r.type == "total")
                ma <- bquote(paste(.(ma), "(", tilde(R)[M]^{2}, "," , tilde(R)[Y]^2, "), sgn", 
                                    (lambda[2]*lambda[3])==.(sgnlab)))
            else ma <- main
            
            contour(R2M, R2Y, effmat, levels = levels, main = ma, 
                    xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim, lwd = lwd, ...)
            contour(R2M, R2Y, effmat, levels = 0, lwd = lwd + 1, add = TRUE, drawlabels = FALSE)
                        
	        if(is.null(xlab) & r.type == "residual")
    	        title(xlab=expression(paste(R[M]^{2},"*")), line=2.5, cex.lab=.9)
	        else if(is.null(xlab) & r.type == "total")
    	        title(xlab=expression(paste(tilde(R)[M]^{2})), line=2.5, cex.lab=.9)
    	    if(is.null(ylab) & r.type == "residual")
         	    title(ylab=expression(paste(R[Y]^2,"*")), line=2.5, cex.lab=.9)
	        else if(is.null(ylab) & r.type == "total")
    	        title(ylab=expression(paste(tilde(R)[Y]^2)), line=2.5, cex.lab=.9)
        	axis(2, at = seq(0, 1, by = .1))
        	axis(1, at = seq(0, 1, by = .1))
		}
    }
}
