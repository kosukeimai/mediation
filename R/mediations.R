mediations <- function(datasets, treatment, mediators, outcome, 
                    covariates=NULL, families=c("gaussian", "gaussian"),
                    tau.m=.5, tau.y=.5, LowerY=NULL, UpperY=NULL, interaction=FALSE,
                    conf.level=.95, sims=500, boot=FALSE, ...) {
    data <- names(datasets)
    labels <- c()
    out <- list()
    count <- 1
    for (i in 1:length(treatment)) {
        d1 <- sprintf("datasets$%s", data[i])
        dataarg <- eval(parse(text=d1))
        for (o in 1:length(outcome)) {
            for (j in 1:length(mediators)) {
                # create model formulas
                if(is.null(covariates)) {
                    f1 <- sprintf("%s ~ %s + %s", mediators[j], treatment[i])
                    f2 <- sprintf("%s ~ %s + %s ", outcome[o], treatment[i], mediators[j])
                    if (interaction) {
                        f2 <- sprintf("%s ~ %s*%s ",outcome[o], treatment[i], mediators[j])
                    }
                } else {
                    f1 <- sprintf("%s ~ %s + %s", mediators[j], treatment[i], covariates)
                    f2 <- sprintf("%s ~ %s + %s + %s", outcome[o], treatment[i], 
                                                                mediators[j], covariates)
                    if (interaction) {
                        f2 <- sprintf("%s ~ %s*%s + %s", outcome[o], treatment[i], 
                                                                mediators[j], covariates)
                    }
                }
                
                fmla.m <- paste(f1)
                fmla.y <- paste(f2)
                
                # Mediations does not currently support the use of GAM's. Logits 
                # are not supported bc. users should use probits so they can do 
                # sensitivity analyses.
                
                if(families[1] == "binomial") {  # run Mediator model using new data/specification
                    result1 <- glm(fmla.m, family=binomial("probit"), data=dataarg)
                } else if(families[1] == "quantile") {
                    result1 <- rq(fmla.m, data=dataarg, tau=tau.m)
                } else if(families[1] == "oprobit") {
                    result1 <- polr(fmla.m, method = "probit", data=dataarg)
                } else if (families[1] == "gaussian") {
                    result1 <- glm(fmla.m, family="gaussian", data=dataarg)
                } else {
                    print("mediations does not support this model for the mediator")
                }
                
                if(families[2] == "binomial") {  # run Outcome model using new data/specification
                    result2 <- glm(fmla.y, family=binomial("probit"), data=dataarg)
                } else if(families[2] == "quantile") {
                    result2 <- rq(fmla.y, data=dataarg, tau=tau.y)
                } else if(families[2] == "tobit") {
                    result2 <- vglm(fmla.y, tobit(Lower=LowerY,Upper=UpperY), data=dataarg)
                } else if(families[2]== "oprobit"){
                    result2 <- polr(fmla.y, method = "probit", data=dataarg)
                } else if(families[2]== "gaussian"){
                    result2 <- glm(fmla.y, family="gaussian", data=dataarg)
                } else {
                    print("mediations does not support this model for the outcome")
                }
                
                out[[(count)]] <- mediate(result1, result2, sims=sims, 
                                        treat=treatment[i], mediator=mediators[j],
                                        conf.level=conf.level, boot=boot, ...)
                rm(result1, result2)
                labels[(count)] <- sprintf("%s.%s.%s", outcome[o],data[i], mediators[j])
                count <- count + 1
            }
        }
    }
    names(out) <- labels
    class(out) <- "mediations"
    out
}



plot.mediations <- function(x, which = names(x),
            ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for(i in 1:length(which)){
        plot.mediate(x[[i]], xlab = which[i], ...)
    }
}


plot.mediations.order <- function(x, which = names(x),
            ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    for(i in 1:length(which)){
        plot.order.mediate(x[[i]], xlab = which[i], ...)
    }
}


summary.mediations <- function(object, ...){
    structure(object, class = c("summary.mediations", class(object)))
}



print.summary.mediations <- function(x, ...){
    clp <- 100 * x[[1]]$conf.level
    name.list <- names(x)
    for(i in 1:length(name.list)){
        cat("\n Causal Mediation Analysis \n\n")
        cat("Specification",name.list[i], "\n\n")
        if(x[[i]]$boot==TRUE){
            cat("Confidence Intervals Based on Nonparametric Bootstrap\n\n")
        } else {
            cat("Quasi-Bayesian Confidence Intervals\n\n")
        }
#    if (isS4(x$model.y)){
#    printone<-FALSE
#    } else {
        printone <- x[[i]]$INT == FALSE && (class(x[[i]]$model.y)[1] %in% c("lm", "rq") ||
            (inherits(x[[i]]$model.y, "glm") && x[[i]]$model.y$family$family == "gaussian"
             && x[[i]]$model.y$family$link == "identity"))
#        }
        if (printone){
            # Print only one set of values if lmY/quanY without interaction
            cat("Mediation Effect: ", format(x[[i]]$d1, digits=4), clp, "% CI ", 
                    format(x[[i]]$d1.ci, digits=4), "\n")
            cat("Direct Effect: ", format(x[[i]]$z0, digits=4), clp, "% CI ", 
                    format(x[[i]]$z0.ci, digits=4), "\n")
            cat("Total Effect: ", format(x[[i]]$tau.coef, digits=4), clp, "% CI ", 
                    format(x[[i]]$tau.ci, digits=4), "\n")
            cat("Proportion of Total Effect via Mediation: ", format(x[[i]]$n0, digits=4), 
                    clp, "% CI ", format(x[[i]]$n0.ci, digits=4),"\n\n")
        } else {
            cat("Mediation Effect_0: ", format(x[[i]]$d0, digits=4), clp, "% CI ", 
                    format(x[[i]]$d0.ci, digits=4), "\n")
            cat("Mediation Effect_1: ", format(x[[i]]$d1, digits=4), clp, "% CI ", 
                    format(x[[i]]$d1.ci, digits=4), "\n")
            cat("Direct Effect_0: ", format(x[[i]]$z0, digits=4), clp, "% CI ", 
                    format(x[[i]]$z0.ci, digits=4), "\n")
            cat("Direct Effect_1: ", format(x[[i]]$z1, digits=4), clp, "% CI ", 
                    format(x[[i]]$z1.ci, digits=4), "\n")
            cat("Total Effect: ", format(x[[i]]$tau.coef, digits=4), clp, "% CI ", 
                    format(x[[i]]$tau.ci, digits=4), "\n")
            cat("Proportion of Total Effect via Mediation_0: ", format(x[[i]]$n0, digits=4), 
                    clp, "% CI ", format(x[[i]]$n0.ci, digits=4),"\n")
            cat("Proportion of Total Effect via Mediation_1: ", format(x[[i]]$n1, digits=4), 
                    clp, "% CI ", format(x[[i]]$n1.ci, digits=4),"\n\n")
            cat("Mediation Effect (Average): ", format(x[[i]]$d.avg, digits=4), clp, "% CI ", 
                    format(x[[i]]$d.avg.ci, digits=4), "\n")
            cat("Direct Effect (Average): ", format(x[[i]]$z.avg, digits=4), clp, "% CI ", 
                    format(x[[i]]$z.avg.ci, digits=4), "\n")
            cat("Proportion of Total Effect via Mediation (Average): ", format(x[[i]]$n.avg, digits=4), 
                    clp, "% CI ", format(x[[i]]$n.avg.ci, digits=4),"\n\n")
        }
    }
    invisible(x)
}


