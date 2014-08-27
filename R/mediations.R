mediations <- function(datasets, treatment, mediators, outcome, 
                    covariates=NULL, families=c("gaussian", "gaussian"),
                    tau.m=.5, tau.y=.5, LowerY=NULL, UpperY=NULL, interaction=FALSE,
                    conf.level=.95, sims=500, boot=FALSE, weights=NULL, ...) {
    data <- names(datasets)
    labels <- c()
    out <- list()
    count <- 1
    weight.storage <- weights
    
    for (i in 1:length(treatment)) {
        d1 <- sprintf("datasets$%s", data[i])
        dataarg <- eval(parse(text=d1))
        for (o in 1:length(outcome)) {
            for (j in 1:length(mediators)) {
                # create model formulas
                if(is.null(covariates)) {
                    f1 <- sprintf("%s ~ %s ", mediators[j], treatment[i])
                    if (interaction) {
                        f2 <- sprintf("%s ~ %s * %s", outcome[o], treatment[i], mediators[j])
                    } else {
                        f2 <- sprintf("%s ~ %s + %s", outcome[o], treatment[i], mediators[j])
                    }
                } else {
                    f1 <- sprintf("%s ~ %s + %s", mediators[j], treatment[i], covariates)
                    if (interaction) {
                        f2 <- sprintf("%s ~ %s * %s + %s", outcome[o], treatment[i],
                                        mediators[j], covariates)
                    } else {
                        f2 <- sprintf("%s ~ %s + %s + %s", outcome[o], treatment[i], 
                                        mediators[j], covariates)
                    }
                }

                if(!is.null(weights)) {
                    weight1 <- sprintf("dataarg$%s", weights)
                    weight <- as.data.frame(eval(parse(text=weight1)))
                } else {
                    dataarg$weight <- weight <- rep(1,nrow(dataarg))
                }
                
                # run Mediator model using new data/specification
                if(families[1] == "binomial") {  
                    result1 <- glm(f1, family=binomial("probit"), weights=weight,
                                    data=dataarg)
                } else if(families[1] == "quantile") {
                    if(!is.null(weights)) {
                        stop("Weights not supported with quantile regression")
                    } else {
                        result1 <- quantreg::rq(f1, data=dataarg, tau=tau.m)
                    }
                } else if(families[1] == "oprobit") {
                    result1 <- polr(f1, method = "probit", weights=weight, data=dataarg, Hess=TRUE)
                } else if (families[1] == "gaussian") {
                    result1 <- glm(f1, family="gaussian", weights=weight, data=dataarg)
                } else {
                    stop("mediations does not support this model for the mediator")
                }
                
                # run Outcome model using new data/specification
                if(families[2] == "binomial") {  
                    result2 <- glm(f2, family=binomial("probit"), weights=weight,
                                    data=dataarg)
                } else if(families[2] == "quantile") {
                    if(!is.null(weights)) {
                        stop("Weights not supported with quantile regression")
                    } else {
                        result2 <- quantreg::rq(f2, data=dataarg, tau=tau.y)
                    }
                } else if(families[2] == "tobit") {
                    result2 <- VGAM::vglm(f2, VGAM::tobit(Lower=LowerY,Upper=UpperY), weights=weight,
                                    data=dataarg, model=TRUE)
                } else if(families[2]== "oprobit"){
                    result2 <- polr(f2, method = "probit", weights=weight, data=dataarg, Hess=TRUE)
                } else if(families[2]== "gaussian"){
                    result2 <- glm(f2, family="gaussian", weights=weight, data=dataarg)
                } else {
                    print("mediations does not support this model for the outcome")
                }
                
                if(is.null(weight.storage)){
                    out[[(count)]] <- mediate(result1, result2, sims=sims, 
                                        treat=treatment[i], mediator=mediators[j],
                                        conf.level=conf.level, boot=boot, ...)
                } else {
                    out[[(count)]] <- mediate(result1, result2, sims=sims, 
                                        treat=treatment[i], mediator=mediators[j],
                                        conf.level=conf.level, boot=boot,  ...)
                    weights <- weight.storage
                }
                
                rm(result1, result2)
                labels[(count)] <- sprintf("%s.%s.%s", outcome[o],treatment[i], mediators[j])
                count <- count + 1
            }
        }
        if(!is.null(weight.storage)){
            weights <- weight.storage
        }
    }
    names(out) <- labels
    if(families[2]== "oprobit") {
        class(out) <- "mediations.order"
    } else {
        class(out) <- "mediations"
    }
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
        plot.mediate.order(x[[i]], xlab = which[i], ...)
    }
}



summary.mediations <- function(object, ...){
    structure(object, class = c("summary.mediations", class(object)))
}



print.summary.mediations <- function(x, ...){
    name.list <- names(x)
    for(i in 1:length(name.list)){
        cat("Specification", name.list[i], "\n") 
        print(summary.mediate(x[[i]]))
    }
}



summary.mediations.order <- function(object, ...){
    structure(object, class = c("summary.mediations.order", class(object)))
}



print.summary.mediations.order <- function(x, ...){
    name.list <- names(x)
    for(i in 1:length(name.list)){
        cat("Specification", name.list[i], "\n") 
        print(summary.mediate.order(x[[i]])  )
    }
}
