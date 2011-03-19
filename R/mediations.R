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
                
                # Mediations does not currently support the use of GAM's. Logits 
                # are not supported bc. users should use probits so they can do 
                # sensitivity analyses.
                
                if(families[1] == "binomial") {  # run Mediator model using new data/specification
                    result1 <- glm(f1, family=binomial("probit"), data=dataarg)
                } else if(families[1] == "quantile") {
                    result1 <- rq(f1, data=dataarg, tau=tau.m)
                } else if(families[1] == "oprobit") {
                    result1 <- polr(f1, method = "probit", data=dataarg)
                } else if (families[1] == "gaussian") {
                    result1 <- glm(f1, family="gaussian", data=dataarg)
                } else {
                    print("mediations does not support this model for the mediator")
                }
                
                if(families[2] == "binomial") {  # run Outcome model using new data/specification
                    result2 <- glm(f2, family=binomial("probit"), data=dataarg)
                } else if(families[2] == "quantile") {
                    result2 <- rq(f2, data=dataarg, tau=tau.y)
                } else if(families[2] == "tobit") {
                    result2 <- vglm(f2, tobit(Lower=LowerY,Upper=UpperY), data=dataarg)
                } else if(families[2]== "oprobit"){
                    result2 <- polr(f2, method = "probit", data=dataarg)
                } else if(families[2]== "gaussian"){
                    result2 <- glm(f2, family="gaussian", data=dataarg)
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
        plot.order.mediate(x[[i]], xlab = which[i], ...)
    }
}

summary.mediations <- function(object, ...){
    structure(object, class = c("summary.mediations", class(object)))
}

print.summary.mediations <- function(x, ...){
    name.list <- names(x)
    for(i in 1:length(name.list)){
        cat("Specification",name.list[i], "\n") 
        print(summary.mediate(x[[i]])  )
    }
    }
    


    
summary.mediations.order <- function(object, ...){
    structure(object, class = c("summary.mediations.order", class(object)))
}

print.summary.mediations.order <- function(x, ...){
    name.list <- names(x)
    for(i in 1:length(name.list)){
        cat("Specification",name.list[i], "\n") 
        print(summary.mediate.order(x[[i]])  )
    }
    }
    
