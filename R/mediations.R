mediations <- function(datasets, treatment, mediators, outcome, 
                    covariates=NULL, families=c("gaussian", "gaussian"),
                    tau.m=.5, tau.y=.5, LowerY=NULL, UpperY=NULL, interaction=FALSE,
                    conf.level=.95, sims=500, boot=FALSE, weight=NULL, ...) {
    data <- names(datasets)
    labels <- c()
    out <- list()
    count <- 1
    weight.storage<-weight
    
    
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
                        f2 <- sprintf("%s ~ %s + %s + %s*%s ",outcome[o], treatment[i], mediators[j], treatment[i], mediators[j])
                    }
                } else {
                    f1 <- sprintf("%s ~ %s + %s", mediators[j], treatment[i], covariates)
                    f2 <- sprintf("%s ~ %s + %s + %s", outcome[o], treatment[i], 
                                                                mediators[j], covariates)
                    if (interaction) {
                        f2 <- sprintf("%s ~ %s + %s + %s*%s + %s", outcome[o], treatment[i], 
                                                                mediators[j], , treatment[i], mediators[j], covariates)
                    }
                }
                

                # Mediations does not currently support the use of GAM's. Logits 
                # are not supported bc. users should use probits so they can do 
                # sensitivity analyses.
                if(is.null(weight)) {
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
                } else {
                    weight1 <- sprintf("dataarg$%s", weight)
                    weight <- eval(parse(text=weight1))
                    print(length(weight))
                    weight<-as.data.frame(weight)
                    if(families[1] == "binomial") {  # run Mediator model using new data/specification
                        result1 <- glm(f1, family=binomial("probit"), weights=weight, data=dataarg)
                    } else if(families[1] == "quantile") {
                        stop("Weights not supported with quantile regression")
                    } else if(families[1] == "oprobit") {
                        result1 <- polr(f1, method = "probit", weights=weight, data=dataarg)
                    } else if (families[1] == "gaussian") {
                        print("test6")
                        print(f2)
                        print(d1)
                        result1 <- glm(f2, family="gaussian", data=dataarg)#use this so only has obs from m and y models            
                        odata.m <- model.frame(result1)
                        newdata <- merge(odata.m, weight, sort=FALSE,
                                    by=c("row.names"))
                        rownames(newdata) <- newdata$Row.names
                        newdata <- newdata[,-1L]
                        #rm(odata.m,weight)
                        weight.temp<-weight
                        newdata<-as.data.frame(newdata)
                        result1 <- glm(f1, family="gaussian", weights=weight, data=newdata)
                    } else {
                        print("mediations does not support this model for the mediator")
                    }
                    
                    if(families[2] == "binomial") {  # run Outcome model using new data/specification
                        result2 <- glm(f2, family=binomial("probit"), weights=weight,  data=dataarg)
                    } else if(families[2] == "quantile") {
                        stop("Weights not supported with quantile regression")
                    } else if(families[2] == "tobit") {
                        result2 <- vglm(f2, tobit(Lower=LowerY,Upper=UpperY), weights=weight, data=dataarg)
                    } else if(families[2]== "oprobit"){
                        result2 <- polr(f2, method = "probit", weights=weight, data=dataarg)
                    } else if(families[2]== "gaussian"){
                        weight<-weight.temp
                        result2 <- glm(f2, family="gaussian", data=dataarg)
            
                        odata.y <- model.frame(result2)
                        newdata <- merge(odata.y, weight, sort=FALSE,
                                    by=c("row.names"))
                        rownames(newdata) <- newdata$Row.names
                        newdata <- newdata[,-1L]
                        print(nrow(newdata))
                        #print(newdata$weight)
                        rm(odata.y, weight)
                        print(length(newdata$weight))
                        #print(is.vector(newdata$weight))
                        #weights<-newdata$weight
                        #print(length(weight))
                        print(names(newdata))
                        newdata<-as.data.frame(newdata)
                        print(names(newdata))
                        result2 <- glm(f2, family="gaussian", weights=weight, data=newdata)
                    
                    print("test2")
                        #result2 <- glm(f2, family="gaussian", weights=weight, data=newdata)
                   
                        #result2 <- glm(f2, family="gaussian", weights=weight, data=dataarg)
                        print(names(result2))
                        print("test3")
                        #print(result2$weights)
                    } else {
                        print("mediations does not support this model for the outcome")
                    }
                }
                
                if(is.null(weight.storage)){
                out[[(count)]] <- mediate(result1, result2, sims=sims, 
                                        treat=treatment[i], mediator=mediators[j],
                                        conf.level=conf.level, boot=boot, ...)
                } else {
                out[[(count)]] <- mediate(result1, result2, sims=sims, 
                                        treat=treatment[i], mediator=mediators[j],
                                        conf.level=conf.level, boot=boot, dropobs=FALSE, ...)
                weight<-weight.storage
                }
                
                rm(result1, result2)
                labels[(count)] <- sprintf("%s.%s.%s", outcome[o],treatment[i], mediators[j])
                count <- count + 1
            }
        }
        if(!is.null(weight.storage)){
        weight<-weight.storage
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
        cat("Specification", name.list[i], "\n") 
        print(summary.mediate(x[[i]])  )
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
