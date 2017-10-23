#' Causal Mediation Analysis for Multiple Outcome/Treatment/Mediator 
#' Combinations
#' 
#' 'mediations' can be used to process a set of outcome/treatment/mediator 
#' combinations through the \code{\link{mediate}} function to produce a series 
#' of causal mediation analysis results.
#' 
#' @details This function processes multiple treatment/mediators/outcome
#'   variable combinations to produce a collected set of output ready for
#'   analysis or graphing. In principle, this is a function designed to
#'   facilitate running causal mediation analyses on multiple models that share
#'   the same basic specification (i.e. the types of parametric models and the
#'   set of pre-treatment covariates) except the treatment, mediator and outcome
#'   variables can differ across specifications. The function works by looping 
#'   over a set of data frames that are pre-loaded into the workspace. Each one 
#'   of these data frames has a specific treatment variable that is used for 
#'   analysis with that data frame. Then the code runs causal mediation analysis
#'   via \code{\link{mediate}} on every combination of the treatment, mediator, 
#'   and outcomes specified in these arguments. This allows the users to explore
#'   whether different mediators transmit the effect of the treatment variable
#'   on a variety of outcome variables. A single set of pre-treatment control 
#'   variables can be specified in 'covariates', which will be used throughout.
#'   
#'   The 'mediations' function can be used with either multiple mediators and a 
#'   single outcome, a single mediator and multiple outcomes, or multiple 
#'   mediators and outcomes. For example, with three different treatments, user 
#'   will create three different data frames, each containing a treatment 
#'   variable. In addition, if there are also four different mediators, each of 
#'   these will be contained in each data frame, along with the outcome
#'   variable. The function will estimate all of the combinations of treatment
#'   variables and mediators instead of separate lines of code being written for
#'   each one.
#'   
#'   Individual elements of the output list (see "Value") may be passed through 
#'   \code{\link[=summary.mediate]{summary}} and 
#'   \code{\link[=plot.mediate]{plot}} for tabular and graphical summaries of
#'   the results. Alternatively, the entire output may be directly passed to 
#'   \code{\link[=summary.mediations]{summary}} or 
#'   \code{\link[=plot.mediations]{plot}} for all results to be inspected.
#'   
#'   The default value of 'covariates' is 'NULL' and no covariate will be 
#'   included in either mediator or outcome models without a custom value. It 
#'   should be noted that users typically should have pre-treatment covariates
#'   to make the sequential ignorability assumption more plausible.
#'   
#'   There are several limitations to the code. First, it works only with a 
#'   subset of the model types that will be accommodated if 'mediate' is used 
#'   individually (see the 'families' argument above for details). Second, one 
#'   cannot specify separate sets of covariates for different 
#'   treatment/mediator/outcome combinations. Users should use 'mediate' 
#'   separately for individual models if more flexibility is required in their 
#'   specific applications.
#'   
#' @param datasets a named list of data frames. Each data frame has a separate 
#'   treatment variable. The names of each data frame must begin with the exact 
#'   name of the treatment variable that is contained in that dataset (see 
#'   example below).
#' @param treatment a vector of character strings indicating the names of the 
#'   treatment variables, with length equal to the length of 'datasets'. Each 
#'   treatment variable must be included in the data frame listed in the same 
#'   position of list 'datasets' and its name must match the first part of the 
#'   corresponding data frame.
#' @param mediators a vector of character strings indicating the names of the 
#'   mediators contained within each data frame. All of the mediators will be 
#'   used with each treatment variable and hence must be included in each data 
#'   frame of 'datasets'.
#' @param outcome a vector of character strings indicating the names of the 
#'   outcome variables contained within each data frame. All of the outcomes
#'   will be used with each treatment variable and must be in each data frame.
#' @param covariates a character string representing the set of pre-treatment 
#'   covariate names (as they appear in the data frame) to be included in each 
#'   model. The value must take the form of standard model formula, with each 
#'   additive component separated by "+", etc. (see example below). All 
#'   covariates must be in each data frame. Default is 'NULL'.
#' @param families a vector of length two specifying the types of the mediator 
#'   and outcome models. Currently only supports "gaussian" (for linear 
#'   regression), "binomial" (for binary probit), "oprobit" (for ordered probit)
#'   and "quantile" (for quantile regression, see 'tau'). For the outcome the 
#'   tobit model ("tobit") is also available in addition to the mediator model 
#'   options.
#' @param tau.m a numeric value specifying the quantile to be used for a 
#'   quantile regression for the mediator model. Only relevant if the first 
#'   element of 'families' is "quantile". See \code{rq}.
#' @param tau.y a numeric value specifying the quantile to be used for a 
#'   quantile regression for the outcome model. Only relevant if the second 
#'   element of 'families' is "quantile". See \code{rq}.
#' @param LowerY a numeric value indicating the lower bound for the tobit 
#'   outcome model. See \code{tobit}.
#' @param UpperY a numeric value indicating the upper bound for the tobit 
#'   outcome model. See \code{tobit}.
#' @param interaction a logical value indicating whether the treatment and 
#'   mediator variables should be interacted. This will apply to applications of
#'   \code{\link{mediate}} to all the treatment/mediator/outcome combinations.
#' @param conf.level confidence level used in each application of the 
#'   \code{\link{mediate}} function.
#' @param sims an integer indicating the desired number of simulations for 
#'   inference. This will apply to all applications of 'mediate' to all the 
#'   treatment/mediator/outcome combinations.
#' @param boot a logical value, indicating whether or not nonparametric 
#'   bootstrap should be used in each \code{\link{mediate}} application.
#' @param weights a single valued vector of a character string indicating a 
#'   weight variable to be used in all model fitting.
#' @param ...  other arguments passed to \code{\link{mediate}}, such as 
#'   'robustSE', 'dropobs', etc.
#'   
#' @return An object of class "mediations" (or "mediations.order" if the outcome
#'   model is ordered probit), a list of "mediate" ("mediate.order") objects
#'   produced by applications of \code{\link{mediate}} for the specified 
#'   treatment/mediator/outcome combinations.  The elements are named based on 
#'   the names of the outcome, treatment, and mediator variables, each separated
#'   by a "." (see example below).
#'   
#' @author Dustin Tingley, Harvard University, 
#'   \email{dtingley@@gov.harvard.edu}; Teppei Yamamoto, Massachusetts Institute
#'   of Technology, \email{teppei@@mit.edu}.
#'   
#' @seealso \code{\link{mediate}}, \code{\link{summary.mediations}}, 
#'   \code{\link{plot.mediations}}, \code{rq}, \code{tobit}.
#'   
#' @export
#' @examples
#'  
#' \dontrun{
#' # Hypothetical example
#' 
#' datasets <- list(T1 = T1, T2 = T2)
#'     # List of data frames corresponding to the two different treatment variables 
#'     #"T1vsCont" and "T2vsCont". 
#'     # Each data set has its respective treatment variable.
#'     
#' mediators <- c("M1", "M2") 
#'     # Vector of mediator names, all included in each data frame.
#' 
#' outcome <- c("Ycont1","Ycont2")
#'     # Vector of outcome variable names, again all included in each data frame.
#'     
#' treatment <- c("T1vsCont", "T2vsCont")
#'     # Vector of treatment variables names; must begin with identical strings with dataset 
#'     # names in 'datasets'.
#'     
#' covariates <- c("X1 + X2")
#'     # Set of covariates (in each data set), entered using the standard model formula format.
#' 
#' x <- mediations(datasets, treatment, mediators, outcome, covariates,
#'     families=c("gaussian","gaussian"), interaction=FALSE, 
#'     conf.level=.90, sims=50) 
#'     # Runs 'mediate' iteratively for each variable combinations, with 'lm' on both mediator 
#'     # and outcome model.
#' 
#' summary(x)  # tabular summary of results for all model combinations
#' plot(x)  # graphical summary of results for all model combinations at once
#' 
#' plot(x$Ycont1.T1vsCont.M1) 
#'     # Individual 'mediate' outputs are stored as list elements and can be 
#'     # accessed using the usual "$" operator.
#' }
#' 
mediations <- function(datasets, treatment, mediators, outcome, 
                       covariates=NULL, families=c("gaussian", "gaussian"),
                       tau.m=.5, tau.y=.5, LowerY=NULL, UpperY=NULL, 
                       interaction=FALSE, conf.level=.95, sims=500, 
                       boot=FALSE, weights=NULL, ...) {
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


#' Plotting Indirect, Direct, and Total Effects from Multiple Mediation Analyses
#' 
#' Function to plot results from multiple causal mediation analyses conducted 
#' via the \code{\link{mediations}} funciton. Output is a series of plots 
#' generated via \code{\link{plot.mediate}} for each treatment/mediator/outcome 
#' combination specified in the input 'mediations' object.
#' 
#' @aliases plot.mediations plot.mediations.order
#' 
#' @param x output from the mediations function.
#' @param which subset of names(x), indicating which model combinations to be 
#'   plotted. Default is to plot all.
#' @param ask logical. If 'TRUE', the user is asked for input before a new 
#'   figure is plotted.  Default is to ask only if the number of plots on
#'   current screen is fewer the number implied by 'which'.
#' @param ...  arguments passed to the \code{\link{plot.mediate}} function for 
#'   individual plots.
#'   
#' @return \code{mediations} returns an object of class \code{mediations}.  The 
#'   function \code{summary} is used to obtain a table of the results. The plot 
#'   function instead plots these quantities. All additional parameters desired 
#'   for the plotting of an output from \code{mediate} can be passed through.
#'   
#' @author Dustin Tingley, Harvard University, 
#'   \email{dtingley@@gov.harvard.edu}; Teppei Yamamoto, Massachusetts Institute
#'   of Technology, \email{teppei@@mit.edu}.
#'   
#' @seealso \code{\link{mediations}}, \code{\link{plot.mediate}}, 
#'   \code{\link{plot}}.
#'   
#' @export
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


#' @export
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


#' Summarizing Output from Multiple Mediation Analyses
#' 
#' The 'summary.mediations' function produces a summary of results from multiple
#' causal analyses conducted via \code{\link{mediations}}.  Output is a series
#' of \code{\link{summary.mediate}} outputs for all the 
#' treatment/mediator/outcome combinations used in the input 'mediations' 
#' object.
#' 
#' @aliases summary.mediations summary.mediations.order print.summary.mediations
#'   print.summary.mediations.order
#'   
#' @param object output from mediations function.
#' @param x output from summary.mediations function.
#' @param ...  additional arguments affecting the summary produced.
#' 
#' @author Dustin Tingley, Harvard University, 
#'   \email{dtingley@@gov.harvard.edu}; Teppei Yamamoto, Massachusetts Institute
#'   of Technology, \email{teppei@@mit.edu}.
#'   
#' @seealso \code{\link{mediations}}, \code{\link{summary.mediate}}, 
#'   \code{\link{summary}}.
#'   
#' @export
summary.mediations <- function(object, ...){
    structure(object, class = c("summary.mediations", class(object)))
}

#' @rdname summary.mediations
#' @export
print.summary.mediations <- function(x, ...){
    name.list <- names(x)
    for(i in 1:length(name.list)){
        cat("Specification", name.list[i], "\n") 
        print(summary.mediate(x[[i]]))
    }
}


#' @export
summary.mediations.order <- function(object, ...){
    structure(object, class = c("summary.mediations.order", class(object)))
}


#' @export
print.summary.mediations.order <- function(x, ...){
    name.list <- names(x)
    for(i in 1:length(name.list)){
        cat("Specification", name.list[i], "\n") 
        print(summary.mediate.order(x[[i]])  )
    }
}
