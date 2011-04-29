\name{mediations} 
\alias{mediations} 
\title{Causal Mediation Analysis for Multiple Outcome/Treatment/Mediator Combinations} 
\description{ 
'mediations' can be used to process a set of outcome/treatment/mediator combinations through the \code{\link{mediate}} function to produce a series of causal mediation analysis results.} 
\usage{

mediations(datasets, treatment, mediators, outcome, 
        covariates = NULL, families = c("gaussian", "gaussian"),
        tau.m = 0.5, tau.y = 0.5, LowerY = NULL, UpperY = NULL, interaction = FALSE, 
        conf.level = .95, sims = 500, boot = FALSE, weights=NULL, ...)

}

\arguments{
\item{datasets}{ a named list of data frames. Each data frame has a separate treatment variable. The names of each data frame must begin with the exact name of the treatment variable that is contained in that dataset (see example below).}
\item{treatment}{ a vector of character strings indicating the names of the treatment variables, with length equal to the length of 'datasets'. Each treatment variable must be included in the data frame listed in the same position of list 'datasets' and its name must match the first part of the corresponding data frame.}
\item{mediators}{ a vector of character strings indicating the names of the mediators contained within each data frame. All of the mediators will be used with each treatment variable and hence must be included in each data frame of 'datasets'.}
\item{outcome}{ a vector of character strings indicating the names of the outcome variables contained within each data frame. All of the outcomes will be used with each treatment variable and must be in each data frame.}
\item{covariates}{ a character string representing the set of pre-treatment covariate names (as they appear in the data frame) to be included in each model. The value must take the form of standard model formula, with each additive component separated by "+", etc. (see example below). All covariates must be in each data frame. Default is 'NULL'.}
\item{families}{ a vector of length two specifying the types of the mediator and outcome models. Currently only supports "gaussian" (for linear regression), "binomial" (for binary probit), "oprobit" (for ordered probit) and "quantile" (for quantile regression, see 'tau'). For the outcome the tobit model ("tobit") is also available in addition to the mediator model options.}
\item{tau.m}{ a numeric value specifying the quantile to be used for a quantile regression for the mediator model. Only relevant if the first element of 'families' is "quantile". See \code{rq}.}
\item{tau.y}{ a numeric value specifying the quantile to be used for a quantile regression for the outcome model. Only relevant if the second element of 'families' is "quantile". See \code{rq}.}
\item{LowerY}{ a numeric value indicating the lower bound for the tobit outcome model. See \code{tobit}.}
\item{UpperY}{ a numeric value indicating the upper bound for the tobit outcome model. See \code{tobit}.}
\item{interaction}{ a logical value indicating whether the treatment and mediator variables should be interacted. This will apply to applications of \code{\link{mediate}} to all the treatment/mediator/outcome combinations.}
\item{conf.level}{ confidence level used in each application of the \code{\link{mediate}} function.}
\item{sims}{ an integer indicating the desired number of simulations to run \code{\link{mediate}} for. This will apply to all applications of 'mediate' to all the treatment/mediator/outcome combinations.}
\item{boot}{ a logical value, indicating whether or not nonparametric bootstrap should be used in each \code{\link{mediate}} application.}
\item{weights}{a single valued vector of a character string indicating a weight variable to be used in all model fitting.}
\item{...}{ other arguments passed to \code{\link{mediate}}, such as 'robustSE', 'dropobs', etc.}

}

\details{This function processes multiple treatment/mediators/outcome variable combinations to produce a collected set of  output ready for analysis or graphing. In principle, this is a function designed to facilitate running causal  mediation analyses on multiple models that share the same basic specification (i.e. the types of parametric models and the set of pre-treatment covariates) except the treatment, mediator and outcome variables can differ across specifications. The function works by looping over a set of data frames that are pre-loaded into the workspace. Each one of these data frames has a specific treatment variable that is used for analysis with that data frame. Then the code runs causal mediation analysis via \code{\link{mediate}} on every combination of the treatment, mediator, and outcomes specified in these arguments. This allows the users to explore whether different mediators transmit the effect of the treatment variable on a variety of outcome variables. A single set of pre-treatment control variables can be specified in 'covariates', which will be used throughout. 

The 'mediations' function can be used with either multiple mediators and a single outcome, a single mediator and multiple outcomes, or multiple mediators and outcomes. For example, with three different treatments, user will create three different data frames, each containing a treatment variable. In addition, if there are also four different mediators, each of these will be contained in each data frame, along with the outcome variable. The function will estimate all of the combinations of treatment variables and mediators instead of separate lines of code being written for each one.

Individual elements of the output list (see "Value") may be passed through \code{\link[=summary.mediate]{summary}} and  \code{\link[=plot.mediate]{plot}} for tabular and graphical summaries of the results. Alternatively, the entire output may be directly passed to \code{\link[=summary.mediations]{summary}} or \code{\link[=plot.mediations]{plot}} for all results to be inspected. 

The default value of 'covariates' is 'NULL' and no covariate will be included in either mediator or outcome models without a custom value. It should be noted that users typically should have pre-treatment covariates to make the sequential ignorability assumption more plausible.

There are several limitations to the code. First, it works only with a subset of the model types that will be accommodated if 'mediate' is used individually (see the 'families' argument above for details). Second, one cannot specify separate sets of covariates for different treatment/mediator/outcome combinations. Users should use 'mediate' separately for individual models if more flexibility is required in their specific applications.} 

\value{ An object of class "mediations" (or "mediations.order" if the outcome model is ordered probit), a list of "mediate" ("mediate.order") objects produced by applications of \code{\link{mediate}} for the specified treatment/mediator/outcome combinations.  The elements are named based on the names of the outcome, treatment, and mediator variables, each separated by a "." (see example below). } 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu},
Teppei Yamamoto, Princeton University, \email{tyamamot@princeton.edu}.}
 
\seealso{\code{\link{mediate}}, \code{\link{summary.mediations}}, \code{\link{plot.mediations}}, \code{rq}, \code{tobit}.}

\examples{ 
\dontrun{
# Hypothetical example

datasets <- list(T1 = T1, T2 = T2)
    # List of data frames corresponding to the two different treatment variables "T1vsCont" and "T2vsCont". 
    # Each data set has its respective treatment variable.
    
mediators <- c("M1", "M2") 
    # Vector of mediator names, all included in each data frame.

outcome <- c("Ycont1","Ycont2")
    # Vector of outcome variable names, again all included in each data frame.
    
treatment <- c("T1vsCont", "T2vsCont")
    # Vector of treatment variables names; must begin with identical strings with dataset names in 'datasets'.
    
covariates <- c("X1 + X2")
    # Set of covariates (in each data set), entered using the standard model formula format.

x <- mediations(datasets, treatment, mediators, outcome, covariates,
    families=c("gaussian","gaussian"), interaction=FALSE, 
    conf.level=.90, sims=50) 
    # Runs 'mediate' iteratively for each variable combinations, with 'lm' on both mediator and outcome model.

summary(x)  # tabular summary of results for all model combinations
plot(x)  # graphical summary of results for all model combinations at once

plot(x$Ycont1.T1vsCont.M1) 
    # Individual 'mediate' outputs are stored as list elements and can be accessed using the usual "$" operator.
}
} 
