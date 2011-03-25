\name{mediations} 
\alias{mediations} 
\title{Causal Mediation Analysis for Multiple Treatment/Mediator/Outcome Combinations} 
\description{ 
'mediations' can be used to process a set of treatment/mediator/outcome combinations through the \code{\link{mediate}} function to produce a series of causal mediation analysis results.} 
\usage{

mediations(datasets, treatment, mediators, outcome, 
        covariates=NULL, family=c("gaussian", "gaussian"),
        tau.m=.5,tau.y=.5, LowerY=NULL, UpperY=NULL, interaction=FALSE,
        conf.level=.95, sims=500, boot=FALSE, ...)



\details{This function processes multiple sets of data/treatment variables/mediators to produce a collected set of output ready for analysis or graphing. For example, if you have 3 different treatment variables and a control condition, as well as 4 different mediators that you want to see if they mediate the effect of the different treatment conditions, this function estimates all of the combinations instead of separate lines of code being written for each one.
Individual elements of the output list may be passed through \code{\link[=summary.mediate]{summary}} and \code{\link[=plot.mediate]{plot}} for tabular and graphical summaries of the results. Alternatively, the entire output may be directly passed to \code{\link[=summary.mediations]{summary}} or \code{\link[plot.mediations]{plot}} to inspect all results. 
}


\value{
\item{datasets}{ a list of data frames. Each dataset has a separate treatment/control contrast. Names of datasets must begin the the exact name of the treatment variable (see example below).}
\item{treatment}{Vector of treatment names.}
\item{mediators}{Vector of mediators names to be used with each treatment variable.}
\item{outcome}{Single outcome variable to estimate all treatment/mediation combinations for.}
\item{covariates}{Set of pre-treatment covariates to be included in each model. Default is NULL but users typically should have pre-treatment variables to control for to the make the sequential ignorability assumption more plausible.}
\item{family}{Two element list specifying the link function to be used in the mediator and outcome models, respectively. Currently only supports "gaussian", "binomial" (for a probit link), and "quantile" (for quantile regression, see tau) for the mediator. For the outcome equation tobit ("tobit") is also available in addition to the mediator model options.}
\item{tau_m}{Specifies the quantile to be used for a quantile regression for the mediator model.}
\item{tau_y}{Specifies the quantile to be used for a quantile regression for the outcome model.}
\item{LowerY}{Lower bound for tobit model.}
\item{UpperY}{Upper bound for tobit model.}
\item{interaction}{Whether the treatment and mediator variables should be interacted.}
\item{conf.level}{Desired confidence interval calculations from mediate.}
\item{sims}{Desired number of simulations to run mediate for.}
\item{boot}{Whether or not bootstrap should be used.}
\item{...}{Other arguments passed to \code{\link{mediate}}, such as 'robustSE', 'dropobs', etc.}

}

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu} }
 
\seealso{\code{\link{mediate}}, \code{link{summary.mediations}}, \code{\link{plot.mediations}} }

\examples{ 
#Hypothetical example

#datasets <- list(T1=T1, T2=T2) List of dataset names corresponding to the two different treatment variables T1 and T2. Each data set has it's respective treatment/contrast variable.
#mediators <- c("M1", "M2" ) List of mediators (in each data set)
#outcome<-c("Ycont1","Ycont2") List of outcome variables (in each data set).
#treatment <- c("T1", "T2") List of treatment variables
#covariates <- c("X1+X2") Set of conditioning pre-treatment covariates (in each data set)
#x <- mediations(datasets, treatment, mediators,outcome, covariates,family=c("gaussian","gaussian"), interaction=FALSE,conf.level=.90, sims=50) #runs mediator and outcome models with OLS
#plot(x) plots the results of each mediation analysis; names(x) will return each mediation result, and plot(x$A-Name) will plot that result.
} 
