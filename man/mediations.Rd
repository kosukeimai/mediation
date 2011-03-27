\name{mediations} 
\alias{mediations} 
\title{Causal Mediation Analysis for Multiple Treatment/Mediator Combinations} 
\description{ 
'mediations' can be used to process a set of treatment/mediator combinations through the \code{\link{mediate}} function to produce a series of causal mediation analysis results.} 
\usage{

mediations(datasets, treatment, mediators, outcome, 
        covariates=NULL, family=c("gaussian", "gaussian"),
        tau.m=.5,tau.y=.5, LowerY=NULL, UpperY=NULL, interaction=FALSE,
        conf.level=.95, sims=500, boot=FALSE, ...)



\details{This function processes multiple treatment/mediators variables combinations to produce a collected set of 
output ready for analysis or graphing. In principle, this is a function designed to facilitate looping across many 
different models instead of having to write out syntax for each separate model. The function works by having an 
outside loop over a set of data frames that are pre-loaded into memory. Each one of these data frames has a specific 
treatment variable that is used for analysis with that data frame. Then the function loops over a set of separate 
outcome variables. Within this loop the function then loops over a set of separate mediator variables. This allows 
the users to explore whether different mediators transmit the effect of the treatment variable on a variety of 
outcome variables. A single set of pre-treatment control variables is throughout.  After completing mediation 
analysis for each mediator, for each outcome, the function goes to the outerloop and loads a new data frame that has 
a separate treatment variables and repeats the process. Some times users will not have multiple different treatment 
or outcome variables, and instead will use the function to explore several separate mediators. For example, if you 
have 3 different treatments and a control condition, the user will create three different data frames, each 
containing a treatment variable equal to 0 if a control condition. In addition, if there are also 4 different 
mediators, each of these will be contained in each data frame, and an outcome variable.. The function will estimate 
all of the combinations of treatment variables and mediators instead of separate lines of code being written for 
each one. The output of mediations is an object of type "mediations". Individual elements of the object using a 
naming convention based on the variables names of the output, the treatment, and the mediator, each separated by a 
".". Individual elements of the output list may be passed through \code{\link[=summary.mediate]{summary}} and 
\code{\link[=plot.mediate]{plot}} for tabular and graphical summaries of the results. Alternatively, the entire 
output may be directly passed to \code{\link[=summary.mediations]{summary}} or \code{\link[plot.mediations]{plot}} 
to inspect all results. There are several limitations to the code. First, it works only with a somewhat limited set 
of models compared to the extremely flexible mediate function. Second, one cannot specify separate sets of 
covariates. In this case, separate uses of mediations can easily be used. For example, across different outcome 
variables, the same set of pre-treatment covariates might not be used to help establish sequential ignorability. } 


\value{
\item{datasets}{A list of character strings corresponding to the names of the data frames. Each dataset has a separate treatment/control contrast. The listed names of each dataset must begin the the exact name of the treatment variable that is contained in that dataset (see example below).}
\item{treatment}{Vector of treatment variable names contained in one of the data frames. Each element is a character expression and part of a vector within c(). Each treatment name must be the first part of each dataset name.}
\item{mediators}{Vector of mediators names contained within each data frame. Each element is a character expression and part of a vector within c(). Each mediator will be used with each treatment variable and hence must be in each dataframe of datasets.}
\item{outcome}{Single outcome variable to estimate all treatment/mediation combinations for.}
\item{covariates}{Single valued character vector that is a set of all pre-treatment covariates names (as they appear in the data frame) to be included in each model. Between each variable is a "+" sign indicating their additive entry to each model. All covariates must be in each model. Default is NULL but users typically should have pre-treatment variables to control for to the make the sequential ignorability assumption more plausible.}
\item{family}{Two element list specifying the link function to be used in the mediator and outcome models, respectively. Currently only supports "gaussian", "binomial" (for a probit link), "oprobit" (for ordered probit) and "quantile" (for quantile regression, see tau). For the outcome equation tobit ("tobit") is also available in addition to the mediator model options.}
\item{tau.m}{Specifies the quantile to be used for a quantile regression for the mediator model.}
\item{tau.y}{Specifies the quantile to be used for a quantile regression for the outcome model.}
\item{LowerY}{Lower bound for tobit model.}
\item{UpperY}{Upper bound for tobit model.}
\item{interaction}{Whether the treatment and mediator variables should be interacted. This will apply to all estimations.}
\item{conf.level}{Desired confidence interval calculations from mediate. This will apply to all estimations.}
\item{sims}{Desired number of simulations to run mediate for. This will apply to all estimations.}
\item{boot}{Whether or not bootstrap should be used. This will apply to all estimations.}
\item{...}{Other arguments passed to \code{\link{mediate}}, such as 'robustSE', 'dropobs', etc.}

}

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu} }
 
\seealso{\code{\link{mediate}}, \code{link{summary.mediations}}, \code{\link{plot.mediations}} }

\examples{ 
#Hypothetical example

#datasets <- list(T1=T1, T2=T2) List of dataset names corresponding to the two different treatment variables T1vsCont and T2vsCont. Each data set has it's respective treatment/contrast variable.
#mediators <- c("M1", "M2" ) List of different mediators, M1 and M2 (in each data set)
#outcome<-c("Ycont1","Ycont2") List of different outcome variables (in each data set).
#treatment <- c("T1vsCont", "T2vsCont") List of treatment variables. Could also just be denoted c("T1","T2")
#covariates <- c("X1+X2") Set of conditioning pre-treatment covariates (in each data set)
#x <- mediations(datasets, treatment, mediators,outcome, covariates,family=c("gaussian","gaussian"), interaction=FALSE,conf.level=.90, sims=50) #runs mediator and outcome models with OLS on both mediator and outcome model.
#plot(x$Ycont1.T1vsCont.M1) plots the result of a model with Ycont1 outcome, T1vsCont treatment, and the mediator M1; names(x) will return the name of each mediation result. summary displays the summary statistics.
} 
