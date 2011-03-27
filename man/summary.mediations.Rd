\name{summary.mediations} 
\alias{summary.mediations} 
\alias{print.summary.mediations} 

\title{Summarizing Output from Multiple Mediation Analyses} 
\description{ 
Function to report results from mediation analysis. Reported categories are mediation effect, direct effect, total effect, and proportion of total effect mediated. 
All quantities reported with confidence intervals specified in mediations. Prints summaries for all models estimated in mediations.
} 
\usage{

  \method{summary}{mediations}(object, ...) 

  \method{print}{summary.mediations}(x, ...)

} 

\arguments{ 
\item{object}{Output from mediations function.} 
\item{x}{Output from mediations function.} 
\item{...}{Additional arguments to be passed.}
}

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}}
 
\seealso{\code{\link{mediations}}, \code{\link[=summary.mediate]{summary}}  }

\examples{ 

#Hypothetical example

#datasets <- list(T1=T1, T2=T2) List of dataset names corresponding to the two different treatment variables T1vsCont and T2vsCont. Each data set has it's respective treatment/contrast variable.
#mediators <- c("M1", "M2" ) List of different mediators, M1 and M2 (in each data set)
#outcome<-c("Ycont1","Ycont2") List of different outcome variables (in each data set).
#treatment <- c("T1vsCont", "T2vsCont") List of treatment variables. Could also just be denoted c("T1","T2")
#covariates <- c("X1+X2") Set of conditioning pre-treatment covariates (in each data set)
#x <- mediations(datasets, treatment, mediators,outcome, covariates,family=c("gaussian","gaussian"), interaction=FALSE,conf.level=.90, sims=50) #runs mediator and outcome models with OLS on both mediator and outcome model.
#summary(x$Ycont1.T1vsCont.M1) plots the result of a model with Ycont1 outcome, T1vsCont treatment, and the mediator M1
#summary(x) provides summary results for all of the separate mediation analyses in x.
} 
