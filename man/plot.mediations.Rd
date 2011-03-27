\name{plot.mediations} \alias{plot.mediations} \title{Plotting Indirect, Direct, and Total Effects from Multiple 
Mediation Analyses} \description{ Plots results from mediations function, which calls the mediatie function. Y axis 
lists indirect, direct, and total effects and x-axis plots the respective magnitudes. Most standard options for plot 
function available. } \usage{ 

\method{plot}{mediations}(model,ask=FALSE,...) {

} 

\arguments{ 
\item{model}{Output from function mediations.} 

}

\details{This function passes results from the mediations function through \code{\link[=plot.mediate]{plot}}. Either 
separate plots can be generated (and par'ed together) or all results from mediations may be plotted one at a time, 
with the user scrolling through the separate results. All limitations to \code{\link[=plot.mediate]{plot}} apply 
here.} 

\arguments{ 
\item{x}{object of class \code{mediate} as produced by \code{mediate} or \code{mediations} function.}
\item{ask}{Whether to cycle through the plots one at a time. If TRUE then user controls advancing through each step.}
}


\value{\code{mediations} returns an object of class \code{mediations}.  The function \code{summary} is used to 
obtain a table of the results. The plot function instead plots these quantities. All additional parameters desired 
for the plotting of an output from \code{mediate} can be passed through.

}

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu} }
 
\seealso{\code{\link{mediations}}, \code{\link[=plot.mediate]{plot}}  }

\examples{ 
#datasets <- list(T1=T1, T2=T2) List of dataset names corresponding to the two different treatment variables T1vsCont and T2vsCont. Each data set has it's respective treatment/contrast variable.
#mediators <- c("M1", "M2" ) List of different mediators, M1 and M2 (in each data set)
#outcome<-c("Ycont1","Ycont2") List of different outcome variables (in each data set).
#treatment <- c("T1vsCont", "T2vsCont") List of treatment variables. Could also just be denoted c("T1","T2")
#covariates <- c("X1+X2") Set of conditioning pre-treatment covariates (in each data set)
#x <- mediations(datasets, treatment, mediators,outcome, covariates,family=c("gaussian","gaussian"), interaction=FALSE,conf.level=.90, sims=50) #runs mediator and outcome models with OLS on both mediator and outcome model.
#plot(x$Ycont1.T1vsCont.M1) plots the result of a model with Ycont1 outcome, T1vsCont treatment, and the mediator M1; names(x) will return the name of each mediation result. summary displays the summary statistics.



} 
