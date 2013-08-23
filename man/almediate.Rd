\name{almediate}
\alias{almediate}
\title{Function for combining outputs from mediations function and calculating quantities of interest. For use with multiple imputation procedures.}
\description{
'almediate' takes the output from \code{\link{mediations}} and stacks the different vectors. Next it outputs these stacked vectors in the format of a \code{\link{mediate}} object.}
\usage{

almediate(g)

}

\arguments{
\item{g}{Output from mediations that used the same models and variables but run on different datasets}
}

\details{almediate is designed to help users process multiple datasets where missing values have been imputed. First create multiple datasets using your preferred imputation software.
Next pass the data sets, as shown below, through the mediations function. Finally pass the output of mediations through the almediate function. This will output an object that can then be passed through the standard summary and plot commands.  }
This function is not completely developed. It does not support models for ordered outcomes, inherits the limitations of the mediations function, and does not pass the information required for calculation of p-values.
\value{ An object of class "mediate". }

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}}

\seealso{\code{\link{mediate}}, \code{\link{mediations}}.}

\examples{
\dontrun{
# Hypothetical example

##To use mediations, must make list of multiple datasets. Then, must also repeat the treatment assignment list as many times as you have data sets.
#T12<-T1
#datasets <- list(T1=T1, T12=T12) #T1 and T12 are separate datasets produced by the users multiple imputation procedure.
#mediators <- c("M1")
#outcome<-c("Ycont1")
#treatment <- c("T1","T1") #note how the treatment indicator is repeared
#covariates <- c("X1+X2")
#olsols <- mediations(datasets, treatment, mediators,outcome, covariates,families=c("gaussian","gaussian"),interaction=FALSE,conf.level=.90, sims=50)
#output<-almediate(olsols)
#plot(output)
}
}
