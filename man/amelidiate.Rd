\name{amelidiate}
\alias{amelidiate}
\title{Function for combining outputs from mediations function and calculating quantities of interest. For use with multiple imputation procedures.}
\description{
'amelidiate' takes the output from \code{\link{mediations}} and stacks the different vectors. Next it outputs these stacked vectors in the format of a \code{\link{mediate}} object.}
\usage{

amelidiate(g)

}

\arguments{
\item{g}{ output from mediations that used the same models and variables but run on different datasets.}
}

\details{\code{amelidiate} is designed to help users process multiple datasets where missing values have been imputed. First create multiple datasets using your preferred imputation software.
Next pass the data sets, as shown in the example below, to the \code{\link{mediations}} function. Finally pass the output of mediations through the \code{amelidiate} function. This will output an object that can then be passed through the standard summary and plot commands.

This function is not completely developed. It does not support models for ordered outcomes, inherits the limitations of the mediations function, and does not pass the information required for calculation of p-values.
}

\value{ An object of class "mediate". }

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu}}

\seealso{\code{\link{mediate}}, \code{\link{mediations}}.}

\examples{
\dontrun{
# Hypothetical example

## To use mediations, must make list of multiple datasets. Then, 
## must also repeat the treatment assignment list as many times 
## as you have data sets.
# datasets <- list(D1=D1, D2=D2) # list of multiply imputed data sets
# mediators <- c("M1")
# outcome <- c("Ycont1")
# treatment <- c("T1","T1") # note how the treatment indicator is repeated
# covariates <- c("X1+X2")
# olsols <- mediations(datasets, treatment, mediators, outcome, covariates, 
#                      families=c("gaussian","gaussian"), interaction=FALSE,
#                      conf.level=.90, sims=1000)
# output <- amelidiate(olsols)
# summary(output)
# plot(output)
}
}
