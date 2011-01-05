\name{mediations} 
\alias{mediations} 
\title{Processing of multiple treatment variable/mediator variable combinations by mediate function.} 
\description{ 
Processes a set of treatment variable/mediator combinations through the mediate function for several standard models.} 
\usage{

mediations(datasets, treatment, mediators, covariates, family=c("gaussian", "gaussian"),LowerY=NULL,UpperY=NULL, interaction=FALSE, conf.level=.95, sims=1000) {

} 

\arguments{ 
\item{model}{Output from function mediate or function mediations.} 

}

\details{This function processes multiple sets of data/treatment variables/mediators to produce a collected set of output ready for analysis or graphing. For example, if you have 3 different treatment variables and a control condition, as well as 4 different mediators that you want to see if they mediate the effect of the different treatment conditions, this function estimates all of the combinations instead of separate lines of code being written for each one. Calling the $plot object can be passed through plot.mediate
}

\section{Warning}{Users wishing additional functionalities can easily update the function.}

\value{
\item{datasets}{List of dataset names, each using a naming convention that includes the treatment variable name.}
\item{treatment}{Vector of treatment names.}
\item{mediators}{Vector of mediators names to be used with each treatment variable.}
\item{covariates}{Set of pre-treatment covariates to be included in each model.}
\item{family}{Two element list specifying the link function to be used in the mediator and outcome models, respectively. Currently only supports gaussian and binomial (for a probit link) for the mediator. For the outcome equation tobit is also available.}
\item{LowerY}{Lower bound for tobit model.}
\item{UpperY}{Upper bound for tobit model.}
\item{interaction}{Whether the treatment and mediator variables should be interacted.}
\item{conf.level}{Desired confidence interval calculations from mediate.}
\item{sims}{Desired number of simulations to run mediate for.}
}

\references{Imai, Kosuke, Luke Keele and Dustin Tingley (2010) A General Approach to Causal Mediation Analysis, Psychological Methods, Vol. 15, No. 4 (December), pp. 309-334.
Imai, Kosuke, Luke Keele and Teppei Yamamoto (2010) Identification, Inference, and Sensitivity Analysis for Causal Mediation Effects, Statistical Science, Vol. 25, No. 1 (February), pp. 51-71.
Imai, Kosuke, Luke Keele, Dustin Tingley, and Teppei Yamamoto. (2009) "Causal Mediation Analysis Using R" in Advances in Social Science Research Using R, ed. H. D. Vinod New York: Springer. Forthcoming.} 

\author{Dustin Tingley, Harvard University, \email{dtingley@gov.harvard.edu} }
 
\seealso{See also \code{\link{mediate}} and \code{\link{plot.mediate}} }

\examples{ 
#Hypothetical example

#datasets <- list(Anxiety=Anxiety, Guilt=Guilt, Anger=Anger)  This gives a list of the names of different data sets that contain a particular treatment/control condition as well as any mediators (can be multiple)
#mediators <- c("BasNeg_Hostility", "BasNeg_Guilt", "GD_NegAff", "GD_PosAff","BasNeg_Fear","BasPos_Joviality") A set of mediator that are used with each of the treatment variables. This permits investigation of multiple treatment variables for a single treatment variable.
#treatment <- c("AnxietyVsCont", "GuiltVsCont", "AngerVsCont")#A list of the treatment variables, each of which belongs to one of the data sets. Further, note that the first word corresponds to the name of each data set. Preprocess your datasets to take advantage of this program. 
#covariates <- c("gender_p+Trust_GSS_pre+ideo7") #Pre-treatment covariates contained in each model.
#plot.int <- mediations(datasets, treatment, mediators, covariates,family=c("gaussian","tobit"),interaction=TRUE,UpperY=5,LowerY=0,conf.level=.90, sims=5000)$plot
#above model uses lm for mediator model and tobit for outcome model. Interaction between T and M is included.
} 
