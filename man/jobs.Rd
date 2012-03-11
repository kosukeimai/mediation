\name{jobs}
\docType{data}
\alias{jobs}
\title{JOBS II data}
\description{
Job Search Intervention Study (JOBS
II). JOBS II is a randomized field experiment that investigates the efficacy of a job training
intervention on unemployed workers. The program is designed to not only increase reemployment
among the unemployed but also enhance the mental health of the job seekers. In the JOBS II
field experiment, 1,801 unemployed workers received a pre-screening questionnaire and were then
randomly assigned to treatment and control groups. Those in the treatment group participated
in job-skills workshops. In the workshops, respondents learned job-search skills and coping strategies
for dealing with setbacks in the job-search process. Those in the control condition received
a booklet describing job-search tips. In follow-up interviews, the two key outcome variables were
measured; a continuous measure of depressive symptoms based on the Hopkins Symptom Checklist,
and a binary variable, representing whether the respondent had become employed.
}
\usage{data}
\format{A data matrix with 899 rows and 16 columns, containing no missing values. The data are provided only for illustrative purposes and not for inference about program efficacy, for which the original data source should be consulted.

\describe{
	\item{econ_hard:}{ Level of economic hardship pre-treatment with values from 1 to 5.}
	\item{depress1:}{ Measure of depressive symptoms pre-treatment.}
	\item{sex:}{ Indicator variable for sex. 1 = female}
	\item{age:}{ Age in years.}
	\item{occp:}{ Factor with seven categories for various occupations.}
	\item{marital:}{ Factor with five categories for marital status.}
	\item{nonwhite:}{ Indicator variable for race. 1 = nonwhite.}
	\item{educ:}{ Factor with five categories for educational attainment.}
	\item{income:}{ Factor with five categories for level of income.}
	\item{job_seek:}{ A continuous scale measuring the level of job-search self-efficacy with values from 1 to 5. The mediator variable.}
	\item{depress2:}{ Measure of depressive symptoms post-treatment.}
	\item{work1:}{ Indicator variable for employment. 1 = employed.}
	\item{job_dich:}{ The job_seek measure recoded into two categories of high and low. 1 = high job search self-efficacy.}
	\item{job_disc:}{ The job_seek measure recoded into four categories from lowest to highest.}
	\item{treat:}{ Indicator variable for whether participant was randomly selected for the JOBS II training program. 1 = participation.}
	\item{control:}{ Indicator variable for whether participant was randomly selected to not participate in the JOBS II training program. 1 = non-participation.}
	}

}

\source{The complete JOBS II data is available from the data archives at www.icpsr.umich.edu/}
\references{
Vinokur, A. and Schul, Y. (1997). Mastery and inoculation against setbacks as active ingredients
in the jobs intervention for the unemployed. Journal of Consulting and Clinical Psychology 65,
5.
}
\keyword{datasets}
