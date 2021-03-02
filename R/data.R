#' Example Data for the Design Functions
#' 
#' A random subsample of the simulated data used in Imai, Tingley, and
#' Yamamoto (2012). The data contains 1000 rows and 7 columns with no missing 
#' values.
#' 
#' @usage boundsdata
#' 
#' @format A data frame containing the following variables, which are 
#'   interpreted as results from a hypothetical randomized trial. See the source 
#'   for a full description.
#' \describe{
#'   \item{out:}{ The binary outcome variable under the parallel design.}
#'   \item{out.enc:}{ The binary outcome variable under the parallel 
#'   encouragement design.}
#'   \item{med:}{ The binary mediator under the parallel design.}
#'   \item{med.enc:}{ The binary mediator under the parallel encouragement 
#'   design.}
#'   \item{ttt:}{ The binary treatment variable.}
#'   \item{manip:}{ The design indicator, or the variable indicating whether the 
#'   mediator is manipulated under the parallel design.}
#'   \item{enc:}{ The trichotomous encouragement variable under the parallel 
#'   encouragement design. Equals 0 if subject received no encouragement; 1 if 
#'   encouraged for the mediator value of 1; and -1 if encouraged for the 
#'   mediator value of 0.}
#' }
#' 
#' @details Conditioning on 'manip' = 0 will simulate a randomized trial under 
#'   the single experiment design, where 'out' and 'med' equal observed outcome 
#'   and mediator values, respectively.
#'   
#'   Unconditionally, using 'out', 'med', 'ttt' and 'manip' will simulate an 
#'   experiment under the parallel design.
#'   
#'   The 'out.enc' and 'med.enc' variables represent the outcome and mediator 
#'   values observed when subjects received the encouragement indicated in 
#'   'enc'. Therefore, using 'out.enc', 'med.enc', 'ttt' and 'enc' will simulate 
#'   an experiment under the parallel encouragement design.
#'   
#'   Note that all the observed responses are generated from an underlying 
#'   distribution of potential outcomes and mediators (not shown in this 
#'   dataset) satisfying the assumptions described in Imai, Tingley and 
#'   Yamamoto (2012). The full simulation code is available as a companion 
#'   replication archive for the article.
#' 
#' @source Imai, K., Tingley, D. and Yamamoto, T. (2012) Experimental Designs 
#'   for Identifying Causal Mechanisms. Journal of the Royal Statistical 
#'   Society, Series A (Statistics in Society).
#' 
#' @keywords datasets
"boundsdata"

#' Example Data for the Crossover Encouragement Design
#' 
#' A randomly generated dataset containing 2000 rows and 7 columns with no 
#' missing values.
#' 
#' @usage CEDdata
#' 
#' @format A data frame containing the following variables, which are 
#'   interpreted as results from a hypothetical randomized trial employing the 
#'   crossover encouragement design. 
#' \describe{
#'   \item{T1:}{ The binary treatment indicator in the first stage.}
#'   \item{M1:}{ The binary mediator variable recorded in the first stage.}
#'   \item{Y1:}{ The binary outcome variable recorded in the first stage.}
#'   \item{T2:}{ The binary treatment in the second stage. Equal to 1 - T1 by 
#'   design.}
#'   \item{Z:}{ The binary encouragement indicator for the second stage.}
#'   \item{M2:}{ The binary mediator recorded in the second stage.}
#'   \item{Y2:}{ The binary outcome recorded in the second stage.}
#' }
#' 
#' @details Note that all the observed responses are generated from an 
#'   underlying distribution of potential outcomes and mediators (not shown in 
#'   this dataset) satisfying the assumptions described in Imai, Tingley and 
#'   Yamamoto (2012).
#'   
#' @source Imai, K., Tingley, D. and Yamamoto, T. (2012) Experimental Designs 
#'   for Identifying Causal Mechanisms. Journal of the Royal Statistical 
#'   Society, Series A (Statistics in Society).
#'   
#' @keywords datasets
"CEDdata"

#' Brader, Valentino and Suhay (2008) Framing Experiment Data
#' 
#' The \code{framing} data contains 265 rows and 15 columns of data from a 
#' framing experiment conducted by Brader, Valentino and Suhay (2008).
#' 
#' @usage framing
#'   
#' @format A data frame containing the following variables: 
#' \describe{ 
#'   \item{immigr:}{ A four-point scale measuring subjects' attitudes toward
#'   increased immigration. Larger values indicate more negative attitudes.} 
#'   \item{english:}{ A four-point scale indicating whether subjects favor or 
#'   oppose a law making English the official language of the U.S.} 
#'   \item{cong_mesg:}{ Whether subjects requested sending an anti-immigration 
#'   message to Congress on their behalf.} \item{anti_info:}{ Whether subjects
#'   wanted to receive information from anti-immigration organizations.} 
#'   \item{tone:}{ 1st treatment; whether the news story is framed positively or
#'   negatively.} \item{eth:}{ 2nd treatment; whether the news story features a
#'   Latino or European immigrant.} \item{cond:}{ Four level measure recording
#'   joint treatment status of tone and eth.} \item{treat:}{ Product of the two
#'   treatment variables. In the original study the authors only find this cell
#'   to be significant.} \item{emo:}{ Measure of subjects' negative feeling
#'   during the experiment. A numeric scale ranging between 3 and 12 where 3
#'   indicates the most negative feeling.} \item{anx:}{ A four-point scale
#'   measuring subjects' anxiety about increased immigration.} \item{p_harm:}{
#'   Subjects' perceived harm caused by increased immigration. A numeric scale
#'   between 2 and 8.} \item{age:}{ Subjects' age.} \item{educ:}{ Subjects'
#'   highest educational attainments.} \item{gender:}{ Subjects' gender.} 
#'   \item{income:}{ Subjects' income, measured as a 19-point scale.} 
#' }
#'   
#' @source Brader, T., Valentino, N. and Suhay, E. (2008). What triggers public
#'   opposition to immigration? Anxiety, group cues, and immigration threat.
#'   American Journal of Political Science 52, 4, 959--978.
#'   
#' @keywords datasets
"framing"

#' JOBS II data
#' 
#' Job Search Intervention Study (JOBS II). JOBS II is a randomized field 
#' experiment that investigates the efficacy of a job trainingintervention on 
#' unemployed workers. The program is designed to not only increase 
#' reemploymentamong the unemployed but also enhance the mental health of the 
#' job seekers. In the JOBS IIfield experiment, 1,801 unemployed workers 
#' received a pre-screening questionnaire and were thenrandomly assigned to 
#' treatment and control groups. Those in the treatment group participatedin 
#' job-skills workshops. In the workshops, respondents learned job-search skills
#' and coping strategiesfor dealing with setbacks in the job-search process. 
#' Those in the control condition receiveda booklet describing job-search tips. 
#' In follow-up interviews, the two key outcome variables weremeasured; a 
#' continuous measure of depressive symptoms based on the Hopkins Symptom 
#' Checklist,and a binary variable, representing whether the respondent had 
#' become employed.
#' 
#' @usage jobs
#'   
#' @format A data matrix with 899 rows and 17 columns, containing no missing 
#'   values. The data are provided only for illustrative purposes and not for 
#'   inference about program efficacy, for which the original data source should
#'   be consulted. \describe{ \item{econ_hard:}{ Level of economic hardship 
#'   pre-treatment with values from 1 to 5.} \item{depress1:}{ Measure of 
#'   depressive symptoms pre-treatment.} \item{sex:}{ Indicator variable for 
#'   sex. 1 = female} \item{age:}{ Age in years.} \item{occp:}{ Factor with 
#'   seven categories for various occupations.} \item{marital:}{ Factor with 
#'   five categories for marital status.} \item{nonwhite:}{ Indicator variable 
#'   for race. 1 = nonwhite.} \item{educ:}{ Factor with five categories for 
#'   educational attainment.} \item{income:}{ Factor with five categories for 
#'   level of income.} \item{job_seek:}{ A continuous scale measuring the level 
#'   of job-search self-efficacy with values from 1 to 5. The mediator 
#'   variable.} \item{depress2:}{ Measure of depressive symptoms 
#'   post-treatment.} \item{work1:}{ Indicator variable for employment. 1 = 
#'   employed.} \item{job_dich:}{ The job_seek measure recoded into two 
#'   categories of high and low. 1 = high job search self-efficacy.} 
#'   \item{job_disc:}{ The job_seek measure recoded into four categories from 
#'   lowest to highest.} \item{treat:}{ Indicator variable for whether 
#'   participant was randomly selected for the JOBS II training program. 1 = 
#'   assignment to participation.} \item{comply:}{ Indicator variable for 
#'   whether participant actually participated in the JOBS II program. 1 = 
#'   participation.} \item{control:}{ Indicator variable for whether participant
#'   was randomly selected to not participate in the JOBS II training program. 1
#'   = non-participation.} }
#'   
#' @source The complete JOBS II data is available from the data archives at 
#'   www.icpsr.umich.edu/
#'   
#' @references Vinokur, A. and Schul, Y. (1997). Mastery and inoculation against
#'   setbacks as active ingredients in the jobs intervention for the
#'   unemployed. Journal of Consulting and Clinical Psychology 65(5):867-77.
#'   
#' @keywords datasets
"jobs"

#' School-level data
#' 
#' The original data source is the Education Longitudinal Study of 2002. To deal
#' with the issue on individually identifiable information, we generated 
#' hypothetical student-level data using a multiple imputation method. The 
#' Education Longitudinal Study of 2002 used a two-stage sample selection 
#' process. First, a national sample of schools was selected using stratified 
#' probability proportional to size (PPS), and school contacting resulted in 
#' 1,221 eligible public, Catholic, and other private schools from a population 
#' of approximately 27,000 schools containing 10th grade students. Of the 
#' eligible schools, 752 participated in the study. In the second stage of 
#' sample selection, a sample of approximately 26 sophomores, from within each 
#' of the participating public and private schools was selected. Each school was
#' asked to provide a list of 10th grade students, and quality assurance (QA) 
#' checks were performed on each list that was received.
#' 
#' @usage school
#'   
#' @format A data matrix with 568 rows and 5 columns, containing no missing 
#'   values. The data are provided only for illustrative purposes and not for 
#'   inference about education effectiveness, for which the original data source
#'   should be consulted. \describe{ \item{SCH_ID:}{ School indicator.} 
#'   \item{coed:}{ Indicator variable for coeducation. 1 = coeducation.} 
#'   \item{smorale:}{ Measure of student morale in the school. 4 levels.} 
#'   \item{free:}{ Percent of 10th grade students receiving free lunch. 1 to 7 
#'   levels.} \item{catholic:}{ Indicator variable for catholic school. 1 = 
#'   catholic school.} }
#'   
#' @source The complete student-level data is available from the data archives 
#'   at www.icpsr.umich.edu/
#'   
#' @references United States Department of Education. National Center for
#'   Education Statistics
#'   
#' @keywords datasets
"school"

#' Hypothetical student-level data
#' 
#' The original data source is the Education Longitudinal Study of 2002. To deal
#' with the issue on individually identifiable information, we generated 
#' hypothetical student-level data using a multiple imputation method. The 
#' Education Longitudinal Study of 2002 used a two-stage sample selection 
#' process. First, a national sample of schools was selected using stratified 
#' probability proportional to size (PPS), and school contacting resulted in 
#' 1,221 eligible public, Catholic, and other private schools from a population 
#' of approximately 27,000 schools containing 10th grade students. Of the 
#' eligible schools, 752 participated in the study. In the second stage of 
#' sample selection, a sample of approximately 26 sophomores, from within each 
#' of the participating public and private schools was selected. Each school was
#' asked to provide a list of 10th grade students, and quality assurance (QA) 
#' checks were performed on each list that was received.
#' 
#' @usage student
#'   
#' @format A data matrix with 9,679 rows and 17 columns, containing no missing 
#'   values. The data are provided only for illustrative purposes and not for 
#'   inference about education effectiveness, for which the original data source
#'   should be consulted. \describe{ \item{SCH_ID:}{ School indicator.} 
#'   \item{fight:}{ Indicator variable for fight at school. 1 = fight.} 
#'   \item{attachment:}{ Indicator variable for attachment to school. 1 = like.}
#'   \item{work:}{ Indicator variable for part-time job. 1 = work.} 
#'   \item{score:}{ Measure of math score.} \item{late:}{ Frequency in which the
#'   student was late for school. 5 levels.} \item{coed:}{ Indicator variable 
#'   for coeducation. 1 = coeducation.} \item{smorale:}{ Measure of student 
#'   morale in the school. 4 levels.} \item{gender:}{ Indicator variable for 
#'   gender. 1 = female.} \item{income:}{ Total family income. 16 levels.} 
#'   \item{free:}{ Percent of 10th grade students receiving free lunch. 1 to 7 
#'   levels.} \item{pared:}{ Parents highest level of education. 8 levels} 
#'   \item{catholic:}{ Indicator variable for catholic school. 1 = catholic 
#'   school.} }
#'   
#' @source The complete student-level data is available from the data archives 
#'   at www.icpsr.umich.edu/
#'   
#' @references United States Department of Education. National Center for
#'   Education Statistics
#'   
#' @keywords datasets
"student"

