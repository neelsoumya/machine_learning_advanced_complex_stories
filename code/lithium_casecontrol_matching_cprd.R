####################################################################################################################################
# Analysis for lithium CPRD control data
#   case control matching
#   
# reference britt Killan ~ca. April 2020
#   email:
#  I have downloaded the Gold data for your 23,070 cases, which are defined by the 10 creatinine and 31 lithium codes. 
# For the control population, defined by the creatinine codes and excluding ones with a lithium code, 
# there are nearly 5,975,651 potential patients (please find attached). From this list you would like data 
# for an age- and gender matched subset in a 10:1 ratio. CPRD advised that creating this list of controls should 
# be done by the research team. CPRD is providing denominator files which are needed for this, please find these here:
#  https://www.cprdrepos.phpc.cam.ac.uk/CODEBROWSERCPRD/Denominators_2020_03.zip
# If you can let me know when you have downloaded this file that would be great, as I will remove it again from this location.
# In the "acceptable patients" file you will find for each patient the gender and yob (year of birth) which can be used to
# create the list of controls. I have attached a presentation on the denominator files but please get in touch if you have
# any questions. Also attached is the index file that gives the patient list for your potential controls.
# 
###################################################################################################################################

##############################
# Load packages and settings
##############################
library(sqldf)
library(ggplot2)
# library(knitr)
# library(rmarkdown)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(png)
library(grid)
library(gridExtra)
library(lme4)
library(lmerTest)
# library(pheatmap)
LIBRARY_PREFIX <- "https://egret.psychol.cam.ac.uk/rlib/"
source(paste0(LIBRARY_PREFIX, "cris_common.R"))
# library(rstan)
# library(rstanarm)
# library(shiny)
# library(shinystan)
# library(parallel)
# library(bridgesampling)
# library(bayesplot)
# library(rpart)
# library(partykit)
library(stringr)
library(dplyr)
library(data.table)
library(lubridate)
# library(changepoint)

# rstan::rstan_options(auto_write=TRUE)
# options(mc.cores = parallel::detectCores())

##########################################################################################
# Functions by Shanquan Chen and Rudolf Cardinal to recalculate eGFR
##########################################################################################

#if NA exists, gfr = NA, else gfr is calculated based on formula
egfr_f <- function(scr, age, gender, ethnic){
  if(length(scr) == 1){
    gfr <- egfr_f_default(scr, age, gender, ethnic)
  } else {
    gfr <- egfr_f_df(scr, age, gender, ethnic)
  }
  return(gfr)
}

egfr_f_default <- function(scr, age, gender, ethnic){
  if(any(is.na(c(scr, age, gender, ethnic)))){
    gfr <- NA
  } else {
    if(str_sub(str_to_lower(gender), 1, 1) == "f"){
      k <- 0.7
      alpha <- -0.329
      female <- 1
    } else {
      k <- 0.9
      alpha <- -0.411
      female <- 1/1.018
    }
    
    black <- ifelse(str_sub(str_to_lower(ethnic), 1, 1) == "b", 1, 1/1.159)
    
    gfr <- 141 * min(scr/k, 1)^alpha * max(scr/k, 1)^(-1.209) * 0.993^age * 1.018 * female * 1.159 * black
    
  }
  return(gfr)
}

egfr_f_df <- function(scr, age, gender, ethnic){
  dat <- data.frame(scr, age, gender, ethnic)
  gfr <- dat %>% rowwise() %>% mutate(egfr = egfr_f_default(scr, age, gender, ethnic))
  gfr <- gfr$egfr
  return(gfr)
}


# example
# 
# mdat <- data.table(
#   scr = runif(100000, 30, 1000),
#   # ... normal range in micromolar is 62-115; plausible e.g. 30-1000
#   age = runif(100000, 0, 100),
#   gender = ifelse(runif(100000, 0, 1) >= 0.5, "male", "female"),
#   ethnic = ifelse(runif(100000, 0, 1) >= 0.5, "whi or oth", "bla")
# )

# mdat$egfr <- egfr_f(scr = mdat$scr, age = mdat$age, gender = mdat$gender, ethnic = mdat$ethnic)


# RNC:

millimolar_from_mg_per_dl <- function(mg_per_dl, molecular_mass_g_per_mol)
{
  # see crate_anon/nlp_manager/regex_units.py for working
  return(mg_per_dl * 10 / molecular_mass_g_per_mol)
}

mg_per_dl_from_millimolar <- function(millimolar, molecular_mass_g_per_mol)
{
  return(millimolar * molecular_mass_g_per_mol / 10)
}

CREATININE_MOLECULAR_MASS_G_PER_MOL = 113.12  # https://pubchem.ncbi.nlm.nih.gov/compound/creatinine

creatinine_micromolar_from_mg_per_dl <- function(creatinine_mg_per_dl)
{
  return(1000 * millimolar_from_mg_per_dl(creatinine_mg_per_dl,
                                          CREATININE_MOLECULAR_MASS_G_PER_MOL))
}

creatinine_mg_per_dl_from_micromolar <- function(creatinine_micromolar)
{
  return(mg_per_dl_from_millimolar(creatinine_micromolar / 1000,
                                   CREATININE_MOLECULAR_MASS_G_PER_MOL))
}

# Unit tests
# Check:
# creatinine_micromolar_from_mg_per_dl(1)  # 88.4017, correct as per https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470344484.app9
# creatinine_mg_per_dl_from_micromolar(1)  # 0.011312, correct as per https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470344484.app9

rnc_egfr_mdrd <- function(creatinine_micromolar, age_y, is_female, is_black)
{
  # https://www.niddk.nih.gov/health-information/communication-programs/nkdep/laboratory-evaluation/glomerular-filtration-rate/estimating
  # GFR (mL/min/1.73 m2) = 175 × (Scr)-1.154 × (Age)-0.203 × (0.742 if female) × (1.212 if African American)
  
  creatinine_mg_per_dl <- creatinine_mg_per_dl_from_micromolar(creatinine_micromolar)
  egfr_ml_per_min <- 
    175 *
    (creatinine_mg_per_dl ^ -1.154) *
    (age_y ^ -0.203) *
    ifelse(is_female, 0.742, 1) *
    ifelse(is_black, 1.212, 1)
  return(egfr_ml_per_min)
}

rnc_egfr_ckd_epi <- function(creatinine_micromolar, age_y, is_female, is_black)
{
  # https://www.niddk.nih.gov/health-information/communication-programs/nkdep/laboratory-evaluation/glomerular-filtration-rate/estimating
  # GFR = 141 × min (Scr /κ, 1)α × max(Scr /κ, 1)-1.209 × 0.993Age × 1.018 [if female] × 1.159 [if black]
  
  creatinine_mg_per_dl <- creatinine_mg_per_dl_from_micromolar(creatinine_micromolar)
  kappa <- ifelse(is_female, 0.7, 0.9)  # 0.7 female, 0.9 male
  alpha <- ifelse(is_female, -0.329, -0.411)
  egfr_ml_per_min <- 
    141 *
    pmin(creatinine_mg_per_dl / kappa, 1) ^ alpha *
    pmax(creatinine_mg_per_dl / kappa, 1) ^ -1.209 *
    0.993 ^ age_y *
    ifelse(is_female, 1.018, 1) *
    ifelse(is_black, 1.159, 1)
  return(egfr_ml_per_min)
}

# mdat$egfr_rnc_mdrd <- rnc_egfr_mdrd(creatinine_micromolar = mdat$scr, 
#                                     age_y = mdat$age,
#                                     is_female = mdat$gender == "female", 
#                                     is_black = mdat$ethnic == "bla")
# mdat$egfr_rnc_ckd_epi <- rnc_egfr_ckd_epi(creatinine_micromolar = mdat$scr,
#                                           age_y = mdat$age,
#                                           is_female = mdat$gender == "female", 
#                                           is_black = mdat$ethnic == "bla")


##############################
# Constants
##############################
DAYS_PER_YEAR <- 365.25


##############################
# Load data from Britt
##############################

# Load all data information
# this includes all data that is research acceptable
#   control + treatment group
dt_demographics <- data.table(
  read.csv('/Users/sbanerjee/Downloads/Denominators_2020_03/acceptable_pats_from_utspracts_2020_03.txt', #/Users/soumya/Downloads/Denominators_2020_03_1/acceptable_pats_from_utspracts_2020_03.txt',#  '  //cip-rg-bi01/banerjes/data_lithium/data_March2019/DR_B1523963675_lbm28_1050_NONPID_CONTROL_DEMOGS_NONPID_V2019-11-12.csv', 
           sep = '\t', header = TRUE, 
           stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
)

#########################
# Computational efficiency considerations
#########################
# randomly sample a few patients
set.seed(12)
#i_NUM_PATIENTS_SAMPLE = 10000
#temp_list_row_numbers = 1:nrow(dt_demographics)

# now sample some patients from them
#i_temp_sample_patient_index = sample(temp_list_row_numbers, i_NUM_PATIENTS_SAMPLE)
# CAUTION: overwrite demographics table
#dt_demographics <- dt_demographics[i_temp_sample_patient_index, ]

# then set key
setkey(dt_demographics, patid)


# load treatment group
dt_treatment_group <- data.table(
  read.csv('/Users/sbanerjee/Downloads/326861_define_cases_Define_results/define_cases_Define_results.txt', 
           sep = '\t', header = TRUE, 
           stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
)

#dt_treatment_group <- dt_treatment_group[
#  STUDY_SUBJECT_DIGEST %in% dt_demographics$STUDY_SUBJECT_DIGEST]

# quote="" # for dealing with ' etc in files
#dt_treatment_group[,
#                     year := lubridate::year(indexdate) ]

setkey(dt_treatment_group, patid)


# selected control
# load treatment group
dt_control_group <- data.table(
  read.csv('/Users/sbanerjee/Downloads/326898_define_controls_Define_results/define_controls_Define_results.txt', 
           sep = '\t', header = TRUE, 
           stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
)

#dt_treatment_group <- dt_treatment_group[
#  STUDY_SUBJECT_DIGEST %in% dt_demographics$STUDY_SUBJECT_DIGEST]

# quote="" # for dealing with ' etc in files
#dt_treatment_group[,
#                     year := lubridate::year(indexdate) ]

setkey(dt_control_group, patid)


#################################
# sanity checks on data
#################################
#sqldf::sqldf('select count(distinct(patid)) from dt_treatment_group')
#sqldf::sqldf('select count(distinct(patid)) from dt_demographics')
uniqueN(dt_demographics)
uniqueN(dt_treatment_group)
uniqueN(dt_control_group)

# merge on yob and gender and region

# Algorithm
# TODO: 
# 1. remove treatment group patients from dt_demographics
# 2. for every patient 
#      age gender region match control 
#      pick 10 at random without replacement
#    end for   
# 3. store in another table and save to disk
  
  
# index date crartinien date

# set seed


# now merge based on year

# Creatinine only:
#dt_meditech_labtests <- dt_meditech_labtests[
#  TestNameAbbreviated == 'CREA' | TestNameAbbreviated == 'CCREA' # *** change to the actual Epic value only
#  ]











##################################



# read in Epic
dt_epic_labtests <- data.table(
  read.csv('//cip-rg-bi01/banerjes/data_lithium/data_March2019/DR_B1523963675_lbm28_1050_NONPID_CONTROL_TEST_EPIC_v2019-11-12.csv', 
           sep = ',', header = TRUE, 
           stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
)

dt_epic_labtests <- dt_epic_labtests[
  STUDY_SUBJECT_DIGEST %in% dt_demographics$STUDY_SUBJECT_DIGEST]

# quote="" # for dealing with ' etc in files
setkey(dt_epic_labtests, STUDY_SUBJECT_DIGEST)
# Creatinine only:
dt_epic_labtests <- dt_epic_labtests[
  TestNameAbbreviated == 'CREA' | TestNameAbbreviated == 'CCREA' # *** change to the actual Epic value only
  ]



# read in prescription
dt_prescription <- data.table(
  read.csv('//cip-rg-bi01/banerjes/data_lithium/data_March2019/DR_B1523963675_lbm28_1050_NONPID_CONTROL_MEDICATION_v2019-11-12.csv', 
           sep = ',', header = TRUE, 
           stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
)
dt_prescription <- dt_prescription[
  STUDY_SUBJECT_DIGEST %in% dt_demographics$STUDY_SUBJECT_DIGEST]

setkey(dt_prescription, STUDY_SUBJECT_DIGEST)




# read in diagnosis
dt_diagnosis <- data.table(
  read.csv('//cip-rg-bi01/banerjes/data_lithium/data_March2019/DR_B1523963675_lbm28_1050_NONPID_CONTROL_DIAGNOSIS_PL_v2019-11-12.csv', 
           sep = ',', header = TRUE, 
           stringsAsFactors=FALSE, na.strings="..") # ,strip.white = TRUE)
)
dt_diagnosis <- dt_diagnosis[
  STUDY_SUBJECT_DIGEST %in% dt_demographics$STUDY_SUBJECT_DIGEST]

setkey(dt_diagnosis, STUDY_SUBJECT_DIGEST)






#########################
# Summary statistics
#########################

n_patients <- uniqueN(dt_demographics$STUDY_SUBJECT_DIGEST)


#########################
# Data munging
#########################

# Combine Meditech and Epic
list_tables  <- list(dt_epic_labtests, dt_meditech_labtests)
dt_all_tests <- data.table::rbindlist(list_tables, use.names = TRUE, fill = TRUE)

# TODO: pick columns that are needed and drop unncecessary columns
# dt_all_tests <- rbind(dt_epic_labtests, dt_meditech_labtests)


# FREE MEMORY On account of having a tiny machine:
rm(dt_epic_labtests)
rm(dt_meditech_labtests)

dt_all_tests[, creatinine := ResultNumeric]

# DUMMY DATA
#dt_meditech_labtests <- data.table(
#  STUDY_SUBJECT_DIGEST = c("s1", "s2", "s2", "s3", "s3", "s3", "s3"),
#  TestNameAbbreviated = c(rep("CREA", 6), "JUNK"),
#  result = 1:7
#)

# Create sum of tests for each patient
dt_patients_with_two_creatinines <- dt_all_tests[
  ,
  .(n_creatinine_tests = .N),
  by = 'STUDY_SUBJECT_DIGEST'  # by subject
  ]

# Patients who have at least 2 serum creatinine measurements
dt_patients_with_two_creatinines <- dt_patients_with_two_creatinines[n_creatinine_tests >= 2]
n_patient_with_two_creat_tests   <- nrow(dt_patients_with_two_creatinines)


# final list of patinets who 
#  1. with at least one Lithium prescripton event
#     OR
#     with at least one serum measurment of Lithium
# AND
# 2. Patients who have at least 2 serum creatinine measurements

# Code up gender
dt_demographics[, is_male   := GENDER_DESC == "Male"]
dt_demographics[, is_female := GENDER_DESC == "Female"]
# Restrict to males/females in age range 12-110 inclusive
dt_demographics <- dt_demographics[
  CURRENT_AGE >= 12 & CURRENT_AGE <= 110 & (is_male | is_female)
  ]
# Restrict to patients with >=2 creatinine measurements
dt_demographics <- merge(dt_demographics, dt_patients_with_two_creatinines, by="STUDY_SUBJECT_DIGEST")

n_patients_final <- nrow(dt_demographics)


##############################
# Summary of data
##############################

#df_summary_of_data <- data.frame(stringsAsFactors = FALSE )
#df_summary_of_data <- rbind ( df_summary_of_data, 
#                              cbind('Total number of patients in clinical data', n_patients )
#)
#df_summary_of_data <- rbind ( df_summary_of_data, 
#                              cbind('Total number of patients who have had Lithium at any point and at least two creatinine results', n_patient_with_two_creat_tests )
#)
#df_summary_of_data <- rbind ( df_summary_of_data, 
#                              cbind('Total number of patients who have had Lithium at any point and have had at least two creatinine results (EPIC only and are aged between 12 and 110 years and have known gender)', n_patients_final )
#                              )

#colnames(df_summary_of_data)[1] <- 'Feature'
#colnames(df_summary_of_data)[2] <- 'Number of patients'




##############################
# Data munging
##############################
# for control only look at creatinine

# Restrict all tests to those with two creatinines:
dt_all_tests <- merge(dt_all_tests, dt_patients_with_two_creatinines, by="STUDY_SUBJECT_DIGEST")


# date data munging
#dt_all_tests[
#  ,
#  ResultDate := lubridate::dmy_hms(ResultDate)
#]

dt_all_tests$ResultDate <- lubridate::make_date(year  = lubridate::year(dt_all_tests$ResultDate),
                                                month = lubridate::month(dt_all_tests$ResultDate),
                                                day   = lubridate::day(dt_all_tests$ResultDate)
)

#dt_all_tests[
#  ,
#  ResultDate := lubridate::make_date( year = lubridate::year(ResultDate), month = lubridate::month(ResultDate), day = lubridate::day(ResultDate) )
#  ]

dt_all_tests[
  ,
  ResultDate_asdate := as.Date(ResultDate)
  ]

# Calculate DOB assuming 1st of each month
dt_demographics$date_of_birth_approx <- lubridate::make_date(year  = dt_demographics$YEAR_OF_BIRTH,
                                                             month = dt_demographics$MONTH_OF_BIRTH,
                                                             day   = 1
)
#dt_demographics[
#  ,
#  approx_dob := lubridate::make_date(year=YEAR_OF_BIRTH, month=MONTH_OF_BIRTH, day=1)
#]

dt_tests_with_age <- merge(dt_all_tests, dt_demographics, by="STUDY_SUBJECT_DIGEST")


# FREE up more memory
rm(dt_all_tests)

# check what this is calculated as
# dt_tests_with_age$age_calculated_asof_test <- as.numeric( dt_tests_with_age$ResultDate - dt_tests_with_age$date_of_birth_approx )   

dt_tests_with_age[
  ,
  age_calculated_asof_test := as.numeric(ResultDate, date_of_birth_approx)
  ]

#dt_tests_with_age[
#  ,
#  age_at_test_y := as.numeric(difftime(ResultDate, approx_dob, unit="days")) / DAYS_PER_YEAR
#]

#dt_tests_with_age[
#  ,
#  egfr := SOME_EGFR_FUNCTION(creatinine, age_at_test, is_male, ...)
#]

# sort by date and get first test date
setkeyv(dt_tests_with_age, c("STUDY_SUBJECT_DIGEST", "ResultDate"))
dt_first_test_per_patient <- dt_tests_with_age[, .SD[1], by=c("STUDY_SUBJECT_DIGEST")] # *** check sorting
dt_first_test_per_patient[
  ,
  first_creatinine_datetime := ResultDate
  ]

# dt_first_test_per_patient$first_creatinine_datetime <- dt_first_test_per_patient$ResultDate

# *** merge dt_tests_with_age with dt_first_test_per_patient but only with the ResultDate column of the second,
# and rename it to "first_creatinine_datetime"

#  *** calculate time-since-first-creatinine as (basically) ResultDate - first_creatinine_datetime
dt_tests_with_age <- merge(dt_tests_with_age, dt_first_test_per_patient, by=c("STUDY_SUBJECT_DIGEST"))

dt_tests_with_age[
  ,
  time_off_computed := as.numeric( ResultDate.x - first_creatinine_datetime )
  ]

#dt_tests_with_age$time_off_computed <- as.numeric( dt_tests_with_age$ResultDate - dt_tests_with_age$first_creatinine_datetime ) #/ DAYS_PER_YEAR
#dt_tests_with_age$time_off_computed <- as.numeric(difftime(dt_tests_with_age$ResultDate, dt_tests_with_age$first_creatinine_datetime, unit="days")) #/ DAYS_PER_YEAR



# fix HACK
# df_epic_labtests_withage_egfr$time_off_computed <- as.numeric( df_epic_labtests_withage_egfr$ResultDate - df_epic_labtests_withage_egfr$date_of_birth_approx )   


##############################
# Calculate eGFR
##############################

# rename columns from GENDER_DESC.x to GENDER_DESC
setnames(dt_tests_with_age
         ,c("ResultNumeric.x", "age_calculated_asof_test.x", "GENDER_DESC.x" , "ETHNIC_GROUP_DESC.x")
         ,c("ResultNumeric"  , "age_calculated_asof_test"  , "GENDER_DESC"   , "ETHNIC_GROUP_DESC"  )
)

dt_tests_with_age[
  ,
  egfr_recalc := rnc_egfr_ckd_epi(creatinine_micromolar = ResultNumeric,
                                  age_y = age_calculated_asof_test/365, #dt_tests_with_age$age_at_test_y * DAYS_PER_YEAR, # ,
                                  is_female = GENDER_DESC == "Female", 
                                  is_black = ETHNIC_GROUP_DESC == "Black Caribbean")
  ]

# TODO: remove outliers or set to baseline at 200 ?


#dt_tests_with_age$egfr_recalc <- 
#  rnc_egfr_ckd_epi(creatinine_micromolar = dt_tests_with_age$ResultNumeric,
#                   age_y = dt_tests_with_age$age_calculated_asof_test/365, #dt_tests_with_age$age_at_test_y * DAYS_PER_YEAR, # ,
#                   is_female = dt_tests_with_age$GENDER_DESC == "Female", 
#                   is_black = dt_tests_with_age$ETHNIC_GROUP_DESC == "Black Caribbean")


##############################
# Plots and data visualization
##############################

theme_set(theme_classic())
gp <- ggplot(dt_tests_with_age, aes(x=time_off_computed, y=egfr_recalc) )
gp <- gp + geom_point()
gp <- gp + geom_smooth(span=0.3)
# gp <- gp + facet_wrap(GENDER_DESC ~ .)
gp <- gp + xlab('Days after Lithium withdrawal (days)')
gp <- gp + ylab('eGFR')
# gp <- gp + title('All patients who have taken Lithium at some point')
gp

#hist(df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime$egfr_recalc)
theme_set(theme_classic())
gp <- ggplot(dt_tests_with_age, aes(x=egfr_recalc))
gp <- gp + geom_histogram()
gp <- gp + xlab('eGFR')
gp <- gp + ylab('Count')
gp


hist(dt_tests_with_age$age_calculated_asof_test/DAYS_PER_YEAR)
hist(dt_tests_with_age$egfr_recalc)




##############################
# Free up more memory
##############################
rm(dt_first_test_per_patient)
rm(dt_patients_with_two_creatinines)
#rm(dt_prescription)
#rm(dt_demographics)
memory.limit()
memory.limit(size = 9000) # in MB
gc()

# remove unnecessary columns
dt_tests_with_age[,SpecimenId.x:=NULL]
dt_tests_with_age[,TestGroupName.x:=NULL]
dt_tests_with_age[,TestName.x:=NULL]
dt_tests_with_age[,TestNameAbbreviated.x:=NULL]
dt_tests_with_age[,ReferenceLow.x:=NULL]
dt_tests_with_age[,ReferenceHigh.x:=NULL]
dt_tests_with_age[,ResultUnit.x:=NULL]
dt_tests_with_age[,SUBMITTER_NAME.x:=NULL]
dt_tests_with_age[,SPECIMEN_NUMBER.x:=NULL]
dt_tests_with_age[,CURRENT_AGE.x:=NULL]
dt_tests_with_age[,CURRENT_AGE.y:=NULL]
dt_tests_with_age[,n_creatinine_tests.x.x:=NULL]
dt_tests_with_age[,n_creatinine_tests.y.x:=NULL]
dt_tests_with_age[,date_of_birth_approx.x:=NULL]
dt_tests_with_age[,SpecimenId.y:=NULL]
dt_tests_with_age[,TestGroupName.y:=NULL]
dt_tests_with_age[,TestName.y:=NULL]
dt_tests_with_age[,TestNameAbbreviated.y:=NULL]
dt_tests_with_age[,ResultValue.y:=NULL]
dt_tests_with_age[,ReferenceLow.y:=NULL]
dt_tests_with_age[,ReferenceHigh.y:=NULL]
dt_tests_with_age[,ResultUnit.y:=NULL]
dt_tests_with_age[,ResultDate.y:=NULL]
dt_tests_with_age[,LAST_UPDATED.y:=NULL]
dt_tests_with_age[,ResultUnit.y:=NULL]
dt_tests_with_age[,Method.y:=NULL]
dt_tests_with_age[,ResultFlag.y:=NULL]
dt_tests_with_age[,ResultNumeric.y:=NULL]
dt_tests_with_age[,ResultInRange.y:=NULL]
dt_tests_with_age[,ORDERING_DEPARTMENT_NAME.y:=NULL]
dt_tests_with_age[,COLLECTED_DATETIME.y:=NULL]
dt_tests_with_age[,ORDERED_DATETIME.y:=NULL]
dt_tests_with_age[,ResultValueCorrected.y:=NULL]
dt_tests_with_age[,RESULTING_LAB_NAME.y:=NULL]
dt_tests_with_age[,RESULTING_SECTION_NAME.y:=NULL]
dt_tests_with_age[,SUBMITTER_NAME.y:=NULL]
dt_tests_with_age[,SPECIMEN_NUMBER.y:=NULL]
dt_tests_with_age[,creatinine.y:=NULL]
dt_tests_with_age[,n_creatinine_tests.x.y:=NULL]
dt_tests_with_age[,ResultDate_asdate.y:=NULL]
dt_tests_with_age[,GENDER_DESC.y:=NULL]
dt_tests_with_age[,CURRENT_AGE.y:=NULL]
dt_tests_with_age[,YEAR_OF_BIRTH.y:=NULL]
dt_tests_with_age[,MONTH_OF_BIRTH.y:=NULL]
dt_tests_with_age[,DATE_OF_DEATH.y:=NULL]
dt_tests_with_age[,ETHNIC_GROUP_DESC.y:=NULL]
dt_tests_with_age[,is_male.y:=NULL]
dt_tests_with_age[,is_female.y:=NULL]
dt_tests_with_age[,n_creatinine_tests.y.y:=NULL]
dt_tests_with_age[,date_of_birth_approx.y:=NULL]
dt_tests_with_age[,age_calculated_asof_test.y:=NULL]



# Remove any rows with NAs 
#   or Inf in eGFR
# https://stackoverflow.com/questions/28878005/remove-rows-with-na-from-data-table-in-r
dt_tests_with_age <- dt_tests_with_age[complete.cases(dt_tests_with_age$egfr_recalc * 0), ]








##############################
# TODO: remove instances of AKI
##############################
## Alternative method of detecting AKI

b_detect_aki_later = FALSE

if (b_detect_aki_later == TRUE)
{
  
  
  i_change_percentage_egfr_aki_negative = -30
  i_change_percentage_egfr_aki_positive =  30
  i_acute_duration = 7 # days to detect acute incident in 
  
  
  # As another computational mechanism to remove instances of AKI, we utilize another approach to detect acute cases.
  # We use a percentage change threshold of `r toString(i_change_percentage_egfr_aki_negative)` and `r toString(i_change_percentage_egfr_aki_positive)` for eGFR values that occur within an acute duration of `r toString(i_acute_duration)` days.
  
  
  
  
  df_some_subject_responder = sqldf::sqldf( "select distinct(STUDY_SUBJECT_DIGEST) as subject_id 
                                                    from dt_tests_with_age")
  
  # df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_responder
  
  
  # creta e a new data frame
  df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki <- data.frame(NULL, stringsAsFactors = FALSE)
  
  
  for (temp_subject_id in df_some_subject_responder$subject_id)
  {
    
    # store this subject_id in a data frame
    temp_data_frame = data.frame(temp_subject_id, stringsAsFactors = FALSE)
    colnames(temp_data_frame)[1] <- "subject_id"
    
    #}
    
    # randomly se;ect some patients from a group (responder here) 
    i_num_subjects_temp = dim(df_some_subject_responder)[1]
    i_pick_subjects_rndom_index <- sample(0:i_num_subjects_temp, 1, replace = FALSE )
    # df_some_subject_responder_curated = data.frame( head(df_some_subject_responder$subject_id , n = 10) , stringsAsFactors = FALSE )
    df_some_subject_responder_curated = data.frame( df_some_subject_responder$subject_id[i_pick_subjects_rndom_index], stringsAsFactors = FALSE )
    
    colnames(df_some_subject_responder_curated)[1] <- "subject_id"
    
    
    df_some_subject_responder_curated_egfr <- 
      sqldf::sqldf("select *
                                   from dt_tests_with_age
                                   where STUDY_SUBJECT_DIGEST in 
                                   (select subject_id 
                                   from temp_data_frame)
                                   ")
    
    
    
    xxx <- diff( as.ts( df_some_subject_responder_curated_egfr$egfr_recalc  ) )
    
    
    # detect a greater than or less than 30% change 
    ts_change_detector <- diff(as.ts( df_some_subject_responder_curated_egfr$egfr_recalc  )) / as.ts( df_some_subject_responder_curated_egfr$egfr_recalc )  * 100
    
    # fraction remove threhdols
    # if less than 30% of previous value then remove that rows
    idx_to_remove_temp1 <- which(ts_change_detector < i_change_percentage_egfr_aki_negative)
    
    idx_to_remove_temp2 <- which(ts_change_detector > i_change_percentage_egfr_aki_positive)
    
    
    
    
    idx_temp_for_this_subject <- which(dt_tests_with_age$STUDY_SUBJECT_DIGEST %in% temp_subject_id)
    
    # assign this subject's rows to a new temorary data frame
    new_df <- dt_tests_with_age[idx_temp_for_this_subject,]
    
    # remove rows idx_to_remove_temp
    # check first if index is indeed greater than 0
    if( !is.null( dim(idx_to_remove_temp1)) || length(idx_to_remove_temp1)  )
    {
      if ( !is.null( dim(idx_to_remove_temp2)) || length(idx_to_remove_temp2) )
      {
        # additional logic to check if duration of interval is 1 week
        interval_temp <- as.Date( new_df$ResultDate[c(idx_to_remove_temp1,idx_to_remove_temp2)] )
        duration_interval_temp <- as.numeric( sum(diff(interval_temp)))
        
        # save data for this patient with AKI to plot later
        saved_patient_temp_withaki <- new_df
        
        # is this duraation less than 1 week
        if (duration_interval_temp <= i_acute_duration)
        {
          new_df <- new_df[-c(idx_to_remove_temp1,idx_to_remove_temp2),]
        }
      }
    }
    
    # finally at end of for loiop
    # create a new data frame
    # append to a new data frame
    df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki <- rbind(df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki, 
                                                                                new_df)
    
  } 
  
  
  # TODO: do regression with 
  df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki
  
  dim(df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime)
  
  dim(df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki)
  
  
  
  
  i_num_points_before_aki = dim(df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime)[1]
  i_num_points_after_aki  =  dim(df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki)[1]
  
  
  
  
  ### Characteristics of some patients after removing instances of AKI
  
  #Here are the characteristics of some random patients after removing instances of AKI.
  #In total, `r toString(i_num_points_before_aki - i_num_points_after_aki)` points were removed using this procedure.
  
  
  
  
  
  # ```{r, echo=FALSE, fig.cap="eGFR time course for some patients."}
  
  # TESTING 
  # how many points removed due to AKI
  
  df_some_subject_responder = sqldf::sqldf( "select distinct(STUDY_SUBJECT_DIGEST) as subject_id 
                                                    from df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki")# df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_responder
  
  # randomly se;ect some patients from a group (responder here) 
  i_num_subjects_temp = dim(df_some_subject_responder)[1]
  i_pick_subjects_rndom_index <- sample(0:i_num_subjects_temp, 4, replace = FALSE )
  # df_some_subject_responder_curated = data.frame( head(df_some_subject_responder$subject_id , n = 10) , stringsAsFactors = FALSE )
  df_some_subject_responder_curated = data.frame( df_some_subject_responder$subject_id[i_pick_subjects_rndom_index], stringsAsFactors = FALSE )
  
  colnames(df_some_subject_responder_curated)[1] <- "subject_id"
  
  # get  all their eGFR values
  df_some_subject_responder_curated_egfr_noaki <- 
    sqldf::sqldf("select *
                         from df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki
                         where STUDY_SUBJECT_DIGEST in 
                         (select subject_id 
                         from df_some_subject_responder_curated)
                         ")
  
  
  
  # pot them over time shopwing when they ceased Lithium
  theme_set(theme_classic())
  gp <- ggplot(df_some_subject_responder_curated_egfr_noaki, aes(x=time_off_negative_also, y=egfr_recalc) )
  gp <- gp + geom_point()
  gp <- gp + geom_smooth() # span=0.3
  gp <- gp + facet_wrap(STUDY_SUBJECT_DIGEST ~ .)
  gp <- gp + xlab('Time after lithium cessation (days)')
  gp <- gp + ylab('eGFR')
  # gp <- gp + title('Responders')
  gp
  
  
  
  
  ### Characteristics of a patient with AKI
  
  #An example patient profile with an instance of AKI is shown below.
  
  
  #```{r, echo=FALSE, fig.cap="eGFR time course of an example patient with AKI"}
  
  
  theme_set(theme_classic())
  gp <- ggplot(saved_patient_temp_withaki, aes(x=time_off_negative_also, y=egfr_recalc) )
  gp <- gp + geom_point()
  gp <- gp + geom_smooth() # span=0.3
  # gp <- gp + facet_wrap(STUDY_SUBJECT_DIGEST ~ .)
  gp <- gp + xlab('Time after lithium cessation (days)')
  gp <- gp + ylab('eGFR')
  # gp <- gp + title('Responders')
  gp
  
  
  
}







## Summary of clinical data after accounting for AKI


### Patients who have never been on Lithium 


# how many patients have died
#sqldf::sqldf("select count(distinct(STUDY_SUBJECT_DIGEST))
#             from df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime 
# where DATE_OF_DEATH in ('') ")

# df_summary_of_data <- data.frame(stringsAsFactors = FALSE )
# df_summary_of_data <- rbind ( df_summary_of_data, 
#                               cbind('Total number of patients in clinical data', df_total_patienst$count )
# )
# df_summary_of_data <- rbind ( df_summary_of_data, 
#                               cbind('Total number of patients who have had Lithium at any point', dim(df_all_patients_temp)[1] )
# )
# df_summary_of_data <- rbind ( df_summary_of_data, 
#                               cbind('Total number of patients who have had Lithium at any point and at least two creatinine results', dim(df_temp)[1] )
# )
# df_summary_of_data <- rbind ( df_summary_of_data, 
#                               cbind('Total number of patients who have had Lithium and come off Lithium', length(unique(df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime$STUDY_SUBJECT_DIGEST   )) )
# )
# df_summary_of_data <- rbind ( df_summary_of_data, 
#                               cbind('Total number of patients who have had Lithium, come off Lithium and had no AKI', length(unique(df_epic_labtests_withage_egfr_on_lithium_anytime_off_anytime_noaki$STUDY_SUBJECT_DIGEST   )) )
# )
# 
# 
# #df_summary_of_data <- rbind ( df_summary_of_data, 
# #                              cbind('Total number of patients who have had Lithium at any point and have had at least two creatinine results (EPIC only and are aged between 12 and 110 years)', dim(df_patients_lithium_creat_test_final)[1] )
# #                              )
# 
# colnames(df_summary_of_data)[1] <- 'Feature'
# colnames(df_summary_of_data)[2] <- 'Number of patients'
# 
# knitr::kable(df_summary_of_data, caption = "Summary of clinical data")



# In total, `r toString(i_num_points_before_aki - i_num_points_after_aki)` points were removed using this procedure.











##############################
# Linear models
##############################

# protoyping code
# TODO: unload and rm all and then load rds
# rm(list=ls())
# dt_tests_with_age <- readRDS(file="source_data/dt_tests_with_age.rds")

lm_model_ALLPATIENTS_real_aki_diab_hyper <- lmerTest::lmer(egfr_recalc ~  time_off_computed + (1 + time_off_computed | STUDY_SUBJECT_DIGEST), # + age + (1 + time_after |+  (1 + time_on | STUDY_SUBJECT_DIGEST) + (1 + time_off_computed | STUDY_SUBJECT_DIGEST)  STUDY_SUBJECT_DIGEST) 
                                                           data = dt_tests_with_age#, 
                                                           # control = lmerControl('bobyqa')#,
                                                           # nAGQ = 10
)
# TODO: try lme4
#lm_model_ALLPATIENTS_real_aki_diab_hyper <- lme4::lmer(egfr_recalc ~  time_off_computed + (1 + time_off_computed | STUDY_SUBJECT_DIGEST), # + age + (1 + time_after |+  (1 + time_on | STUDY_SUBJECT_DIGEST) + (1 + time_off_computed | STUDY_SUBJECT_DIGEST)  STUDY_SUBJECT_DIGEST) 
#                                                           data = dt_tests_with_age#, 
#                                                           # control = lmerControl('bobyqa')#,
#                                                           # nAGQ = 10
#                                                          )

# TODO: try REML i.e set REML = FALSE
#lm_model_ALLPATIENTS_real_aki_diab_hyper <- lme4::lmer(egfr_recalc ~  time_off_computed + (1 + time_off_computed | STUDY_SUBJECT_DIGEST), # + age + (1 + time_after |+  (1 + time_on | STUDY_SUBJECT_DIGEST) + (1 + time_off_computed | STUDY_SUBJECT_DIGEST)  STUDY_SUBJECT_DIGEST) 
#                                                       data = dt_tests_with_age,
#                                                       REML = FALSE
#                                                       # control = lmerControl('bobyqa')#,
#                                                       # nAGQ = 10
#                                                      )

# TODO: try REML i.e set REML = FALSE
#lm_model_ALLPATIENTS_real_aki_diab_hyper <- lmerTest::lmer(egfr_recalc ~  time_off_computed + (1 + time_off_computed | STUDY_SUBJECT_DIGEST), # + age + (1 + time_after |+  (1 + time_on | STUDY_SUBJECT_DIGEST) + (1 + time_off_computed | STUDY_SUBJECT_DIGEST)  STUDY_SUBJECT_DIGEST) 
#                                                           data = dt_tests_with_age, 
#                                                           REML = FALSE
#                                                           # control = lmerControl('bobyqa')#,
#                                                           # nAGQ = 10
#                                                          )


###################################
# Visualize parameters
###################################
cris$visualize_fixed_effects_from_lmer(lmer_result = lm_model_ALLPATIENTS_real_aki_diab_hyper)
cris$fixed_effects_from_lmer(lmer_result = lm_model_ALLPATIENTS_real_aki_diab_hyper)

df_coef_ALLPATIENTS_aki_diab_hyper <- as.data.frame( coef(summary(lm_model_ALLPATIENTS_real_aki_diab_hyper)) )
df_coef_ALLPATIENTS_aki_diab_hyper$variable_name   <- rownames(df_coef_ALLPATIENTS_aki_diab_hyper)
#zz4["time_on",]$Estimate
#zz4["(Intercept)",] <- 0

df_coef_ALLPATIENTS_filtered <- df_coef_ALLPATIENTS_aki_diab_hyper[2,] # ignore intercept and gender

# COMMENt by Shanquan Chen on reordering
df_coef_ALLPATIENTS_filtered$variable_name <- factor(df_coef_ALLPATIENTS_filtered$variable_name, 
                                                     levels = c("time_off_computed")
                                                     # levels = c("time_on", "time_off_computed", "GENDER_DESCMale", "age_calculated_asof_test")
)

#theme_set(theme_classic())
theme_set(theme_gray())
gp <- ggplot(data=df_coef_ALLPATIENTS_filtered, aes(x=variable_name, y=Estimate))
gp <- gp + geom_bar(stat = "identity")
gp <- gp + xlab("Feature Name")
gp <- gp + ylab("Estimated Coefficient")
gp <- gp + coord_flip()
gp



# TODO: subsample patients start with 10,000

# TODO: run from cmd line in Windows

# TODO: save fit model as rds object and potentialy move to CPFT

# TODO: as Vince for more RAM, more HD and more machines (chunk data into multiple machines)

# TODO: use rlibb::debugfunc::wtf() to see how much space is taken by each object

#############################################
# TODO: correct for diabetes hypertension
#############################################

# TODO: RECODE using %like%
# https://stackoverflow.com/questions/19149663/r-data-table-i-myvar-like-somethingsomethingelsesomethingmore

df_diabetes <- dt_diagnosis[ICD10_1 %like% 'E14% | E13% | E12% | E11% | E10%', ]

df_diabetes <- sqldf::sqldf("select distinct(STUDY_SUBJECT_DIGEST) as subject_id 
                            from dt_diagnosis 
                            where ICD10_1 like 'E14%' or ICD10_1 like 'E13%' or ICD10_1 like 'E12%' or ICD10_1 like 'E11%' or ICD10_1 like 'E10%' ")

# now create a diabetes flag
dt_tests_with_age$diabetes <- 0

# whoch patients have diabetes?
idx_to_change <- which(dt_tests_with_age$STUDY_SUBJECT_DIGEST %in% df_diabetes$subject_id)

# set their flags to 1
dt_tests_with_age$diabetes[idx_to_change] <- 1



# detect hypertension and renal function affected
df_hypertension <-
  sqldf::sqldf("select distinct(STUDY_SUBJECT_DIGEST) as subject_id
                from dt_prescription  
               where PHARM_CLASS in ('Diuretics') or PHARM_SUBCLASS in ('ACE inhibitors')  ")

# sqldf::sqldf("select distinct(STUDY_SUBJECT_DIGEST) from df_prescription  where   ")

# generate flag for hypertension
dt_tests_with_age$hypertension <- 0

# which subjects have hypertension with drugs which affect kidney fuunction
idx_to_change <- which( dt_tests_with_age$STUDY_SUBJECT_DIGEST %in%  df_hypertension$subject_id )

# change their hypertension flag to 1
dt_tests_with_age$hypertension[idx_to_change] <- 1



lm_model_ALLPATIENTS_real_aki_diab_hyper <- lmerTest::lmer(egfr_recalc ~  time_off_computed + diabetes + hypertension + (1 + time_off_computed | STUDY_SUBJECT_DIGEST), # + age + (1 + time_after |+  (1 + time_on | STUDY_SUBJECT_DIGEST) + (1 + time_off_computed | STUDY_SUBJECT_DIGEST)  STUDY_SUBJECT_DIGEST) 
                                                           data = dt_tests_with_age#, 
                                                           # control = lmerControl('bobyqa')#,
                                                           # nAGQ = 10
)

###################################
# Visualize parameters
###################################
cris$visualize_fixed_effects_from_lmer(lmer_result = lm_model_ALLPATIENTS_real_aki_diab_hyper)
cris$fixed_effects_from_lmer(lmer_result = lm_model_ALLPATIENTS_real_aki_diab_hyper)

df_coef_ALLPATIENTS_aki_diab_hyper <- as.data.frame( coef(summary(lm_model_ALLPATIENTS_real_aki_diab_hyper)) )
df_coef_ALLPATIENTS_aki_diab_hyper$variable_name   <- rownames(df_coef_ALLPATIENTS_aki_diab_hyper)
#zz4["time_on",]$Estimate
#zz4["(Intercept)",] <- 0

df_coef_ALLPATIENTS_filtered <- df_coef_ALLPATIENTS_aki_diab_hyper[2,] # ignore intercept and gender

# COMMENt by Shanquan Chen on reordering
df_coef_ALLPATIENTS_filtered$variable_name <- factor(df_coef_ALLPATIENTS_filtered$variable_name, 
                                                     levels = c("time_off_computed")
                                                     # levels = c("time_on", "time_off_computed", "GENDER_DESCMale", "age_calculated_asof_test")
)

#theme_set(theme_classic())
theme_set(theme_gray())
gp <- ggplot(data=df_coef_ALLPATIENTS_filtered, aes(x=variable_name, y=Estimate))
gp <- gp + geom_bar(stat = "identity")
gp <- gp + xlab("Feature Name")
gp <- gp + ylab("Estimated Coefficient")
gp <- gp + coord_flip()
gp


##################################
# TODO: correct for AKI
##################################

##################################
# save to disk
#################################

# TODO: later append df_epic_labtests_withage_egfr to treatment group data frame 
#       with time_on = 0, all other param = 0

# TODO: other fields related to lithium
df_epic_labtests_withage_egfr_recomputed_CONTROL                           = dt_tests_with_age
df_epic_labtests_withage_egfr_recomputed_CONTROL$time_on                   <- 0
df_epic_labtests_withage_egfr_recomputed_CONTROL$time_before               <- 0
df_epic_labtests_withage_egfr_recomputed_CONTROL$time_after                <- 0
df_epic_labtests_withage_egfr_recomputed_CONTROL$control                   <- 1
df_epic_labtests_withage_egfr_recomputed_CONTROL$time_off_computed_control <- df_epic_labtests_withage_egfr_recomputed_CONTROL$time_off_computed
df_epic_labtests_withage_egfr_recomputed_CONTROL$I_on_reversible           <- 0


# saveRDS(dt_tests_with_age, file = "source_data/dt_tests_with_age.rds")
saveRDS(df_epic_labtests_withage_egfr_recomputed_CONTROL, file = "source_data/df_epic_labtests_withage_egfr_recomputed_CONTROL.rds")



# TODO: use propensity score or other matching criterion to get case matched control