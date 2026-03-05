# -----------------------------------------------------------------------------
# PROGRAM NAME:  timeConditionalPS_Hybrid
# AUTHOR:        John Tazare
# NOTES:        Estimate the time-conditional propensity score for PNU analysis
#
# REQUIRES: Hybrid exSets to be run and exSetsInfo (dataset of characteristics 
# for inidividuals in the exposure sets)
# -----------------------------------------------------------------------------
# Flag to enable debugging; if TRUE, script will only be run on 1,000 patients
# selected at random
debugMode = FALSE

# Load libraries --------------------------------------------------------------
list <- c("survival")
new.packages <- list[!(list %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

if("package:plyr" %in% search()) detach("package:plyr", unload=TRUE)
library(dbplyr)
library(tidyverse)
library(lubridate)
library(scales)
library(DBI)
library(sparklyr)
library(survival)

select <- dplyr::select

# Load data -------------------------------------------------------------------

# Exposure set covariates
if (debugMode == TRUE) {
  if (file.exists("data/debug_exSetsInfo.csv")) {
    # Load in patient summary
    exSetsInfo <- read_csv("data/debug_exSetsInfo.csv") %>%
      mutate_at(c("patid"), as.character)
  } else {
    message("Error: debug_exSetsInfo.csv, not found.")
  }
} else{
  if (file.exists("data/exSetsInfo.csv")) {
    # Load in patient summary
    exSetsInfo <- read_csv("data/exSetsInfo.csv") %>%
      mutate_at(c("patid"), as.character)
  } else {
    message("Error: exSetsInfo.csv, not found.")
  }
}

# Fit conditional logistic regression model for each of the hybrid exposure sets 

  exSet = "hybridExSets"
  # Extract type for file names
  type = str_replace(exSet, "ExSets", "")
  
  # Exposure set data
  if (debugMode == TRUE) {
    if (file.exists(paste0("data/debug_", exSet, ".csv"))) {
      # Load in patient summary
      exSets <- read_csv(paste0("data/debug_", exSet, ".csv")) %>%
        mutate_at(c("patid"), as.character)
    } else {
      message(paste0("Error: debug_", exSet, ".csv", " not found."))
    }
  } else{
    if (file.exists(paste0("data/", exSet, ".csv"))) {
      # Load in patient summary
      exSets <- read_csv(paste0("data/", exSet, ".csv")) %>%
        mutate_at(c("patid"), as.character) 
    } else {
      message(paste0("Error: ", exSet, ".csv", " not found."))
    }
  }
  
  # merge covariate information to pooled exposure sets
  tempPnuCohort <- exSets %>%
    rename(indexdate = eventdate) %>%
    left_join(exSetsInfo, by = c("patid", "indexdate")) %>%
    mutate(patid = as.factor(patid)) %>%
    mutate(sampleWts = 1/actSampleFrac)
  
  # Vector of variables to summarize
  vars <- colnames(tempPnuCohort[,!names(tempPnuCohort) %in% c("patid",
                                                               "exposureSet",
                                                               "exSetSubject",
                                                               "indexdate",
                                                               "exposed",
                                                               "bmi",
                                                               "calendarYear",
                                                               "sampleWts", 
                                                               "actSampleFrac")])
  
  # Fit time-conditional propensity score
  # Specify model
  outcome <- "exposed"
  strata = "strata(exposureSet) + "
  vars = paste(vars, collapse = " + ")
  rhs = paste(strata, vars)
  
  psModelFunction <- as.formula(paste(outcome, rhs,
                                      sep = " ~ "))
  
  print(psModelFunction)
  
  # Conditional logistic regression stratified on exposure set to
  # estimate the time conditional propensity scores
  message(paste0("Fitting conditional logistic regression model for ", exSet, " at ", format(Sys.time(),"%H:%M:%S")))
  psModel <- clogit(psModelFunction , data = tempPnuCohort, method = c("approximate"), weights = tempPnuCohort$sampleWts)
  message(paste0("Fitted conditional logistic regression model for ", exSet, " at ", format(Sys.time(),"%H:%M:%S")))  
  
  #}
  
  summary(psModel)
  
#}
