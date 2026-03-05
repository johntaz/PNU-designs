# -----------------------------------------------------------------------------
# PROGRAM NAME:  PNU_Analysis
# AUTHOR:        John Tazare
# NOTES:        Take the matched cohort and merge on the follow up
#               perform -> perform the cox analysis
#
# REQUIRES: pnuFollowUp pnuMatchedCohort
# -----------------------------------------------------------------------------
# Flag to enable debugging; if TRUE, script will only be run on 1,000 patients
# selected at random
debugMode = FALSE

# Load libraries --------------------------------------------------------------
list <- c("survival", "tableone")
new.packages <- list[!(list %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

library(dbplyr)
library(tidyverse)
library(lubridate)
library(scales)
library(DBI)
library(sparklyr)
library(survival)
require(devtools)
library(tableone)

select <- dplyr::select

# Load data -------------------------------------------------------------------

# Load ISD and switcher infoif(debugMode == TRUE) {
if(debugMode == TRUE) {
  # Load in patient summary 
  studyDrugUsersInfo <- read_csv("data/debug_ studyDrugUsersInfo.csv") %>% 
    mutate_at(c("patid"), as.character) %>% 
    rename(indexdate = eventdate)
} else {
  studyDrugUsersInfo <- read_csv("data/studyDrugUsersInfo.csv") %>% 
    mutate_at(c("patid"), as.character) %>% 
    rename(indexdate = eventdate)
}

# Loop over exposure sets (time, hybrid, rx)
exSetType <- c("rxExSets", "timeExSets", "hybridExSets")

for (val in exSetType) {
  exSet = val
  # Extract type for file names
  type = str_replace(exSet, "ExSets", "")
  
  # Matched data
  if (debugMode == TRUE) {
    if (file.exists(paste0("data/", type, "DebugPnuMatched.csv"))) {
      # Load in data
      pnuMatched <-
        read_csv(paste0("data/", type, "DebugPnuMatched.csv")) %>%
        mutate_at(c("patid"), as.character)
    } else {
      message(paste0("Error: ", type, "DebugPnuMatched.csv", " not found."))
    }
  } else{
    if (file.exists(paste0("data/", type, "PnuMatched.csv"))) {
      # Load in data
      pnuMatched <-
        read_csv(paste0("data/", type, "PnuMatched.csv")) %>%
        mutate_at(c("patid"), as.character)
    } else {
      message(paste0("Error: ", type, "PnuMatched.csv", " not found."))
    }
  }
  
  # Follow up data
  if (debugMode == TRUE) {
    if (file.exists(paste0("data/", type, "DebugFollowUpPNU.csv.csv"))) {
      # Load in data
      pnuFollowUp <-
        read_csv(paste0("data/", type, "DebugFollowUpPNU.csv")) %>%
        mutate_at(c("patid"), as.character)
    } else {
      message(paste0("Error: ", type, "DebugFollowUpPNU.csv", " not found."))
    }
  } else{
    if (file.exists(paste0("data/", type, "FollowUpPNU.csv"))) {
      # Load in data
      pnuFollowUp <-
        read_csv(paste0("data/", type, "FollowUpPNU.csv")) %>%
        mutate_at(c("patid"), as.character)
    } else {
      message(paste0("Error: ", type, "FollowUpPNU.csv", " not found."))
    }
  }
  # Join
  pnuMatchedCohort <- pnuFollowUp %>%
    left_join(pnuMatched) %>%
    rename(exposed = exSetSubject)
  
  # Estimate Treatment Effect -----------------------------------------------
  
  # Matched
  dataMatched <- pnuMatchedCohort %>%
    mutate(time = (censorDate - indexdate) + 1) # to quantify time at risk
  
  dataMatched6mth <- pnuMatchedCohort %>% 
    mutate(time = (censorDate - indexdate) + 1) %>%
    mutate(outcome6mth = ifelse(time > 183, outcome*0, outcome)) %>% 
    mutate(time6mth = ifelse(time> 183 , 183, time)) 
  
  
  # Matched Analysis
  anPnuMatched <- coxph(Surv(time, outcome) ~ exposed , dataMatched)
  summary(anPnuMatched)
  
  anPnuMatched6mth <- coxph(Surv(time6mth, outcome6mth) ~ exposed , dataMatched6mth)
  summary(anPnuMatched)
  
  # Subgroup analysis - keeping only incident/prevalent matched sets
  dataSubgroup <- dataMatched %>%
    left_join(select(studyDrugUsersInfo, patid, indexdate, switcher)) %>%
    group_by(exposureSet) %>% 
    mutate(switcher = replace_na(switcher, 0)) %>% 
    mutate(subgroupFlag = max(switcher)) %>%
    ungroup()
  
  dataPrevSubgroup <- dataSubgroup %>% 
    filter(subgroupFlag==1)
  
  pnuPrevSubgroup <-coxph(Surv(time, outcome) ~ exposed, dataPrevSubgroup)
  summary(pnuPrevSubgroup)
  
  dataInciSubgroup <- dataSubgroup %>% 
    filter(subgroupFlag!=1)
  
  pnuInciSubgroup <-coxph(Surv(time, outcome) ~ exposed, dataInciSubgroup)
  summary(pnuInciSubgroup)
  
  
  dataSubgroup6mth <- dataMatched6mth %>%
    left_join(select(studyDrugUsersInfo, patid, indexdate, switcher)) %>%
    group_by(exposureSet) %>% 
    mutate(switcher = replace_na(switcher, 0)) %>% 
    mutate(subgroupFlag = max(switcher)) %>%
    ungroup()
  
  dataPrevSubgroup6mth <- dataSubgroup6mth %>% 
    filter(subgroupFlag==1)
  
  pnuPrevSubgroup6mth <-coxph(Surv(time6mth, outcome6mth) ~ exposed, dataPrevSubgroup6mth)
  summary(pnuPrevSubgroup6mth)
  
  dataInciSubgroup6mth <- dataSubgroup6mth %>% 
    filter(subgroupFlag!=1)
  
  pnuInciSubgroup6mth <-coxph(Surv(time6mth, outcome6mth) ~ exposed, dataInciSubgroup6mth)
  summary(pnuInciSubgroup6mth)
  
  # Results
  pnuResults <- list(Matched  = ShowRegTable(anPnuMatched, printToggle = FALSE),
                     SubgroupPrev = ShowRegTable(pnuPrevSubgroup, printToggle = FALSE),
                     SubgroupInci = ShowRegTable(pnuInciSubgroup, printToggle = FALSE),
                     Matched6mth  = ShowRegTable(anPnuMatched6mth, printToggle = FALSE),
                     SubgroupPrev6mth = ShowRegTable(pnuPrevSubgroup6mth, printToggle = FALSE),
                     SubgroupInci6mth = ShowRegTable(pnuInciSubgroup6mth, printToggle = FALSE))
  print(pnuResults, quote = FALSE)
  
  # Output data -------------------------------------------------------------
  # Improve appearance of these
  if (debugMode == TRUE) {
    write_csv(data.frame(pnuResults),
              paste0("output/", type, "DebugPnuResults.csv"))
  } else{
    write_csv(data.frame(pnuResults),
              paste0("output/", type, "PnuResults.csv"))
  }
  
}
