# -----------------------------------------------------------------------------
# PROGRAM NAME:  RxExSets
# AUTHOR:        John Tazare
# NOTES:        Create prescription based exposure sets
# REQUIRES: Generation of a base cohort
# -----------------------------------------------------------------------------
debugMode = FALSE

set.seed(26005)
samplingFrac =  # for sampling comparators from the exposure sets 

library(plyr)
library(dplyr)

select <- dplyr::select
# Define function -------------------------------------------------------------

# Rx-based exposure sets for incident study drug users 
exSets <- function(a) {
  if (a %% 500 == 0) {
print(paste0("Started generation of exposure set ", a, ", out of ", nStudyDrugUsers ))
  }
tempSdUser <- studyDrugUsersInfo %>%
  arrange(eventdate) %>%
  slice(a) %>%
  mutate(exSetSubject=1) %>% 
  mutate(actSampleFrac = 1)

if (tempSdUser$switcher == 0) {
  # Exposure Set Logic: ISDusers #

  # Each exposure set has all patients with 1 comparator prescription 
  # This is regardless of future no. of Rxs or if becomes swithcer 
  
tempComparator1 <- comparatorsPool %>% 
  filter(!patid %in% tempSdUser$patid) %>%
  filter(numberCompRxs == 1 )

tempComparator2 <- tempComparator1 %>% 
  sample_frac(samplingFrac) %>% # apply sampling fraction
  select(patid, eventdate) 

}

if (tempSdUser$switcher == 1) {
  # Exposure Set Logic: Switchers 
  
  # Each exposure set has all patients with the same number of 
  # comparator prescriptions as the switcher 
  
  tempComparator1 <- comparatorsPool %>% 
    filter(!patid %in% tempSdUser$patid) %>%
    filter(numberPriorCompRxs == tempSdUser$numberPriorCompRxs[1]) 
  
  tempComparator2 <- tempComparator1 %>%
    sample_frac(samplingFrac) %>%
    select(patid, eventdate)
}
  
  tempSdUser <- tempSdUser %>%
    select(-c("rxYear", "switcher"))
  # Bind the two 
  tempExSet<- bind_rows(tempSdUser, tempComparator2) %>%
    mutate(exSetSubject = replace_na(exSetSubject, 0)) %>%
    mutate(actSampleFrac = replace_na(actSampleFrac, nrow(tempComparator2)/nrow(tempComparator1))) %>% # calculate actual sampling fractions 
    mutate(exposureSet = a) %>% 
    select(patid, eventdate, exSetSubject, exposureSet, actSampleFrac)
  
  return(tempExSet)
  
}

# Load Data ------ ------------------------------------------------------------
if(debugMode == TRUE) {
    # Load in patient summary 
    baseCohort <- read_csv("data/debug_baseCohort.csv") %>% 
      mutate_at(c("patid"), as.character)
  } else {
    baseCohort <- read_csv("data/baseCohort.csv") %>% 
      mutate_at(c("patid"), as.character)
}

if(debugMode == TRUE) {
  # Load in patient summary 
  studyDrugUsersInfo <- read_csv("data/debug_studyDrugUsersInfo.csv") %>% 
    mutate_at(c("patid"), as.character)
} else {
  studyDrugUsersInfo <- read_csv("data/studyDrugUsersInfo.csv") %>% 
    mutate_at(c("patid"), as.character)
}

# Form the exposure sets -----------------------------------------------------

comparatorsPool <- baseCohort %>% 
  filter(!drug=="COX2i")
 
# study drug users
nStudyDrugUsers <- nrow(studyDrugUsersInfo) 
# Individual exposure sets for incident study drug users
indExSets <- lapply(X=1:nStudyDrugUsers, FUN=exSets) 
# Combine the individual exposure sets for incident study drug users
tempExSets <- rbind.fill(indExSets) 

# Number of individuals in each exposure set
tempExSetPtSummary <- tempExSets %>%
  dplyr::group_by(exposureSet) %>% 
  dplyr::summarise(nExposureSet = n()) 

# Identify any exposure sets where there's only one person i.e. no Comparators available
cleanExSets <- tempExSetPtSummary %>%
  filter(nExposureSet == 1) %>% 
  select(exposureSet)

# Exposure set subjects with no comparators 
excludedStudyDrugUsers <- tempExSets %>% 
  filter(exposureSet %in% cleanExSets$exposureSet & exSetSubject == 1) %>%
  select(patid)
  
# Drop these exposure sets from final exposure sets data set
rxExSets <- tempExSets %>%
  filter(!exposureSet %in% cleanExSets$exposureSet)

# Summary of number of individuals in each exposure set
rxExSetsPtSummary <- tempExSetPtSummary %>%
  filter(!exposureSet %in% cleanExSets$exposureSet)

# Clean Up
rm(list=ls()[! ls() %in% c("rxExSets", "debugMode", "rxExSetsPtSummary", "excludedStudyDrugUsers")])

# Output data ------------------------------------------------------------- 
if(debugMode == TRUE) {
  write_csv(excludedStudyDrugUsers, "data/debug_rxExcludedStudyDrugUsers.csv") 
} else{
  write_csv(excludedStudyDrugUsers, "data/rxExcludedStudyDrugUsers.csv")
}

if(debugMode == TRUE) {
  write_csv(rxExSets, "data/debug_rxExSets.csv") 
} else{
  write_csv(rxExSets, "data/rxExSets.csv")
}

if(debugMode == TRUE) {
  write_csv(rxExSetsPtSummary, "data/debug_rxExSetsPtSummary.csv") 
} else{
  write_csv(rxExSetsPtSummary, "data/rxExSetsPtSummary.csv")
}

