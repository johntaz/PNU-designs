# -----------------------------------------------------------------------------
# PROGRAM NAME:  HybridExSets
# AUTHOR:        John Tazare
# NOTES:        Create Hybrid based exposure sets
# REQUIRES: Generation of a base cohort
# -----------------------------------------------------------------------------

debugMode = FALSE

set.seed(26005)
samplingFrac =  # for sampling comparators from the exposure sets 

library(plyr)
library(dplyr)

select <- dplyr::select

# Time window in Days 
timeWindow = 15

# Define functions  -------------------------------------------------------

# Hybrid exposure sets for incident study drug users 
exSets <- function(a) {
  if (a %% 500 == 0) {
    print(paste0("Started generation of exposure set ", a, ", out of ", nStudyDrugUsers))
  }
  
  tempSdUser <- studyDrugUsersInfo %>%
    arrange(eventdate) %>%
    slice(a) %>%
    mutate(exSetSubject=1) %>% 
    mutate(actSampleFrac =1)
    
if (tempSdUser$switcher == 0) {
  # Exposure Set Logic: ISDUsers # 
  
  # Each exposure set has all patients with the same no. of prior prescriptions 
  # of the comparator, occuring during the same calendar year and the time since cohort entry is +/- 15 days 
  
  # Find all the relevant comp users
  tempComparator1 <- comparatorsPool %>%
    filter(!patid %in% tempSdUser$patid) %>%
    filter(compRx==1) %>%
    filter(timeSinceIndex <= timeWindow & rxYear==studyDrugUsersInfo$rxYear[a] & numberCompRxs==1)
  
  tempComparator2 <- tempComparator1 %>% 
    sample_frac(samplingFrac) %>% # apply sampling fraction
    select(patid, eventdate)
  
}

if (tempSdUser$switcher == 1) {
    # Find all the relevant comp users
    tempComparator1 <- comparatorsPool %>%
      filter(!patid %in% tempSdUser$patid) %>%
      filter(compRx == 1) %>%
      filter(numberPriorCompRxs == studyDrugUsersInfo$numberPriorCompRxs[a]) %>%
      filter(rxYear== studyDrugUsersInfo$rxYear[a]) %>%
      mutate(maxTime = studyDrugUsersInfo$timeSinceIndex[a] + timeWindow) %>%
      mutate(minTime = studyDrugUsersInfo$timeSinceIndex[a] - timeWindow) %>% 
      filter(timeSinceIndex >= minTime  & timeSinceIndex <= maxTime) 
    
    tempComparator2 <- tempComparator1 %>% 
      sample_frac(samplingFrac) %>% # apply sampling fraction
      select(patid, eventdate)
  }
  
  tempSdUser <- tempSdUser %>%
    select(-c("rxYear", "switcher"))
  # Bind the two 
  tempExSet<- bind_rows(tempSdUser, tempComparator2) %>%
    mutate(exSetSubject = replace_na(exSetSubject, 0)) %>%
    mutate(actSampleFrac = replace_na(actSampleFrac, nrow(tempComparator2)/nrow(tempComparator1))) %>%
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
hybridExSets <- tempExSets %>%
  filter(!exposureSet %in% cleanExSets$exposureSet)

# Summary of number of individuals in each exposure set
hybridExSetsPtSummary <- tempExSetPtSummary %>%
  filter(!exposureSet %in% cleanExSets$exposureSet)

# Clean Up
rm(list=ls()[! ls() %in% c("hybridExSets", "debugMode", "hybridExSetsPtSummary", "excludedStudyDrugUsers")])

# Output data -------------------------------------------------------------
if(debugMode == TRUE) {
  write_csv(excludedStudyDrugUsers, "data/debug_hybridExcludedStudyDrugUsers.csv") 
} else{
  write_csv(excludedStudyDrugUsers, "data/hybridExcludedStudyDrugUsers.csv")
}

if(debugMode == TRUE) {
  write_csv(hybridExSets, "data/debug_hybridExSets.csv") 
} else{
  write_csv(hybridExSets, "data/hybridExSets.csv")
}

if(debugMode == TRUE) {
  write_csv(hybridExSetsPtSummary, "data/debug_hybridExSetsPtSummary.csv") 
} else{
  write_csv(hybridExSetsPtSummary, "data/hybridExSetsPtSummary.csv")
}
















