# -----------------------------------------------------------------------------
# PROGRAM NAME:  matching_Time
# AUTHOR:        John Tazare
# NOTES:        Perform the matching process for PNU analysis and return
#               the matched cohort with cohort index dates.
#
# REQUIRES: timeConditionalPS and pnuExclusions (i.e. data for establishing
# timing of exclusions/outcomes)
# -----------------------------------------------------------------------------

# Flag to enable debugging; if TRUE, script will only be run on 1,000 patients
# selected at random
debugMode = FALSE

# Functions ---------------------------------------------------------------------
# Time window in days
timeWindow = floor(365.25/12)

timeExSets <- function(a) {
  tempSdUser <- matchingOrder %>%
    slice(a) %>%
    mutate(exSetSubject=1)
  

  if (tempSdUser$switcher == 0) {
    # Exposure Set Logic: ISDUsers #
    
    # Each exposure set has all comparator prescriptions within 1 month of
    # cohort entry
    
    tempComparator <- patientPool %>% 
      filter(!patid %in% tempSdUser$patid) %>%
      filter(compRx==1 & timeSinceIndex <= timeWindow) %>% # 1month
      select(patid, eventdate)
    
  }
  
  if (tempSdUser$switcher == 1) {
    # Exposure Set Logic: Switchers #
    
    # Given the time interval betweeen the 'switch' prescription and index date  
    # for the switcher
    # Each exposure set finds all comparator prescriptions given to a patient with at most 1 month either side 
    # either of this duration
    
    tempComparator <- patientPool %>%
      filter(compRx==1) %>%
      filter(!patid %in% tempSdUser$patid) %>%
      mutate(maxTime = tempSdUser$timeSinceIndex[1] + timeWindow) %>%
      mutate(minTime = tempSdUser$timeSinceIndex[1] - timeWindow) %>% 
      filter(timeSinceIndex >= minTime  & timeSinceIndex <= maxTime) %>% 
      select(patid, eventdate)
  }
  
  tempSdUser <- tempSdUser %>%
    select(-c("rxYear", "switcher"))
  # Bind the two 
  tempExSet <- bind_rows(tempSdUser, tempComparator) %>%
    mutate(exSetSubject = replace_na(exSetSubject, 0)) %>%
    mutate(exposureSet = replace_na(tempSdUser$exposureSet[1])) %>% 
    select(patid, eventdate, exSetSubject, exposureSet)
  
  return(tempExSet)
}
  

# Load libraries --------------------------------------------------------------
list <- c("ggplot2", "gridExtra")
new.packages <- list[!(list %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

library(dbplyr)
library(tidyverse)
library(lubridate)
library(scales)
library(gridExtra)

# Load data -------------------------------------------------------------------
# Base cohort 
if(debugMode == TRUE) {
  # Load in patient summary 
  baseCohort <- read_csv("data/debug_baseCohort.csv") %>% 
    mutate_at(c("patid"), as.character)
} else {
  baseCohort <- read_csv("data/baseCohort.csv") %>% 
    mutate_at(c("patid"), as.character)
}

# studyDrug users
if(debugMode == TRUE) {
  # Load in patient summary 
  studyDrugUsersInfo <- read_csv("data/debug_studyDrugUsersInfo.csv") %>% 
    mutate_at(c("patid"), as.character)
} else {
  studyDrugUsersInfo <- read_csv("data/studyDrugUsersInfo.csv") %>% 
    mutate_at(c("patid"), as.character)
}

# Exclusion data
if (debugMode == TRUE) {
  if (file.exists("data/debug_flagCirrPts.csv")) {
    # Load in patient summary
    exclusionEvents <- read_csv("data/debug_flagCirrPts.csv") %>%
      mutate_at(c("patid"), as.character)
  } else {
    message("Error: debug_flagCirrPts.csv, not found.")
  }
} else{
  if (file.exists("data/flagCirrPts.csv")) {
    # Load in patient summary
    exclusionEvents <- read_csv("data/flagCirrPts.csv") %>%
      mutate_at(c("patid"), as.character)
  } else {
    message("Error: flagCirrPts.csv, not found.")
  }
}

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

# Exclude study drug users
if (debugMode == TRUE) {
  if (file.exists("data/debug_timeExcludedStudyDrugUsers.csv")) {
    # Load in patient summary
    excludedStudyDrugUsers <- read_csv("data/debug_timeExcludedStudyDrugUsers.csv") %>%
      mutate_at(c("patid"), as.character)
  } else {
    message("Error: debug_timeExcludedStudyDrugUsers.csv, not found.")
  }
} else{
  if (file.exists("data/timeExcludedStudyDrugUsers.csv")) {
    
    excludedStudyDrugUsers <- read_csv("data/timeExcludedStudyDrugUsers.csv") %>%
      mutate_at(c("patid"), as.character)
  } else {
    message("Error: timeExcludedStudyDrugUsers.csv, not found.")
  }
}

# load _exSets generated - for exSet Nos
if (debugMode == TRUE) {
  if (file.exists("data/debug_timeExSets.csv")) {
    # Load in patient summary
    exSetNums <- read_csv("data/debug_timeExSets.csv") %>%
      mutate_at(c("patid"), as.character) %>% 
      filter(exSetSubject == 1)
  } else {
    message("Error: debug_timeExSets.csv, not found.")
  }
} else{
  if (file.exists("data/timeExSets.csv")) {
    
    exSetsNums <- read_csv("data/timeExSets.csv") %>%
      mutate_at(c("patid"), as.character) %>% 
      filter(exSetSubject == 1)
  } else {
    message("Error: timeExSets.csv, not found.")
  }
}
# Matching procedure  ---------------------------------------------------------

# Patient pool
patientPool <- baseCohort %>% 
  filter(!drug=="COX2i")

# Identify order 
matchingOrder <- studyDrugUsersInfo %>%
  filter(!patid %in% excludedStudyDrugUsers$patid) %>%
  left_join(exSetsNums, by=c("patid", "eventdate")) %>%
  select(-c(exSetSubject)) %>% 
  arrange(eventdate) %>% 
  mutate(order = 1:n()) %>%
  mutate(matchingExclusionFlag = 0)  

# Fit for rxExSets exposure sets
  exSet = "timeExSets"
  # Extract type for file names
  type = str_replace(exSet, "ExSets", "")
  
  # NOTE: These steps must be carried out in chronological order
  #       sorted by the exposure set subject with the earliest index date.

  # Matching Loop
  for (row in 1:nrow(matchingOrder)) {

    # Check if study drug user is excluded during matching process
    if (nrow(matchingOrder %>% slice(row) %>% filter(matchingExclusionFlag == 1)) == 1) {
      next
    }
    
    # Since patients will be removed from matchingOrder need to check this and
    # move onto next set if already dropped
    if (nrow(matchingOrder %>% filter(order == row)) == 0) {
      next
    }
    
    if (row %% 50 == 0) {
      print(paste(
        "Matching exposure set ",
        row,
        " out of ",
        nrow(matchingOrder)
      ))
    } 

    # Form Exposure set
    tempExSet <-  timeExSets(row)
    
    # Check there is at least one possible match and if not go to next (important for near the end of matching cycle)
    if (nrow(tempExSet) <= 1){
      next
    }
    
    # left_join exclusion criteria and if study drug user has exclusion criteria n
    tempExSet1 <- tempExSet %>%
      group_by(patid) %>% 
      left_join(exclusionEvents) %>% 
      mutate(exclusionFlag = 
               ifelse(possExclusionEvent == 1 & exclusionDate <= eventdate, TRUE, FALSE))
  
    # If exposure set subject meets exclusion criteria:
    if (tempExSet1$exclusionFlag[1] == TRUE) {
      # Exclude patid out right from ANY selection
      excludePatid <- tempExSet1 %>%
        filter(exSetSubject == 1) %>%
        slice(1)
      
      # Output message
      message(
        exSet,
        ": Exposure set subject meets exclusion criteria. Exposure set excluded from analysis. Patient ",
        excludePatid$patid,
        " is excluded from all further inclusion."
      )
      
      # Moves to next element of the loop
      patientPool <- patientPool %>% 
        filter(!patid %in% excludePatid$patid)
      }
  
# Estimate ps for all patients in the exposure set
   else{
     # Join covariate information
     tempExSet2 <- tempExSet1 %>% 
       rename(indexdate = eventdate) %>%
       group_by(patid, indexdate) %>%
       left_join(exSetsInfo) %>%
       ungroup()
    
     # Estimate propensity score
     tempExSet2$odds <- exp(predict(psModel, type = "lp", newdata = tempExSet2))
     tempExSet2$pscore <- tempExSet2$odds / (1 + tempExSet2$odds)
     
     
     
     # save exset subject info
     exSetSubjectInfo <-  tempExSet2 %>%
       filter(exSetSubject == 1) %>%
       select(patid,
              indexdate,
              exposureSet,
              exSetSubject,
              exclusionFlag,
              pscore)
     
     # Step 1: The exposure set subject's PScore must lie within the range of
     # PScores for that exposure set
     
     # Summarise PScore for all comparators
     exSetSummaryPS <- tempExSet2 %>%
       filter(exSetSubject == 0) %>% # Only include comparators
       summarise(minPscore = min(pscore),
                 maxPscore = max(pscore))
     
     # Identify Exposure Sets where positivity assumption is violated
     exSetsPositivityViolation <- select(tempExSet2, patid, exSetSubject, pscore) %>%
       filter(exSetSubject ==1) %>%
       left_join(exSetSummaryPS, by=character()) %>%
       mutate(positivityViolation =
                ifelse(pscore < minPscore |
                         pscore > maxPscore , TRUE, FALSE)) 
     
     
     # Remove exposureSet(s) where positivity assumption is violated
     if (exSetsPositivityViolation$positivityViolation[1] == TRUE) {
       next
    
        }
     
     # Generate matched cohort
     if (row == 1) {
       pnuMatched <- data.frame(
         patid = as.character(),
         indexdate = as.character(),
         exposureSet = as.numeric(),
         exSetSubject = as.numeric(),
         exclusionFlag = as.logical(),
         pscore = as.numeric()
       ) %>%
         mutate_at(c("indexdate"), dmy)
     }
     
     # Generate relevant info for comparators in the exposure set
     comparatorInfo <- tempExSet2 %>%
       arrange(-exSetSubject) %>%
       mutate(exSetSubjectPscore = first(pscore)) %>%
       # Generate measure of distance from comparator Pscore to ex set subject pscore
       filter(exSetSubject == 0) %>%
       mutate(psDistance = abs(exSetSubjectPscore - pscore)) %>%
       select(patid, indexdate, exposureSet, exSetSubject, exclusionFlag, pscore, psDistance) %>%
       # Arrange by closest pscore
       arrange(psDistance)
     
     comparatorCandidate <- comparatorInfo %>%
       # Identify comparator with closet pscore
       slice(1) %>%
       select(patid,
              indexdate,
              exposureSet,
              exSetSubject,
              exclusionFlag,
              pscore)
     
     
     # If comparator meets exclusion criteria
     if (comparatorCandidate$exclusionFlag[1] == TRUE) {
       # 1. They are excluded from any further selection - as comparator and/or exSetSubject
       patientPool <- patientPool %>%
         filter(!patid %in% comparatorCandidate$patid)
       
       # Update matchingOrder to reflect that this patient can no longer be matched as a study drug user
       matchingOrder <- matchingOrder %>%
         mutate(matchingExclusionFlag = replace(matchingExclusionFlag, patid ==  comparatorCandidate$patid[1], 1))
       
       # 2. Attempt to find an alternative match
       for (rowComp in 2:nrow(comparatorInfo)) {
         comparatorCandidate <- comparatorInfo %>%
           slice(rowComp) %>%
           select(patid,
                  indexdate,
                  exposureSet,
                  exSetSubject,
                  exclusionFlag,
                  pscore)
         
         # If 2:nth person meets exclusion criteria
         if (comparatorCandidate$exclusionFlag[1] == TRUE) {
           # 1. They are excluded from any further selection - as comparator and/or exSetSubject
           patientPool <- patientPool %>%
             filter(!patid %in% comparatorCandidate$patid)
          
           # Update matchingOrder to reflect that this patient can no longer be matched as a study drug user
           matchingOrder <- matchingOrder %>%
             mutate(matchingExclusionFlag = replace(matchingExclusionFlag, patid ==  comparatorCandidate$patid[1], 1))
           
         } else{
         # Append matched pair
         pnuMatched <- pnuMatched %>%
           add_row(
             patid = exSetSubjectInfo$patid ,
             indexdate = exSetSubjectInfo$indexdate ,
             exposureSet = exSetSubjectInfo$exposureSet ,
             exSetSubject = exSetSubjectInfo$exSetSubject,
             exclusionFlag = exSetSubjectInfo$exclusionFlag,
             pscore = exSetSubjectInfo$pscore
           ) %>%
           
           add_row(
             patid = comparatorCandidate$patid ,
             indexdate = comparatorCandidate$indexdate ,
             exposureSet = comparatorCandidate$exposureSet ,
             exSetSubject = comparatorCandidate$exSetSubject,
             exclusionFlag = comparatorCandidate$exclusionFlag,
             pscore = comparatorCandidate$pscore
           )
         
         # Exclude comparator from being matched again (can still be a exSetSubject)
         patientPool <- patientPool %>%
           filter(!patid %in% comparatorCandidate$patid) 
         # move out of loop
         break
         }
         # does compRow = nrow of comparator info if yes, next
     }
   }
     # Else comparator has no reason to be excluded so match and remove from further exposure sets
     else {
       # Append matched pair
       pnuMatched <- pnuMatched %>%
         add_row(
           patid = exSetSubjectInfo$patid ,
           indexdate = exSetSubjectInfo$indexdate ,
           exposureSet = exSetSubjectInfo$exposureSet ,
           exSetSubject = exSetSubjectInfo$exSetSubject,
           exclusionFlag = exSetSubjectInfo$exclusionFlag,
           pscore = exSetSubjectInfo$pscore
         ) %>%
         
         add_row(
           patid = comparatorCandidate$patid ,
           indexdate = comparatorCandidate$indexdate ,
           exposureSet = comparatorCandidate$exposureSet ,
           exSetSubject = comparatorCandidate$exSetSubject,
           exclusionFlag = comparatorCandidate$exclusionFlag,
           pscore = comparatorCandidate$pscore
         )
       
       # Exclude comparator from being matched again (can still be a exSetSubject)
       patientPool <- patientPool %>%
         filter(!patid %in% comparatorCandidate$patid)
      
     }
     
     
   } 
    
}
  
  # Plot Overlap
  matched <- pnuMatched %>%
    mutate(trtlabel = ifelse(exSetSubject == 1,
                             yes = 'COX-2 Inhibitor',
                             no = 'NSAID')) %>%
    ggplot(aes(x = pscore, linetype = trtlabel)) +
    scale_linetype_manual(values = c("longdash", "dotted")) +
    geom_density(alpha = 0.5) +
    xlab('Probability of receiving COX-2 inihibitor') +
    ylab('Density') +
    scale_fill_discrete('') +
    scale_color_discrete('') +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    theme(strip.text = element_text(colour = 'black')) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    theme(
      legend.position = c(0.2, .9),
      legend.direction = 'vertical',
      panel.background = element_rect(fill = "white", colour = "white"),
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  matched
  
  # Add ExSet type paste
  if (debugMode == TRUE) {
    ggsave(paste0("output/debug_", type, "OverlapMatchedPNU.tiff"),
           device = tiff())
  } else{
    ggsave(paste0("output/", type, "OverlapMatchedPNU.tiff"),
           device = tiff())
  }
  
  # Clean up
  
  # Output data -------------------------------------------------------------
  if (debugMode == TRUE) {
    write_csv(data.frame(pnuMatched),
              paste0("data/", type, "DebugPnuMatched.csv"))
  } else{
    write_csv(data.frame(pnuMatched),
              paste0("data/", type, "PnuMatched.csv"))
  }
  
