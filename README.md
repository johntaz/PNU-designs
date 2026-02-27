# Prevalent new-user designs in R

This repository is for implementing the prevalent new-user design in the following study, under the time-conditional propensity score matching approach:

> Tazare J, Gibbons DC, Smeeth L, Ali MS, Gillespie IA, Cunnington M, Logie J, Williamson EJ, Douglas IJ,. Empirical Comparison of Exposure Set Definitions in the Prevalent New-User Design. **Pharmacoepidemiol Drug Saf**. 2026. doi: 10.1002/pds.70339.

## Repository Structure

The following scripts assume construction of a base cohort including incident new-users of Drug A and Drug B and switchers (i.e. prevalent new-users) from Drug B to Drug A. 

### Exposure Set Definitions

These scripts construct exposure sets based on the three definitions presented in the paper.

| Script | Exposure Set Type | Description |
|---|---|---|
| `RxExSets.R` | Prescription-based | Based on prior prescription counts |
| `TimeExSets.R` | Time-based | Based on time since base cohort entry (within a caliper) |
| `HybridExSets.R` | Hybrid | Based on prior prescriptions, time since cohort entry, and calendar year |

### Time-Conditional Propensity Score Models

These scripts fit time conditonal propensity score (TCPS) models across the relevant exposure sets. Scripts follow the naming convention `timeConditionalPS_<exposureset>.R`.

| Script | Exposure Set |
|---|---|
| `timeConditionalPS_Rx.R` | Prescription-based |
| `timeConditionalPS_Time.R` | Time-based |
| `timeConditionalPS_Hybrid.R`  | Hybrid |

### Matching

These scripts perform the TCPS matching procedure — matching users of Drug A to Drug B within exposure sets. Scripts follow the naming convention `matching_<exposureset>.R`.

| Script | Exposure Set |
|---|---|
| `matching_Rx.R` | Prescription-based |
| `matching_Time.R` | Time-based |
| `matching_Hybrid.R`  | Hybrid |

### Analysis

| Script | Description |
|---|---|
| `PNU_Analysis.R` | Example outcome analysis using one of the matched cohorts generated |

---

## Workflow Overview

For each exposure set definition, the intended workflow is:

```
Define exposure sets  →  Fit TCPS model  →  Perform matching  →  Analyse matched cohort
  (RxExSets.R etc.)     (_timeConditionalPS_*.R)  (_Matching_*.R)    (PNU_Analysis.R)
```

---
