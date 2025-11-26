# Imputing Boat-Based Fishing Effort (Phone Diary Survey) from Vessel Activity (Remote Camera Data)

Primary author: Stephen Bradshaw

Contact: stephen.bradshaw@dpird.wa.gov.au

Date: 2025-11-11

This repository supports a publication (titled above) to estimate boat-based fishing effort within a Bayesian framework using STAN (executed using the statistical programming language R).

## Introduction
This repo provides a series of scripts:

(i) demo_Code.R 
  - Loads simulated data
  - Establishes structure for STAN execution
  - Execute STAN code to estimate boat-based fishing
  - Report on metrics and plot

(ii) demo_rSTAN.R
  - Imputation model

## Files
```plaintext
Parent Folder
├── 00_src
│   └── functions.R
├── 01_data
│   └── demoData.rds
├── 02_scripts
│   ├── demo_Code.R
│   └── demo_rSTAN.R
└── outputs
    ├── plot_imputed_DEMO BOAT RAMP.png
    ├── stanOuts_20251111.rds
    └── tab_summStats.csv
```
