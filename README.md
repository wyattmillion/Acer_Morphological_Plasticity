# Acer_Morphological_Plasticity

This repository contains all raw data and R scripts used in "Evidence for adaptive morphological plasticity in the Caribbean coral, Acropora cervicornis." 


Description of files:

1) MorphologicalPlasticityAnalysis.R 
  - contains all of the data organization and statistical analysis for all trait related tests including cox proportional hazard models, linear mixed effect models for absolute size and growth, joint regression for plasticity, and various tests on fragmentation
  - This script uses Acer_Growth_Data.csv or AcerMorphologyData_Imputed2.csv

2) EnvironmentalAnalysis_Complement.R 
  - contains code for generating maps, calculating temperature characteristics, extracting water quality metrics from SERC dataset, associated statistical analysis of water quality, and all Bayesian-related steps
  - this script used AcerMorphologyData_Imputed2.csv, OutplantTempData.csv, serc.csv, SERC_WQ_AUG2020_uM.csv, MCMCSupportHighStatV4.R

3) serc.csv 
  - Coordinates of SERC water quality monitoring stations associated with outplant reef sites 

4) Acer_Growth_Data.csv 
  - raw trait data for all ramets throughout the entire experimental period. TLE, surface area, volume, volume and surface area of the convex hull, status (alive,dead, or missing) and fragmentation based of visual assessments are included for each survey point (T0, T3 months, T6 months, T9 months, T12 months)
  - also includes survival status and time of death (TOD) in relation to T12 status to be used for cox proportional hazard models

5) AcerMorphologyData_Imputed2.csv
  - includes all data from Acer_Growth_Data.csv but also has data for ramets at T0 that did not have measures for 3D traits because model building was unsuccessful and for T3 ramets at EDR whose photographs were lost

6) MCMCSupportHighStatV4.R
  - support file for Bayesian negative binomial generalized linear model
 
7) OutplantTempData.csv
  - hourly temperature data over the experimental period (April 2018 to April 2019) for each reef site used here

8) SERC_WQ_AUG2020_uM.csv
  - biogeochemical microMolar concentrations for water quality montioring stations associated with outplanted reef sites. Data was downloaded from http://serc.fiu.edu/wqmnetwork/FKNMS-CD/index.htm in August 2020

9) outplant site coordinates.csv
  - coordinates for outplant reef sites used in this experiment
