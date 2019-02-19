


##################################################################################
#### MASTER FILE FOR ANALYSIS Primary outcome efficacy
##################################################################################

#load needed libraries
library(meta)
library(metafor)
library(netmeta)
library(readxl)
library(devtools)
install_github("esm-ispm-unibe-ch/NMAJags")
library(NMAJags)
library(R2jags)

#get the data and 
Schizo <- read_excel("Schizo.xlsx",na = "NA")

#exclude the open label studies
outopen=(Schizo$Blinding_type<2)|is.na(Schizo$Blinding_type)
SchizoOBJ=Schizo[outopen==0,]
write.table(SchizoOBJ,"SchizoWITHOUTopenlabel.csv", sep="\t")

#Analysis of the primary outcome
source("Schizo primary outcome analysis.R")
SchizoEFFall=SchizoEFF

#Meta-regressions for placebo response
#The script excludes drugs with less than 100 patients randomised from the SchizoEFF database
source("modelNMRPlacChange.R")#load the JAGS model
source("modelNMRPlacChangeExchB.R")#load the JAGS model
source("Placebo response on efficacy meta-regression.R")


#Meta-regressions and sensitivity analyses
#The script excludes drugs with less than 100 patients randomised from the SchizoEFF database
source("modelNMRContinuousCommonB.R")#load the JAGS model
source("Placebo response on efficacy meta-regression.R")



rm(list=ls())
