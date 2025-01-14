library(TMB)
library(mvtnorm)

#source("simRR.R")   # Simulate run reconstruction
source(here("analysis/R/run-reconstructions/fitRR.R"))   # Fit run reconstruction
#source("fitSim.R")  # Fit run reconstruction to simmed data
#source("calcRunTiming.R")
source(here("analysis/R/run-reconstructions/tools.R"))   # Background functions
source(here("analysis/R/run-reconstructions/plot.R"))    # Plotting functions

# Compile TMB objective function
compile(here("analysis/R/run-reconstructions/yukonChinookRunRecon.cpp"))
