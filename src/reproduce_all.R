#################################
# REPRODUCE NUMERICAL SIMULATIONS
#################################

# clear workspace
rm(list = ls())

# source files containing functions
source("CRT.R")
source("utils_CRT.R")
source("utils_Gaussian.R")
source("run_one_experiment.R")

# load libraries
library(glmnet)
library(reshape2)
library(MASS)
library(knockoff)
library(R.utils)
library(tidyverse)

# run numerical simulations
run_one_experiment("main_simulation")
run_one_experiment("knockoffs_variability")

# plot results
source("plot_results.R")