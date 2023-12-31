# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
#library(tidyverse)
#library(survival)
#library(timereg)
#library(dplyr)
#library(data.table)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble", "targets","survival","timereg","dplyr","data.table") # packages that your targets need to run
)


# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/functions.R")

# Replace the target list below with your own:
list(
  tar_target(
    name = n_sim,
    command = 10000
  ),
  tar_target(
    name = x_vals,
    command = (0:10)
  ),
  tar_target(
    name = n_reps,
    command = 1000
  ),
  
  tar_target(
    name = data,
    command = generate_survival_data(n = n_sim, x_vals = c(-1,1)),
     format = "feather" # efficient storage for large data frames0
  ),
  
  #model diagnostics
  tar_target(
    name = dat_observed,
    command = proportion_observed(data)
  ),
  
  tar_target(
    name = data_plots,
    command = visualize_data(data)
  ),
  
  #building models
  tar_target(
    name = oracle_T_model,
    command = oracle_model(data)
    # format = "feather" # efficient storage for large data frames
  ),
  tar_target(
    name = non_oracle_T_model,
    command = non_oracle_model(data)
  ),
    
  tar_target(
    name = non_oracle_T_and_C_model,
    command = non_oracle_model(data)
  ),
  
  
  
    
)

