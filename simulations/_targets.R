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
  
  #tar_target(
  #  name = surv_cox_cens_cox,
  #  command = generate_survival_data(n = n_sim, x_vals = x_vals,surv_is_cox = T, cens_is_cox = T)
  #  # format = "feather" # efficient storage for large data frames
  #),
  #tar_target(
  #  name = surv_cox_cens_not_cox,
  #  command = generate_survival_data(n = n_sim, x_vals = x_vals,surv_is_cox = T, cens_is_cox = F)
  #  # format = "feather" # efficient storage for large data frames
  #),
  #tar_target(
  #  name = surv_not_cox_cens_cox,
  #  command = generate_survival_data(n = n_sim, x_vals = x_vals,surv_is_cox = F, cens_is_cox = T)
  #  # format = "feather" # efficient storage for large data frames
  #),
  #tar_target(
  #  name = surv_not_cox_cens_not_cox,
  #  command = generate_survival_data(n = n_sim, x_vals = x_vals,surv_is_cox = F, cens_is_cox = F)
  #  # format = "feather" # efficient storage for large data frames
  #)
  tar_target(name = test_of_data_generating, command = replicate(n = n_reps, 
                                                                 check_mod(n = n_sim, 
                                                                           surv_is_cox = F, 
                                                                           treatment_effect = 0, 
                                                                           x_vals = c(0,10))
  )
  )
)

