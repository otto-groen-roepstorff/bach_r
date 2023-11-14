# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble") # packages that your targets need to run
)


# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/functions.R")

# Replace the target list below with your own:
list(
  tar_target(
    name = data,
    command = generate_survival_data(n = 100,x_vals = (0:5),is_cox = T)
    # format = "feather" # efficient storage for large data frames
  ),
  tar_target(
    name = model,
    command = coefficients(lm(y ~ x, data = data))
  )
)

