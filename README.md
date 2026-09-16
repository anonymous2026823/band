# Supplementary Code for AISTATS 2027 Submission

This repository contains the code required to reproduce the experiments and results reported in our paper submitted to AISTATS 2027. This material is provided strictly for anonymous peer review.

## Setup Instructions

The R scripts utilize the `here` package to manage project-relative file paths. To ensure all required source files are located correctly, please open and activate the provided `band.Rproj` project in RStudio before running any scripts.

## Repository Contents

* **`loglikelihood.R`**: Reproduces the log-likelihood comparisons reported in **Section 5.5**.
* **`u_shaped_bias.R`**: Reproduces the U-shaped log-likelihood profile experiment reported in **Section 5.4**.
* **`real_data_experiments.R`**: Reproduces the real-data experiments reported in **Section 6**.

## Dependencies

This repository requires **R (version 4.1.0 or higher)**. 

To ensure reproducibility and easy evaluation, you can install all required dependencies by running the following snippet in your R console. This script will check your environment and install only the packages you are currently missing.

```R
# Define the list of required packages
required_packages <- c(
  "here",       # Project-relative paths
  "Matrix",     # Sparse matrix operations
  "igraph",     # Graph and network structures
  "ggplot2",    # Plotting and visualization
  "mclust",     # Gaussian mixture modeling baseline
  "Rcpp"        # C++ integration for computational speed
)

# Identify which packages are missing
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# Install the missing packages
if(length(new_packages)) {
  install.packages(new_packages, repos = "[http://cran.us.r-project.org](http://cran.us.r-project.org)")
} else {
  message("All dependencies are already installed!")
}
