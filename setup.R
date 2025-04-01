# setup.R - Script to set up the R environment
# This script installs the renv package if not available, initializes a new
# environment if no lockfile is present, and installs the packages listed in
# the script if a lockfile is not present.

# Install renv if not available
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize or restore the environment
if (file.exists("renv.lock")) {
  renv::restore()  # Restore packages from lockfile
} else {
  renv::init()  # Initialize a new environment
}

# Install required packages
if (!requireNamespace("pak", quietly = TRUE)) {
  renv::install("pak")
}

pak::pak(c("reprex", "Rcpp"))
pak::pak(c("tidyverse", "ggplot2", "igraph", "tidygraph", "ggraph", "remotes"))
# remotes::install_github("PixelgenTechnologies/pixelatorR@v0.12.0")
pak::pak("PixelgenTechnologies/pixelatorR@v0.12.0")

# Save the environment state
renv::snapshot()
