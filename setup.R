# setup.R - Script to set up the R environment
# This script installs the renv package if not available, initializes a new
# environment if no lockfile is present, and installs the packages listed
# here or collected by `renv::dependencies()`

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

renv::install("PixelgenTechnologies/pixelatorR@pna-636-migrate-code-base")

# Save the environment state
renv::snapshot()
