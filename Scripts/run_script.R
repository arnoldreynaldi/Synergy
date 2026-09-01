file_path <- "C:/Projects/Synergy/Scripts"
setwd(file_path)

# First thing is to run the setup file
source("setup.R")

# Source all the helper files
source("./R/process_file.R")

# Source all the processing files
source("./processing/01_combine_files.R")

# Source all the plotting and analysis files
source("./analysis/01_plotting_with_summary.R")
source("./analysis/02_normalised_remaining_quantity.R")