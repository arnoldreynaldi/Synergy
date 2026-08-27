file_path <- "C:/Users/z3443623/OneDrive - UNSW/Synergy/Scripts"
setwd(file_path)

# First thing is to run the setup file
source("setup.R")

# Source all the helper files
source("./R/process_file.R")

# Source all the processing files
source("./processing/01_combine_files.R")

# Source all the plotting and analysis files
source("./analysis/01_plotting.R")
