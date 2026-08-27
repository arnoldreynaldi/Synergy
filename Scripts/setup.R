# Load packages
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(glue)
library(ggplot2)

file_path <- "C:/Users/z3443623/OneDrive - UNSW/Synergy"

onedrivefolder = glue("C:/Users/z3443623/OneDrive - UNSW/Synergy")

figure_folder = glue("{onedrivefolder}/output/plots/")

table_folder = glue("{onedrivefolder}/output/tables/")

# Folders to read from
data_folders <- c(
  glue("{onedrivefolder}/data/T cell data/CD8 T cell data"),
  glue("{onedrivefolder}/data/T cell data/CD4 T cell data")
)


# Get all xls/xlsx files from both folders
file_list <- list.files(
  path       = data_folders,
  pattern    = "\\.xlsx?$",
  full.names = TRUE
)

# Define the output folder
output_folder_Tcell_master <- glue("{onedrivefolder}/data/T cell data")

# Define the columns that are NOT measurements
non_measurement_cols <- c("Exp label", "Vaccine", "Sex", "Mouse #", "Exclude (Y/N)")