# Load packages
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(glue)


file_path <- "C:/Users/z3443623/OneDrive - UNSW/Synergy"
setwd(file_path)

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


# ===================================================================
# STEP 1: Collect all possible measurement column names from all files
#         by excluding known non‑measurement columns by name
# ===================================================================
# Define the columns that are NOT measurements
non_measurement_cols <- c("Exp label", "Vaccine", "Sex", "Mouse #", "Exclude (Y/N)")

measurement_cols <- map(file_list, function(f) {
  # Read only the header row of Sheet 3
  headers <- suppressMessages(read_excel(f, sheet = 3, n_max = 0) %>% colnames())
  
  # Keep everything except the non‑measurement columns
  setdiff(headers, non_measurement_cols)
}) %>% unlist() %>% unique()

if (length(measurement_cols) == 0) stop("No measurement columns found.")

# ===================================================================
# STEP 2: Process one file
# ===================================================================
process_file <- function(file_path, master_measurements) {
  
  # ---- 2a. Metadata from Sheet 1 (cells C3:C6) ----
  day_cell       <- read_excel(file_path, sheet = 1, range = "C4", col_names = FALSE)
  exp_date_cell  <- read_excel(file_path, sheet = 1, range = "C3", col_names = FALSE)
  rep_num_cell   <- read_excel(file_path, sheet = 1, range = "C5", col_names = FALSE)
  cell_type_cell <- read_excel(file_path, sheet = 1, range = "C6", col_names = FALSE)
  
  Day       <- day_cell[[1]]
  Exp_Date  <- exp_date_cell[[1]]
  if (inherits(Exp_Date, "POSIXct")) Exp_Date <- as.character(Exp_Date)
  Rep_Num   <- rep_num_cell[[1]]
  Cell_Type <- cell_type_cell[[1]]
  
  # ---- 2b. Mouse mapping table (B8:E31) ----
  mapping <- suppressMessages(
    read_excel(file_path, sheet = 1, range = "B8:E31", col_names = TRUE)
  )
  # Rename using actual column names (no positional indexing)
  mapping <- mapping %>%
    rename(Exp_label = `Exp label`,
           Mouse_number = `Mouse #`,
           Vaccine_mapping = Vaccine,
           Sex = Sex) %>%
    mutate(across(c(Exp_label, Mouse_number, Sex), as.character)) %>%
    filter(!is.na(Exp_label))
  
  # ---- 2c. Read Sheet 3 ----
  df <- suppressMessages(read_excel(file_path, sheet = 3))
  
  # Check that required columns exist by name
  if (!"Exp label" %in% names(df)) stop("Sheet 3 in ", basename(file_path), " does not have 'Exp label' column.")
  if (!"Vaccine" %in% names(df)) stop("Sheet 3 in ", basename(file_path), " does not have 'Vaccine' column.")
  
  # Rename first column to Exp_label; Vaccine already correct
  df <- df %>%
    rename(Exp_label = `Exp label`) %>%
    mutate(Exp_label = as.character(Exp_label))
  
  # ---- 2d. Ensure all master measurement columns exist ----
  existing_measurement_cols <- intersect(master_measurements, colnames(df))
  missing_measurement_cols  <- setdiff(master_measurements, colnames(df))
  
  df <- df %>% select(Exp_label, Vaccine, all_of(existing_measurement_cols))
  
  for (mcol in missing_measurement_cols) {
    df[[mcol]] <- NA_real_
  }
  
  # ---- 2e. Clean measurement columns ----
  df <- df %>%
    mutate(across(all_of(master_measurements),
                  ~ as.numeric(gsub("%", "", .)) %>% suppressWarnings()))
  
  # ---- 2f. Join with mouse mapping ----
  df <- df %>%
    left_join(mapping %>% select(Exp_label, Mouse_number, Sex), by = "Exp_label")
  
  # ---- 2g. Pivot to long format ----
  df_long <- df %>%
    pivot_longer(
      cols      = all_of(master_measurements),
      names_to  = "Measurement_type",
      values_to = "Value"
    )
  
  # ---- 2h. Add metadata and Unit ----
  df_long <- df_long %>%
    mutate(
      Day       = Day,
      Exp_Date  = Exp_Date,
      Rep_Num   = Rep_Num,
      Cell_Type = Cell_Type,
      Unit = if_else(grepl("%", Measurement_type), "percentage", "number")
    )
  
  # ---- 2i. Final column order ----
  df_long <- df_long %>%
    select(Exp_label, Mouse_number, Vaccine, Sex,
           Day, Exp_Date, Rep_Num, Cell_Type,
           Unit, Measurement_type, Value)
  
  return(df_long)
}

# ===================================================================
# STEP 3: Process all files and combine
# ===================================================================
combined_data <- map_dfr(file_list, process_file, master_measurements = measurement_cols)

head(combined_data)

# Write the CSV
write.csv(combined_data, file.path(output_folder_Tcell_master, "Data_Master_T_Cell.csv"), row.names = FALSE)