measurement_cols <- map(file_list, function(f) {
  # Read only the header row of Sheet 3
  headers <- suppressMessages(read_excel(f, sheet = 3, n_max = 0) %>% colnames())
  
  # Keep everything except the non‑measurement columns
  setdiff(headers, non_measurement_cols)
}) %>% unlist() %>% unique()

if (length(measurement_cols) == 0) stop("No measurement columns found.")

combined_data <- map_dfr(file_list, process_file, master_measurements = measurement_cols)

# Rename Vaccine into four categories
combined_data <- combined_data %>%
  mutate(
    Vaccine = case_when(
      grepl("TMd21", Vaccine) ~ "TMd21",
      grepl("SOL",   Vaccine) ~ "SOL",
      grepl("IC",    Vaccine) ~ "IC",
      grepl("TM",    Vaccine) ~ "TM"
    )
  )

# Write the CSV
write.csv(combined_data, file.path(output_folder_Tcell_master, "Data_Master_T_Cell.csv"), row.names = FALSE)





# 
# 
# library(ggplot2)
# library(dplyr)
# 
# # -------------------------------------------------------------------
# # 1. Clean data: rename Vaccine, converting to NA, convert types
# # -------------------------------------------------------------------
# 
# combined_data <- combined_data %>%
#   mutate(
#     Vaccine = case_when(
#       grepl("TMd21", Vaccine) ~ "TMd21",
#       grepl("SOL",   Vaccine) ~ "SOL",
#       grepl("IC",    Vaccine) ~ "IC",
#       grepl("TM",    Vaccine) ~ "TM",
#       TRUE ~ Vaccine
#     )
#   )
# 
# plot_data <- combined_data %>%
#   mutate(
#     Day     = as.numeric(Day),
#     Rep_Num = as.character(Rep_Num),
#     Vaccine = as.character(Vaccine)
#   ) %>%
#   filter(
#     !is.na(Day),
#     !is.na(Rep_Num), Rep_Num != "NA", Rep_Num != "NaN",
#     !is.na(Vaccine), Vaccine != "NA", Vaccine != "NaN"
#   ) %>%
#   mutate(Rep_Num = factor(Rep_Num))
# 
# # -------------------------------------------------------------------
# # 2. Unique combinations of Measurement_type + Cell_Type,
# #    ordered: all CD4 before CD8
# # -------------------------------------------------------------------
# plot_combinations <- plot_data %>%
#   distinct(Measurement_type, Cell_Type) %>%
#   arrange(factor(Cell_Type, levels = c("CD4", "CD8")), Measurement_type)
# 
# # -------------------------------------------------------------------
# # 3. Save one PNG per combination, skipping empty ones
# # -------------------------------------------------------------------
# output_folder <- "measurement_plots"
# dir.create(output_folder, showWarnings = FALSE)
# 
# for (i in seq_len(nrow(plot_combinations))) {
#   
#   meas      <- plot_combinations$Measurement_type[i]
#   cell_type <- plot_combinations$Cell_Type[i]
#   
#   # Filter data for this specific measurement AND cell type
#   sub_data <- plot_data %>%
#     filter(Measurement_type == meas, Cell_Type == cell_type)
#   
#   # Skip if no rows or if all Value are NA
#   if (nrow(sub_data) == 0 || sum(!is.na(sub_data$Value)) == 0) {
#     next
#   }
#   
#   # Unit (percentage or number)
#   unit_label <- unique(sub_data$Unit)[1]
#   
#   # Build plot (title without cell type, as you requested)
#   p <- ggplot(sub_data, aes(x = Day, y = Value, color = Vaccine, shape = Rep_Num)) +
#     geom_point(size = 2.5, alpha = 0.7, na.rm = TRUE) +
#     stat_summary(
#       aes(group = Vaccine, color = Vaccine),
#       fun = mean,
#       geom = "line",
#       size = 1,
#       na.rm = TRUE
#     ) +
#     labs(
#       title = paste0(meas, "  [", unit_label, "]"),
#       x = "Day",
#       y = "Value",
#       color = "Vaccine",
#       shape = "Replicate"
#     ) +
#     theme_bw(base_size = 12) +
#     theme(
#       plot.title = element_text(size = 11, face = "bold"),
#       legend.position = "right"
#     ) +
#     scale_color_discrete(na.translate = FALSE) +
#     scale_shape_discrete(na.translate = FALSE)
#   
#   # Create safe filename: "<CellType>_<measurement>.png"
#   safe_meas <- gsub("[^A-Za-z0-9]", "_", meas)
#   file_name <- paste0(cell_type, "_", safe_meas, ".png")
#   file_path <- file.path(output_folder, file_name)
#   
#   ggsave(
#     filename = file_path,
#     plot = p,
#     width = 8,
#     height = 6,
#     dpi = 300
#   )
# }