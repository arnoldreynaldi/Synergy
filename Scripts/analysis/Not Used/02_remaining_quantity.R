# -------------------------------------------------------------------
# Compute medians for all combinations 
# -------------------------------------------------------------------
median_list <- list()

for (i in seq_len(nrow(plot_combinations))) {
  
  meas      <- plot_combinations$Measurement_type[i]
  cell_type <- plot_combinations$Cell_Type[i]
  
  sub_data <- plot_data %>%
    filter(Measurement_type == meas, Cell_Type == cell_type)
  
  if (nrow(sub_data) == 0 || sum(!is.na(sub_data$Value)) == 0) next
  
  unit_label <- unique(sub_data$Unit)[1]
  
  # For count data, keep only positive values so log scale works
  if (unit_label == "number") {
    sub_data <- sub_data %>% filter(Value > 0)
  }
  
  # Compute median per Vaccine and Day
  med_summary <- sub_data %>%
    group_by(Vaccine, Day) %>%
    summarise(Median = median(Value, na.rm = TRUE), .groups = "drop")
  
  # Save this summary 
  med_summary <- med_summary %>%
    mutate(Measurement_type = meas, Cell_Type = cell_type, Unit = unit_label)
  
  median_list[[length(median_list) + 1]] <- med_summary
}

# Combine all median dataframes into one
all_medians <- bind_rows(median_list)

# -------------------------------------------------------------------
# Normalize medians by baseline (minimum Day) within each Vaccine/Measurement/Cell_Type
# -------------------------------------------------------------------
normalised_medians <- all_medians %>%
  group_by(Measurement_type, Cell_Type, Vaccine) %>%
  mutate(
    baseline_Day = min(Day),
    baseline_Median = Median[Day == baseline_Day][1],
    normalised = Median / baseline_Median
  ) %>%
  ungroup() %>%
  filter(!is.na(normalised), is.finite(normalised))

# -------------------------------------------------------------------
# Plot and save normalised median lines
# -------------------------------------------------------------------
norm_plot_folder <- file.path(figure_folder, "Normalised_median")
if (!dir.exists(norm_plot_folder)) {
  dir.create(norm_plot_folder, recursive = TRUE)
}

norm_table_folder <- file.path(table_folder, "Normalised_median")
if (!dir.exists(norm_table_folder)) {
  dir.create(norm_table_folder, recursive = TRUE)
}

for (i in seq_len(nrow(plot_combinations))) {
  
  meas      <- plot_combinations$Measurement_type[i]
  cell_type <- plot_combinations$Cell_Type[i]
  
  sub_norm <- normalised_medians %>%
    filter(Measurement_type == meas, Cell_Type == cell_type)
  
  if (nrow(sub_norm) == 0) next
  
  unit_label <- unique(sub_norm$Unit)[1]
  
  p <- ggplot(sub_norm, aes(x = Day, y = normalised, color = Vaccine)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = paste0(meas, "  [", unit_label, "]  (normalised)"),
      x = "Day",
      y = "normalised median (fold change)",
      color = "Vaccine"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "right"
    ) +
    scale_y_log10()   
  
  safe_meas <- gsub("[^A-Za-z0-9]", "_", meas)
  file_name <- paste0(cell_type, "_", safe_meas, "_normalised_median.png")
  file_path <- file.path(norm_plot_folder, file_name)
  
  ggsave(
    filename = file_path,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# save normalised CSV in table folder
write.csv(normalised_medians, file.path(norm_table_folder, "normalised_medians.csv"), row.names = FALSE)