plot_data <- combined_data %>%
  mutate(
    Day     = as.numeric(Day),
    Rep_Num = as.character(Rep_Num),
    Vaccine = as.character(Vaccine)
  ) %>%
  filter(
    !is.na(Day),
    !is.na(Rep_Num), Rep_Num != "NA", Rep_Num != "NaN",
    !is.na(Vaccine), Vaccine != "NA", Vaccine != "NaN"
  ) %>%
  mutate(Rep_Num = factor(Rep_Num))

# -------------------------------------------------------------------
# Unique combinations of Measurement_type + Cell_Type,
# ordered: all CD4 before CD8
# -------------------------------------------------------------------
plot_combinations <- plot_data %>%
  distinct(Measurement_type, Cell_Type) %>%
  arrange(factor(Cell_Type, levels = c("CD4", "CD8")), Measurement_type)

# -------------------------------------------------------------------
# Save one PNG per combination, skipping empty ones
# -------------------------------------------------------------------
for (i in seq_len(nrow(plot_combinations))) {
  
  meas      <- plot_combinations$Measurement_type[i]
  cell_type <- plot_combinations$Cell_Type[i]
  
  # Filter data for this specific measurement AND cell type
  sub_data <- plot_data %>%
    filter(Measurement_type == meas, Cell_Type == cell_type)
  
  # Skip if all Value are NA
  if (nrow(sub_data) == 0 || sum(!is.na(sub_data$Value)) == 0) {
    next
  }
  
  # Unit (percentage or number)
  unit_label <- unique(sub_data$Unit)[1]
  
  p <- ggplot(sub_data, aes(x = Day, y = Value, color = Vaccine, shape = Rep_Num)) +
    geom_point(size = 2.5, alpha = 0.7, na.rm = TRUE) +
    stat_summary(
      aes(group = Vaccine, color = Vaccine),
      fun = mean,
      geom = "line",
      size = 1,
      na.rm = TRUE
    ) +
    labs(
      title = paste0(meas, "  [", unit_label, "]"),
      x = "Day",
      y = "Value",
      color = "Vaccine",
      shape = "Replicate"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "right"
    ) +
    scale_color_discrete(na.translate = FALSE) +
    scale_shape_discrete(na.translate = FALSE)
  
  # Create filename
  safe_meas <- gsub("[^A-Za-z0-9]", "_", meas)
  file_name <- paste0(cell_type, "_", safe_meas, ".png")
  file_path <- file.path(figure_folder, file_name)
  
  ggsave(
    filename = file_path,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}