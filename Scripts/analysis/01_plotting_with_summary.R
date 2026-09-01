plot_data <- combined_data %>%
  mutate(
    Day = as.numeric(Day),
    Rep_Num = as.character(Rep_Num),
    Vaccine = as.character(Vaccine)
  ) %>%
  filter(
    !is.na(Day),
    !is.na(Rep_Num), Rep_Num != "NA", Rep_Num != "NaN",
    !is.na(Vaccine), Vaccine != "NA", Vaccine != "NaN"
  ) %>%
  mutate(Rep_Num = factor(Rep_Num))

plot_combinations <- plot_data %>%
  distinct(Measurement_type, Cell_Type) %>%
  arrange(factor(Cell_Type, levels = c("CD4", "CD8")), Measurement_type)

# Function to plot and save raw data for a given summary statistic
save_raw_plots <- function(stat = c("median", "mean")) {
  stat <- match.arg(stat)
  stat_fun <- if (stat == "median") median else mean
  
  out_folder <- file.path(figure_folder, paste0("Raw_data_", stat))
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  
  for (i in seq_len(nrow(plot_combinations))) {
    meas <- plot_combinations$Measurement_type[i]
    cell_type <- plot_combinations$Cell_Type[i]
    
    sub_data <- plot_data %>%
      filter(Measurement_type == meas, Cell_Type == cell_type)
    
    if (nrow(sub_data) == 0 || sum(!is.na(sub_data$Value)) == 0) next
    
    unit_label <- unique(sub_data$Unit)[1]
    if (unit_label == "number") sub_data <- filter(sub_data, Value > 0)
    
    p <- ggplot(sub_data, aes(x = Day, y = Value, color = Vaccine, shape = Rep_Num)) +
      geom_boxplot(
        aes(group = interaction(Vaccine, Day), fill = Vaccine),
        alpha = 0.3, outlier.shape = NA,
        position = position_dodge(width = 0.8),
        show.legend = FALSE, na.rm = TRUE
      ) +
      geom_point(size = 2.5, alpha = 0.7, na.rm = TRUE) +
      stat_summary(
        aes(group = Vaccine, color = Vaccine),
        fun = stat_fun, geom = "line", size = 1, na.rm = TRUE
      ) +
      labs(
        title = paste0(meas, "  [", unit_label, "]"),
        x = "Day", y = "Value",
        color = "Vaccine", shape = "Replicate"
      ) +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(size = 11, face = "bold"),
            legend.position = "right") +
      scale_color_discrete(na.translate = FALSE) +
      scale_shape_discrete(na.translate = FALSE)
    
    if (unit_label == "number") p <- p + scale_y_log10()
    
    safe_meas <- gsub("[^A-Za-z0-9]", "_", meas)
    file_name <- paste0(cell_type, "_", safe_meas, ".png")
    ggsave(file.path(out_folder, file_name), p, width = 8, height = 6, dpi = 300)
  }
}

# Generate both versions
save_raw_plots("median")
save_raw_plots("mean")