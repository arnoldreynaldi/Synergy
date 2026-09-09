# =====================================================
# Master script: Decay modelling for all CD4/CD8 number measurements
# - White background plots
# - Remove Mouse_number 11171 for CD4
# - Save plots and tables in separate directories per cell type and measurement type
# - Consistent vaccine colour scheme
# =====================================================

library(dplyr)
library(ggplot2)

# Set global theme to white background
theme_set(theme_bw())

# Base directories
figure_folder <- "C:/Projects/Synergy/output/plots/"
table_folder  <- "C:/Projects/Synergy/output/tables/"

# Vaccine order and colours
vaccine_levels <- c("TM", "TMd21", "SOL", "IC")
formulation_cols <- c(
  TM    = "#2F75B5",
  TMd21 = "#D99A2B",
  SOL   = "#4A9E7D",
  IC    = "#B86691"
)

# Function to sanitize folder names
sanitize <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

# Main analysis function
run_analysis <- function(data, cell_type, measurement_type) {
  
  cat("\n====================\n")
  cat("Running analysis for:", cell_type, "-", measurement_type, "\n")
  
  # Filter data
  subset_df <- data %>%
    filter(Cell_Type == cell_type,
           Measurement_type == measurement_type,
           Unit == "number") %>%
    mutate(log10_Value = log10(Value))
  
  # Remove rows with non-finite log values
  subset_df <- subset_df %>% filter(is.finite(log10_Value))
  
  # Remove Mouse_number 11171 for CD4
  if (cell_type == "CD4") {
    subset_df <- subset_df %>% filter(Mouse_number != 11171)
  }
  
  # Check minimum observations per vaccine
  n_per_vaccine <- subset_df %>% group_by(Vaccine) %>% summarise(n = n())
  if (any(n_per_vaccine$n < 3)) {
    cat("Skipping: insufficient data (some vaccine has less than 3 observations).\n")
    return(NULL)
  }
  
  # Ensure vaccine factor levels match the desired order
  subset_df$Vaccine <- factor(subset_df$Vaccine, levels = vaccine_levels)
  vaccines <- levels(droplevels(subset_df$Vaccine))
  if (length(vaccines) < 2) {
    cat("Skipping: fewer than 2 vaccine groups present.\n")
    return(NULL)
  }
  
  # ---------------------------
  # Simple linear model (per vaccine)
  # ---------------------------
  lm_fits <- subset_df %>%
    group_by(Vaccine) %>%
    do(mod = lm(log10_Value ~ Day, data = .))
  
  lm_params <- bind_rows(lapply(vaccines, function(v) {
    mod <- lm_fits$mod[[which(lm_fits$Vaccine == v)]]
    ci <- confint(mod)
    data.frame(Vaccine = v,
               Parameter = c("Intercept", "Slope"),
               Estimate = coef(mod),
               Lower = ci[,1], Upper = ci[,2])
  }))
  
  lm_fitted <- bind_rows(lapply(vaccines, function(v) {
    dat <- subset_df %>% filter(Vaccine == v)
    day_seq <- seq(min(dat$Day), max(dat$Day), length.out = 50)
    pred <- predict(lm_fits$mod[[which(lm_fits$Vaccine == v)]],
                    newdata = data.frame(Day = day_seq))
    data.frame(Vaccine = v, Day = day_seq, Fitted = pred)
  }))
  
  # ---------------------------
  # Best breakpoint per vaccine (simple piecewise)
  # ---------------------------
  bp_candidates <- c(14, 28)
  best_bp <- sapply(vaccines, function(v) {
    dat <- subset_df %>% filter(Vaccine == v)
    aics <- sapply(bp_candidates, function(bp) {
      mod <- lm(log10_Value ~ Day + pmax(0, Day - bp), data = dat)
      AIC(mod)
    })
    bp_candidates[which.min(aics)]
  })
  
  subset_df <- subset_df %>%
    mutate(Day_plus = pmax(0, Day - best_bp[Vaccine]))
  
  # ---------------------------
  # Piecewise model (per vaccine)
  # ---------------------------
  piecewise_models <- list()
  piecewise_params <- bind_rows(lapply(vaccines, function(v) {
    dat <- subset_df %>% filter(Vaccine == v)
    mod <- lm(log10_Value ~ Day + Day_plus, data = dat)
    piecewise_models[[v]] <<- mod
    coefs <- coef(mod); vc <- vcov(mod); ci <- confint(mod)
    slope_after <- coefs[2] + coefs[3]
    se_after <- sqrt(vc[2,2] + vc[3,3] + 2*vc[2,3])
    data.frame(Vaccine = v,
               Parameter = c("Intercept", "Slope_before", "Slope_after"),
               Estimate = c(coefs[1], coefs[2], slope_after),
               Lower = c(ci[1,1], ci[2,1], slope_after - 1.96*se_after),
               Upper = c(ci[1,2], ci[2,2], slope_after + 1.96*se_after))
  }))
  
  piecewise_fitted <- bind_rows(lapply(vaccines, function(v) {
    dat <- subset_df %>% filter(Vaccine == v)
    day_seq <- seq(min(dat$Day), max(dat$Day), length.out = 50)
    newdata <- data.frame(Day = day_seq, Day_plus = pmax(0, day_seq - best_bp[v]))
    pred <- predict(piecewise_models[[v]], newdata = newdata)
    data.frame(Vaccine = v, Day = day_seq, Fitted = pred)
  }))
  
  # ---------------------------
  # Coefficient plots (coloured by vaccine)
  # ---------------------------
  plot_coef <- function(data, param, title) {
    data %>% filter(Parameter == param) %>%
      ggplot(aes(x = Vaccine, y = Estimate, color = Vaccine)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
      scale_color_manual(values = formulation_cols, drop = FALSE) +
      labs(title = title, y = "Estimate") +
      theme_bw() +
      theme(legend.position = "none")
  }
  
  p1 <- plot_coef(lm_params, "Intercept", "Intercept (Linear)")
  p2 <- plot_coef(lm_params, "Slope", "Slope (Linear)")
  p3 <- plot_coef(piecewise_params, "Intercept", "Intercept (Piecewise)")
  p4 <- plot_coef(piecewise_params, "Slope_before", "Slope before breakpoint")
  p5 <- plot_coef(piecewise_params, "Slope_after", "Slope after breakpoint")
  
  # Model vs data plots with vaccine colour scheme
  plot_lm_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
    geom_point(alpha = 0.6) +
    geom_line(data = lm_fitted, aes(y = Fitted), size = 1) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Linear model fit (log10 scale)", y = "log10(Value)", x = "Day") +
    facet_wrap(~Vaccine) +
    theme_bw()
  
  plot_piecewise_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
    geom_point(alpha = 0.6) +
    geom_line(data = piecewise_fitted, aes(y = Fitted), size = 1) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Best piecewise model fit (log10 scale)", y = "log10(Value)", x = "Day") +
    facet_wrap(~Vaccine) +
    theme_bw()
  
  # ---------------------------
  # Slope grouping analysis
  # ---------------------------
  generate_partitions <- function(n) {
    if (n == 1) return(list(1L))
    result <- list()
    prev <- generate_partitions(n - 1)
    for (p in prev) {
      max_cluster <- max(p)
      for (c in 1:max_cluster) result <- c(result, list(c(p, c)))
      result <- c(result, list(c(p, max_cluster + 1L)))
    }
    result
  }
  
  n_vaccines <- length(vaccines)
  partitions <- generate_partitions(n_vaccines)
  
  fit_lm_grouping <- function(cl) {
    slope_cluster <- factor(cl[match(subset_df$Vaccine, vaccines)])
    if (length(unique(cl)) == 1) {
      mod <- lm(log10_Value ~ Vaccine + Day, data = subset_df)
    } else {
      mod <- lm(log10_Value ~ Vaccine + Day:slope_cluster, data = subset_df)
    }
    list(aic = AIC(mod), cluster = cl)
  }
  lm_group_results <- lapply(partitions, fit_lm_grouping)
  
  fit_piecewise_grouping <- function(cl_before, cl_after) {
    before_cluster <- factor(cl_before[match(subset_df$Vaccine, vaccines)])
    after_cluster  <- factor(cl_after[match(subset_df$Vaccine, vaccines)])
    formula_str <- "log10_Value ~ Vaccine"
    if (length(unique(cl_before)) == 1) formula_str <- paste(formula_str, "+ Day")
    else formula_str <- paste(formula_str, "+ Day:before_cluster")
    if (length(unique(cl_after)) == 1) formula_str <- paste(formula_str, "+ Day_plus")
    else formula_str <- paste(formula_str, "+ Day_plus:after_cluster")
    mod <- lm(as.formula(formula_str), data = subset_df)
    list(aic = AIC(mod), cl_before = cl_before, cl_after = cl_after)
  }
  
  group_combos <- expand.grid(idx_before = 1:length(partitions),
                              idx_after  = 1:length(partitions))
  piecewise_group_results <- vector("list", nrow(group_combos))
  for (i in 1:nrow(group_combos)) {
    cl_before <- partitions[[group_combos$idx_before[i]]]
    cl_after  <- partitions[[group_combos$idx_after[i]]]
    piecewise_group_results[[i]] <- fit_piecewise_grouping(cl_before, cl_after)
  }
  
  # Readable AIC table
  format_group <- function(cl) {
    groups <- split(vaccines, cl)
    paste(sapply(groups, paste, collapse = " & "), collapse = " | ")
  }
  
  linear_table <- data.frame(
    Model = "Linear",
    Slope_group = sapply(lm_group_results, function(x) format_group(x$cluster)),
    k = n_vaccines + sapply(lm_group_results, function(x) length(unique(x$cluster))),
    AIC = sapply(lm_group_results, function(x) x$aic)
  )
  
  piecewise_table <- data.frame(
    Model = "Piecewise",
    Slope_group = paste0("Early: ", sapply(piecewise_group_results, function(x) format_group(x$cl_before)),
                         "  |  Late: ", sapply(piecewise_group_results, function(x) format_group(x$cl_after))),
    k = n_vaccines + sapply(piecewise_group_results, function(x) length(unique(x$cl_before))) +
      sapply(piecewise_group_results, function(x) length(unique(x$cl_after))),
    AIC = sapply(piecewise_group_results, function(x) x$aic)
  )
  
  full_table <- rbind(linear_table, piecewise_table)
  full_table$Delta_AIC <- full_table$AIC - min(full_table$AIC)
  full_table <- full_table[order(full_table$AIC), ]
  
  # ---------------------------
  # Save outputs
  # ---------------------------
  plot_base <- file.path(figure_folder, paste0(cell_type, "_analysis"), sanitize(measurement_type))
  table_base <- file.path(table_folder, paste0(cell_type, "_analysis"), sanitize(measurement_type))
  dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_base, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(plot_base, "1_Intercept_Linear.png"), p1, width = 6, height = 4)
  ggsave(file.path(plot_base, "2_Slope_Linear.png"), p2, width = 6, height = 4)
  ggsave(file.path(plot_base, "3_Intercept_Piecewise.png"), p3, width = 6, height = 4)
  ggsave(file.path(plot_base, "4_Slope_before_Piecewise.png"), p4, width = 6, height = 4)
  ggsave(file.path(plot_base, "5_Slope_after_Piecewise.png"), p5, width = 6, height = 4)
  ggsave(file.path(plot_base, "6_LM_fit.png"), plot_lm_fit, width = 8, height = 5)
  ggsave(file.path(plot_base, "7_Piecewise_fit.png"), plot_piecewise_fit, width = 8, height = 5)
  
  write.csv(lm_params, file.path(table_base, "linear_params.csv"), row.names = FALSE)
  write.csv(piecewise_params, file.path(table_base, "piecewise_params.csv"), row.names = FALSE)
  write.csv(full_table, file.path(table_base, "AIC_table.csv"), row.names = FALSE)
  
  cat("Saved outputs to:\n  Plots: ", plot_base, "\n  Tables: ", table_base, "\n")
  return(full_table)
}

# =====================================================
# Main loop
# =====================================================
combos <- combined_data %>%
  filter(Unit == "number", Cell_Type %in% c("CD4", "CD8")) %>%
  distinct(Cell_Type, Measurement_type) %>%
  arrange(Cell_Type, Measurement_type)

cat("Total analyses to run:", nrow(combos), "\n")

for (i in 1:nrow(combos)) {
  ct <- combos$Cell_Type[i]
  mt <- combos$Measurement_type[i]
  tryCatch(
    run_analysis(combined_data, ct, mt),
    error = function(e) {
      cat("ERROR for", ct, "-", mt, ":", e$message, "\n")
    }
  )
}

cat("\nAll done.\n")