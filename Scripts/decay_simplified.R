sanitize <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

safe_predict <- function(model, newdata) {
  if (is.null(model)) return(rep(NA, nrow(newdata)))
  tryCatch(predict(model, newdata), error = function(e) rep(NA, nrow(newdata)))
}

wrap_title <- function(txt, width = 22) {
  paste(strwrap(txt, width = width), collapse = "\n")
}

theme_pub <- function(base_size = 18) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle    = element_text(size = base_size - 3),
      axis.title       = element_text(size = base_size, face = "bold"),
      axis.text        = element_text(size = base_size - 2, color = "black"),
      strip.text       = element_text(size = base_size - 1, face = "bold"),
      legend.title     = element_text(size = base_size, face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      legend.position  = "none",
      plot.margin      = margin(6, 6, 6, 6)
    )
}

run_analysis <- function(data, cell_type, measurement_type, forced_bp = NULL) {
  subset_df <- data %>%
    filter(Cell_Type == cell_type,
           Measurement_type == measurement_type,
           Unit == "number") %>%
    mutate(log10_Value = log10(Value))
  
  subset_df <- subset_df %>% filter(is.finite(log10_Value))
  subset_df <- subset_df %>% filter(!Mouse_number %in% subject_to_exclude)
  
  n_per_vaccine <- subset_df %>% group_by(Vaccine) %>% summarise(n = n())
  if (any(n_per_vaccine$n < 1)) return(NULL)
  
  subset_df$Vaccine <- factor(subset_df$Vaccine, levels = vaccine_levels)
  vaccines <- levels(droplevels(subset_df$Vaccine))
  if (length(vaccines) < 1) return(NULL)
  
  # Per-vaccine linear model
  lm_fits <- subset_df %>%
    group_by(Vaccine) %>%
    do(mod = lm(log10_Value ~ Day, data = .))
  
  lm_params <- bind_rows(lapply(vaccines, function(v) {
    mod <- lm_fits$mod[[which(lm_fits$Vaccine == v)]]
    if (is.null(mod)) return(NULL)
    ci <- confint(mod)
    data.frame(Vaccine = v,
               Parameter = c("Intercept", "Slope"),
               Estimate = coef(mod),
               Lower = ci[,1], Upper = ci[,2])
  }))
  
  lm_fitted <- bind_rows(lapply(vaccines, function(v) {
    dat <- subset_df %>% filter(Vaccine == v)
    day_seq <- seq(min(dat$Day), max(dat$Day), length.out = 50)
    pred <- safe_predict(lm_fits$mod[[which(lm_fits$Vaccine == v)]],
                         data.frame(Day = day_seq))
    data.frame(Vaccine = v, Day = day_seq, Fitted = pred)
  }))
  
  # ---- Breakpoint candidate selection ----
  if (is.null(forced_bp)) {
    bp_candidates <- c(14, 28)
  } else if (identical(forced_bp, "second_timepoint")) {
    all_days <- sort(unique(subset_df$Day))
    bp_candidates <- if (length(all_days) >= 2) all_days[2] else all_days[1]
    message("Forcing breakpoint at SECOND timepoint = Day ", bp_candidates,
            " for Cell_Type = ", cell_type, ", measurement_type = ", measurement_type)
  } else if (identical(forced_bp, "third_timepoint")) {
    all_days <- sort(unique(subset_df$Day))
    bp_candidates <- if (length(all_days) >= 3) all_days[3] else all_days[length(all_days)]
    message("Forcing breakpoint at THIRD timepoint = Day ", bp_candidates,
            " for Cell_Type = ", cell_type, ", measurement_type = ", measurement_type)
  } else if (is.numeric(forced_bp)) {
    bp_candidates <- forced_bp
    message("Forcing breakpoint at Day ", forced_bp,
            " for Cell_Type = ", cell_type, ", measurement_type = ", measurement_type)
  } else {
    stop("`forced_bp` must be NULL, numeric, \"second_timepoint\", or \"third_timepoint\"")
  }
  
  best_bp_per_vaccine <- sapply(vaccines, function(v) {
    dat <- subset_df %>% filter(Vaccine == v)
    aics <- sapply(bp_candidates, function(bp) {
      mod <- tryCatch(lm(log10_Value ~ Day + pmax(0, Day - bp), data = dat),
                      error = function(e) NULL)
      if (is.null(mod)) return(Inf)
      AIC(mod)
    })
    if (all(is.infinite(aics))) return(bp_candidates[1])
    bp_candidates[which.min(aics)]
  })
  
  subset_df <- subset_df %>%
    mutate(Day_plus = pmax(0, Day - best_bp_per_vaccine[Vaccine]))
  
  piecewise_models <- list()
  piecewise_params <- bind_rows(lapply(vaccines, function(v) {
    dat <- subset_df %>% filter(Vaccine == v)
    mod <- tryCatch(lm(log10_Value ~ Day + Day_plus, data = dat),
                    error = function(e) NULL)
    piecewise_models[[v]] <<- mod
    if (is.null(mod)) return(NULL)
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
    if (is.null(piecewise_models[[v]])) return(NULL)
    dat <- subset_df %>% filter(Vaccine == v)
    day_seq <- seq(min(dat$Day), max(dat$Day), length.out = 50)
    newdata <- data.frame(Day = day_seq,
                          Day_plus = pmax(0, day_seq - best_bp_per_vaccine[v]))
    pred <- safe_predict(piecewise_models[[v]], newdata)
    data.frame(Vaccine = v, Day = day_seq, Fitted = pred)
  }))
  
  plot_coef <- function(data, param, title) {
    data %>% filter(Parameter == param) %>%
      ggplot(aes(x = Vaccine, y = Estimate, color = Vaccine)) +
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, linewidth = 1) +
      scale_color_manual(values = formulation_cols, drop = FALSE) +
      labs(title = wrap_title(title), y = "Estimate") +
      theme_pub()
  }
  
  p1 <- plot_coef(lm_params, "Intercept", "Intercept (Linear)")
  p2 <- plot_coef(lm_params, "Slope", "Slope (Linear)")
  p3 <- plot_coef(piecewise_params, "Intercept", "Intercept (Piecewise)")
  p4 <- plot_coef(piecewise_params, "Slope_before", "Slope before breakpoint")
  p5 <- plot_coef(piecewise_params, "Slope_after", "Slope after breakpoint")
  
  # ---- half-life transformation ----
  # No cap. Lower <- slope Lower, Upper <- slope Upper (monotone transform,
  # no swap). Positive slope -> NA for that bound.
  # Then: if Upper half-life is NA, replace with 100 (plotting convenience).
  slope_to_halflife <- function(params, slope_param, hl_param) {
    d <- params[params$Parameter == slope_param, , drop = FALSE]
    out <- data.frame(
      Vaccine   = d$Vaccine,
      Parameter = hl_param,
      Estimate  = ifelse(d$Estimate < 0, -log10(2) / d$Estimate, NA_real_),
      Lower     = ifelse(d$Lower    < 0, -log10(2) / d$Lower,    NA_real_),
      Upper     = ifelse(d$Upper    < 0, -log10(2) / d$Upper,    NA_real_)
    )
    out$Upper[is.na(out$Upper)] <- 100
    out
  }
  
  hl_linear <- slope_to_halflife(lm_params,        "Slope",        "Half_life_linear")
  hl_before <- slope_to_halflife(piecewise_params, "Slope_before", "Half_life_before")
  hl_after  <- slope_to_halflife(piecewise_params, "Slope_after",  "Half_life_after")
  
  plot_halflife <- function(data, param, title) {
    data %>%
      filter(Parameter == param) %>%
      ggplot(aes(x = Vaccine, y = Estimate, color = Vaccine)) +
      geom_point(size = 4, na.rm = TRUE) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2,
                    linewidth = 1, na.rm = TRUE) +
      scale_color_manual(values = formulation_cols, drop = FALSE) +
      labs(title = wrap_title(title), y = "Half-life (days)") +
      theme_pub()
  }
  
  p_hl_linear <- plot_halflife(hl_linear, "Half_life_linear", "Half-life from linear slope")
  p_hl_before <- plot_halflife(hl_before, "Half_life_before", "Half-life before breakpoint")
  p_hl_after  <- plot_halflife(hl_after,  "Half_life_after",  "Half-life after breakpoint")
  
  plot_lm_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(data = lm_fitted, aes(y = Fitted), linewidth = 1.1) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = wrap_title("Linear model fit (log10 scale)", width = 40),
         y = "log10(Cell count)", x = "Day") +
    facet_wrap(~Vaccine) + theme_pub()
  
  plot_piecewise_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(data = piecewise_fitted, aes(y = Fitted), linewidth = 1.1) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = wrap_title("Best piecewise model fit (log10 scale)", width = 40),
         y = "log10(Cell count)", x = "Day") +
    facet_wrap(~Vaccine) + theme_pub()
  
  per_vaccine_aic <- data.frame(
    Vaccine = vaccines,
    AIC_lm = sapply(lm_fits$mod, function(m) if(is.null(m)) NA else AIC(m)),
    AIC_piecewise = sapply(piecewise_models, function(m) if(is.null(m)) NA else AIC(m)),
    Breakpoint = best_bp_per_vaccine
  )
  
  plot_base <- file.path(figure_folder, paste0(cell_type, "_analysis"), sanitize(measurement_type))
  table_base <- file.path(table_folder, paste0(cell_type, "_analysis"), sanitize(measurement_type))
  dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_base, recursive = TRUE, showWarnings = FALSE)
  
  w_single <- 3.5; h_single <- 3.2; dpi_use <- 300
  w_facet <- 6.5; h_facet <- 4.0
  
  ggsave(file.path(plot_base, "1_Intercept_Linear.png"), p1,
         width = w_single, height = h_single, dpi = dpi_use)
  ggsave(file.path(plot_base, "2_Slope_Linear.png"), p2,
         width = w_single, height = h_single, dpi = dpi_use)
  ggsave(file.path(plot_base, "3_Intercept_Piecewise.png"), p3,
         width = w_single, height = h_single, dpi = dpi_use)
  ggsave(file.path(plot_base, "4_Slope_before_Piecewise.png"), p4,
         width = w_single, height = h_single, dpi = dpi_use)
  ggsave(file.path(plot_base, "5_Slope_after_Piecewise.png"), p5,
         width = w_single, height = h_single, dpi = dpi_use)
  ggsave(file.path(plot_base, "6_LM_fit.png"), plot_lm_fit,
         width = w_facet, height = h_facet, dpi = dpi_use)
  ggsave(file.path(plot_base, "7_Piecewise_fit.png"), plot_piecewise_fit,
         width = w_facet, height = h_facet, dpi = dpi_use)
  ggsave(file.path(plot_base, "10_Half_life_Linear.png"), p_hl_linear,
         width = w_single, height = h_single, dpi = dpi_use)
  ggsave(file.path(plot_base, "11_Half_life_before_Piecewise.png"), p_hl_before,
         width = w_single, height = h_single, dpi = dpi_use)
  ggsave(file.path(plot_base, "12_Half_life_after_Piecewise.png"), p_hl_after,
         width = w_single, height = h_single, dpi = dpi_use)
  
  write.csv(lm_params, file.path(table_base, "linear_params.csv"), row.names = FALSE)
  write.csv(piecewise_params, file.path(table_base, "piecewise_params.csv"), row.names = FALSE)
  write.csv(hl_linear, file.path(table_base, "half_life_linear.csv"), row.names = FALSE)
  write.csv(hl_before, file.path(table_base, "half_life_before_piecewise.csv"), row.names = FALSE)
  write.csv(hl_after,  file.path(table_base, "half_life_after_piecewise.csv"), row.names = FALSE)
  write.csv(per_vaccine_aic, file.path(table_base, "AIC_per_vaccine.csv"), row.names = FALSE)
  
  return(per_vaccine_aic)
}

forced_bp_config <- list(
  CD4 = list(
    "Number of Ly6C+ FOXP3- Teff cells" = "second_timepoint"
  ),
  CD8 = list(
    ".default" = "third_timepoint"
  )
)

get_forced_bp <- function(cell_type, measurement_type) {
  ct_cfg <- forced_bp_config[[cell_type]]
  if (is.null(ct_cfg)) return(NULL)
  if (!is.null(ct_cfg[[measurement_type]])) return(ct_cfg[[measurement_type]])
  if (!is.null(ct_cfg[[".default"]])) return(ct_cfg[[".default"]])
  NULL
}

combos <- combined_data %>%
  filter(Unit == "number", Cell_Type %in% c("CD4", "CD8")) %>%
  distinct(Cell_Type, Measurement_type) %>%
  arrange(Cell_Type, Measurement_type)

for (i in 1:nrow(combos)) {
  ct <- combos$Cell_Type[i]
  mt <- combos$Measurement_type[i]
  fb <- get_forced_bp(ct, mt)
  tryCatch(
    run_analysis(combined_data, ct, mt, forced_bp = fb),
    error = function(e) {
    }
  )
}