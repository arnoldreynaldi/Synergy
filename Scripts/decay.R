sanitize <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

safe_predict <- function(model, newdata) {
  if (is.null(model)) return(rep(NA, nrow(newdata)))
  tryCatch(predict(model, newdata), error = function(e) rep(NA, nrow(newdata)))
}

# ---- Global plot theme: smaller canvas, bigger text/labels ----
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
  if (any(n_per_vaccine$n < 1)) {
    return(NULL)
  }
  
  subset_df$Vaccine <- factor(subset_df$Vaccine, levels = vaccine_levels)
  vaccines <- levels(droplevels(subset_df$Vaccine))
  if (length(vaccines) < 1) {
    return(NULL)
  }
  
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
    message("Forcing breakpoint at second timepoint = Day ", bp_candidates,
            " for Cell_Type = ", cell_type,
            ", measurement_type = ", measurement_type)
  } else if (is.numeric(forced_bp)) {
    bp_candidates <- forced_bp
    message("Forcing breakpoint at Day ", forced_bp,
            " for Cell_Type = ", cell_type,
            ", measurement_type = ", measurement_type)
  } else {
    stop("`forced_bp` must be NULL, a numeric value, or \"second_timepoint\"")
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
  
  # Coefficient plots
  plot_coef <- function(data, param, title) {
    data %>% filter(Parameter == param) %>%
      ggplot(aes(x = Vaccine, y = Estimate, color = Vaccine)) +
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, linewidth = 1) +
      scale_color_manual(values = formulation_cols, drop = FALSE) +
      labs(title = title, y = "Estimate") +
      theme_pub()
  }
  
  p1 <- plot_coef(lm_params, "Intercept", "Intercept (Linear)")
  p2 <- plot_coef(lm_params, "Slope", "Slope (Linear)")
  p3 <- plot_coef(piecewise_params, "Intercept", "Intercept (Piecewise)")
  p4 <- plot_coef(piecewise_params, "Slope_before", "Slope before breakpoint")
  p5 <- plot_coef(piecewise_params, "Slope_after", "Slope after breakpoint")
  
  # --- half-life transformation from slope parameters ---
  HL_CAP <- 1000  # days
  
  slope_to_halflife <- function(params, slope_param, hl_param) {
    params %>%
      filter(Parameter == slope_param) %>%
      mutate(
        Parameter = hl_param,
        Estimate = ifelse(Estimate < 0, -log10(2) / Estimate, NA_real_),
        Lower    = ifelse(Lower < 0, -log10(2) / Lower, NA_real_),
        Upper    = ifelse(Upper < 0, -log10(2) / Upper,
                          ifelse(Lower < 0, Inf, NA_real_))
      ) %>%
      mutate(
        Estimate = ifelse(is.na(Estimate) | !is.finite(Estimate) | Lower > 0,
                          HL_CAP, Estimate),
        Lower    = ifelse(is.na(Lower) | !is.finite(Lower),
                          NA_real_, Lower),
        Upper    = ifelse(is.infinite(Upper) | is.na(Upper) | Upper > HL_CAP,
                          HL_CAP, Upper),
        Lower    = ifelse(!is.na(Lower) & Lower > HL_CAP, HL_CAP, Lower)
      )
  }
  
  hl_linear <- slope_to_halflife(lm_params, "Slope", "Half_life_linear")
  hl_before <- slope_to_halflife(piecewise_params, "Slope_before", "Half_life_before")
  hl_after  <- slope_to_halflife(piecewise_params, "Slope_after",  "Half_life_after")
  
  plot_halflife <- function(data, param, title) {
    data %>%
      filter(Parameter == param) %>%
      ggplot(aes(x = Vaccine, y = Estimate, color = Vaccine)) +
      geom_point(size = 4, na.rm = TRUE) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2,
                    linewidth = 1, na.rm = TRUE) +
      geom_hline(yintercept = HL_CAP, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = formulation_cols, drop = FALSE) +
      coord_cartesian(ylim = c(0, HL_CAP * 1.05)) +
      labs(title = title, y = "Half-life (days, capped at 1000)") +
      theme_pub()
  }
  
  p_hl_linear <- plot_halflife(hl_linear, "Half_life_linear", "Half-life from linear slope")
  p_hl_before <- plot_halflife(hl_before, "Half_life_before", "Half-life before breakpoint")
  p_hl_after  <- plot_halflife(hl_after,  "Half_life_after",  "Half-life after breakpoint")
  
  plot_lm_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(data = lm_fitted, aes(y = Fitted), linewidth = 1.1) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Linear model fit (log10 scale)",
         y = "log10(Cell count)", x = "Day") +
    facet_wrap(~Vaccine) + theme_pub()
  
  plot_piecewise_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(data = piecewise_fitted, aes(y = Fitted), linewidth = 1.1) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Best piecewise model fit (log10 scale)",
         y = "log10(Cell count)", x = "Day") +
    facet_wrap(~Vaccine) + theme_pub()
  
  # Full search
  subset_df <- subset_df %>% select(-Day_plus)
  
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
  partitions <- generate_partitions(length(vaccines))
  
  bp_combos <- expand.grid(rep(list(bp_candidates), length(vaccines)))
  names(bp_combos) <- vaccines
  bp_combos_list <- asplit(bp_combos, 1)
  
  fit_linear_grouping <- function(cl) {
    slope_cluster <- factor(cl[match(subset_df$Vaccine, vaccines)])
    if (length(unique(cl)) == 1) {
      mod <- tryCatch(lm(log10_Value ~ Vaccine + Day, data = subset_df),
                      error = function(e) NULL)
    } else {
      mod <- tryCatch(lm(log10_Value ~ Vaccine + Day:slope_cluster, data = subset_df),
                      error = function(e) NULL)
    }
    if (is.null(mod)) return(NULL)
    list(aic = AIC(mod), k = length(vaccines) + length(unique(cl)),
         bp_vec = NULL, cl_before = cl, cl_after = NULL, model = mod)
  }
  
  fit_piecewise_full <- function(bp_vec, cl_before, cl_after) {
    data_temp <- subset_df %>%
      mutate(Day_plus = pmax(0, Day - bp_vec[Vaccine]))
    before_cluster <- cl_before[match(data_temp$Vaccine, vaccines)]
    after_cluster  <- cl_after[match(data_temp$Vaccine, vaccines)]
    formula_str <- "log10_Value ~ Vaccine"
    if (length(unique(cl_before)) == 1) {
      formula_str <- paste(formula_str, "+ Day")
    } else {
      data_temp$group_before <- factor(before_cluster)
      formula_str <- paste(formula_str, "+ Day:group_before")
    }
    if (length(unique(cl_after)) == 1) {
      formula_str <- paste(formula_str, "+ Day_plus")
    } else {
      data_temp$group_after <- factor(after_cluster)
      formula_str <- paste(formula_str, "+ Day_plus:group_after")
    }
    mod <- tryCatch(lm(as.formula(formula_str), data = data_temp),
                    error = function(e) NULL)
    if (is.null(mod)) return(NULL)
    list(aic = AIC(mod),
         k = length(vaccines) + length(unique(cl_before)) + length(unique(cl_after)),
         bp_vec = bp_vec, cl_before = cl_before, cl_after = cl_after, model = mod)
  }
  
  all_models <- list()
  for (i in seq_along(partitions)) {
    res <- fit_linear_grouping(partitions[[i]])
    if (!is.null(res)) all_models[[length(all_models)+1]] <- res
  }
  for (bp_vec in bp_combos_list) {
    for (i in seq_along(partitions)) {
      for (j in seq_along(partitions)) {
        res <- fit_piecewise_full(bp_vec, partitions[[i]], partitions[[j]])
        if (!is.null(res)) all_models[[length(all_models)+1]] <- res
      }
    }
  }
  
  all_models <- Filter(function(x) !is.null(x) && !is.null(x$model), all_models)
  if (length(all_models) == 0) {
    return(NULL)
  }
  
  aic_values <- sapply(all_models, function(x) x$aic)
  sorted_idx <- order(aic_values)
  all_models_sorted <- all_models[sorted_idx]
  aic_sorted <- aic_values[sorted_idx]
  
  best_overall <- all_models_sorted[[1]]
  best_aic <- aic_sorted[1]
  
  delta_aic <- aic_sorted - best_aic
  eligible <- which(delta_aic < 10)
  if (length(eligible) == 0) {
    parsimonious <- best_overall
    parsimonious_delta_aic <- 0
  } else {
    eligible_models <- all_models_sorted[eligible]
    eligible_k <- sapply(eligible_models, function(x) x$k)
    min_k <- min(eligible_k)
    candidates <- which(eligible_k == min_k)
    best_candidate <- candidates[which.min(aic_sorted[eligible][candidates])]
    parsimonious <- eligible_models[[best_candidate]]
    parsimonious_delta_aic <- aic_sorted[eligible][best_candidate] - best_aic
  }
  parsimonious_k <- parsimonious$k
  
  format_group <- function(cl) {
    if (is.null(cl)) return("")
    groups <- split(vaccines, cl)
    paste(sapply(groups, paste, collapse = " & "), collapse = " | ")
  }
  
  combined_table <- data.frame(
    Model = ifelse(sapply(all_models_sorted, function(x) is.null(x$bp_vec)), "Linear", "Piecewise"),
    Breakpoints = sapply(all_models_sorted, function(x) {
      if (is.null(x$bp_vec)) return("")
      paste(paste0(names(x$bp_vec), "=", x$bp_vec), collapse = ", ")
    }),
    Early_group = sapply(all_models_sorted, function(x) format_group(x$cl_before)),
    Late_group = sapply(all_models_sorted, function(x) format_group(x$cl_after)),
    k = sapply(all_models_sorted, function(x) x$k),
    AIC = aic_sorted
  )
  combined_table$Delta_AIC <- combined_table$AIC - best_aic
  
  make_plot <- function(model_obj, title_prefix, extra_info = "") {
    if (is.null(model_obj) || is.null(model_obj$model)) {
      return(ggplot() + labs(title = paste(title_prefix, "- unavailable")) + theme_pub())
    }
    day_grid <- seq(min(subset_df$Day), max(subset_df$Day), length.out = 100)
    pred_data <- expand.grid(Vaccine = vaccines, Day = day_grid)
    pred_data$Vaccine <- factor(pred_data$Vaccine, levels = vaccine_levels)
    
    if (is.null(model_obj$bp_vec)) {
      cl <- model_obj$cl_before
      if (length(unique(cl)) > 1) {
        pred_data$slope_cluster <- factor(cl[match(pred_data$Vaccine, vaccines)])
      }
      pred_data$Fitted <- safe_predict(model_obj$model, pred_data)
      plot_title <- paste(title_prefix, "(Linear)")
      subtitle <- paste("Slope grouping:", format_group(cl))
    } else {
      bp_vec <- model_obj$bp_vec
      pred_data$Day_plus <- pmax(0, pred_data$Day - bp_vec[as.character(pred_data$Vaccine)])
      if (length(unique(model_obj$cl_before)) > 1) {
        pred_data$group_before <- factor(model_obj$cl_before[match(pred_data$Vaccine, vaccines)])
      }
      if (length(unique(model_obj$cl_after)) > 1) {
        pred_data$group_after <- factor(model_obj$cl_after[match(pred_data$Vaccine, vaccines)])
      }
      pred_data$Fitted <- safe_predict(model_obj$model, pred_data)
      plot_title <- paste(title_prefix, "(Piecewise)")
      subtitle <- paste("Breakpoints:", paste(paste0(names(bp_vec), "=", bp_vec), collapse = ", "),
                        "\nEarly:", format_group(model_obj$cl_before),
                        " | Late:", format_group(model_obj$cl_after))
    }
    if (nzchar(extra_info)) {
      subtitle <- paste(subtitle, extra_info, sep = "\n")
    }
    ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_line(data = pred_data, aes(y = Fitted), linewidth = 1.1) +
      scale_color_manual(values = formulation_cols, drop = FALSE) +
      labs(title = plot_title, subtitle = subtitle,
           y = "log10(Cell count)", x = "Day") +
      facet_wrap(~Vaccine) +
      theme_pub()
  }
  
  plot_best_overall <- make_plot(best_overall, "Best overall model")
  plot_parsimonious <- make_plot(
    parsimonious,
    "Parsimonious model (ΔAIC<10, min k)",
    extra_info = paste0("ΔAIC = ", round(parsimonious_delta_aic, 2),
                        ", k = ", parsimonious_k)
  )
  
  plot_base <- file.path(figure_folder, paste0(cell_type, "_analysis"), sanitize(measurement_type))
  table_base <- file.path(table_folder, paste0(cell_type, "_analysis"), sanitize(measurement_type))
  dir.create(plot_base, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_base, recursive = TRUE, showWarnings = FALSE)
  
  # Smaller canvas + higher dpi so text stays readable when figure is shrunk
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
  ggsave(file.path(plot_base, "8_Best_overall_model_fit.png"), plot_best_overall,
         width = w_facet, height = h_facet, dpi = dpi_use)
  ggsave(file.path(plot_base, "9_Parsimonious_model_fit.png"), plot_parsimonious,
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
  
  per_vaccine_aic <- data.frame(
    Vaccine = vaccines,
    AIC_lm = sapply(lm_fits$mod, function(m) if(is.null(m)) NA else AIC(m)),
    AIC_piecewise = sapply(piecewise_models, function(m) if(is.null(m)) NA else AIC(m)),
    Breakpoint = best_bp_per_vaccine
  )
  write.csv(per_vaccine_aic, file.path(table_base, "AIC_per_vaccine.csv"), row.names = FALSE)
  
  write.csv(combined_table, file.path(table_base, "AIC_table_full_search.csv"), row.names = FALSE)
  
  return(combined_table)
}

# ---- Forced breakpoint config ----
# Names MUST match combined_data$Measurement_type EXACTLY.
forced_bp_config <- list(
  CD4 = list(
    "Number of Ly6C+ FOXP3- Teff cells" = "second_timepoint"
  )
)

get_forced_bp <- function(cell_type, measurement_type) {
  ct_cfg <- forced_bp_config[[cell_type]]
  if (is.null(ct_cfg)) return(NULL)
  ct_cfg[[measurement_type]]
}

# Main loop
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