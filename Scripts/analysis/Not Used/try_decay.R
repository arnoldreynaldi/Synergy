
sanitize <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

safe_predict <- function(model, newdata) {
  if (is.null(model)) return(rep(NA, nrow(newdata)))
  tryCatch(predict(model, newdata), error = function(e) rep(NA, nrow(newdata)))
}

run_analysis <- function(data, cell_type, measurement_type) {
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
  
  # Per‑vaccine linear model
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
  
  # Per‑vaccine piecewise model
  bp_candidates <- c(14, 28)
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
  
  plot_lm_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
    geom_point(alpha = 0.6) +
    geom_line(data = lm_fitted, aes(y = Fitted), size = 1) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Linear model fit (log10 scale)", y = "log10(Cell count)", x = "Day") +
    facet_wrap(~Vaccine) + theme_bw()
  
  plot_piecewise_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
    geom_point(alpha = 0.6) +
    geom_line(data = piecewise_fitted, aes(y = Fitted), size = 1) +
    scale_color_manual(values = formulation_cols, drop = FALSE) +
    labs(title = "Best piecewise model fit (log10 scale)", y = "log10(Cell count)", x = "Day") +
    facet_wrap(~Vaccine) + theme_bw()
  
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
  
  # Filter out any NULLs
  all_models <- Filter(function(x) !is.null(x) && !is.null(x$model), all_models)
  if (length(all_models) == 0) {
    return(NULL)
  }
  
  # Sort by AIC
  aic_values <- sapply(all_models, function(x) x$aic)
  sorted_idx <- order(aic_values)
  all_models_sorted <- all_models[sorted_idx]
  aic_sorted <- aic_values[sorted_idx]
  
  best_overall <- all_models_sorted[[1]]
  best_aic <- aic_sorted[1]
  
  # Parsimonious selection: ΔAIC < 10 and minimum k
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
  
  # Build combined AIC table
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
  
  # Helper to make a fit plot safely, with optional extra info
  make_plot <- function(model_obj, title_prefix, extra_info = "") {
    if (is.null(model_obj) || is.null(model_obj$model)) {
      return(ggplot() + labs(title = paste(title_prefix, "- unavailable")))
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
      geom_point(alpha = 0.6) +
      geom_line(data = pred_data, aes(y = Fitted), size = 1) +
      scale_color_manual(values = formulation_cols, drop = FALSE) +
      labs(title = plot_title, subtitle = subtitle,
           y = "log10(Cell count)", x = "Day") +
      facet_wrap(~Vaccine) +
      theme_bw()
  }
  
  # Generate plots
  plot_best_overall <- make_plot(best_overall, "Best overall model")
  plot_parsimonious <- make_plot(
    parsimonious,
    "Parsimonious model (ΔAIC<10, min k)",
    extra_info = paste0("ΔAIC = ", round(parsimonious_delta_aic, 2),
                        ", k = ", parsimonious_k)
  )
  
  # Save outputs
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
  ggsave(file.path(plot_base, "8_Best_overall_model_fit.png"), plot_best_overall, width = 8, height = 5)
  ggsave(file.path(plot_base, "9_Parsimonious_model_fit.png"), plot_parsimonious, width = 8, height = 5)
  
  write.csv(lm_params, file.path(table_base, "linear_params.csv"), row.names = FALSE)
  write.csv(piecewise_params, file.path(table_base, "piecewise_params.csv"), row.names = FALSE)
  
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

# Main loop
combos <- combined_data %>%
  filter(Unit == "number", Cell_Type %in% c("CD4", "CD8")) %>%
  distinct(Cell_Type, Measurement_type) %>%
  arrange(Cell_Type, Measurement_type)

for (i in 1:nrow(combos)) {
  ct <- combos$Cell_Type[i]
  mt <- combos$Measurement_type[i]
  tryCatch(
    run_analysis(combined_data, ct, mt),
    error = function(e) {
    }
  )
}