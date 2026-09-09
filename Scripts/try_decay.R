get_subset <- function(data, cell_types, measurement_types, unit = "number") {
  data %>%
    filter(Cell_Type %in% cell_types,
           Unit == unit,
           Measurement_type %in% measurement_types)
}

# Use the function
subset_df <- get_subset(
  data = combined_data,
  cell_types = "CD4",
  measurement_types = "Number of tetramer+ cells",
  unit = "number"
)
#subset_df <- subset_df %>% filter(Mouse_number != 11171)

# ---------------------------
# 1. Prepare data
# ---------------------------
subset_df <- subset_df %>% mutate(log10_Value = log10(Value))
vaccines <- c("TM", "TMd21", "SOL", "IC")

# ---------------------------
# 2. Simple linear model (per vaccine)
# ---------------------------
lm_fits <- subset_df %>%
  group_by(Vaccine) %>%
  do(mod = lm(log10_Value ~ Day, data = .))

# Extract parameters
lm_params <- bind_rows(lapply(vaccines, function(v) {
  mod <- lm_fits$mod[[which(lm_fits$Vaccine == v)]]
  ci <- confint(mod)
  data.frame(Vaccine = v,
             Parameter = c("Intercept", "Slope"),
             Estimate = coef(mod),
             Lower = ci[,1], Upper = ci[,2])
}))

# Generate fitted lines for linear model (mean only)
lm_fitted <- bind_rows(lapply(vaccines, function(v) {
  dat <- subset_df %>% filter(Vaccine == v)
  day_seq <- seq(min(dat$Day), max(dat$Day), length.out = 50)
  pred <- predict(lm_fits$mod[[which(lm_fits$Vaccine == v)]],
                  newdata = data.frame(Day = day_seq))
  data.frame(Vaccine = v, Day = day_seq, Fitted = pred)
}))

# ---------------------------
# 3. Best breakpoint per vaccine (simple piecewise)
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

# Add Day_plus using each vaccine's best breakpoint
subset_df <- subset_df %>%
  mutate(Day_plus = pmax(0, Day - best_bp[Vaccine]))

# ---------------------------
# 4. Piecewise model (per vaccine)
# ---------------------------
piecewise_models <- list()
piecewise_params <- bind_rows(lapply(vaccines, function(v) {
  dat <- subset_df %>% filter(Vaccine == v)
  mod <- lm(log10_Value ~ Day + Day_plus, data = dat)
  piecewise_models[[v]] <<- mod   # store model for later predictions
  coefs <- coef(mod); vc <- vcov(mod); ci <- confint(mod)
  slope_after <- coefs[2] + coefs[3]
  se_after <- sqrt(vc[2,2] + vc[3,3] + 2*vc[2,3])
  data.frame(Vaccine = v,
             Parameter = c("Intercept", "Slope_before", "Slope_after"),
             Estimate = c(coefs[1], coefs[2], slope_after),
             Lower = c(ci[1,1], ci[2,1], slope_after - 1.96*se_after),
             Upper = c(ci[1,2], ci[2,2], slope_after + 1.96*se_after))
}))

# Generate fitted lines for piecewise model (mean only)
piecewise_fitted <- bind_rows(lapply(vaccines, function(v) {
  dat <- subset_df %>% filter(Vaccine == v)
  day_seq <- seq(min(dat$Day), max(dat$Day), length.out = 50)
  newdata <- data.frame(Day = day_seq,
                        Day_plus = pmax(0, day_seq - best_bp[v]))
  pred <- predict(piecewise_models[[v]], newdata = newdata)
  data.frame(Vaccine = v, Day = day_seq, Fitted = pred)
}))

# ---------------------------
# 5. Coefficient plots (five separate)
# ---------------------------
plot_coef <- function(data, param, title, color) {
  data %>% filter(Parameter == param) %>%
    ggplot(aes(x = Vaccine, y = Estimate)) +
    geom_point(size = 3, color = color) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = color) +
    labs(title = title, y = "Estimate") +
    theme_minimal()
}

p1 <- plot_coef(lm_params, "Intercept", "Intercept (Linear)", "steelblue")
p2 <- plot_coef(lm_params, "Slope", "Slope (Linear)", "darkred")
p3 <- plot_coef(piecewise_params, "Intercept", "Intercept (Piecewise)", "forestgreen")
p4 <- plot_coef(piecewise_params, "Slope_before", "Slope before breakpoint", "purple")
p5 <- plot_coef(piecewise_params, "Slope_after", "Slope after breakpoint", "orange")

# ---------------------------
# 6. Model vs data plots (mean lines)
# ---------------------------
plot_lm_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
  geom_point(alpha = 0.6) +
  geom_line(data = lm_fitted, aes(y = Fitted), size = 1) +
  labs(title = "Linear model fit (log10 scale)",
       y = "log10(Value)", x = "Day") +
  facet_wrap(~Vaccine) +
  theme_minimal()

plot_piecewise_fit <- ggplot(subset_df, aes(x = Day, y = log10_Value, color = Vaccine)) +
  geom_point(alpha = 0.6) +
  geom_line(data = piecewise_fitted, aes(y = Fitted), size = 1) +
  labs(title = "Best piecewise model fit (log10 scale)",
       y = "log10(Value)", x = "Day") +
  facet_wrap(~Vaccine) +
  theme_minimal()

# Display all plots
print(p1); print(p2); print(p3); print(p4); print(p5)
print(plot_lm_fit)
print(plot_piecewise_fit)

# ---------------------------
# 7. Slope grouping analysis
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
partitions <- generate_partitions(4)

# Fit linear models for each grouping
fit_lm_grouping <- function(cl) {
  slope_cluster <- factor(cl[match(subset_df$Vaccine, vaccines)])
  if (length(unique(cl)) == 1) {
    mod <- lm(log10_Value ~ Vaccine + Day, data = subset_df)
  } else {
    mod <- lm(log10_Value ~ Vaccine + Day:slope_cluster, data = subset_df)
  }
  list(aic = AIC(mod), cluster = cl)
}
lm_results <- lapply(partitions, fit_lm_grouping)

# Fit piecewise models for each combination of before/after groupings
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

group_combos <- expand.grid(idx_before = 1:15, idx_after = 1:15)
piecewise_results <- vector("list", nrow(group_combos))
for (i in 1:nrow(group_combos)) {
  cl_before <- partitions[[group_combos$idx_before[i]]]
  cl_after  <- partitions[[group_combos$idx_after[i]]]
  piecewise_results[[i]] <- fit_piecewise_grouping(cl_before, cl_after)
}

# ---------------------------
# 8. AIC table
# ---------------------------
format_group <- function(cl) {
  groups <- split(vaccines, cl)
  paste(sapply(groups, paste, collapse = " & "), collapse = " | ")
}

linear_table <- data.frame(
  Model = "Linear",
  Slope_group = sapply(lm_results, function(x) format_group(x$cluster)),
  k = 4 + sapply(lm_results, function(x) length(unique(x$cluster))),
  AIC = sapply(lm_results, function(x) x$aic)
)

piecewise_table <- data.frame(
  Model = "Piecewise",
  Slope_group = paste0("Early: ", sapply(piecewise_results, function(x) format_group(x$cl_before)),
                       "  |  Late: ", sapply(piecewise_results, function(x) format_group(x$cl_after))),
  k = 4 + sapply(piecewise_results, function(x) length(unique(x$cl_before))) +
    sapply(piecewise_results, function(x) length(unique(x$cl_after))),
  AIC = sapply(piecewise_results, function(x) x$aic)
)

full_table <- rbind(linear_table, piecewise_table)
full_table$Delta_AIC <- full_table$AIC - min(full_table$AIC)
full_table <- full_table[order(full_table$AIC), ]

cat("\n===== TOP MODELS (by AIC) =====\n")
print(head(full_table, 10), row.names = FALSE)

# Optional: save full table
write.csv(full_table, "all_models_AIC.csv", row.names = FALSE)