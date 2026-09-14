# rm(list = ls())

library(MASS)
library(graph)
library(here)
library(dplyr)
library(knitr)
library(kableExtra)
library(ggplot2)
library(tidyr)
library(scales) # Required for the trans_new() function

# Activate the project via band.Rproj
source(here("histogram_wrappers.R"))
source(here("density.R"))
source(here("dag.R"))
source(here("dgp.R"))
source(here("compute_H.R"))
source(here("evaluation_pipeline.R"))

cube_root <- function(n, bin_factor = 1) {
  ceiling(bin_factor * n**(1/3))
}

# ==========================================
# 1. EXPERIMENT CONFIGURATION & GRIDS
# ==========================================
R_trials <- 50
n_test <- 4000
mix_probability <- 0.5
baseline_bin_factors <- c(0.5)

# The specific AR lags we want to test
k_values <- c(0, 1, 2, 3, 4)

# Create GMM parameter grid 
gmm_grid <- expand.grid(
  p = c(20),
  n_total = c(10000),
  active_dgp = "GMM",
  r_dist = c(0, 1, 2, 3, 4),
  prob_edge = NA,
  stringsAsFactors = FALSE
)

experiment_grid <- gmm_grid
total_runs <- nrow(experiment_grid) * R_trials

# ==========================================
# 2. INITIALIZE MASTER STORAGE
# ==========================================
results_ll <- data.frame(
  Run_ID = 1:total_runs,
  Scenario_ID = integer(total_runs),
  active_dgp = character(total_runs),
  p = integer(total_runs),
  n_total = integer(total_runs),
  r_dist = numeric(total_runs),
  prob_edge = numeric(total_runs),
  Trial = integer(total_runs),
  Oracle = numeric(total_runs),
  Baseline_DAG_k0 = numeric(total_runs),
  Baseline_DAG_k1 = numeric(total_runs),
  Baseline_DAG_k2 = numeric(total_runs),
  Baseline_DAG_k3 = numeric(total_runs),
  Baseline_DAG_k4 = numeric(total_runs),
  stringsAsFactors = FALSE
)

results_time <- results_ll 

cat(sprintf("\nTotal Scenarios: %d\n", nrow(experiment_grid)))
cat(sprintf("Total Independent Runs: %d\n\n", total_runs))

# ==========================================
# 3. MAIN EXPERIMENT NESTED LOOP
# ==========================================
run_index <- 1 

for (s in 1:nrow(experiment_grid)) {
  
  p          <- experiment_grid$p[s]
  n_total    <- experiment_grid$n_total[s]
  active_dgp <- experiment_grid$active_dgp[s]
  r_dist     <- experiment_grid$r_dist[s]
  prob_edge  <- experiment_grid$prob_edge[s]
  
  cat("\n**************************************************\n")
  cat(sprintf("SCENARIO %d OF %d\n", s, nrow(experiment_grid)))
  cat(sprintf("DGP: %s | p: %d | n_total: %d | r_dist: %g\n", active_dgp, p, n_total, r_dist))
  cat("**************************************************\n")
  
  for (trial in 1:R_trials) {
    
    # A. Generate Independent Data
    full_train_data <- generate_gmm_data(n = n_total, p = p, r = r_dist, mix_prob = mix_probability)
    x_test <- generate_gmm_data(n = n_test, p = p, r = r_dist, mix_prob = mix_probability)
    
    # B. Pre-process Data 
    train_idx <- sample(seq_len(n_total), size = floor(0.7 * n_total))
    train_data <- full_train_data[train_idx, , drop = FALSE]
    val_data <- full_train_data[-train_idx, , drop = FALSE]
    
    # Storage for this specific trial's lag results
    trial_ll_k <- numeric(length(k_values))
    trial_time_k <- numeric(length(k_values))
    names(trial_ll_k) <- paste0("k", k_values)
    names(trial_time_k) <- paste0("k", k_values)
    
    # C. Loop over all target AR Lags
    for (i in seq_along(k_values)) {
      k <- k_values[i]
      initial_dag_path <- get_manual_ar_dag(p, k)
      
      exec_time <- system.time({
        b_res <- tune_baseline_dag(train_data, 
                                   val_data, 
                                   x_test, 
                                   initial_dag_path, 
                                   baseline_bin_factors,
                                   epsilons = c(0.05))
      })
      
      trial_ll_k[i] <- b_res$test_log_likelihood
      trial_time_k[i] <- exec_time["elapsed"]
    }
    
    # Run Baseline comparators
    time_oracle <- system.time({ oracle_ll <- evaluate_oracle_baseline(x_test, active_dgp, p, r_dist, mix_probability) })
    
    # D. Save Results to Storage
    results_ll[run_index, c("Scenario_ID", "active_dgp", "p", "n_total", "r_dist", "prob_edge", "Trial")] <- 
      list(s, active_dgp, p, n_total, r_dist, prob_edge, trial)
    results_ll$Oracle[run_index] <- oracle_ll
    results_ll$Baseline_DAG_k0[run_index]  <- trial_ll_k["k0"]
    results_ll$Baseline_DAG_k1[run_index]  <- trial_ll_k["k1"]
    results_ll$Baseline_DAG_k2[run_index]  <- trial_ll_k["k2"]
    results_ll$Baseline_DAG_k3[run_index]  <- trial_ll_k["k3"]
    results_ll$Baseline_DAG_k4[run_index]  <- trial_ll_k["k4"]
    
    results_time[run_index, c("Scenario_ID", "active_dgp", "p", "n_total", "r_dist", "prob_edge", "Trial")] <- 
      list(s, active_dgp, p, n_total, r_dist, prob_edge, trial)
    results_time$Oracle[run_index] <- time_oracle["elapsed"]
    results_time$Baseline_DAG_k0[run_index]  <- trial_time_k["k0"]
    results_time$Baseline_DAG_k1[run_index]  <- trial_time_k["k1"]
    results_time$Baseline_DAG_k2[run_index]  <- trial_time_k["k2"]
    results_time$Baseline_DAG_k3[run_index]  <- trial_time_k["k3"]
    results_time$Baseline_DAG_k4[run_index]  <- trial_time_k["k4"]
    
    # E. Display Trial Results
    cat("\n==================================================\n")
    cat(sprintf("SCENARIO %d | TRIAL %d: HELD-OUT TEST LOG-LIKELIHOODS\n", s, trial))
    cat("==================================================\n")
    cat("True Oracle:       ", oracle_ll, "\n")
    for(k in k_values) {
      cat(sprintf("Baseline DAG (k=%-2d): %f\n", k, trial_ll_k[paste0("k", k)]))
    }
    cat("==================================================\n\n")
    
    run_index <- run_index + 1
  }
}

# ==========================================
# 4. FINAL EXPORT 
# ==========================================
cat("\nSaving experiment results to disk...\n")
saveRDS(results_ll, "experiment_gmm_ar_lags_loglikelihoods.rds")



# ==========================================
# 5. AGGREGATE AND FORMAT TABLE
# ==========================================

# Function for Oracle (mean only)
format_mean_only <- function(x) {
  sprintf("%.1f", mean(x, na.rm = TRUE))
}

# Function for Baselines (mean and standard deviation)
format_mean_sd <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  sprintf("%.1f (%.1f)", m, s)
}

results_summary <- results_ll %>%
  group_by(Scenario_ID, active_dgp, p, n_total, r_dist, prob_edge) %>%
  summarize(
    Oracle = format_mean_only(Oracle),
    across(c(Baseline_DAG_k0, Baseline_DAG_k1, Baseline_DAG_k2, Baseline_DAG_k3, Baseline_DAG_k4), format_mean_sd),
    .groups = "drop"
  )

gmm_data_clean <- results_summary %>% 
  filter(active_dgp == "GMM") %>% 
  arrange(n_total, p, r_dist) %>%
  select(
    `$r$` = r_dist, 
    Oracle, 
    `$\\text{n\\_lag}=0$` = Baseline_DAG_k0, 
    `$\\text{n\\_lag}=1$` = Baseline_DAG_k1, 
    `$\\text{n\\_lag}=2$` = Baseline_DAG_k2, 
    `$\\text{n\\_lag}=3$` = Baseline_DAG_k3, 
    `$\\text{n\\_lag}=4$` = Baseline_DAG_k4
  )

# 1. Generate ONLY the tabular core (no caption means no \begin{table} wrapper)
tabular_core <- kable(gmm_data_clean, 
                      format = "latex", 
                      booktabs = TRUE, 
                      escape = FALSE, 
                      linesep = "")

# Inject the [t] into the tabular environment
tabular_core <- sub("\\begin{tabular}", "\\begin{tabular}[t]", tabular_core, fixed = TRUE)

# 2. Define your exact custom LaTeX header (remember to escape backslashes in R strings)
latex_header <- "\\begin{table}
\\centering
\\caption{Data are simulated from Example~\\ref{example2} with $(n, p, K) = (10000, 20, 2)$ over 50 independent replications, corresponding to the data used in Figure~\\ref{fig:1}. Standard deviations are reported in parentheses. All results are based on the $\\text{BAND}_0$ estimator with autoregressive DAG with $\\text{n\\_lag}$ lags.}\\label{tab:1}
\\setlength{\\tabcolsep}{1pt}
\\centering
\\fontsize{8}{8}\\selectfont\n"

# 3. Concatenate the header, the tabular core, and the closing tag
final_latex_table <- paste0(
  latex_header,
  as.character(tabular_core),
  "\n\\end{table}\n"
)

# 4. Print the final LaTeX code to the console
cat(final_latex_table)




###
###
###
# Output Figure:
# 1. Define a custom piecewise transformation
compress_trans <- function(compression_factor = 4) {
  trans <- function(x) {
    ifelse(is.na(x), NA, ifelse(x >= -36, x, -36 + (x + 36) / compression_factor))
  }
  inv <- function(x) {
    ifelse(is.na(x), NA, ifelse(x >= -36, x, -36 + (x + 36) * compression_factor))
  }
  trans_new("compress", trans, inv)
}

# 2. Prepare the numeric data for plotting
plot_data <- results_ll %>%
  filter(active_dgp == "GMM") %>%
  group_by(r_dist) %>%
  summarize(
    k0 = mean(Baseline_DAG_k0, na.rm = TRUE),
    k1 = mean(Baseline_DAG_k1, na.rm = TRUE),
    k2 = mean(Baseline_DAG_k2, na.rm = TRUE),
    k3 = mean(Baseline_DAG_k3, na.rm = TRUE),
    k4 = mean(Baseline_DAG_k4, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Reshape from wide to long format
  pivot_longer(
    cols = starts_with("k"),
    names_to = "n_lag",
    values_to = "Mean_LogLikelihood"
  ) %>%
  # Clean up the n_lag column so it is purely numeric
  mutate(n_lag = as.numeric(gsub("k", "", n_lag)))

# 3. Plot: Apply transformations and new legend styling
plot_B <- ggplot(plot_data, aes(
  x = as.factor(r_dist), 
  y = Mean_LogLikelihood, 
  color = as.factor(n_lag),
  group = n_lag 
)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  
  # Trim the lightest colors by setting end = 0.85
  scale_color_viridis_d(option = "viridis", end = 0.85) + 
  
  # Apply the custom y-axis scale
  scale_y_continuous(
    trans = compress_trans(compression_factor = 4),
    breaks = c(-50, -45, -40, -36, -34, -32, -30, -28) 
  ) +
  
  labs(
    x = "r",
    y = "Mean Log-Likelihood",
    color = "n_lag"
  ) +
  theme_bw() +
  theme(
    # Anchor legend to the top-right corner
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1), 
    legend.direction = "horizontal",
    
    # Make the font and spacing even smaller
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm"),   # Shrinks the colored line segments
    legend.margin = margin(t = 2, r = 4, b = 2, l = 4), # Tighter padding inside the box
    
    # White background box
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3)
  )

print(plot_B)