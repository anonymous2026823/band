# rm(list = ls())

library(MASS)
library(graph)
library(here)
# Load the dplyr package
library(dplyr)
library(knitr) # for printing latex tables
library(kableExtra) # Required for kable_styling to adjust font size

# Activate the project via band.Rproj
# The working directory should be set to the project root: .../band
source(here("histogram_wrappers.R")) # histogram wraper
source(here("density.R")) # estimate_g
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

# Algorithm Hyperparameters
search_space <- expand.grid(
  w0 = c(0.001),
  gamma = c(1),
  iota = c(3),
  bin_factor = c(0.5, 1, 2, 3),
  epsilon = c(0.05)
)
baseline_bin_factors <- unique(search_space$bin_factor)

# Create GMM parameter grid (p: 2 x n: 2 x r_dist: 3 = 12 scenarios)
gmm_grid <- expand.grid(
  p = c(5, 100),
  n_total = c(500, 4000),
  active_dgp = "GMM",
  r_dist = c(0, 1, 10),
  prob_edge = NA,
  stringsAsFactors = FALSE
)

# Create Linear DAG parameter grid (p: 2 x n: 2 x prob_edge: 2 = 8 scenarios)
linear_grid <- expand.grid(
  p = c(5, 100),
  n_total = c(500, 4000),
  active_dgp = "LINEAR",
  r_dist = NA,
  prob_edge = c(0.05, 0.1), # (0.1, 0.2)
  stringsAsFactors = FALSE
)

# Master list of all 20 scenarios
experiment_grid <- rbind(gmm_grid, linear_grid)
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
  Baseline_DAG = numeric(total_runs),
  Enhanced_DAG = numeric(total_runs),
  mclust = numeric(total_runs),
  rvinecopulib = numeric(total_runs),
  stringsAsFactors = FALSE
)

results_time <- results_ll # Duplicate the exact structure for times

cat(sprintf("\nTotal Scenarios: %d\n", nrow(experiment_grid)))
cat(sprintf("Total Independent Runs: %d\n\n", total_runs))

# ==========================================
# 3. MAIN EXPERIMENT NESTED LOOP
# ==========================================
run_index <- 1 # Tracks the global row in our storage data frames

for (s in 1:nrow(experiment_grid)) {
  
  # Extract parameters for the current scenario
  p          <- experiment_grid$p[s]
  n_total    <- experiment_grid$n_total[s]
  active_dgp <- experiment_grid$active_dgp[s]
  r_dist     <- experiment_grid$r_dist[s]
  prob_edge  <- experiment_grid$prob_edge[s]
  
  cat("\n**************************************************\n")
  cat(sprintf("SCENARIO %d OF %d\n", s, nrow(experiment_grid)))
  cat(sprintf("DGP: %s | p: %d | n_total: %d\n", active_dgp, p, n_total))
  if (active_dgp == "GMM") cat(sprintf("r_dist: %g\n", r_dist))
  if (active_dgp == "LINEAR") cat(sprintf("prob_edge: %.g\n", prob_edge))
  cat("**************************************************\n")
  
  for (trial in 1:R_trials) {
    
    # ------------------------------------------------
    # A. Generate Independent Data
    # ------------------------------------------------
    if (active_dgp == "GMM") {
      full_train_data <- generate_gmm_data(n = n_total, p = p, r = r_dist, mix_prob = mix_probability)
      x_test <- generate_gmm_data(n = n_test, p = p, r = r_dist, mix_prob = mix_probability)
    } else if (active_dgp == "LINEAR") {
      full_train_data <- generate_linear_dag_data(n = n_total, p = p, prob_edge = prob_edge)
      x_test <- generate_linear_dag_data(n = n_test, p = p, prob_edge = prob_edge)
    }
    
    # ------------------------------------------------
    # B. Pre-process Data 
    # ------------------------------------------------
    train_idx <- sample(seq_len(n_total), size = floor(0.7 * n_total))
    train_data <- full_train_data[train_idx, , drop = FALSE]
    val_data <- full_train_data[-train_idx, , drop = FALSE]
    
    initial_dag_path <- get_dag_structure_ges(full_train_data)
    
    # ------------------------------------------------
    # C. Run Evaluations & Record Time (Silent execution)
    # ------------------------------------------------
    time_baseline <- system.time({
      baseline_results <- tune_baseline_dag(train_data, 
                                            val_data, 
                                            x_test, 
                                            initial_dag_path, 
                                            baseline_bin_factors, 
                                            epsilons = 0.05)
    })
    
    time_enhanced <- system.time({
      enhanced_results <- tune_enhanced_dag(train_data, 
                                            val_data, 
                                            x_test, 
                                            initial_dag_path, 
                                            search_space)
    })
    
    time_mclust <- system.time({
      mclust_ll <- evaluate_mclust_model(train_data, x_test)
    })
    
    time_rvinecopulib <- system.time({
      rvinecopulib_ll <- evaluate_rvinecopulib_model(train_data, x_test)
    })
    
    time_oracle <- system.time({
      if (active_dgp == "GMM") {
        oracle_ll <- evaluate_oracle_baseline(x_test, active_dgp, p, r_dist, mix_probability)
      } else {
        oracle_ll <- evaluate_oracle_baseline(x_test, active_dgp, p)
      }
    })
    
    # ------------------------------------------------
    # D. Save Results to Storage
    # ------------------------------------------------
    # Log-Likelihoods
    results_ll[run_index, c("Scenario_ID", "active_dgp", "p", "n_total", "r_dist", "prob_edge", "Trial")] <- 
      list(s, active_dgp, p, n_total, r_dist, prob_edge, trial)
    results_ll$Oracle[run_index]       <- oracle_ll
    results_ll$Baseline_DAG[run_index] <- baseline_results$test_log_likelihood
    results_ll$Enhanced_DAG[run_index] <- enhanced_results$test_log_likelihood
    results_ll$mclust[run_index]       <- mclust_ll
    results_ll$rvinecopulib[run_index] <- rvinecopulib_ll
    
    # Execution Times
    results_time[run_index, c("Scenario_ID", "active_dgp", "p", "n_total", "r_dist", "prob_edge", "Trial")] <- 
      list(s, active_dgp, p, n_total, r_dist, prob_edge, trial)
    results_time$Oracle[run_index]       <- time_oracle["elapsed"]
    results_time$Baseline_DAG[run_index] <- time_baseline["elapsed"]
    results_time$Enhanced_DAG[run_index] <- time_enhanced["elapsed"]
    results_time$mclust[run_index]       <- time_mclust["elapsed"]
    results_time$rvinecopulib[run_index] <- time_rvinecopulib["elapsed"]
    
    # ------------------------------------------------
    # E. Display Trial Results
    # ------------------------------------------------
    cat("\n==================================================\n")
    cat(sprintf("SCENARIO %d | TRIAL %d: HELD-OUT TEST LOG-LIKELIHOODS\n", s, trial))
    cat("==================================================\n")
    cat("True Oracle:       ", oracle_ll, "\n")
    cat("Baseline DAG:      ", baseline_results$test_log_likelihood, "\n")
    cat("Enhanced DAG:      ", enhanced_results$test_log_likelihood, "\n")
    cat("mclust:            ", mclust_ll, "\n")
    cat("rvinecopulib:      ", rvinecopulib_ll, "\n")
    
    cat("\n==================================================\n")
    cat(sprintf("SCENARIO %d | TRIAL %d: EXECUTION TIMES (Seconds)\n", s, trial))
    cat("==================================================\n")
    cat("True Oracle:       ", time_oracle["elapsed"], "\n")
    cat("Baseline DAG:      ", time_baseline["elapsed"], "\n")
    cat("Enhanced DAG:      ", time_enhanced["elapsed"], "\n")
    cat("mclust:            ", time_mclust["elapsed"], "\n")
    cat("rvinecopulib:      ", time_rvinecopulib["elapsed"], "\n")
    cat("==================================================\n\n")
    
    run_index <- run_index + 1
  }
}


# cat("\n not saving experiment results \n")
# ==========================================
# 4. FINAL EXPORT 
# ==========================================
cat("\nSaving experiment results to disk...\n")
saveRDS(results_ll, "experiment_loglikelihoods.rds")
saveRDS(results_time, "experiment_execution_times.rds")


# ==========================================
# 5. AGGREGATE AND FORMAT TABLE
# ==========================================
# ==============================================================================
# --- Helper Functions ---
# ==============================================================================

# Define a helper function to format "Mean (SD)" for the estimators
format_mean_sd <- function(x) {
  m <- mean(x, na.rm = TRUE)
  
  # Check for artificial penalty floor and center the en-dash in the cell
  if (!is.na(m) && m < -1000) {
    return("\\multicolumn{1}{c}{--}")
  }
  
  s <- sd(x, na.rm = TRUE)
  sprintf("%.1f (%.1f)", m, s)
}

# Define a helper function to format Oracle (mean only)
format_mean_only <- function(x) {
  m <- mean(x, na.rm = TRUE)
  
  if (!is.na(m) && m < -1000) {
    return("\\multicolumn{1}{c}{--}")
  }
  
  sprintf("%.1f", m)
}

# ==============================================================================
# --- Data Aggregation ---
# ==============================================================================

# Aggregate the results
results_summary <- results_ll %>%
  group_by(Scenario_ID, active_dgp, p, n_total, r_dist, prob_edge) %>%
  summarize(
    Oracle = format_mean_only(Oracle),
    across(c(Baseline_DAG, Enhanced_DAG, mclust, rvinecopulib), format_mean_sd),
    .groups = "drop"
  )

# Sort by n, then p, then r_dist before merging columns for GMM
gmm_data_clean <- results_summary %>% 
  filter(active_dgp == "GMM") %>% 
  arrange(n_total, p, r_dist) %>%
  mutate(`p/n/r` = paste(p, n_total, r_dist, sep = "/")) %>% 
  select(`p/n/r`, Oracle, `BAND$_{0}$` = Baseline_DAG, BAND = Enhanced_DAG, mclust, rvine = rvinecopulib)

# Sort by n, then p, then prob_edge before merging columns for LINEAR
linear_data_clean <- results_summary %>% 
  filter(active_dgp == "LINEAR") %>% 
  arrange(n_total, p, prob_edge) %>%
  mutate(`p/n/prob` = paste(p, n_total, prob_edge, sep = "/")) %>% 
  select(`p/n/prob`, Oracle, `BAND$_{0}$` = Baseline_DAG, BAND = Enhanced_DAG, mclust, rvine = rvinecopulib)


# ==============================================================================
# --- Table Generation & Formatting ---
# ==============================================================================

library(knitr)

# 1. Generate ONLY the tabular cores (no caption means no \begin{table} wrapper)
gmm_tabular_core <- kable(gmm_data_clean, 
                          format = "latex", 
                          booktabs = TRUE, 
                          escape = FALSE, 
                          linesep = "")

linear_tabular_core <- kable(linear_data_clean, 
                             format = "latex", 
                             booktabs = TRUE, 
                             escape = FALSE, 
                             linesep = "")

# 2. Inject the [t] alignment into both tabular environments
gmm_tabular_core <- sub("\\begin{tabular}", "\\begin{tabular}[t]", gmm_tabular_core, fixed = TRUE)
linear_tabular_core <- sub("\\begin{tabular}", "\\begin{tabular}[t]", linear_tabular_core, fixed = TRUE)

# 3. Define the exact custom LaTeX headers
gmm_header <- "\\begin{table}
\\centering
\\caption{Mean log-likelihoods for the Gaussian Mixture Model (Setting 1) across varying $(p, n)$ and cluster separations ($r$) over 50 replications. Standard deviations are reported in parentheses.}\\label{tab:gmm}
\\setlength{\\tabcolsep}{1pt}
\\centering
\\fontsize{6.5}{8}\\selectfont\n"

linear_header <- "\\begin{table}
\\centering
\\caption{Mean log-likelihoods for the Linear Gaussian DAG (Setting 2) across varying $(p, n)$ and edge probabilities (prob) over 50 replications. Standard deviations are reported in parentheses. The symbol `-' indicates severe model failure where log-likelihoods $\\le -1000$.}\\label{tab:linear}
\\setlength{\\tabcolsep}{1pt}
\\centering
\\fontsize{6}{8}\\selectfont\n"

# 4. Concatenate the headers, the modified tabular cores, and the closing tags
final_gmm_table <- paste0(
  gmm_header,
  as.character(gmm_tabular_core),
  "\n\\end{table}\n"
)

final_linear_table <- paste0(
  linear_header,
  as.character(linear_tabular_core),
  "\n\\end{table}\n"
)

# 5. Print the final LaTeX tables to the console
cat(final_gmm_table)
cat("\n\n")
cat(final_linear_table)


