# rm(list = ls())
library(R.utils)
library(MASS)
library(graph)
library(here)
library(dplyr)
library(knitr)
library(kableExtra)
library(readxl) 

# Activate the project via band.Rproj
source(here("histogram_wrappers.R")) 
source(here("density.R")) 
source(here("dag.R"))
source(here("compute_H.R"))
source(here("evaluation_pipeline.R"))

cube_root <- function(n, bin_factor = 1) {
  ceiling(bin_factor * n**(1/3))
}

# ==========================================
# 1. EXPERIMENT CONFIGURATION
# ==========================================
R_trials <- 20
dataset_list <- c("SONAR", "IONOSPHERE", "PIMA", 
                  "OIL_SPILL", "PHONEME", "HABERMAN", "CONCRETE", "ONLINE_NEWS")
# dataset_list <- c("ONLINE_NEWS") # "CONCRETE", "MAGIC"
total_runs <- length(dataset_list) * R_trials

search_space <- expand.grid(
  w0 = c(0.001), 
  gamma = c(1),
  iota = c(3),
  bin_factor = c(0.5, 1, 2, 3),
  epsilon = c(0.05)
)
baseline_bin_factors <- unique(search_space$bin_factor)

# Initialize master storage arrays
results_ll <- data.frame(
  Dataset = character(total_runs),
  Trial = integer(total_runs),
  Baseline_DAG = numeric(total_runs),
  Enhanced_DAG = numeric(total_runs),
  mclust = numeric(total_runs),
  rvinecopulib = numeric(total_runs),
  stringsAsFactors = FALSE
)

results_time <- results_ll

results_hyper <- data.frame(
  Dataset = character(total_runs),
  Trial = integer(total_runs),
  Base_bin = numeric(total_runs),
  Enh_bin = numeric(total_runs),
  Enh_eps = numeric(total_runs),
  Enh_iota = numeric(total_runs),
  stringsAsFactors = FALSE
)

row_idx <- 1

# ==========================================
# 2. MASTER LOOP OVER DATASETS
# ==========================================
for (d_name in dataset_list) {
  
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("LOADING DATASET: %s\n", d_name))
  cat(sprintf("==================================================\n"))
  
  if (d_name == "SONAR") {
    url <- "https://raw.githubusercontent.com/jbrownlee/Datasets/master/sonar.csv"
    raw_data <- read.csv(url, header = FALSE)
    continuous_data <- raw_data[, 1:60]
    
  } else if (d_name == "IONOSPHERE") {
    url <- "https://raw.githubusercontent.com/jbrownlee/Datasets/master/ionosphere.csv"
    raw_data <- read.csv(url, header = FALSE)
    continuous_data <- raw_data[, c(1, 3:34)]
    
  } else if (d_name == "PIMA") {
    url <- "https://raw.githubusercontent.com/jbrownlee/Datasets/master/pima-indians-diabetes.csv"
    raw_data <- read.csv(url, header = FALSE)
    continuous_data <- raw_data[, 1:8]
    
  } else if (d_name == "OIL_SPILL") {
    url <- "https://raw.githubusercontent.com/jbrownlee/Datasets/master/oil-spill.csv"
    raw_data <- read.csv(url, header = FALSE)
    continuous_data <- raw_data[, 1:49]
    
  } else if (d_name == "PHONEME") {
    url <- "https://raw.githubusercontent.com/jbrownlee/Datasets/master/phoneme.csv"
    raw_data <- read.csv(url, header = FALSE)
    continuous_data <- raw_data[, 1:5]
    
  } else if (d_name == "HABERMAN") {
    url <- "https://raw.githubusercontent.com/jbrownlee/Datasets/master/haberman.csv"
    raw_data <- read.csv(url, header = FALSE)
    continuous_data <- raw_data[, 1:3]
    
  } else if (d_name == "ONLINE_NEWS") {
    local_zip <- "OnlineNewsPopularity.zip"
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00332/OnlineNewsPopularity.zip", 
                  destfile = local_zip, mode = "wb", method = "curl")
    unzip(local_zip)
    raw_data <- read.csv("OnlineNewsPopularity/OnlineNewsPopularity.csv")
    unlink(local_zip)
    unlink("OnlineNewsPopularity", recursive = TRUE)
    continuous_data <- raw_data[, 3:60]
    rm(raw_data)
    gc()
    
  } else if (d_name == "CONCRETE") {
    local_file <- "Concrete_Data.xls"
    download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/concrete/compressive/Concrete_Data.xls", 
                  destfile = local_file, mode = "wb", method = "curl")
    raw_data <- read_excel(local_file)
    unlink(local_file)
    continuous_data <- raw_data[, 1:8]
  }
  
  # Clean and store the full dataset safely
  # 1. Drop missing values
  continuous_data <- na.omit(continuous_data)
  
  # 2. Automatically detect and remove zero-variance columns
  col_vars <- apply(continuous_data, 2, var)
  continuous_data <- continuous_data[, col_vars > 0, drop = FALSE]
  
  # 3. Scale safely
  X_full <- scale(continuous_data)
  N_full <- nrow(X_full)
  p <- ncol(X_full)
  
  cat(sprintf("Dataset %s ready. Total Available Samples: %d | Features: %d\n", d_name, N_full, p))
  
  
  # ==========================================
  # 3. INNER LOOP OVER TRIALS (WITH DETAILED TRACKING)
  # ==========================================
  successful_trials <- 1
  
  while (successful_trials <= R_trials) {
    
    cat(sprintf("\n**************************************************\n"))
    cat(sprintf("%s | SUCCESSFUL TRIAL %d OF %d\n", d_name, successful_trials, R_trials))
    cat(sprintf("**************************************************\n"))
    
    X_real <- X_full
    n_total <- N_full
    
    # 2. Perform the 60/20/20 split on the newly sampled data
    shuffled_idx <- sample(seq_len(n_total))
    
    n_train <- floor(0.60 * n_total)
    n_val   <- floor(0.20 * n_total)
    
    train_idx <- shuffled_idx[1:n_train]
    val_idx   <- shuffled_idx[(n_train + 1):(n_train + n_val)]
    test_idx  <- shuffled_idx[(n_train + n_val + 1):n_total]
    
    train_data <- X_real[train_idx, , drop = FALSE]
    val_data   <- X_real[val_idx, , drop = FALSE]
    x_test     <- X_real[test_idx, , drop = FALSE]
    
    run_failed <- FALSE
    failed_component <- ""
    
    # 1. GES Structure Learning
    initial_dag_path <- tryCatch({
      get_dag_structure_ges(rbind(train_data, val_data))
    }, error = function(e) {
      failed_component <<- "GES structure learning"
      run_failed <<- TRUE
      NULL
    })
    if (run_failed) {
      cat(sprintf("--> FAILED at: %s. Discarding this trial run.\n", failed_component))
      next
    }
   
    # 2. Baseline DAG tuning
    time_baseline <- system.time({
      baseline_results <- tune_baseline_dag(train_data, val_data, x_test, initial_dag_path, 
                                            baseline_bin_factors, epsilons = 0.01)
    })
    
    # 3. Enhanced DAG tuning
    time_enhanced <- system.time({
      enhanced_results <- tune_enhanced_dag(train_data, val_data, x_test, initial_dag_path, 
                                            search_space)
    })
    
    # 4. mclust evaluation
    time_mclust <- system.time({
      mclust_ll <- evaluate_mclust_model(train_data, x_test)
    })
    
    # 5. rvinecopulib evaluation
    time_rvinecopulib <- system.time({
      rvinecopulib_ll <- evaluate_rvinecopulib_model(train_data, x_test)
    })
    
    # If all steps pass successfully:
    results_ll[row_idx, c("Dataset", "Trial")] <- list(d_name, successful_trials)
    results_ll$Baseline_DAG[row_idx] <- baseline_results$test_log_likelihood
    results_ll$Enhanced_DAG[row_idx] <- enhanced_results$test_log_likelihood
    results_ll$mclust[row_idx]        <- mclust_ll
    results_ll$rvinecopulib[row_idx] <- rvinecopulib_ll
    
    results_time[row_idx, c("Dataset", "Trial")] <- list(d_name, successful_trials)
    results_time$Baseline_DAG[row_idx] <- time_baseline["elapsed"]
    results_time$Enhanced_DAG[row_idx] <- time_enhanced["elapsed"]
    results_time$mclust[row_idx]       <- time_mclust["elapsed"]
    results_time$rvinecopulib[row_idx] <- time_rvinecopulib["elapsed"]
    
    results_hyper[row_idx, c("Dataset", "Trial")] <- list(d_name, successful_trials)
    results_hyper$Base_bin[row_idx] <- baseline_results$best_bin_factor
    results_hyper$Enh_bin[row_idx]  <- enhanced_results$best_params$bin_factor
    results_hyper$Enh_eps[row_idx]  <- enhanced_results$best_params$epsilon
    results_hyper$Enh_iota[row_idx] <- enhanced_results$best_params$iota
    
    cat("HELD-OUT TEST LOG-LIKELIHOODS\n")
    cat("Baseline DAG:      ", baseline_results$test_log_likelihood, "\n")
    cat("Enhanced DAG:      ", enhanced_results$test_log_likelihood, "\n")
    cat("mclust:            ", mclust_ll, "\n")
    cat("rvinecopulib:      ", rvinecopulib_ll, "\n")
    
    row_idx <- row_idx + 1
    successful_trials <- successful_trials + 1 
  }
} 

# ==========================================
# SAVE RESULTS TO RDS
# ==========================================
cat("\nSaving results to working directory...\n")
saveRDS(results_ll, file = "experiment_results_ll.rds")
saveRDS(results_time, file = "experiment_results_time.rds")
saveRDS(results_hyper, file = "experiment_results_hyper.rds")
cat("Done. Files saved successfully.\n")



# ==========================================
# AGGREGATE AND PRINT SUMMARY (LATEX)
# ==========================================
cat("\n==================================================\n")
cat("FINAL AVERAGE LOG-LIKELIHOODS (SD)\n")
cat("==================================================\n")

# Calculate means, find the winner, and format row-by-row
summary_table <- results_ll %>%
  group_by(Dataset) %>%
  summarize(
    m_base = mean(Baseline_DAG, na.rm = TRUE), s_base = sd(Baseline_DAG, na.rm = TRUE),
    m_enh  = mean(Enhanced_DAG, na.rm = TRUE), s_enh  = sd(Enhanced_DAG, na.rm = TRUE),
    m_mc   = mean(mclust, na.rm = TRUE),         s_mc   = sd(mclust, na.rm = TRUE),
    m_rv   = mean(rvinecopulib, na.rm = TRUE),   s_rv   = sd(rvinecopulib, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    # Values <= -1000 or NA are treated as unavailable and cannot be winners.
    v_base = ifelse(!is.na(m_base) & m_base > -1000, round(m_base, 1), -Inf),
    v_enh  = ifelse(!is.na(m_enh)  & m_enh  > -1000, round(m_enh, 1), -Inf),
    v_mc   = ifelse(!is.na(m_mc)   & m_mc   > -1000, round(m_mc, 1), -Inf),
    v_rv   = ifelse(!is.na(m_rv)   & m_rv   > -1000, round(m_rv, 1), -Inf),
    
    # Identify the maximum valid value for the current row
    row_max = max(c(v_base, v_enh, v_mc, v_rv), na.rm = TRUE),
    
    # Format values; use -- for NA or values <= -1000
    Baseline_DAG = if(v_base == -Inf) {
      "\\multicolumn{1}{c}{--}"
    } else {
      val <- sprintf("%.1f (%.1f)", m_base, s_base)
      if(v_base == row_max) sprintf("\\textbf{%s}", val) else val
    },
    
    Enhanced_DAG = if(v_enh == -Inf) {
      "\\multicolumn{1}{c}{--}"
    } else {
      val <- sprintf("%.1f (%.1f)", m_enh, s_enh)
      if(v_enh == row_max) sprintf("\\textbf{%s}", val) else val
    },
    
    mclust = if(v_mc == -Inf) {
      "\\multicolumn{1}{c}{--}"
    } else {
      val <- sprintf("%.1f (%.1f)", m_mc, s_mc)
      if(v_mc == row_max) sprintf("\\textbf{%s}", val) else val
    },
    
    rvinecopulib = if(v_rv == -Inf) {
      "\\multicolumn{1}{c}{--}"
    } else {
      val <- sprintf("%.1f (%.1f)", m_rv, s_rv)
      if(v_rv == row_max) sprintf("\\textbf{%s}", val) else val
    }
  ) %>%
  ungroup() %>%
  select(Dataset, Baseline_DAG, Enhanced_DAG, mclust, rvinecopulib) %>%
  mutate(Dataset = gsub("_", "\\\\_", Dataset))


# Generate the tabular core using kable
tabular_core <- kable(
  summary_table,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  linesep = '',
  col.names = c("Dataset", "BAND$_{0}$", "BAND", "mclust", "rvine"),
  align = c("l", "l", "l", "l", "l")
)

# Convert to character and replace the default tabular with tabular[t]
tabular_core <- as.character(tabular_core)
tabular_core <- sub(
  "\\\\begin\\{tabular\\}\\{lllll\\}",
  "\\\\begin{tabular}[t]{lllll}",
  tabular_core
)

# Construct final LaTeX table
latex_out <- paste0(
  "\\begin{table}\n",
  "\\centering\n",
  "\\caption{Mean log-likelihoods for the real datasets over 20 replications. ",
  "Standard deviations are reported in parentheses. ",
  "Bold values indicate the best-performing method for each dataset. ",
  "The symbol `--' indicates severe model failure, with mean log-likelihood ",
  "$\\le -1000$.}\\label{tab:real_data}\n",
  "\\setlength{\\tabcolsep}{1pt}\n",
  "\\centering\n",
  "\\fontsize{7}{8}\\selectfont\n",
  tabular_core,
  "\n\\end{table}\n"
)

cat("\n")
cat(latex_out)