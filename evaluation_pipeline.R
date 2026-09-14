# ==============================================================================
# File: evaluation_pipeline.R
# Description: Contains the hyperparameter tuning pipelines for the baseline and 
#              enhanced DAG models, as well as the evaluation wrappers for all 
#              baseline competitors (mclust, rvinecopulib, and Oracle).
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Tuning Pipelines for BAND Estimators
# ------------------------------------------------------------------------------

#' Tune and Evaluate Baseline DAG
#' 
#' @param train_data The training data split.
#' @param val_data The validation data split.
#' @param x_test The unseen test data.
#' @param initial_dag_path The baseline DAG structure.
#' @param bin_factors A numeric vector of bin_factors to test.
#' @param epsilons A numeric vector of epsilon values to test.
#' 
#' @return A list containing the best parameters and the test log-likelihood.
tune_baseline_dag <- function(train_data, 
                              val_data, 
                              x_test, 
                              initial_dag_path, 
                              bin_factors, 
                              epsilons = c(0.05)) {
  n_val <- nrow(val_data)
  best_val_ll <- -Inf
  best_bf <- NULL
  best_eps <- NULL
  best_models <- NULL
  
  # Create a search space for the baseline parameters
  search_space <- expand.grid(bf = bin_factors, eps = epsilons)
  
  for (i in 1:nrow(search_space)) {
    current_bf <- search_space$bf[i]
    current_eps <- search_space$eps[i]
    
    baseline_breaks <- calculate_global_breaks(
      train_data,
      bins_per_dim = cube_root(nrow(train_data), bin_factor = current_bf)
    )
    
    baseline_models <- pretrain_all_cdes(
      data = train_data,
      dag_structure = initial_dag_path,
      train_func = train_histogram_wrapper,
      global_breaks = baseline_breaks
    )
    
    val_ll_scores <- numeric(n_val)
    for (k in 1:n_val) {
      val_ll_scores[k] <- estimate_g_loglik(
        x_target = as.numeric(val_data[k, ]),
        dag_structure = initial_dag_path,
        trained_models = baseline_models,
        predict_func = predict_histogram_wrapper,
        epsilon = current_eps # Pass current epsilon to validation
      )
    }
    
    mean_val_ll <- mean(val_ll_scores)
    
    if (mean_val_ll > best_val_ll) {
      best_val_ll <- mean_val_ll
      best_bf <- current_bf
      best_eps <- current_eps
      best_models <- baseline_models
    }
  }
  
  # Final Evaluation
  n_test <- nrow(x_test)
  test_ll_scores <- numeric(n_test)
  for (m in 1:n_test) {
    test_ll_scores[m] <- estimate_g_loglik(
      x_target = as.numeric(x_test[m, ]),
      dag_structure = initial_dag_path,
      trained_models = best_models,
      predict_func = predict_histogram_wrapper,
      epsilon = best_eps # Pass best epsilon to test evaluation
    )
  }
  
  return(list(
    best_bin_factor = best_bf,
    best_epsilon = best_eps,
    test_log_likelihood = mean(test_ll_scores)
  ))
}


#' Tune and Evaluate Enhanced DAG
#' 
#' @param train_data The training data split.
#' @param val_data The validation data split.
#' @param x_test The unseen test data.
#' @param initial_dag_path The baseline DAG structure.
#' @param search_space A dataframe of hyperparameter combinations (must include epsilon).
#' 
#' @return A list containing the best parameters, optimized DAG, and test log-likelihood.
tune_enhanced_dag <- function(train_data, 
                              val_data, 
                              x_test, 
                              initial_dag_path, 
                              search_space) {
  
  p <- length(initial_dag_path$b_star)
  n_val <- nrow(val_data)
  
  best_val_ll <- -Inf
  best_params <- NULL
  best_dag <- NULL
  best_models <- NULL
  
  for (i in 1:nrow(search_space)) {
    current_params <- search_space[i, ]
    current_dag <- initial_dag_path
    
    # Extract the epsilon value for this specific grid iteration
    current_epsilon <- current_params$epsilon
    
    current_breaks <- calculate_global_breaks(
      train_data,
      bins_per_dim = cube_root(nrow(train_data), bin_factor = current_params$bin_factor)
    )
    
    # Refinement step
    for (j in 1:p) {
      b_j <- current_dag$b_star[j]
      history_nodes <- if (j == 1) integer(0) else current_dag$b_star[1:(j - 1)]
      
      score_func <- make_H_hat_calculator(
        b_j = b_j,
        history_nodes = history_nodes,
        data = train_data,
        global_breaks = current_breaks,
        epsilon = current_epsilon # Use grid epsilon here
      )
      
      current_dag$S_sets[[j]] <- refine_parent_sets(
        j = j,
        b_star = current_dag$b_star,
        S_hat_j_minus_1 = current_dag$S_sets[[j]],
        n = nrow(train_data),
        w0 = current_params$w0,
        delta_prime = 0.49,
        iota = current_params$iota,
        gamma = current_params$gamma,
        calculate_H_hat = score_func
      )
    }
    
    enhanced_models <- pretrain_all_cdes(
      data = train_data,
      dag_structure = current_dag,
      train_func = train_histogram_wrapper,
      global_breaks = current_breaks
    )
    
    val_ll_scores <- numeric(n_val)
    for (k in 1:n_val) {
      val_ll_scores[k] <- estimate_g_loglik(
        x_target = as.numeric(val_data[k, ]),
        dag_structure = current_dag,
        trained_models = enhanced_models,
        predict_func = predict_histogram_wrapper,
        epsilon = current_epsilon # Use grid epsilon here
      )
    }
    
    mean_val_ll <- mean(val_ll_scores)
    
    if (mean_val_ll > best_val_ll) {
      best_val_ll <- mean_val_ll
      best_params <- current_params
      best_dag <- current_dag
      best_models <- enhanced_models
    }
  }
  
  # Final Evaluation
  n_test <- nrow(x_test)
  test_ll_scores <- numeric(n_test)
  for (m in 1:n_test) {
    test_ll_scores[m] <- estimate_g_loglik(
      x_target = as.numeric(x_test[m, ]),
      dag_structure = best_dag,
      trained_models = best_models,
      predict_func = predict_histogram_wrapper,
      epsilon = best_params$epsilon # Pass the finalized best epsilon to the test evaluation
    )
  }
  
  return(list(
    best_params = best_params,
    optimized_dag = best_dag,
    test_log_likelihood = mean(test_ll_scores)
  ))
}


# ------------------------------------------------------------------------------
# 2. Baseline Evaluators (Competitors and Oracle)
# ------------------------------------------------------------------------------

#' Evaluate mclust Model Baseline
#' 
#' @param train_data The training data split.
#' @param x_test The unseen test data.
#' @param lower_bound The minimum allowed log-likelihood per observation.
#' 
#' @return The test log-likelihood for the mclust model.
evaluate_mclust_model <- function(train_data, x_test, lower_bound = -10000) {
  suppressPackageStartupMessages(require(mclust))
  suppressPackageStartupMessages(require(caret)) 
  
  # 1. Remove constant or zero-variance columns that cause singular matrices
  train_df <- as.data.frame(train_data)
  keep_cols <- sapply(train_df, function(col) sd(col, na.rm = TRUE) > 1e-7)
  
  train_data_clean <- as.matrix(train_df[, keep_cols, drop = FALSE])
  x_test_clean <- as.matrix(as.data.frame(x_test)[, keep_cols, drop = FALSE])
  
  # 2. Wrap Mclust in tryCatch to gracefully handle any residual matrix singularities
  mclust_fit <- tryCatch({
    Mclust(train_data_clean, verbose = FALSE)
  }, error = function(e) {
    NULL
  })
  
  # If mclust fails entirely due to matrix issues, return the lower bound penalty
  if (is.null(mclust_fit)) {
    return(lower_bound)
  }
  
  test_log_densities <- dens(
    data = x_test_clean,
    modelName = mclust_fit$modelName,
    parameters = mclust_fit$parameters,
    logarithm = TRUE
  )
  
  test_log_densities[is.na(test_log_densities) | 
                       is.nan(test_log_densities) | 
                       is.infinite(test_log_densities)] <- lower_bound
  
  bounded_densities <- pmax(test_log_densities, lower_bound)
  
  return(mean(bounded_densities))
}


#' Evaluate rvinecopulib Model Baseline
#' 
#' @param train_data The training data split.
#' @param x_test The unseen test data.
#' @param max_depth The maximum truncation level for the vine copula. Default is 2.
#' 
#' @return The test log-likelihood for the rvinecopulib model.
evaluate_rvinecopulib_model <- function(train_data, x_test, max_depth = 2) {
  suppressPackageStartupMessages(require(rvinecopulib))
  
  cop_controls <- list(
    family_set = "parametric", 
    trunc_lvl = max_depth
  )
  
  # Train the model
  vine_model <- vine(as.matrix(train_data), copula_controls = cop_controls)
  
  # 1. Calculate raw density for the test data
  test_densities <- dvine(as.matrix(x_test), vine_model)
  
  # 2. Take the log FIRST. Any underflowed 0 will safely become -Inf
  log_densities <- log(test_densities)
  
  # 3. Punish catastrophic failures in log-space
  # If the model density underflowed to 0, its true log-likelihood is worse 
  # than -708. We assign a severe penalty so it cannot artificially beat working models.
  log_densities[is.infinite(log_densities) | is.nan(log_densities)] <- -10000
  
  # 4. Calculate the honest mean log-likelihood
  mean_test_loglik <- mean(log_densities)
  
  return(mean_test_loglik)
}


#' Evaluate Oracle Log-Likelihood Baseline
#' 
#' @param x_test The unseen test data generated by your DGP function.
#' @param active_dgp The data generating process flag ("GMM" or "LINEAR").
#' @param p Number of variables.
#' @param r_dist The distance between cluster means (required for GMM).
#' @param mix_prob The mixing probability of the clusters (required for GMM).
#' 
#' @return The exact average true log-likelihood per data point.
evaluate_oracle_baseline <- function(x_test, active_dgp, p, 
                                     r_dist = NULL, 
                                     mix_prob = 0.5) {
  
  n <- nrow(x_test)
  x_matrix <- as.matrix(x_test)
  
  if (active_dgp == "GMM") {
    if (is.null(r_dist)) stop("r_dist must be provided for GMM oracle.")
    
    total_loglik <- 0
    
    for (i in 1:n) {
      x_row <- x_matrix[i, ]
      
      # Log-density of Component 1 
      log_dens1 <- sum(dnorm(x_row, mean = 0, sd = 1, log = TRUE))
      L1 <- log(mix_prob) + log_dens1
      
      # Log-density of Component 2 
      log_dens2 <- sum(dnorm(x_row, mean = r_dist, sd = 1, log = TRUE))
      L2 <- log(1 - mix_prob) + log_dens2
      
      # Calculate combined log density safely
      max_L <- max(L1, L2)
      point_loglik <- max_L + log(exp(L1 - max_L) + exp(L2 - max_L))
      
      total_loglik <- total_loglik + point_loglik
    }
    
    return(total_loglik / n)
    
  } else if (active_dgp == "LINEAR") {
    
    # Extract the true adjacency matrix from the data attributes
    adj_matrix <- attr(x_test, "true_adjacency_matrix")
    
    if (is.null(adj_matrix)) {
      stop("The x_test object is missing the 'true_adjacency_matrix' attribute.")
    }
    
    suppressPackageStartupMessages(require(mvtnorm))
    
    # Calculate the true covariance matrix
    # Based on the data generation: Data = Data %*% Adjacency + Noise
    I_minus_B <- diag(p) - adj_matrix
    W <- solve(I_minus_B)
    
    # The true covariance is the transpose of W multiplied by W
    true_cov_matrix <- t(W) %*% W
    
    # Calculate exact log-likelihood for the test data
    loglik_vals <- dmvnorm(x_matrix, mean = rep(0, p), sigma = true_cov_matrix, log = TRUE)
    
    return(mean(loglik_vals))
    
  } else {
    stop("Invalid active_dgp provided.")
  }
}



#' 
#' 
#' 
#' 
#' #' Tune and Evaluate Baseline DAG
#' #' 
#' #' @param train_data The training data split.
#' #' @param val_data The validation data split.
#' #' @param x_test The unseen test data.
#' #' @param initial_dag_path The baseline DAG structure.
#' #' @param bin_factors A numeric vector of bin_factors to test.
#' #' @param epsilons A numeric vector of epsilon values to test.
#' #' 
#' #' @return A list containing the best parameters and the test log-likelihood.
#' tune_baseline_dag <- function(train_data, 
#'                               val_data, 
#'                               x_test, 
#'                               initial_dag_path, 
#'                               bin_factors, 
#'                               epsilons = c(0.05)) {
#'   n_val <- nrow(val_data)
#'   best_val_ll <- -Inf
#'   best_bf <- NULL
#'   best_eps <- NULL
#'   best_models <- NULL
#'   
#'   # Create a search space for the baseline parameters
#'   search_space <- expand.grid(bf = bin_factors, eps = epsilons)
#'   
#'   for (i in 1:nrow(search_space)) {
#'     current_bf <- search_space$bf[i]
#'     current_eps <- search_space$eps[i]
#'     
#'     baseline_breaks <- calculate_global_breaks(
#'       train_data,
#'       bins_per_dim = cube_root(nrow(train_data), bin_factor = current_bf)
#'     )
#'     
#'     baseline_models <- pretrain_all_cdes(
#'       data = train_data,
#'       dag_structure = initial_dag_path,
#'       train_func = train_histogram_wrapper,
#'       global_breaks = baseline_breaks
#'     )
#'     
#'     val_ll_scores <- numeric(n_val)
#'     for (k in 1:n_val) {
#'       val_ll_scores[k] <- estimate_g_loglik(
#'         x_target = as.numeric(val_data[k, ]),
#'         dag_structure = initial_dag_path,
#'         trained_models = baseline_models,
#'         predict_func = predict_histogram_wrapper,
#'         epsilon = current_eps # Pass current epsilon to validation
#'       )
#'     }
#'     
#'     mean_val_ll <- mean(val_ll_scores)
#'     
#'     if (mean_val_ll > best_val_ll) {
#'       best_val_ll <- mean_val_ll
#'       best_bf <- current_bf
#'       best_eps <- current_eps
#'       best_models <- baseline_models
#'     }
#'   }
#'   
#'   # cat("Best Baseline Parameters - bin_factor:", best_bf, "| epsilon:", best_eps, "\n")
#'   
#'   # Final Evaluation
#'   n_test <- nrow(x_test)
#'   test_ll_scores <- numeric(n_test)
#'   for (m in 1:n_test) {
#'     test_ll_scores[m] <- estimate_g_loglik(
#'       x_target = as.numeric(x_test[m, ]),
#'       dag_structure = initial_dag_path,
#'       trained_models = best_models,
#'       predict_func = predict_histogram_wrapper,
#'       epsilon = best_eps # Pass best epsilon to test evaluation
#'     )
#'   }
#'   
#'   return(list(
#'     best_bin_factor = best_bf,
#'     best_epsilon = best_eps,
#'     test_log_likelihood = mean(test_ll_scores)
#'   ))
#' }
#' 
#' 
#' 
#' 
#' #' Tune and Evaluate Enhanced DAG
#' #' 
#' #' @param train_data The training data split.
#' #' @param val_data The validation data split.
#' #' @param x_test The unseen test data.
#' #' @param initial_dag_path The baseline DAG structure.
#' #' @param search_space A dataframe of hyperparameter combinations (must include epsilon).
#' #' 
#' #' @return A list containing the best parameters, optimized DAG, and test log-likelihood.
#' tune_enhanced_dag <- function(train_data, 
#'                               val_data, 
#'                               x_test, 
#'                               initial_dag_path, 
#'                               search_space) {
#'   
#'   p <- length(initial_dag_path$b_star)
#'   n_val <- nrow(val_data)
#'   
#'   best_val_ll <- -Inf
#'   best_params <- NULL
#'   best_dag <- NULL
#'   best_models <- NULL
#'   
#'   # cat("Starting Enhanced DAG Grid Search...\n")
#'   
#'   for (i in 1:nrow(search_space)) {
#'     current_params <- search_space[i, ]
#'     current_dag <- initial_dag_path
#'     
#'     # Extract the epsilon value for this specific grid iteration
#'     current_epsilon <- current_params$epsilon
#'     
#'     current_breaks <- calculate_global_breaks(
#'       train_data,
#'       bins_per_dim = cube_root(nrow(train_data), bin_factor = current_params$bin_factor)
#'     )
#'     
#'     
#'       # Refinement step
#'     for (j in 1:p) {
#'       b_j <- current_dag$b_star[j]
#'       history_nodes <- if (j == 1) integer(0) else current_dag$b_star[1:(j - 1)]
#'       
#'       score_func <- make_H_hat_calculator(
#'         b_j = b_j,
#'         history_nodes = history_nodes,
#'         data = train_data,
#'         global_breaks = current_breaks,
#'         epsilon = current_epsilon # Use grid epsilon here
#'       )
#'       
#'       current_dag$S_sets[[j]] <- refine_parent_sets(
#'         j = j,
#'         b_star = current_dag$b_star,
#'         S_hat_j_minus_1 = current_dag$S_sets[[j]],
#'         n = nrow(train_data),
#'         w0 = current_params$w0,
#'         delta_prime = 0.49,
#'         iota = current_params$iota,
#'         gamma = current_params$gamma,
#'         calculate_H_hat = score_func
#'       )
#'     }
#'     
#'     enhanced_models <- pretrain_all_cdes(
#'       data = train_data,
#'       dag_structure = current_dag,
#'       train_func = train_histogram_wrapper,
#'       global_breaks = current_breaks
#'     )
#'     
#'     val_ll_scores <- numeric(n_val)
#'     for (k in 1:n_val) {
#'       val_ll_scores[k] <- estimate_g_loglik(
#'         x_target = as.numeric(val_data[k, ]),
#'         dag_structure = current_dag,
#'         trained_models = enhanced_models,
#'         predict_func = predict_histogram_wrapper,
#'         epsilon = current_epsilon # Use grid epsilon here
#'       )
#'     }
#'     
#'     mean_val_ll <- mean(val_ll_scores)
#'     
#'     if (mean_val_ll > best_val_ll) {
#'       best_val_ll <- mean_val_ll
#'       best_params <- current_params
#'       best_dag <- current_dag
#'       best_models <- enhanced_models
#'     }
#'   }
#'   
#'   # cat("Best Enhanced DAG Parameters Found.\n")
#'   
#'   # Final Evaluation
#'   n_test <- nrow(x_test)
#'   test_ll_scores <- numeric(n_test)
#'   for (m in 1:n_test) {
#'     test_ll_scores[m] <- estimate_g_loglik(
#'       x_target = as.numeric(x_test[m, ]),
#'       dag_structure = best_dag,
#'       trained_models = best_models,
#'       predict_func = predict_histogram_wrapper,
#'       epsilon = best_params$epsilon # Pass the finalized best epsilon to the test evaluation
#'     )
#'   }
#'   
#'   return(list(
#'     best_params = best_params,
#'     optimized_dag = best_dag,
#'     test_log_likelihood = mean(test_ll_scores)
#'   ))
#' }
#' 
#' 
#' evaluate_mclust_model <- function(train_data, x_test, lower_bound = -10000) {
#'   suppressPackageStartupMessages(require(mclust))
#'   suppressPackageStartupMessages(require(caret)) # Optional, or use base R below
#'   
#'   # 1. Remove constant or zero-variance columns that cause singular matrices
#'   train_df <- as.data.frame(train_data)
#'   keep_cols <- sapply(train_df, function(col) sd(col, na.rm = TRUE) > 1e-7)
#'   
#'   train_data_clean <- as.matrix(train_df[, keep_cols, drop = FALSE])
#'   x_test_clean <- as.matrix(as.data.frame(x_test)[, keep_cols, drop = FALSE])
#'   
#'   # 2. Wrap Mclust in tryCatch to gracefully handle any residual matrix singularities
#'   mclust_fit <- tryCatch({
#'     Mclust(train_data_clean, verbose = FALSE)
#'   }, error = function(e) {
#'     NULL
#'   })
#'   
#'   # If mclust fails entirely due to matrix issues, return the lower bound penalty
#'   if (is.null(mclust_fit)) {
#'     return(lower_bound)
#'   }
#'   
#'   test_log_densities <- dens(
#'     data = x_test_clean,
#'     modelName = mclust_fit$modelName,
#'     parameters = mclust_fit$parameters,
#'     logarithm = TRUE
#'   )
#'   
#'   test_log_densities[is.na(test_log_densities) | 
#'                        is.nan(test_log_densities) | 
#'                        is.infinite(test_log_densities)] <- lower_bound
#'   
#'   bounded_densities <- pmax(test_log_densities, lower_bound)
#'   
#'   return(mean(bounded_densities))
#' }
#' 
#' 
#' #' #' Evaluate mclust 
#' #' #' 
#' #' #' @param train_data The training data split.
#' #' #' @param x_test The unseen test data.
#' #' #' @param lower_bound The minimum allowed log-likelihood per observation.
#' #' #' 
#' #' #' @return The test log-likelihood for the mclust model.
#' #' evaluate_mclust_model <- function(train_data, x_test, lower_bound = -10000) {
#' #'   # cat("Training and Evaluating mclust...\n")
#' #'   suppressPackageStartupMessages(require(mclust))
#' #'   
#' #'   mclust_model <- Mclust(train_data, verbose = FALSE)
#' #'   
#' #'   test_log_densities <- dens(
#' #'     data = x_test,
#' #'     modelName = mclust_model$modelName,
#' #'     parameters = mclust_model$parameters,
#' #'     logarithm = TRUE
#' #'   )
#' #'   
#' #'   # 1. Catch any NA, NaN, or strict -Inf values and replace them
#' #'   test_log_densities[is.na(test_log_densities) | 
#' #'                        is.nan(test_log_densities) | 
#' #'                        is.infinite(test_log_densities)] <- lower_bound
#' #'   
#' #'   # 2. Cap any extremely large negative numbers at the lower bound
#' #'   bounded_densities <- pmax(test_log_densities, lower_bound)
#' #'   
#' #'   # 3. Return the safe mean log-likelihood
#' #'   return(mean(bounded_densities))
#' #' }
#' #' 
#' 
#' 
#' 
#' 
#' 
#' evaluate_rvinecopulib_model <- function(train_data, x_test, max_depth = 2) {
#'   suppressPackageStartupMessages(require(rvinecopulib))
#'   
#'   cop_controls <- list(
#'     family_set = "parametric", 
#'     trunc_lvl = max_depth
#'   )
#'   
#'   # Train the model
#'   vine_model <- vine(as.matrix(train_data), copula_controls = cop_controls)
#'   
#'   # 1. Calculate raw density for the test data
#'   test_densities <- dvine(as.matrix(x_test), vine_model)
#'   
#'   # 2. Take the log FIRST. Any underflowed 0 will safely become -Inf
#'   log_densities <- log(test_densities)
#'   
#'   # 3. Punish catastrophic failures in log-space
#'   # If the model density underflowed to 0, its true log-likelihood is worse 
#'   # than -708. We assign a severe penalty so it cannot artificially beat working models.
#'   log_densities[is.infinite(log_densities) | is.nan(log_densities)] <- -10000
#'   
#'   # 4. Calculate the honest mean log-likelihood
#'   mean_test_loglik <- mean(log_densities)
#'   
#'   return(mean_test_loglik)
#' }
#' 
#' 
#' #' Evaluate Oracle Log-Likelihood Baseline
#' #' 
#' #' @param x_test The unseen test data generated by your DGP function.
#' #' @param active_dgp The data generating process flag ("GMM" or "LINEAR").
#' #' @param p Number of variables.
#' #' @param r_dist The distance between cluster means (required for GMM).
#' #' @param mix_prob The mixing probability of the clusters (required for GMM).
#' #' @return The exact average true log-likelihood per data point.
#' evaluate_oracle_baseline <- function(x_test, active_dgp, p, 
#'                                      r_dist = NULL, 
#'                                      mix_prob = 0.5) {
#'   
#'   n <- nrow(x_test)
#'   x_matrix <- as.matrix(x_test)
#'   
#'   if (active_dgp == "GMM") {
#'     if (is.null(r_dist)) stop("r_dist must be provided for GMM oracle.")
#'     
#'     total_loglik <- 0
#'     
#'     for (i in 1:n) {
#'       x_row <- x_matrix[i, ]
#'       
#'       # Log-density of Component 1 
#'       log_dens1 <- sum(dnorm(x_row, mean = 0, sd = 1, log = TRUE))
#'       L1 <- log(mix_prob) + log_dens1
#'       
#'       # Log-density of Component 2 
#'       log_dens2 <- sum(dnorm(x_row, mean = r_dist, sd = 1, log = TRUE))
#'       L2 <- log(1 - mix_prob) + log_dens2
#'       
#'       # Calculate combined log density safely
#'       max_L <- max(L1, L2)
#'       point_loglik <- max_L + log(exp(L1 - max_L) + exp(L2 - max_L))
#'       
#'       total_loglik <- total_loglik + point_loglik
#'     }
#'     
#'     return(total_loglik / n)
#'     
#'   } else if (active_dgp == "LINEAR") {
#'     
#'     # Extract the true adjacency matrix from the data attributes
#'     adj_matrix <- attr(x_test, "true_adjacency_matrix")
#'     
#'     if (is.null(adj_matrix)) {
#'       stop("The x_test object is missing the 'true_adjacency_matrix' attribute.")
#'     }
#'     
#'     suppressPackageStartupMessages(require(mvtnorm))
#'     
#'     # Calculate the true covariance matrix
#'     # Based on the data generation: Data = Data %*% Adjacency + Noise
#'     I_minus_B <- diag(p) - adj_matrix
#'     W <- solve(I_minus_B)
#'     
#'     # The true covariance is the transpose of W multiplied by W
#'     true_cov_matrix <- t(W) %*% W
#'     
#'     # Calculate exact log-likelihood for the test data
#'     loglik_vals <- dmvnorm(x_matrix, mean = rep(0, p), sigma = true_cov_matrix, log = TRUE)
#'     
#'     return(mean(loglik_vals))
#'     
#'   } else {
#'     stop("Invalid active_dgp provided.")
#'   }
#' }
