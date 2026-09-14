# ==============================================================================
# File: histogram_wrappers.R
# Description: Core routines for histogram-based conditional density estimation.
#              This module handles global grid calculations, pre-training, and
#              prediction wrappers for the DAG structure.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Grid Calculation
# ------------------------------------------------------------------------------

#' Calculate Global Grid Breaks for the Histogram
#' 
#' @param data A numeric matrix or data frame containing the observations.
#' @param bins_per_dim An integer specifying the number of bins to create along each dimension. Default is 10.
#' 
#' @return A list of length p (number of columns in data), where each element is a numeric vector of break points.
calculate_global_breaks <- function(data, bins_per_dim = 10) {
  p <- ncol(data)
  breaks_list <- vector("list", p)
  
  for (k in 1:p) {
    # Add a slight buffer to the min/max to ensure all data points fall inside
    min_val <- min(data[, k]) - 0.00001
    max_val <- max(data[, k]) + 0.00001
    breaks_list[[k]] <- seq(min_val, max_val, length.out = bins_per_dim + 1)
  }
  return(breaks_list)
}


# ------------------------------------------------------------------------------
# 2. Model Pre-training
# ------------------------------------------------------------------------------

#' Pre-train all Conditional Density Estimators (Universal)
#' 
#' @param data A numeric matrix or data frame containing the training observations.
#' @param dag_structure A list containing the network structure, including the topological ordering (`b_star`) and parent sets (`S_sets`).
#' @param train_func The function used to train individual conditional density estimators (e.g., `train_histogram_wrapper`).
#' @param ... Additional arguments passed to the `train_func`.
#' 
#' @return A list of trained model objects, one for each node in the DAG structure.
pretrain_all_cdes <- function(data, dag_structure, train_func, ...) {
  p <- length(dag_structure$b_star)
  trained_models <- vector("list", p)
  
  for (j in 1:p) {
    b_j <- dag_structure$b_star[j]
    S_j_minus_1 <- dag_structure$S_sets[[j]]
    history_nodes <- if (j == 1) integer(0) else dag_structure$b_star[1:(j - 1)]
    
    # Pass both sparse parents AND full history. The wrapper picks what it wants.
    trained_models[[j]] <- train_func(
      b_j = b_j, 
      S_j_minus_1 = S_j_minus_1, 
      history_nodes = history_nodes,
      data = data, 
      ...
    )
  }
  return(trained_models)
}


# ------------------------------------------------------------------------------
# 3. Histogram Training Wrapper
# ------------------------------------------------------------------------------

#' Train Histogram Wrapper
#' 
#' @param b_j An integer representing the target node index.
#' @param S_j_minus_1 An integer vector containing the indices of the parent nodes.
#' @param history_nodes An integer vector of all nodes preceding the target node in the topological order.
#' @param data A numeric matrix or data frame of the training data.
#' @param global_breaks A list of pre-calculated grid breaks for all dimensions.
#' @param ... Additional arguments (unused in this specific wrapper).
#' 
#' @return A list containing the target index, parent indices, bin assignments, cell lengths, and global breaks required for prediction.
train_histogram_wrapper <- function(b_j, 
                                    S_j_minus_1, 
                                    history_nodes, 
                                    data, 
                                    global_breaks, ...) {
  breaks_bj <- global_breaks[[b_j]]
  lambda_bj_vector <- diff(breaks_bj)
  
  data_bins_bj <- cut(data[, b_j], breaks = breaks_bj, labels = FALSE)
  
  # Extracts sparse parents
  if (length(S_j_minus_1) > 0) {
    data_bins_parents <- lapply(S_j_minus_1, function(p) {
      cut(data[, p], breaks = global_breaks[[p]], labels = FALSE)
    })
  } else {
    data_bins_parents <- NULL
  }
  
  return(list(
    n = nrow(data),
    b_j = b_j,
    S_j_minus_1 = S_j_minus_1, 
    breaks_bj = breaks_bj,
    lambda_bj_vector = lambda_bj_vector,
    data_bins_bj = data_bins_bj,
    data_bins_parents = data_bins_parents,
    global_breaks = global_breaks 
  ))
}


# ------------------------------------------------------------------------------
# 4. Histogram Prediction Wrapper
# ------------------------------------------------------------------------------

#' Predict Conditional Density using Histogram Wrapper
#' 
#' @param model_j The trained model object returned by train_histogram_wrapper.
#' @param x_current_state A numeric vector representing a single test observation.
#' @param x_val A single numeric target value to evaluate the density at.
#' @param epsilon A numeric smoothing parameter. Default is 0.05.
#' @param ... Additional arguments.
#' 
#' @return A single numeric value representing the estimated conditional probability density.
predict_histogram_wrapper <- function(model_j, 
                                      x_current_state, 
                                      x_val, 
                                      epsilon = 0.05, 
                                      ...) {
  n <- model_j$n
  parent_match <- rep(TRUE, n)
  K <- length(model_j$lambda_bj_vector)
  
  # 1. Match the parent state
  if (length(model_j$S_j_minus_1) > 0) {
    for (i in seq_along(model_j$S_j_minus_1)) {
      p_idx <- model_j$S_j_minus_1[i]
      breaks_p <- model_j$global_breaks[[p_idx]]
      
      bin_idx_p <- cut(x_current_state[p_idx], breaks = breaks_p, labels = FALSE, include.lowest = TRUE)
      
      # If the parent value is out of bounds, no training data matches this state
      if (is.na(bin_idx_p)) {
        parent_match <- rep(FALSE, n)
        break
      }
      
      parent_match <- parent_match & (model_j$data_bins_parents[[i]] == bin_idx_p)
    }
  }
  
  count_S <- sum(parent_match, na.rm = TRUE)
  
  # 2. Identify the target bin for the single x_val
  bin_idx_x <- cut(x_val, breaks = model_j$breaks_bj, labels = FALSE, include.lowest = TRUE)
  
  # 3. Handle out-of-bounds (OOB) target values instantly
  if (is.na(bin_idx_x)) {
    return(epsilon / (count_S + epsilon * K))
  }
  
  # 4. Calculate the specific density for this matched bin
  if (count_S == 0) {
    count_S_plus <- 0
  } else {
    match_S_plus <- parent_match & (model_j$data_bins_bj == bin_idx_x)
    count_S_plus <- sum(match_S_plus, na.rm = TRUE)
  }
  
  dens_val <- (count_S_plus + epsilon) / ((count_S + epsilon * K) * model_j$lambda_bj_vector[bin_idx_x])
  
  return(dens_val)
}





