# ==============================================================================
# File: compute_H.R
# Description: Calculates the empirical L2 density score (H_hat) used to 
#              evaluate and select optimal parent sets in the DAG.
# ==============================================================================


# ------------------------------------------------------------------------------
# High-Speed Closure Factory for DAG Search
# ------------------------------------------------------------------------------

#' Optimized Factory for Smoothed H_hat Calculator
#' 
#' @description 
#' This factory function pre-calculates the data binning over the entire history 
#' space. It returns a lightweight, highly optimized function that can instantly 
#' compute the L2 score for any subset of parents. This prevents redundant data 
#' processing during the intensive DAG parent search phase.
#' 
#' @param b_j Integer representing the target node index.
#' @param history_nodes Integer vector of all node indices preceding the target.
#' @param data A numeric matrix or data frame of the continuous observations.
#' @param global_breaks A list of pre-computed break points for each variable.
#' @param epsilon A numeric Laplace smoothing parameter. Default is 0.00001.
#' 
#' @return A function (closure) that takes a candidate parent subset (`candidate_U`)
#'         and returns the numeric L2 integral score (H_hat).
make_H_hat_calculator <- function(b_j, 
                                  history_nodes, 
                                  data, 
                                  global_breaks, 
                                  epsilon = 0.05) {
  n <- nrow(data)
  
  # Step 1: Pre-bin and Pre-calculate Widths once for all data
  target_vals <- data[, b_j]
  target_breaks <- global_breaks[[b_j]]
  target_binned <- findInterval(target_vals, target_breaks, all.inside = TRUE)
  
  K <- length(target_breaks) - 1
  lambda_bj <- diff(target_breaks) 
  
  history_binned <- matrix(0, nrow = n, ncol = length(history_nodes))
  if (length(history_nodes) > 0) {
    for (idx in seq_along(history_nodes)) {
      h_node <- history_nodes[idx]
      history_binned[, idx] <- findInterval(data[, h_node], global_breaks[[h_node]], all.inside = TRUE)
    }
  }
  
  hist_col_map <- setNames(seq_along(history_nodes), history_nodes)
  
  # Step 2: Return high-performance closure utilizing vectorized tabulate
  function(candidate_U) {
    if (length(candidate_U) == 0) {
      count_S_plus <- tabulate(target_binned, nbins = K)
      dens_val <- (count_S_plus + epsilon) / ((n + epsilon * K) * lambda_bj)
      integral <- sum((dens_val^2) * lambda_bj)
      return(integral)
    }
    
    selected_cols <- hist_col_map[as.character(candidate_U)]
    parent_binned <- history_binned[, selected_cols, drop = FALSE]
    
    if (ncol(parent_binned) == 1) {
      parent_keys <- parent_binned[, 1]
    } else {
      parent_keys <- as.integer(interaction(as.data.frame(parent_binned), drop = TRUE))
    }
    
    
    M <- max(parent_keys)
    joint_keys <- (parent_keys - 1L) * K + target_binned
    
    joint_counts <- tabulate(joint_keys, nbins = M * K)
    count_mat <- matrix(joint_counts, nrow = K, ncol = M)
    n_S_counts <- tabulate(parent_keys, nbins = M)
    
    active_idx <- which(n_S_counts > 0)
    C_active <- count_mat[, active_idx, drop = FALSE]
    N_S_active <- n_S_counts[active_idx]
    
    num <- C_active + epsilon
    denom <- outer(lambda_bj, N_S_active + epsilon * K)
    
    dens_mat <- num / denom
    # For each distinct set of parent nodes
    integrals <- colSums(dens_mat^2 * lambda_bj) 
    total_H_hat <- sum(integrals * N_S_active) / n
    
    return(total_H_hat)
  }
}
