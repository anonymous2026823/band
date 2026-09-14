# ==============================================================================
# File: dag.R
# Description: Routines for generating and estimating Directed Acyclic Graph (DAG)
#              structures. Includes manual autoregressive generation, 
#              data-driven structure learning via Greedy Equivalence Search (GES),
#              and iterative graph refinement.
# ==============================================================================

# Installation Helper (Uncomment to install dependencies)
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(c("graph", "RBGL"))
# install.packages("pcalg")

library(pcalg)


# ------------------------------------------------------------------------------
# 1. Manual Autoregressive DAG Generator
# ------------------------------------------------------------------------------

#' Generate a Manual Autoregressive DAG Structure
#' 
#' @description 
#' Creates a deterministic baseline graph where each node is connected to a fixed 
#' number of immediately preceding nodes in the standard topological order.
#' 
#' @param p Integer representing the total number of dimensions (nodes) in the graph.
#' @param k Integer representing the maximum number of preceding nodes to include as parents.
#' 
#' @return A list containing the topological ordering (`b_star`) and parent sets (`S_sets`).
get_manual_ar_dag <- function(p, k) {
  b_star <- 1:p
  S_sets <- vector("list", p)
  
  for (j in 1:p) {
    # If it is the first node or if k is 0, the parent set is strictly empty
    if (j == 1 || k == 0) {
      S_sets[[j]] <- integer(0)
    } else {
      # The parents are the k nodes immediately before j
      # max(1, ...) ensures we do not try to look at node index 0 or negative
      start_idx <- max(1, j - k)
      S_sets[[j]] <- b_star[start_idx:(j - 1)]
    }
  }
  
  return(list(b_star = b_star, S_sets = S_sets))
}


# ------------------------------------------------------------------------------
# 2. Greedy Equivalence Search (GES) Baseline
# ------------------------------------------------------------------------------

#' Get Initial DAG Structure using Greedy Equivalence Search (GES)
#' 
#' @description 
#' Learns an initial directed acyclic graph structure from continuous data using 
#' the GES algorithm with a Gaussian observational score. A slight numerical 
#' jitter is applied to prevent singular covariance matrices during graph search.
#' 
#' @param data A numeric matrix or dataframe containing the observations.
#' @param max_degree The maximum number of edges allowed for any node during the search. 
#'                   Set to a small integer (e.g., 2 or 3) for significant speedups. Default is 10.
#' 
#' @return A list containing the topological ordering (`b_star`) and parent sets (`S_sets`).
get_dag_structure_ges <- function(data, max_degree = 10) {
  suppressPackageStartupMessages(require(pcalg))
  p <- ncol(data)
  
  # Inject Micro-Jitter for Numerical Stability
  # Adds a tiny amount of noise to prevent the C++ backend from encountering 
  # singular covariance matrices and freezing during graph space search.
  jitter_matrix <- matrix(rnorm(nrow(data) * p, mean = 0, sd = 0.00001),
                          nrow = nrow(data), ncol = p)
  safe_data <- data + jitter_matrix
  
  
  # 1. Define the Gaussian observational score (BIC) using the stabilized data
  score <- new("GaussL0penObsScore", data = safe_data)
  
  # 2. Run GES to get the Essential Graph (CPDAG)
  ges_fit <- ges(score, maxDegree = max_degree)
  
  # 3. Convert CPDAG to a consistent DAG
  cpdag_graph <- as(ges_fit$essgraph, "graphNEL")
  dag_graph <- pdag2dag(cpdag_graph)$graph
  dag_matrix <- as(dag_graph, "matrix")
  
  # 4. Extract Topological Ordering (b_star)
  remaining_nodes <- 1:p
  b_star <- integer(0)
  temp_mat <- dag_matrix
  
  while(length(remaining_nodes) > 0) {
    in_degrees <- colSums(temp_mat[remaining_nodes, remaining_nodes, drop = FALSE])
    roots <- remaining_nodes[in_degrees == 0]
    
    if (length(roots) == 0) {
      roots <- remaining_nodes[1] 
    }
    
    b_star <- c(b_star, roots)
    remaining_nodes <- setdiff(remaining_nodes, roots)
  }
  
  # 5. Extract Parent Sets (S_sets)
  S_sets <- vector("list", p)
  for (j in 1:p) {
    target_node <- b_star[j]
    if (j == 1) {
      S_sets[[j]] <- integer(0)
    } else {
      all_parents <- which(dag_matrix[, target_node] == 1)
      history <- b_star[1:(j - 1)]
      S_sets[[j]] <- intersect(all_parents, history)
    }
  }
  
  return(list(b_star = b_star, S_sets = S_sets))
}


# ------------------------------------------------------------------------------
# 3. Iterative Refinement Algorithm
# ------------------------------------------------------------------------------

#' Iterative Refinement for Parent Selection
#' 
#' @param j The index of the current target node in the topological ordering.
#' @param b_star The full topological ordering vector.
#' @param S_hat_j_minus_1 The initial parent set from the baseline DAG.
#' @param n The sample size.
#' @param w0 Tuning parameter (w0 > 0) to control the penalty strength.
#' @param delta_prime Tuning parameter (Delta' < 0.5) to bound maximum subset size.
#' @param iota Integer threshold (iota > 0) to determine when to reset the parent set.
#' @param gamma Lipschitz constant from Assumption 1 (default is 1).
#' @param calculate_H_hat A function that takes a candidate parent set and returns the empirical H_hat.
#' 
#' @return An integer vector representing the refined parent set for node j.
refine_parent_sets <- function(j, 
                               b_star, 
                               S_hat_j_minus_1, 
                               n, 
                               w0 = 0.001, 
                               delta_prime = 0.49, 
                               iota = 3, 
                               gamma = 1,
                               calculate_H_hat) {
  
  # Base case for the first node in the topological ordering
  if (j == 1) {
    return(integer(0))
  }
  
  # B_{j-1} represents the accumulated history
  B_j_minus_1 <- b_star[1:(j - 1)]
  
  # ==========================================
  # Step 1: Initialization
  # ==========================================
  U <- list()
  if (length(S_hat_j_minus_1) > iota) {
    U[[1]] <- integer(0) # Reset to empty set
  } else {
    U[[1]] <- S_hat_j_minus_1
  }
  
  # ==========================================
  # Step 2: Iterative Construction
  # ==========================================
  s_n <- ceiling((log(n))^delta_prime)
  W <- min(s_n, length(B_j_minus_1))
  
  # Store the empirical H_hat scores to avoid recalculating them in the argmax step
  H_hat_scores <- numeric(W + 1)
  H_hat_scores[1] <- calculate_H_hat(U[[1]])
  
  for (t in 0:(W - 1)) {
    U_t <- U[[t + 1]]
    
    # Candidates are nodes in the history that are not yet in U_t
    candidates <- setdiff(B_j_minus_1, U_t)
    
    # Break early if there are no more history nodes to add
    if (length(candidates) == 0) {
      W <- t
      break
    }
    
    best_v <- NA
    best_H <- -Inf
    
    # Find the single best candidate node to add
    for (v in candidates) {
      candidate_U <- c(U_t, v)
      current_H <- calculate_H_hat(candidate_U)
      
      if (current_H > best_H) {
        best_H <- current_H
        best_v <- v
      }
    }
    
    # Update U_{t+1}
    U[[t + 2]] <- c(U_t, best_v)
    H_hat_scores[t + 2] <- best_H
  }
  
  # ==========================================
  # Step 3: Selection via Penalization
  # ==========================================
  
  # Helper function to compute T_n(s)
  calc_T_n <- function(s) {
    (n^(-1 / (s + 2))) * sqrt(s) * (2 * gamma * log(n))^s
  }
  
  final_scores <- numeric(W + 1)
  for (t in 0:W) {
    current_U <- U[[t + 1]]
    subset_size <- length(current_U)
    
    # The penalty term is w0 * T_n(1 + #U_t)
    penalty <- w0 * calc_T_n(1 + subset_size)
    final_scores[t + 1] <- H_hat_scores[t + 1] - penalty
  }
  
  # Find t_0 (which.max returns 1-based index)
  t_0_index <- which.max(final_scores)
  refined_S_j_minus_1 <- U[[t_0_index]]
  
  return(refined_S_j_minus_1)
}

