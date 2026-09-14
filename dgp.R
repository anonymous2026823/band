# ==============================================================================
# File: dgp.R
# Description: Data Generating Process (DGP) functions for synthetic experiments.
#              Includes generators for Gaussian Mixture Models (GMM) and 
#              Linear Gaussian DAGs.
# ==============================================================================

library(MASS)


# ------------------------------------------------------------------------------
# 1. Gaussian Mixture Model (GMM) Data Generator
# ------------------------------------------------------------------------------

#' Generate 2-Component Gaussian Mixture Data with Uncorrelated Features
#' 
#' @param n Sample size.
#' @param p Number of variables.
#' @param r Shift multiplier to separate the component means.
#' @param mix_prob Probability of an observation belonging to the second component.
#' 
#' @return A data frame of generated observations with true component labels attached as an attribute.
generate_gmm_data <- function(n, p, r = 0, mix_prob = 0.5) {
  
  # 1. Create the uncorrelated covariance matrix (Identity matrix)
  Sigma <- diag(p)
  
  # 2. Generate the base data centered at (0, 0, ..., 0)
  mu_base <- rep(0, p)
  data <- MASS::mvrnorm(n, mu = mu_base, Sigma = Sigma)
  
  # 3. Determine component assignments via a binomial draw
  # 0 = Component 1, 1 = Component 2
  is_comp2 <- rbinom(n, size = 1, prob = mix_prob)
  
  # 4. Shift the observations assigned to Component 2 by the vector (r, r, ..., r)
  data[is_comp2 == 1, ] <- data[is_comp2 == 1, ] + r
  
  # 5. Format output
  data_df <- as.data.frame(data)
  colnames(data_df) <- paste0("X", 1:p)
  
  # Attach the true component labels as a hidden attribute
  attr(data_df, "component_labels") <- ifelse(is_comp2 == 0, 1, 2)
  
  return(data_df)
}



# ------------------------------------------------------------------------------
# 2. Linear Gaussian DAG Data Generator
# ------------------------------------------------------------------------------

#' Generate Data from a Standard Linear Gaussian DAG Model with Random Ordering
#' 
#' @param n Sample size.
#' @param p Number of variables.
#' @param prob_edge The probability of an edge existing between any two nodes.
#' 
#' @return A data frame of generated observations with the true adjacency matrix, topological ordering, and parent sets attached as attributes.
generate_linear_dag_data <- function(n, p, prob_edge = 0.2) {
  data <- matrix(0, nrow = n, ncol = p)
  adj_matrix <- matrix(0, nrow = p, ncol = p)
  
  # 1. Randomize the true topological ordering
  true_order <- sample(1:p)
  
  # 2. Generate data sequentially according to the random order
  for (step in 1:p) {
    curr_node <- true_order[step]
    noise <- rnorm(n, mean = 0, sd = 1)
    
    if (step == 1) {
      data[, curr_node] <- noise
    } else {
      potential_parents <- true_order[1:(step - 1)]
      is_parent <- rbinom(length(potential_parents), size = 1, prob = prob_edge)
      actual_parents <- potential_parents[is_parent == 1]
      
      if (length(actual_parents) > 0) {
        signs <- sample(c(-1, 1), size = length(actual_parents), replace = TRUE)
        magnitudes <- runif(length(actual_parents), min = 0.5, max = 1.5)
        weights <- signs * magnitudes
        
        adj_matrix[actual_parents, curr_node] <- weights
        
        parent_data <- data[, actual_parents, drop = FALSE]
        data[, curr_node] <- as.numeric(parent_data %*% weights) + noise
      } else {
        data[, curr_node] <- noise
      }
    }
  }
  
  data_df <- as.data.frame(data)
  colnames(data_df) <- paste0("X", 1:p)
  
  # 3. Automatically build the true_S_sets list
  true_S_sets <- list()
  for (step in 1:p) {
    curr_node <- true_order[step]
    parents <- which(adj_matrix[, curr_node] != 0)
    true_S_sets[[step]] <- as.integer(parents)
  }
  
  # 4. Attach all three components as attributes
  attr(data_df, "true_adjacency_matrix") <- adj_matrix
  attr(data_df, "true_b_star") <- true_order
  attr(data_df, "true_S_sets") <- true_S_sets
  
  return(data_df)
}
