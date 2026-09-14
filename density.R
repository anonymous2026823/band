# ==============================================================================
# File: density.R
# Description: Universal routines for evaluating joint probability densities 
#              and log-likelihoods over a directed acyclic graph (DAG).
# ==============================================================================

#' Estimate Joint Log-Likelihood using Pre-Trained Models (Universal)
#' 
#' @param x_target A numeric vector representing a single multivariate observation.
#' @param dag_structure A list containing the network structure, including the topological ordering (`b_star`).
#' @param trained_models A list of trained conditional density models corresponding to the nodes in the DAG.
#' @param predict_func The specific prediction wrapper function to use (e.g., `predict_histogram_wrapper`).
#' @param epsilon A numeric smoothing parameter to prevent log(0) evaluations. Default is 0.00001.
#' @param ... Additional arguments passed to the `predict_func`.
#' 
#' @return A single numeric value representing the joint log-likelihood of the observation.
estimate_g_loglik <- function(x_target, 
                              dag_structure, 
                              trained_models, 
                              predict_func, 
                              epsilon = 0.05, 
                              ...) {
  b_star <- dag_structure$b_star
  p <- length(b_star)
  
  # Initialize with 0 for log-summation
  joint_log_likelihood <- 0
  
  for (j in 1:p) {
    b_j <- b_star[j]
    
    # The value of the specific variable we are evaluating
    x_val <- x_target[b_j]
    
    # Call the predict wrapper using the pre-trained model and pass epsilon explicitly
    conditional_density <- predict_func(
      model_j = trained_models[[j]], 
      x_current_state = x_target, 
      x_val = x_val, 
      epsilon = epsilon,
      ...
    )
    
    # Accumulate the log-likelihood instead of multiplying
    joint_log_likelihood <- joint_log_likelihood + log(conditional_density)
  }
  
  return(joint_log_likelihood)
}
