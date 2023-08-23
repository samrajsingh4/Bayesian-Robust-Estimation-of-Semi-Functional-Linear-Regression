# Define function to add values
accumulate_scores <- function(accumulator, name, value) {
  accumulator[name] <- accumulator[name] + value
  return(accumulator)
}

# Define function to add values to df
accumulate_data_frame <- function(accumulator, name, df) {
  if (is.null(accumulator[[name]])) {
    accumulator[[name]] <- df
  } else {
    accumulator[[name]] <- accumulator[[name]] + df
  }
  return(accumulator)
}

# Define function to add values for list input
accumulate_list <- function(accumulator, name, value) {
  if (is.null(accumulator[[name]])) {
    accumulator[[name]] <- value
  } else {
    accumulator[[name]] <- rbind(accumulator[[name]], value)
  }
  return(accumulator)
}

# Define function to round values
round_scores <- function(scores, digits) {
  sapply(scores, round, digits = digits)
}

# Define function, used to make df for results for criteria
add_to_results <- function(results, error_type, scores_list) {
  
  # Set model names
  model_names <- names(scores_list$DIC)
  
  # Create a new data frame to store the results
  new_df <- data.frame(Error = rep(error_type, length(model_names)),
                       Model = model_names,
                       DIC = scores_list$DIC,
                       EAIC = scores_list$EAIC,
                       EBIC = scores_list$EBIC,
                       LPML = scores_list$LPML,
                       WAIC1 = scores_list$WAIC1)
  
  # Bind results
  rbind(results, new_df)
}

# Define function, used to make df for results for bias and mse for sim 1
add_to_results_b <- function(results, error_type, scores_list) {
  
  # Set model names
  model_names <- names(scores_list$MSE_b)
  
  # Create a new data frame to store the results
  new_df <- data.frame(Error = rep(error_type, length(model_names)),
                       Model = model_names,
                       MSE_Beta = scores_list$MSE_b,
                       Bias_Beta = scores_list$Bias_b,
                       MSE_Alpha1 = scores_list$MSE_a1,
                       Bias_Alpha1 = scores_list$Bias_a1)
  
  # Bind results
  rbind(results, new_df)
}

# Define function, used to make df for results for bias and mse for sim 2
add_to_results_c <- function(results, error_type, scores_list) {
  
  # Set model names
  model_names <- names(scores_list$MSE_b)
  
  # Create a new data frame to store the results
  new_df <- data.frame(Error = rep(error_type, length(model_names)),
                       Model = model_names,
                       MSE_Beta = scores_list$MSE_b,
                       Bias_Beta = scores_list$Bias_b,
                       MSE_Alpha1 = scores_list$MSE_a1,
                       Bias_Alpha1 = scores_list$Bias_a1,
                       MSE_Alpha2 = scores_list$MSE_a2,
                       Bias_Alpha2 = scores_list$Bias_a2)
  
  # Bind results
  rbind(results, new_df)
}

# Define function to perform FPCA
perform_fpca <- function(basis_expansion, time_grid, num_basis=240, num_scores=9) {
  
  # Transpose basis expansion
  transposed_expansion <- t(basis_expansion)
  
  # Create a B-spline basis
  bspline_basis <- create.bspline.basis(rangeval=range(time_grid), nbasis=num_basis)
  
  # Generate functional data object
  fd_objects <- Data2fd(time_grid, transposed_expansion, bspline_basis)
  
  # Perform FPCA
  pca_results <- pca.fd(fd_objects, nharm=num_basis)
  
  # Evaluate the eigenfunctions
  eigenfunctions_eval <- eval.fd(time_grid, pca_results$harmonics)
  
  # Calculate scores
  scores <- inprod(fd_objects, pca_results$harmonics)
  
  # Keep specified scores
  scores <- scores[, 1:num_scores]
  
  # Get proportion of variance explained by each eigenfunction
  variance_explained <- pca_results$varprop  
  
  # Return list of values
  list(eigenfunctions_eval = eigenfunctions_eval, scores = scores, variance_explained = variance_explained)
}

# Define function that combines MCMC models
combined_MCMC_function <- function(distribution = MCMC_function_Normal, a_0, b_0, omega_0, theta_0, niter, burn_in = 10000, y = y, xi = xi, c = NULL, d = NULL) {
  if(is.null(c) && is.null(d)){
    result <- distribution(a_0, b_0, omega_0, theta_0, niter, burn_in, y, xi)
  } else {
    result <- distribution(a_0, b_0, c, d, omega_0, theta_0, niter, burn_in, y, xi)
  }
  return(result)
}
