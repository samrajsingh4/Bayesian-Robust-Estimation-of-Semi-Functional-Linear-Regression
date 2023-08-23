# Assign error types, and models
error_types <- c("Normal", "t", "Cauchy", "CN", "Mixture")
distributions = list('Normal' = MCMC_function_Normal, 'Student-t' = MCMC_function_T, 'CN' = MCMC_function_CN, 'Slash' = MCMC_function_Slash)

# Initialise result variables
final_results <- list()
final_results_b <- list()
results <- list()
MSE_results <- list()
final_results_DIC <- data.frame()
final_results_MSE <- data.frame()
final_results_EAIC <- data.frame()
final_results_EBIC <- data.frame()
final_results_LPML <- data.frame()

# Set sample size
sample_size <- 100
num_obs <- sample_size + 1
num_bases <- 50
time_grid <- seq(0, 1, length=num_obs)

# Set seed for reproducibility
set.seed(123)
for (error_type in error_types) {
  print(paste("Running simulation for error type:", error_type))
  
  # Initialise accumulators for results
  dic_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  eaic_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  ebic_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  LPML_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  waic1_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  MSE_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  MSE_accum_a1 <- setNames(rep(0, length(names(distributions))), names(distributions))
  MSE_accum_a2 <- setNames(rep(0, length(names(distributions))), names(distributions))
  bias_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  bias_accum_a1 <- setNames(rep(0, length(names(distributions))), names(distributions))
  bias_accum_a2 <- setNames(rep(0, length(names(distributions))), names(distributions))
  tauE_mean_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  tauE_sd_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
  beta_mean_accum <- setNames(vector("list", length(names(distributions))), names(distributions))
  beta_sd_accum <- setNames(vector("list", length(names(distributions))), names(distributions))
  
  # Run MCMC functions 50 times
  for (i in 1:50) {
    
    # Simulate data with specified error and assign required values
    simulated_data <- simulate_functional_data_polynomial(sample_size, num_obs, num_bases, time_grid, error_type)
    basis_expansion <- simulated_data$basis_expansion
    Y <- simulated_data$Y
    true_beta <- simulated_data$true_beta
    true_alpha1 <- simulated_data$true_alpha1
    true_alpha2 <- simulated_data$true_alpha2
    fpca_results <- perform_fpca(basis_expansion, time_grid, num_basis = 100, num_scores = 4)
    eigenfunctions_eval <- fpca_results$eigenfunctions_eval
    scores <- fpca_results$scores
    Z1 <- simulated_data$Z1
    Z2 <- simulated_data$Z2
    scores <- cbind(Z1, Z2, scores)
    
    for(name in names(distributions)) {
      if(name == 'Normal'){
        res <- combined_MCMC_function(distribution = distributions[[name]], a_0 = 0.001, b_0 = 0.001, omega_0 = 1000, theta_0 = -5 + (i*0.1), niter = 50000, burn_in = 10000, y = Y, xi = scores)
      } else if(name == 'CN') {
        res <- combined_MCMC_function(distribution = distributions[[name]], a_0 = 0.001, b_0 = 0.001, c = 1, d = 1, omega_0 = 1000, theta_0 = -5 + (i*0.1), niter = 50000, burn_in = 10000, y = Y, xi = scores)
      } else {
        res <- combined_MCMC_function(distribution = distributions[[name]], a_0 = 0.001, b_0 = 0.001, c = 0.001, d = 1, omega_0 = 1000, theta_0 = -5 + (i*0.1), niter = 50000, burn_in = 10000, y = Y, xi = scores)
      }
      
      # Accumulate DIC, EAIC, EBIC, LPML, MSE, and bias results
      dic_scores_accum <- accumulate_scores(dic_scores_accum, name, res$DIC)
      eaic_scores_accum <- accumulate_scores(eaic_scores_accum, name, res$EAIC)
      ebic_scores_accum <- accumulate_scores(ebic_scores_accum, name, res$EBIC)
      LPML_scores_accum <- accumulate_scores(LPML_scores_accum, name, res$LPML)
      waic1_scores_accum <- accumulate_scores(waic1_scores_accum, name, res$WAIC1)
      
      # Calculate and store MSE and Bias values
      average_beta <- eigenfunctions_eval[,1:4] %*% colMeans(res$beta[,3:6])
      MSE_beta <- sum((true_beta - average_beta)^2) / length(true_beta)
      MSE_accum <- accumulate_scores(MSE_accum, name, MSE_beta)
      MSE_alpha1 <- sum((true_alpha1 - res$beta[,1])^2) / length(res$beta[,1])
      MSE_accum_a1 <- accumulate_scores(MSE_accum_a1, name, MSE_alpha1)
      MSE_alpha2 <- sum((true_alpha2 - res$beta[,2])^2) / length(res$beta[,2])
      MSE_accum_a2 <- accumulate_scores(MSE_accum_a2, name, MSE_alpha2)
      bias_beta <- sum(true_beta - average_beta) / length(true_beta)
      bias_accum <- accumulate_scores(bias_accum, name, bias_beta)
      bias_alpha1 <- sum(true_alpha1 - res$beta[,1]) / length(res$beta[,1])
      bias_accum_a1 <- accumulate_scores(bias_accum_a1, name, bias_alpha1)
      bias_alpha2 <- sum(true_alpha2 - res$beta[,2]) / length(res$beta[,2])
      bias_accum_a2 <- accumulate_scores(bias_accum_a2, name, bias_alpha2)
    }
  }
  
  # Average the accumulated results
  avg_scores <- lapply(list(dic_scores_accum, eaic_scores_accum, ebic_scores_accum, LPML_scores_accum, waic1_scores_accum, MSE_accum, bias_accum, MSE_accum_a1, bias_accum_a1, MSE_accum_a2, bias_accum_a2), function(x) x / 50)
  
  # Round values for selection criteria
  rounded_scores <- list(
    DIC = round_scores(avg_scores[[1]], 2),
    EAIC = round_scores(avg_scores[[2]], 2),
    EBIC = round_scores(avg_scores[[3]], 2),
    LPML = round_scores(avg_scores[[4]], 2),
    WAIC1 = round_scores(avg_scores[[5]], 2)
  )
  
  # Round values for Bias and MSE
  rounded_scores_b <- list(
    MSE_b = round_scores(avg_scores[[6]], 5),
    Bias_b = round_scores(avg_scores[[7]], 5),
    MSE_a1 = round_scores(avg_scores[[8]], 5),
    Bias_a1 = round_scores(avg_scores[[9]], 5),
    MSE_a2 = round_scores(avg_scores[[10]], 5),
    Bias_a2 = round_scores(avg_scores[[11]], 5)
  )
  
  # Put values into final df
  final_results <- add_to_results(final_results, error_type, rounded_scores)
  final_results_b <- add_to_results_c(final_results_b, error_type, rounded_scores_b)
}

# Print tables
print(final_results, digits = 2)
print(final_results_b, digits = 5)

