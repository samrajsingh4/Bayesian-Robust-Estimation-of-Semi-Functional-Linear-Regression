#########################SIMULATION 1###########################################
# Define function to generate random coefficients
generate_random_coeffs <- function(sample_size, num_bases) {
  sapply(1:num_bases, function(j) rnorm(sample_size, 0, 3/j))
}

# Define function to generate fouirer basis expansion X(t)
generate_basis_expansion <- function(sample_size, num_obs, random_coeffs, observation_times) {
  X <- matrix(0, sample_size, num_obs)
  for(i in 1:sample_size) {
    for(j in 1:num_obs) {
      X[i, j] <-  sqrt(3) * sum(random_coeffs[i, ] * sin(1:num_bases * pi * observation_times[j]))
    }
  }
  return(X)
}

# Define function to generate true beta
generate_true_beta <- function(num_obs, w, observation_times) {
  sapply(1:num_obs, function(j) {
    sqrt(3) * sum(1/(w^4) * sin(w * pi * observation_times[j]))
  })
}

# Define function to generate mixture of normal errors
generate_mixture_errors <- function(sample_size) {
  error_term <- numeric(sample_size)
  for (i in 1:sample_size) {
    
    # Generate value from 1 to 3 with specified probabilities
    num <- sample(1:3, 1, prob = c(0.70, 0.20, 0.1))
    
    # Generate error from differing normal distributions
    if (num == 1) {
      error_term[i] <- rnorm(1, 0, 1)  
    } else if (num == 2) {
      error_term[i] <- rnorm(1, 0, 10)  
    } else {
      error_term[i] <- rnorm(1, 0, 100) 
    }
  }
  return(error_term)
}

# Define function that generates simulated data
simulate_functional_data <- function(sample_size, num_obs, num_bases, time_grid, error_type) {
  
  # Set observation times between 0 and 1
  observation_times <- seq(0, 1, length=num_obs)
  
  # Generate simulated data values from functions
  random_coeffs <- generate_random_coeffs(sample_size, num_bases)
  basis_expansion <- generate_basis_expansion(sample_size, num_obs, random_coeffs, observation_times)
  true_beta <- generate_true_beta(num_obs, 1:num_bases, observation_times)
  
  # Calculate the mean response
  mean_response <- (observation_times[2] - observation_times[1]) * c(basis_expansion %*% true_beta)
  
  # Add specified error
  if (error_type == "CN") { # Contaminated Normal
    temp <- runif(sample_size, 0, 1)
    error_term <- ifelse(temp < 0.2, rnorm(sample_size, 0, 20), rnorm(sample_size, 0, 1))
  } else if (error_type == "t") { # t-distribution
    error_term <- rt(sample_size, 2)
  } else if (error_type == "Cauchy") { # Cauchy distribution
    error_term <- rcauchy(sample_size, 0, 1)
  } else if (error_type == "Normal") { # Normal distribution
    error_term <- rnorm(sample_size, 0, 1)
  } else if (error_type == "Mixture") { # Mixture of normal distribution
    error_term <- generate_mixture_errors(sample_size) 
  } else { # Invalid Input error
    stop("Type either 'CN', 't', 'Cauchy', or 'Normal'")
  }
  
  # Set Z values
  Z <- runif(sample_size, -1, 1)
  
  # Calculate non-functional component
  non_functional <- 3 * Z
  
  # Calculate response value
  response <- non_functional + mean_response + error_term
  
  # Return list of values
  list(basis_expansion = basis_expansion, Y = response, true_beta = true_beta, true_alpha = 3, Z = Z)
}

#########################SIMULATION 2###########################################
# Define function to generate random coefficients
generate_random_coeffs <- function(sample_size, num_bases) {
  matrix(rnorm(sample_size * num_bases, 0, 1), ncol = num_bases)
}

# Define function to generate polynomial basis expansion X(t)
generate_polynomial_expansion <- function(sample_size, num_obs, random_coeffs, observation_times) {
  X <- matrix(0, sample_size, num_obs)
  for(i in 1:sample_size) {
    for(j in 1:num_obs) {
      X[i, j] <- sum(random_coeffs[i, ] * observation_times[j]^(0:(num_bases-1)))
    }
  }
  return(X)
}

# Define function to generate true beta using polynomial functions
generate_true_beta_polynomial <- function(num_obs, observation_times) {
  sapply(1:num_obs, function(j) {
    sum(1/(1 + j) * observation_times[j]^(0:(num_bases-1)))
  })
}

# Define function that generates simulated data
simulate_functional_data_polynomial <- function(sample_size, num_obs, num_bases, time_grid, error_type) {
  
  # Set observation times between 0 and 1
  observation_times <- seq(0, 1, length=num_obs)
  
  # Generate simulated data values from functions
  random_coeffs <- generate_random_coeffs(sample_size, num_bases)
  basis_expansion <- generate_polynomial_expansion(sample_size, num_obs, random_coeffs, observation_times)
  true_beta <- generate_true_beta_polynomial(num_obs, observation_times)
  
  # Calculate the mean response
  mean_response <- (observation_times[2] - observation_times[1]) * c(basis_expansion %*% true_beta)
  
  # Add specified error
  if (error_type == "CN") { # Contaminated Normal
    temp <- runif(sample_size, 0, 1)
    error_term <- ifelse(temp < 0.1, rnorm(sample_size, 0, 10), rnorm(sample_size, 0, 1))
  } else if (error_type == "t") { # t-distribution
    error_term <- rt(sample_size, 2)
  } else if (error_type == "Cauchy") { # Cauchy distribution
    error_term <- rcauchy(sample_size, 0, 1)
  } else if (error_type == "Normal") { # Normal distribution
    error_term <- rnorm(sample_size, 0, 1)
  } else if (error_type == "Mixture") { # Mixture of normal distribution
    error_term <- generate_mixture_errors(sample_size) 
  } else { # Invalid Input error
    stop("Type either 'CN', 't', 'Cauchy', or 'Normal'")
  }
  
  # Set Z values
  Z1 <- runif(sample_size, 0, 1)
  Z2 <- rbeta(sample_size, 3, 2)
  
  # Calculate non-functional component
  non_functional <- (2 * Z1) + (3*Z2)
  
  # Calculate response value
  response <- non_functional + mean_response + error_term
  
  # Return list of values
  list(basis_expansion = basis_expansion, Y = response, true_beta = true_beta, true_alpha1 = 2 , true_alpha2 = 3, Z1 = Z1, Z2 = Z2)
}
