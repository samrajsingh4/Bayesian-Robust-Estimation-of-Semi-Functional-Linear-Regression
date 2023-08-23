# Define Normal Model
MCMC_function_Normal <- function(a_0, b_0, omega_0, theta_0, niter, burn_in = 10000, y = y, xi=xi) {
  
  # Set initial values for MCMC variables
  m <- ncol(xi)
  n <- length(y)
  tauE <- 2
  cov_matrix <- matrix(0, m, m)
  mu_star <- rep(0, m)
  theta_0 <- rep(theta_0, m)  
  omega_0 <- diag(omega_0, m)
  
  # Initialise matrices to hold MCMC results
  tauE_MCMC <- rep(NA, niter)
  beta_MCMC <- matrix(NA, nrow = niter, ncol = m)
  loglike_MCMC <- rep(NA, niter) 
  
  # Initialise progress bar for MCMC iterations
  count <- txtProgressBar(min = 1, max = niter, style = 3)
  
  # Begin MCMC
  for (k in 1:niter) {
    
    #------------------beta-----------------------------------------------------
    # Compute xi cross product and xi weighted with y
    xi_cross <- tauE * tcrossprod(t(xi))
    xi_w <- tauE * (t(xi) %*% y)
    
    # Compute inverse of omega_0
    omega_0_inv <- solve(omega_0)
    
    # Update covariance matrix
    cov_matrix <- solve(xi_cross + omega_0_inv)
    
    # Update mean vector
    mu_star <- cov_matrix %*% (omega_0_inv %*% theta_0 + xi_w)
    
    # Generate new value for beta
    beta <- mvrnorm(1, mu_star, cov_matrix)
    
    #-----------------------tauE------------------------------------------------
    # Generate new value for tauE
    tauE <- rgamma(1, a_0 + n/2, b_0 + 0.5*sum((y - xi %*% beta)^2))
    
    #-------------------------DIC-----------------------------------------------
    # Calculate and store the log-likelihood at each iteration
    y1 = (y - xi %*% beta) * sqrt(tauE)
    loglike_MCMC[k] <- sum(dnorm(y1, 0, 1, log = TRUE) + 0.5 * log(tauE))
    
    # Save the accepted values after each iteration
    tauE_MCMC[k] <- tauE
    beta_MCMC[k,] <- beta  
    
    setTxtProgressBar(count, k, title = "Loading", label = "Welcome")
    
  }
  
  # Set burned-in range
  burned_in <- (burn_in + 1):niter
  
  # Compute y1 using matrix
  y1_matrix <- matrix(y, n, length(burned_in)) - xi %*% t(beta_MCMC[burned_in,])
  y1_matrix <- sweep(y1_matrix, 2, sqrt(tauE_MCMC[burned_in]), "*")
  
  # Compute likelihoods 
  likelihoods_matrix <- dnorm(y1_matrix, 0, 1) * sqrt(tauE_MCMC[burned_in])
  
  # Compute CPO values
  CPO <- apply(likelihoods_matrix, 1, function(likelihoods) {
    1 / mean(1 / likelihoods)
  })
  
  # Calculate LPML
  LPML = sum(log(CPO))
  
  # Compute y2 
  y2 = (y - xi %*% colMeans(beta_MCMC[burned_in,])) * sqrt(mean(tauE_MCMC[burned_in]))
  
  # Calculate the average deviance
  D_theta_avg = -2 * mean(loglike_MCMC[burned_in])
  
  # Compute DIC
  DIC = 2 * D_theta_avg - (-2 * sum(dnorm(y2, 0, 1, log = TRUE) + 0.5 * log(mean(tauE_MCMC[burned_in]))))
  
  # Set number of parameters 
  theta = length(beta) + 1
  
  # Compute EAIC and EBIC
  EAIC = D_theta_avg + 2 * theta
  EBIC = D_theta_avg + theta * log(n)
  
  # Calculate p(z)
  log_lik_matrix <- log(likelihoods_matrix)
  p_z <- sum(apply(log_lik_matrix, 1, mean))
  
  # Compute WAIC1
  WAIC1 = 2 * p_z + 2 * D_theta_avg
  
  # Return the results after burn-in
  return(list(tauE = tauE_MCMC[burned_in], beta = beta_MCMC[burned_in,], DIC = DIC, EAIC = EAIC, EBIC = EBIC, LPML = LPML, WAIC1 = WAIC1))
}

# Define Student-t Model
MCMC_function_T <- function(a_0, b_0, c, d, omega_0, theta_0, niter, burn_in = 10000, y = y, xi=xi) {
  
  # Set initial values for MCMC variables
  m <- ncol(xi)
  n <- length(y)
  tauE <- 2
  u.e <- rep(0.1, n)
  lambda.nu <- 0.3
  nu <- 2
  cov_matrix <- matrix(0, m, m)
  mu_star <- rep(0, m)
  theta_0 <- rep(theta_0, m)  
  omega_0 <- diag(omega_0, m)
  
  # Initialise matrices to hold MCMC results
  tauE_MCMC <- rep(NA, niter)
  beta_MCMC <- matrix(NA, nrow = niter, ncol = m)
  u_MCMC <- matrix(NA, nrow = niter, ncol = n)
  nu_MCMC <- rep(NA, niter)
  lambda_MCMC <- rep(NA, niter)
  loglike_MCMC <- rep(NA, niter) 
  
  # Define function to generate truncated gamma values
  truncated_gamma <- function(n, shape, rate, lower = 0, upper = 1) {
    temp <- qgamma(runif(1,pgamma(lower, shape, rate), pgamma(upper, shape, rate)),  shape, rate)
    return(temp)
  }
  
  # Define function to generate truncated exponential values
  truncated_exponential <- function(n, lambda, lower = 2) {
    temp <- qexp(runif(1,pexp(lower, lambda), 1), lambda)
    return(temp)
  }
  
  # Initialise progress bar for MCMC iterations
  count <- txtProgressBar(min = 1, max = niter, style = 3)
  
  # Begin MCMC
  for (k in 1:niter) {
    
    #------------------beta-----------------------------------------------------
    # Compute y and x multiplied by u.e
    x_u <- xi * sqrt(u.e)
    y_u <- y * sqrt(u.e)
    
    # Compute xi cross product and xi weighted
    xi_cross <- tauE * tcrossprod(t(x_u))
    xi_w <- tauE * (t(x_u) %*% y_u)
    
    # Compute inverse of omega_0
    omega_0_inv <- solve(omega_0)
    
    # Update covariance matrix
    cov_matrix <- solve(xi_cross + omega_0_inv)
    
    # Update mean vector
    mu_star <- cov_matrix %*% (omega_0_inv %*% theta_0 + xi_w)
    
    # Generate new value for beta
    beta <- mvrnorm(1, mu_star, cov_matrix)
    
    #---------------------------u.e---------------------------------------------
    # Compute squared residuals
    res_sq <- (y - xi %*% beta)^2
    
    # Generate new value for u.e
    u.e <- rgamma(n, nu/2 + 0.5, nu/2 + (0.5 * tauE * res_sq))
    
    #-----------------------------lambda.nu-------------------------------------
    # Generate new value for lambda.nu  
    lambda.nu <- truncated_gamma(n = 1, shape = 2, rate = nu, lower = c, upper = d)
    
    #-----------------------tauE------------------------------------------------
    # Generate new value for tauE
    tauE <- rgamma(1, a_0 + n/2, b_0 + 0.5*sum(u.e*res_sq))
    
    #---------------------------nu----------------------------------------------
    # Propose a new value for nu
    nu.new <- truncated_exponential(1, lambda.nu)
    
    # Calculate log prior
    lprior.nu <- - nu*lambda.nu
    lprior.nu.new <- - nu.new*lambda.nu
    
    # Calculate log likelihood
    llik.nu <- sum(0.5*nu*log(0.5*nu)-lgamma(0.5*nu)+(0.5*nu-1) * log(u.e) - (u.e*nu)/2)
    llik.nu.new <- sum(0.5*nu.new*log(0.5*nu.new)-lgamma(0.5*nu.new)+(0.5*nu.new-1) * log(u.e) - (u.e*nu.new)/2)
    
    # Calculate acceptance ratio  
    lpost.diff.nu <- (llik.nu.new + lprior.nu.new) - (llik.nu + lprior.nu)
    
    # Accept-reject ratio process
    if (lpost.diff.nu > log(runif(1))) 
    {
      # Set new parameters if accepted
      nu <- nu.new
    } 
    
    #-------------------------DIC-----------------------------------------------
    # Calculate and store the log-likelihood at each iteration
    y1 = (y - xi %*% beta) * sqrt(tauE)
    loglike_MCMC[k] <- sum(dt(y1, nu, log = T) + 0.5 * log(tauE))
    
    #---------------------------------------------------------------------------
    # Save the accepted values after each iteration
    tauE_MCMC[k] <- tauE
    beta_MCMC[k,] <- beta  
    nu_MCMC[k] <- nu
    lambda_MCMC[k] <- lambda.nu
    
    setTxtProgressBar(count, k, title = "Loading", label = "Welcome")
    
  }
  
  # Set burned-in range
  burned_in <- (burn_in + 1):niter
  
  # Compute y1 using matrix
  y1_matrix <- matrix(y, n, length(burned_in)) - xi %*% t(beta_MCMC[burned_in,])
  y1_matrix <- sweep(y1_matrix, 2, sqrt(tauE_MCMC[burned_in]), "*")
  
  # Compute likelihoods 
  likelihoods_matrix <- dt(y1_matrix, mean(nu_MCMC[burned_in])) * sqrt(tauE_MCMC[burned_in])
  
  # Compute CPO values
  CPO <- apply(likelihoods_matrix, 1, function(likelihoods) {
    1 / mean(1 / likelihoods)
  })
  
  # Calculate LPML
  LPML = sum(log(CPO))
  
  # Compute y2 
  y2 = (y - xi %*% colMeans(beta_MCMC[burned_in,])) * sqrt(mean(tauE_MCMC[burned_in]))
  
  # Calculate average deviance
  D_theta_avg = -2 * mean(loglike_MCMC[burned_in])
  
  # Compute DIC
  DIC = 2 * D_theta_avg - (-2 * sum(dt(y2, mean(nu_MCMC[burned_in]), log = TRUE) + 0.5 * log(mean(tauE_MCMC[burned_in]))))
  
  # Set number of parameters 
  theta = length(beta) + 2
  
  # Compute EAIC and EBIC
  EAIC = D_theta_avg + 2 * theta
  EBIC = D_theta_avg + theta * log(n)
  
  # Calculate p(z)
  log_lik_matrix <- log(likelihoods_matrix)
  p_z <- sum(apply(log_lik_matrix, 1, mean))
  
  # Compute WAIC1
  WAIC1 = 2 * p_z + 2 * D_theta_avg
  
  # Return the results
  return(list(tauE = tauE_MCMC[burned_in], beta = beta_MCMC[burned_in,], nu = nu_MCMC[burned_in], lambda = lambda_MCMC[burned_in], DIC = DIC, EAIC = EAIC, EBIC = EBIC, LPML = LPML, WAIC1 = WAIC1))
}

# Define Slash Model
MCMC_function_Slash <- function(a_0, b_0, c, d, omega_0, theta_0, niter, burn_in = 10000, y = y, xi=xi) {
  
  # Set initial values for MCMC variables
  m <- ncol(xi)
  n <- length(y)
  tauE <- 2
  u.e <- rep(0.1, n)
  nu <- 1
  cov_matrix <- matrix(0, m, m)
  mu_star <- rep(0, m)
  theta_0 <- rep(theta_0, m)  
  omega_0 <- diag(omega_0, m)
  
  # Initialise matrices to hold MCMC results
  tauE_MCMC <- rep(NA, niter)
  beta_MCMC <- matrix(NA, nrow = niter, ncol = m)
  u_MCMC <- matrix(NA, nrow = niter, ncol = n)
  nu_MCMC <- rep(NA, niter)
  lambda_MCMC <- rep(NA, niter)
  loglike_MCMC <- rep(NA, niter) 
  
  # Initialise progress bar for MCMC iterations
  count <- txtProgressBar(min = 1, max = niter, style = 3)
  
  # Define function to generate truncated gamma values - different method to avoid NaN values
  truncated_gamma <- function(n, shape, rate, lower = 0, upper = 1) {
    result <- numeric(0)
    while(length(result) < n) {
      temp <- rgamma(2 * n, shape, rate)
      temp <- temp[temp >= lower & temp <= upper]
      result <- c(result, temp)
    }
    return(result[1:n])
  }
  
  # Define function to calculate PDF of Slash
  dSlash = function(y, mu, sigma2, nu, log = F){  
    z = (y - mu) / sqrt(sigma2)
    PDF = log(nu/sqrt(2 * pi * sigma2)) + (nu + 0.5) * log(2/z^2) + pgamma(1, nu + 0.5, rate = 0.5 * (z^2), log.p = T) + lgamma(nu + 0.5)
    ifelse(log == F, return(exp(PDF)), return(PDF))
  }
  
  # Begin MCMC
  for (k in 1:niter) {
    
    #------------------beta-------------------------------------------------------
    # Compute y and x multiplied by u.e
    x_u <- xi * sqrt(u.e)
    y_u <- y * sqrt(u.e)
    
    # Compute xi cross product and xi weighted
    xi_cross <- tauE * tcrossprod(t(x_u))
    xi_w <- tauE * (t(x_u) %*% y_u)
    
    # Compute inverse of omega_0
    omega_0_inv <- solve(omega_0)
    
    # Update covariance matrix
    cov_matrix <- solve(xi_cross + omega_0_inv)
    
    # Update mean vector
    mu_star <- cov_matrix %*% (omega_0_inv %*% theta_0 + xi_w)
    
    # Generate new value for beta
    beta <- mvrnorm(1, mu_star, cov_matrix)
    
    #---------------------------u.e-----------------------------------------------
    # Compute squared residuals
    res_sq <- (y - xi %*% beta)^2
    
    # Generate new value for u.e
    u.e <- truncated_gamma(n = n, shape = nu + 0.5, rate = nu/2 + (0.5 * tauE * res_sq))
    
    #-----------------------tauE-------------------------------------------------
    # Generate new value for tauE
    tauE <- rgamma(1, a_0 + n/2, b_0 + 0.5*sum(u.e*res_sq))
    
    #-----------------------------lambda.nu---------------------------------------
    # Generate new value for lambda.nu  
    lambda.nu <- truncated_gamma(n = 1, shape = 2, rate = nu, lower = c, upper = d)
    
    #---------------------------nu------------------------------------------------
    # Generate new value for nu
    nu <- rgamma(n = 1, shape = n + 1, rate = lambda.nu - sum(log(u.e)))
    
    #-------------------------DIC-------------------------------------------------
    # Calculate and store the log-likelihood at each iteration
    y1 = (y - xi %*% beta) * sqrt(tauE)
    loglike_MCMC[k] <- sum(dSlash(y1, 0, 1, nu, log = T) + 0.5 * log(tauE))
    
    #-----------------------------------------------------------------------------
    # Save the accepted values after each iteration
    tauE_MCMC[k] <- tauE
    beta_MCMC[k,] <- beta  
    nu_MCMC[k] <- nu
    
    setTxtProgressBar(count, k, title = "Loading", label = "Welcome")
    
  }
  
  # Define the burned-in range
  burned_in <- (burn_in+1):niter
  
  # Matrix multiplication for all data points and all burned_in samples
  y1_matrix <- matrix(y, n, length(burned_in)) - xi %*% t(beta_MCMC[burned_in,])
  y1_matrix <- sweep(y1_matrix, 2, sqrt(tauE_MCMC[burned_in]), "*")
  
  # Compute likelihoods matrix
  likelihoods_matrix <- sapply(1:length(burned_in), function(j) {
    dSlash(y1_matrix[, j], 0, 1, nu_MCMC[burned_in[j]], log = FALSE) * sqrt(tauE_MCMC[burned_in[j]])
  })
  
  # Compute CPO values
  CPO <- apply(likelihoods_matrix, 1, function(likelihoods) {
    1 / mean(1 / likelihoods)
  })
  
  # Calculate LPML
  LPML <- sum(log(CPO))
  
  # Compute y2
  y2 <- (y - xi %*% colMeans(beta_MCMC[burned_in,])) * sqrt(mean(tauE_MCMC[burned_in]))
  
  # Calculate the average deviance
  D_theta_avg <- -2 * mean(loglike_MCMC[burned_in])
  
  # Compute DIC
  DIC <- 2 * D_theta_avg - (-2 * sum(dSlash(y2, 0, 1, mean(nu_MCMC[burned_in]), log = TRUE) + 0.5 * log(mean(tauE_MCMC[burned_in]))))
  
  # Set number of parameters 
  theta <- length(beta) + 2
  
  # Compute EAIC and EBIC
  EAIC <- D_theta_avg + 2 * theta
  EBIC <- D_theta_avg + theta * log(n)
  
  # Calculate p(z)
  log_lik_matrix <- log(likelihoods_matrix)
  p_z <- sum(apply(log_lik_matrix, 1, mean))
  
  # Compute WAIC1
  WAIC1 = 2 * p_z + 2 * D_theta_avg
  
  # Return the results after burn-in
  list(tauE = tauE_MCMC[burned_in], beta = beta_MCMC[burned_in,], nu = nu_MCMC[burned_in], DIC = DIC, EAIC = EAIC, EBIC = EBIC, LPML = LPML, WAIC1 = WAIC1)
}

# Define Contaminated Normal Model
MCMC_function_CN <- function(a_0, b_0, c, d, omega_0, theta_0, niter, burn_in = 10000, y = y, xi=xi) {
  
  # Set initial values for MCMC variables
  m <- ncol(xi)
  n <- length(y)
  tauE <- 2
  u.e <- rep(1, n)
  u <- rep(1, n)
  lambda.nu <- 0.3
  nu1 <- 0.4
  nu2 <- 0.4
  cov_matrix <- matrix(0, m, m)
  mu_star <- rep(0, m)
  theta_0 <- rep(theta_0, m)  
  omega_0 <- diag(omega_0, m)
  
  # Initialise matrices to hold MCMC results
  tauE_MCMC <- rep(NA, niter)
  beta_MCMC <- matrix(NA, nrow = niter, ncol = m)
  u_MCMC <- matrix(NA, nrow = niter, ncol = n)
  nu1_MCMC <- rep(NA, niter)
  nu2_MCMC <- rep(NA, niter)
  lambda_MCMC <- rep(NA, niter)
  loglike_MCMC <- rep(NA, niter) 
  
  # Initialise progress bar for MCMC iterations
  count <- txtProgressBar(min = 1, max = niter, style = 3)
  
  # Define function to calculate PDF of CN
  dCN  = function(y, mu, sigma2, nu, gama, log = F){
    PDF = vector(mode = "numeric", length = length(y))
    D1 = dnorm(y, mu, sqrt(sigma2 / gama), log = T)
    D2 = exp(dnorm(y, mu, sqrt(sigma2), log = T) -  D1)
    PDF = D1 + log(nu + (1 - nu) * D2)
    ifelse(log == T, return(PDF), return(exp(PDF)))
  }
  
  # Begin MCMC
  for (k in 1:niter) {
    
    #------------------beta-----------------------------------------------------
    # Compute y and x multiplied by u.e
    x_u <- xi * sqrt(u.e)
    y_u <- y * sqrt(u.e)
    
    # Compute x multiples
    xi_cross <- tauE * tcrossprod(t(x_u))
    xi_w <- tauE * (t(x_u) %*% y_u)
    
    # Compute inverse of omega_0
    omega_0_inv <- solve(omega_0)
    
    # Update covariance matrix
    cov_matrix <- solve(xi_cross + omega_0_inv)
    
    # Update mean vector
    mu_star <- cov_matrix %*% (omega_0_inv %*% theta_0 + xi_w)
    
    # Generate new value for beta
    beta <- mvrnorm(1, mu_star, cov_matrix)
    
    #---------------------------nu1---------------------------------------------
    # Generate new value for nu1
    nu1 <- rbeta(1, 1 + 1/(1-nu2)*sum(1-u), 1 + 1/(1-nu2)*sum(u-nu2))
    
    #---------------------------u.e---------------------------------------------
    # Compute squared residuals
    res_sq <- (y - xi %*% beta)^2
    
    # Set u.e depending on probability
    p1 <- nu1*sqrt(nu2)*exp(-0.5*tauE*nu2*res_sq)
    p2 <- (1-nu1)*exp(-0.5*tauE*res_sq)
    prob_denom <- p1+p2
    prob <- p1/prob_denom
    
    # Check for NaN due to division by 0 and set to 0
    prob[is.nan(prob)] <- 0
    
    # Generate a binary indicator for whether u.e should be nu2 or 1
    indicator <- rbinom(n, size = 1, prob = prob)
    
    # Set u.e to nu2 where the indicator is 1, and 1 where the indicator is 0
    u.e <- nu2 * indicator + (1 - indicator)
    
    #---------------------------nu----------------------------------------------
    # Generate new nu2 and u
    nu2.new <- rbeta(1, c, d)
    u.new <- ifelse(u == 1, 1, ifelse(u == nu2, nu2.new, u))
    
    # Calculate log prior
    lprior.nu2 <- dbeta(nu2, shape1 = c, shape2 = d, log = TRUE)
    lprior.nu2.new <- dbeta(nu2.new, shape1 = c, shape2 = d, log = TRUE)
    
    # Calculate log liklihood
    llik.nu2 <- 0.5 * sum(log(u)) - ((n - sum(u)) / (1 - nu2)) * log(nu1) + ((sum(u) - n * nu2) / (1 - nu2)) * log(1 - nu1)
    llik.nu2.new <- 0.5 * sum(log(u.new)) - ((n - sum(u.new)) / (1 - nu2.new)) * log(nu1) + ((sum(u.new) - n * nu2.new) / (1 - nu2.new)) * log(1 - nu1)
    
    # Calculate acceptance ratio  
    lpost.diff.nu <- (llik.nu2.new + lprior.nu2) - (llik.nu2 + lprior.nu2.new)
    
    # Accept-reject ratio process
    if (lpost.diff.nu > log(runif(1))) 
    {
      # Set new parameters if accepted
      nu2 <- nu2.new
      u <- u.new
    } 
    
    #-----------------------tauE------------------------------------------------
    # Generate new value for tauE
    tauE <- rgamma(1, a_0 + n/2, b_0 + 0.5*sum(u.e*res_sq))
    
    #-------------------------DIC-----------------------------------------------
    # Calculate and store the log-likelihood at each iteration
    y1 = (y - xi %*% beta) * sqrt(tauE)
    loglike_MCMC[k] <- sum(dCN(y1, 0, 1, nu1, nu2, log = T) + 0.5 * log(tauE))
    
    #---------------------------------------------------------------------------
    # Save the accepted values after each iteration
    tauE_MCMC[k] <- tauE
    beta_MCMC[k,] <- beta  
    nu1_MCMC[k] <- nu1
    nu2_MCMC[k] <- nu2
    
    setTxtProgressBar(count, k, title = "Loading", label = "Welcome")
    
  }
  
  # Define the burned-in range
  burned_in <- (burn_in+1):niter
  
  # Compute y1 using matrix
  y1_matrix <- matrix(y, n, length(burned_in)) - xi %*% t(beta_MCMC[burned_in,])
  y1_matrix <- sweep(y1_matrix, 2, sqrt(tauE_MCMC[burned_in]), "*")
  
  # Compute likelihoods matrix
  likelihoods_matrix <- sapply(1:length(burned_in), function(j) {
    dCN(y1_matrix[, j], 0, 1, nu1_MCMC[burned_in[j]], nu2_MCMC[burned_in[j]], log = FALSE) * sqrt(tauE_MCMC[burned_in[j]])
  })
  
  # Compute CPO values
  CPO <- apply(likelihoods_matrix, 1, function(likelihoods) {
    1 / mean(1 / likelihoods)
  })
  
  # If value less than tolerance, set as tolerance - prevents -Inf which often occurred
  CPO[CPO < 1e-10] <- 1e-10
  
  # Calculate LPML
  LPML <- sum(log(CPO))
  
  # Calculate y2
  y2 <- (y - xi %*% colMeans(beta_MCMC[burned_in,])) * sqrt(mean(tauE_MCMC[burned_in]))
  
  # Compute the average deviance
  D_theta_avg <- -2 * mean(loglike_MCMC[burned_in])
  
  # Compute DIC
  DIC <- 2 * D_theta_avg - (-2 * sum(dCN(y2, 0, 1, mean(nu1_MCMC[burned_in]), mean(nu2_MCMC[burned_in]), log = TRUE) + 0.5 * log(mean(tauE_MCMC[burned_in]))))
  
  # Set number of parameters 
  theta <- length(beta) + 3
  
  # Compute EAIC and EBIC
  EAIC <- D_theta_avg + 2 * theta
  EBIC <- D_theta_avg + theta * log(n)
  
  # Calculate p(z)
  likelihoods_matrix[likelihoods_matrix < 1e-10] <- 1e-10
  log_lik_matrix <- log(likelihoods_matrix)
  p_z <- sum(apply(log_lik_matrix, 1, mean))
  
  # Compute WAIC1
  WAIC1 = 2 * p_z + 2 * D_theta_avg
  
  # Return the results after burn-in
  list(tauE = tauE_MCMC[burned_in], beta = beta_MCMC[burned_in,], nu1 = nu1_MCMC[burned_in], nu2 = nu2_MCMC[burned_in], DIC = DIC, EAIC = EAIC, EBIC = EBIC, LPML = LPML, WAIC1 = WAIC1)
}

