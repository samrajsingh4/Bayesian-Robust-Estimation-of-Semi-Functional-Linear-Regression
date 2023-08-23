# load necessary libraries and functions
library(fda)
library(MASS)
library(xtable)
library(fda.usc) 
library(coda)
library(ggplot2)
source('MCMC_Functions.R')
source('Common-Functions.R')

# Load data
data(tecator) 

# Extract t range
t <- tecator$absorp.fdata$argvals

# Perform FPCA on  data
absorp <- tecator$absorp.fdata$data
fpca_results_absorp <- perform_fpca(absorp, t, num_basis=240, num_scores=9)

# Set variable values
scores_absorp <- fpca_results_absorp$scores
Y  <- tecator$y$Fat
protein  <- tecator$y$Protein
moisture <- tecator$y$Water

# Set data matrix
combined_scores <- cbind(scores_absorp, protein, moisture)

# Plot graph of curves
matplot(t, t(absorp), type = "l", lty = 1, col=rainbow(ncol(absorp)), 
        xlab = "Wavelength (nm)", 
        ylab = "Absorbance T(t)",
        main = "Absorbance Spectrum of Meat Samples",
        las=1) 

distributions = list('Normal' = MCMC_function_Normal, 'Student-t' = MCMC_function_T, 'CN' = MCMC_function_CN, 'Slash' = MCMC_function_Slash)

# Initialise accumulators
dic_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
eaic_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
ebic_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
waic1_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
LPML_scores_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
hpd_beta_accum <- setNames(vector("list", length(names(distributions))), names(distributions))
hpd_tauE_accum <- setNames(vector("list", length(names(distributions))), names(distributions))
tauE_mean_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
tauE_sd_accum <- setNames(rep(0, length(names(distributions))), names(distributions))
beta_mean_accum <- setNames(vector("list", length(names(distributions))), names(distributions))
beta_sd_accum <- setNames(vector("list", length(names(distributions))), names(distributions))

# Run MCMC functions 100 times for each distribution
set.seed(123)
for (i in 1:100) {
  print(i)
  for(name in names(distributions)) {
    if(name == 'Normal'){
      res <- combined_MCMC_function(distribution = distributions[[name]], a_0 = 0.001, b_0 = 0.001, omega_0 = 1000, theta_0 = 1 + i, niter = 50000, burn_in = 10000, y = Y, xi = combined_scores)
    } else if(name == 'CN') {
      res <- combined_MCMC_function(distribution = distributions[[name]], a_0 = 0.001, b_0 = 0.001, c = 1, d = 1, omega_0 = 1000, theta_0 = 1 + i, niter = 50000, burn_in = 10000, y = Y, xi = combined_scores)
    } else {
      res <- combined_MCMC_function(distribution = distributions[[name]], a_0 = 0.001, b_0 = 0.001, c = 0.001, d = 1, omega_0 = 1000, theta_0 = 1 + i, niter = 50000, burn_in = 10000, y = Y, xi = combined_scores)
    }
    
    # Accumulate DIC, EAIC, EBIC, and LPML results
    dic_scores_accum <- accumulate_scores(dic_scores_accum, name, res$DIC)
    eaic_scores_accum <- accumulate_scores(eaic_scores_accum, name, res$EAIC)
    ebic_scores_accum <- accumulate_scores(ebic_scores_accum, name, res$EBIC)
    waic1_scores_accum <- accumulate_scores(waic1_scores_accum, name, res$WAIC1)
    LPML_scores_accum <- accumulate_scores(LPML_scores_accum, name, res$LPML)
    
    # Calculate and accumulate HPD
    hpd_beta <- HPDinterval(as.mcmc(res$beta), prob = 0.95)
    hpd_tauE <- HPDinterval(as.mcmc(res$tauE), prob = 0.95)
    hpd_beta_accum <- accumulate_data_frame(hpd_beta_accum, name, as.data.frame(hpd_beta))
    hpd_tauE_accum <- accumulate_data_frame(hpd_tauE_accum, name, as.data.frame(hpd_tauE))
    
    
    # Beta mean and SD accumulations
    beta_mean_accum <- accumulate_list(beta_mean_accum, name, colMeans(res$beta))
    tauE_mean_accum <- accumulate_scores(tauE_mean_accum, name, mean(res$tauE))
    beta_sd_accum <- accumulate_list(beta_sd_accum, name, apply(res$beta, 2, sd))
    tauE_sd_accum <- accumulate_scores(tauE_sd_accum, name, sd(res$tauE))
  }
}

# Calculate average values of selection criteria
avg_scores <- list(
  DIC = dic_scores_accum / 100,
  EAIC = eaic_scores_accum / 100,
  EBIC = ebic_scores_accum / 100,
  LPML = LPML_scores_accum / 100,
  WAIC1 = waic1_scores_accum / 100
)

print(avg_scores)

# Calculate and store the average results
average_hpd_beta <- lapply(hpd_beta_accum, function(x) colMeans(x))
average_hpd_tauE <- lapply(hpd_tauE_accum, function(x) colMeans(x))
average_beta <- lapply(beta_mean_accum, function(x) colMeans(x))
average_sd_beta <- lapply(beta_sd_accum, function(x) colMeans(x))
hpd_beta_accum_divided <- lapply(hpd_beta_accum, function(df) df / 100)
hpd_tauE_accum_divided <- lapply(hpd_tauE_accum, function(df) df / 100)
average_tauE <- tauE_mean_accum / 100
average_sd_tauE <- tauE_sd_accum / 100

# Set models and initialise data
models <- c("Normal", "Student-t", "CN", "Slash")
data_list <- list()

for (model in models) {
  beta_values <- round(average_beta[[model]], 2)
  sd_values <- round(average_sd_beta[[model]], 2)
  
  hpd_beta <- hpd_beta_accum_divided[[model]]
  hpd_beta$lower <- round(hpd_beta$lower, 2)
  hpd_beta$upper <- round(hpd_beta$upper, 2)
  
  for (i in 1:length(beta_values)) {
    data_list[[paste(model, "Beta", i)]] <- c(model, paste("Beta", i), beta_values[i], sd_values[i], hpd_beta$lower[i], hpd_beta$upper[i])
  }
  
  tauE_value <- round(average_tauE[model], 2)
  sd_tauE_value <- round(average_sd_tauE[model], 2)
  
  hpd_tauE <- hpd_tauE_accum_divided[[model]]
  hpd_tauE$lower <- round(hpd_tauE$lower, 2)
  hpd_tauE$upper <- round(hpd_tauE$upper, 2)
  
  data_list[[paste(model, "tauE")]] <- c(model, "tau_E", tauE_value, sd_tauE_value, hpd_tauE$lower, hpd_tauE$upper)
}

# Print results of Mean, SD, HPD
df <- do.call(rbind, data_list)
colnames(df) <- c("Model", "Parameter", "Mean", "SD", "HPD Lower", "HPD Upper")
rownames(df) <- NULL
print(df)

# Set eigenfunctions
eigenfunctions_eval <- fpca_results_absorp$eigenfunctions_eval

# Compute functional coefficient
beta_means_T <- average_beta$`Student-t`[1:9]
result <- eigenfunctions_eval[,1:9] %*% matrix(beta_means_T, ncol = 1)

# Extract the lower and upper values 
lower_values <- hpd_beta_accum_divided$`Student-t`[1:9, "lower"]
upper_values <- hpd_beta_accum_divided$`Student-t`[1:9, "upper"]

# Compute credible interval values
result_lower <- eigenfunctions_eval[,1:9] %*% matrix(lower_values, ncol = 1)
result_upper <- eigenfunctions_eval[,1:9] %*% matrix(upper_values, ncol = 1)

# Plot graph for functional coefficient with CI
plot(t, result, type = "l", lty = 1, col="black", 
     xlab = "Wavelength (nm)", 
     ylab = expression(beta(t)),
     main = "Estimated Functional Coefficient",
     las=1)

# Shade area between upper and lower CI
polygon(c(t, rev(t)), c(result_upper, rev(result_lower)), col = 'grey80', border = NA)

# Add line over shaded area
lines(t, result, col = "black")

# Set variance explained by FPCs
variance_explained <- fpca_results_absorp$variance_explained

# Compute cumulative variance explained
cum_sum <- cumsum(variance_explained)

# Create a dataframe with the data for the plots
explained_data <- data.frame(
  PC = 1:length(variance_explained),
  Var = variance_explained,
  Cum_sum = cum_sum
)

# Calculate how many FPCs explain specified variance
var_80 <- which.max(cum_sum >= .999999)

# Plot the cumulative variance explained 
ggplot(explained_data, aes(x = PC, y = Cum_sum)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0.999999, linetype="dashed", color = "red") +
  geom_vline(xintercept = var_80, linetype="dashed", color = "red") +
  annotate("text", x = var_80, y = 0.999999, label = paste("FPC", var_80), vjust = -2, hjust = -0.5) +
  xlab("Functional Principal Component") +
  ylab("Cumulative Proportion of Variance Explained") +
  ggtitle("Cumulative Variance Explained by FPCs") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )
