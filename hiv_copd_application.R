# -----------------------------------------------------------------------------
# Applying the BPMF model to the HIV-COPD BALF Biocrates and Somascan data. 
# Comparing the model application when there is no sparsity vs. when
# there is sparsity. 
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set-up
# -----------------------------------------------------------------------------

# Packages
library(doParallel)
library(foreach)

# Working directories
data_wd <- "~/BayesianPMFWithGit/data/"
model_wd <- "~/BayesianPMFWithGit/"
results_wd <- "~/BayesianPMFWithGit/data_application_results/"

# Loading in the model functions
source(paste0(model_wd, "bpmf.R"))


# -----------------------------------------------------------------------------
# Preparing the data
# -----------------------------------------------------------------------------

# Loading data in 
load(paste0(data_wd, "BiocratesLavageProcessed.rda"))
load(paste0(data_wd, "HIV_COPD_SomaScan_Normalized_Clean.rda"))

# BALF Biocrates
lavage_processed_no_info <- lavage_processed[, -c(1:2)] # removing meta
bioc_meta_data <- lavage_processed[, 1:2] # saving meta elsewhere
lavage_processed_no_info_log_scale <- scale(log(lavage_processed_no_info + 1), # log/scale
                                            center = TRUE, scale = TRUE)
rownames(lavage_processed_no_info_log_scale) <- bioc_meta_data$Metabolite

# BALF Somascan
somascan_normalized_clean_no_info <- somascan_normalized_clean[, -c(1:29)] # removing meta
soma_meta_data <- somascan_normalized_clean[, 1:29] # saving meta elsewhere
soma_pids <- somascan_normalized_clean$PID
somascan_normalized_clean_no_info_transpose <- t(somascan_normalized_clean_no_info)
somascan_normalized_clean_no_info_transpose_scale <- t(apply(somascan_normalized_clean_no_info_transpose, 1, function(row) {
  scale(log(row + 1), center = TRUE, scale = TRUE)
})) # scaling the rows
colnames(somascan_normalized_clean_no_info_transpose_scale) <- soma_pids

# Reorder the Somascan cols to match the Biocrates cols
somascan_normalized_clean_no_info_transpose_reorder <- 
  somascan_normalized_clean_no_info_transpose_scale[, match(colnames(lavage_processed_no_info_log_scale),
                                                            colnames(somascan_normalized_clean_no_info_transpose_scale))]

# Check
all(colnames(somascan_normalized_clean_no_info_transpose_reorder) == colnames(lavage_processed_no_info_log_scale))
dim(somascan_normalized_clean_no_info_transpose_reorder)
dim(lavage_processed_no_info)

# Combining the sources together
hiv_copd_data <- matrix(list(), nrow = 2, ncol = 1)
hiv_copd_data[[1,1]] <- lavage_processed_no_info_log_scale
hiv_copd_data[[2,1]] <- somascan_normalized_clean_no_info_transpose_reorder

# Preparing the response variable
fev1pp <- matrix(list(), nrow = 1, ncol = 1)
fev1pp[[1,1]] <- matrix(subject_processed$FEV1_percent_predicted, ncol = 1)

# -----------------------------------------------------------------------------
# Model parameters
# -----------------------------------------------------------------------------

# Data dimensions
p1 <- nrow(lavage_processed_no_info_log_scale)
p2 <- nrow(somascan_normalized_clean_no_info_transpose_reorder)
n <- ncol(lavage_processed_no_info_log_scale)

# Error variances
sigma21 <- 1 # error variance for Biocrates
sigma22 <- 1 # error variance for Somascan

# Variance of joint structure
lambda_joint <- sqrt(nrow(lavage_processed_no_info_log_scale) + nrow(somascan_normalized_clean_no_info_transpose_reorder)) + 
  sqrt(ncol(lavage_processed_no_info_log_scale))
sigma2_joint <- 1/(lambda_joint)

# Variance of individual structure for Biocrates
lambda_indiv1 <- sqrt(nrow(lavage_processed_no_info_log_scale)) + sqrt(ncol(lavage_processed_no_info_log_scale))
sigma2_indiv1 <- 1/(lambda_indiv1) 

# Variance of individual structure for Somascan
lambda_indiv2 <- sqrt(nrow(somascan_normalized_clean_no_info_transpose_reorder)) + sqrt(ncol(lavage_processed_no_info_log_scale))
sigma2_indiv2 <- 1/(lambda_indiv2)

# For the regression coefficients, beta
lambda2_intercept <- 1e6
lambda2_joint <- 1
lambda2_indiv1 <- 1
lambda2_indiv2 <- 1
beta_vars <- c(lambda2_intercept, lambda2_joint, lambda2_indiv1, lambda2_indiv2)

# For the response vector
shape <- 0.01
rate <- 0.01

# Putting the model parameters together
model_params <- list(error_vars = c(sigma21, sigma22),
                     joint_var = sigma2_joint,
                     indiv_vars = c(sigma2_indiv1, sigma2_indiv2),
                     beta_vars = beta_vars,
                     response_vars = c(shape = shape, rate = rate))

# Gibbs sampler parameters
nsample <- 5000
burnin <- nsample/2
thinned_iters <- seq(1, nsample, by = 10)
thinned_iters_burnin <- seq(burnin, nsample, by = 10)

# -----------------------------------------------------------------------------
# Training Data Model Fit
# -----------------------------------------------------------------------------

# Fitting the model without sparsity first
fev1pp_training_fit_nonsparse <- bpmf(
  data = hiv_copd_data,
  Y = fev1pp,
  nninit = TRUE,
  model_params = model_params,
  sparsity = FALSE,
  nsample = nsample,
  progress = TRUE
)

# Fitting the model with sparsity
fev1pp_training_fit_sparse <- bpmf(
  data = hiv_copd_data,
  Y = fev1pp,
  nninit = TRUE,
  model_params = model_params,
  sparsity = TRUE,
  nsample = nsample,
  progress = TRUE
)

# Save the results
save(fev1pp_training_fit_nonsparse, fev1pp_training_fit_sparse, file = paste0(results_wd, "training_data_fit.rda"))

# Assessing convergence for both model fits
fev1pp_training_nonsparse_conv <- sapply(thinned_iters, function(sim_iter) {
  # Calculate the log-joint density at each thinned iterations
  log_joint_density(data = hiv_copd_data, 
                    U.iter = fev1pp_training_fit_nonsparse$U.draw[[sim_iter]], 
                    V.iter = fev1pp_training_fit_nonsparse$V.draw[[sim_iter]], 
                    W.iter = fev1pp_training_fit_nonsparse$W.draw[[sim_iter]], 
                    Vs.iter = fev1pp_training_fit_nonsparse$Vs.draw[[sim_iter]],
                    model_params = model_params,
                    ranks = fev1pp_training_fit_nonsparse$ranks,
                    Y = fev1pp,
                    beta.iter = fev1pp_training_fit_nonsparse$beta.draw[[sim_iter]],
                    tau2.iter = fev1pp_training_fit_nonsparse$tau2.draw[[sim_iter]])
})


