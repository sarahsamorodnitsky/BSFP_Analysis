# -----------------------------------------------------------------------------
# Conducting model validation simulations to ensure coverage, MSE, and CI
# width look appropriate. 
#
# This script contains finalized results for validation studies. 
# -----------------------------------------------------------------------------

# Load in the helper functions
source("~/BayesianPMFWithGit/bpmf.R")

# Setting up the data
n <- 50
p.vec <- c(150, 100)
r <- 1
r.vec <- c(1, 1)
ranks <- c(r, r.vec)
q <- 2

# -----------------------------------------------------------------------------
# Coverage simulations - No response, no missingness
# -----------------------------------------------------------------------------

# No response, no missingness

# Setting up the model parameters
model_params <- true_params <- list(error_vars = c(1,1),
                                    joint_var = 1,
                                    indiv_vars = c(1,1),
                                    beta_vars = c(1, 1, rep(1, q)), # Use the same variance for all the effects from each source
                                    response_vars = c(shape = 1, rate = 1))

no_response_no_missing <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks)

# -----------------------------------------------------------------------------
# Coverage simulations - Continuous response, no missingness
# -----------------------------------------------------------------------------

# Continuous response
response_continuous <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous")

# -----------------------------------------------------------------------------
# Coverage simulations - Binary response, no missingness
# -----------------------------------------------------------------------------

# Binary response
response_binary <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary")

# -----------------------------------------------------------------------------
# Coverage simulations - Continuous response with missingness
# -----------------------------------------------------------------------------

# Missing continuous response
response_continuous_missing0.3 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100,  center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.3)
response_continuous_missing0.5 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.5)
response_continuous_missing0.7 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.7)

# -----------------------------------------------------------------------------
# Coverage simulations - Binary response with missingness
# -----------------------------------------------------------------------------

# Missing binary response
response_binary_missing0.3 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.3)
response_binary_missing0.5 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.5)
response_binary_missing0.7 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.7)

# -----------------------------------------------------------------------------
# Coverage simulations - No response, entrywise missing data
# -----------------------------------------------------------------------------

# Entrywise missing data
missing0.3 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.3)
missing0.5 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.5)
missing0.7 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.7)

# -----------------------------------------------------------------------------
# Coverage simulations - No response, columnwise missing data
# -----------------------------------------------------------------------------

# Columnwise missing data
columnwise_missing0.1 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = FALSE, prop_missing = 0.1)
columnwise_missing0.3 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = FALSE, prop_missing = 0.3)
columnwise_missing0.5 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = FALSE, prop_missing = 0.5)

# Considering columnwise missingness with many Gibbs sampling iterations
columnwise_missing0.3 <- validation_simulation(nsample = 20000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = FALSE, prop_missing = 0.3)

# -----------------------------------------------------------------------------
# Coverage simulations - Continuous response, entrywise missing data
# -----------------------------------------------------------------------------

# Continuous response
response_continuous_missing0.1 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.1)
response_continuous_missing0.3 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.3)
response_continuous_missing0.5 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.5)
response_continuous_missing0.7 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.5)
response_continuous_missing0.9 <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.5)

# -----------------------------------------------------------------------------
# Coverage simulation - validate prediction on a test dataset with NO response
# -----------------------------------------------------------------------------

no_response_test_data <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = NULL, missingness = NULL, predict_test_data = TRUE)
continuous_response_test_data <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = NULL, predict_test_data = TRUE)
binary_response_test_data <- validation_simulation(nsample = 2000, n_clust = 10, p.vec = p.vec, n = n, true_params = true_params, model_params = model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary", missingness = NULL, predict_test_data = TRUE)

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------

save(no_response_no_missing, response_continuous, 
     response_continuous_missing0.3, response_continuous_missing0.5, response_continuous_missing0.7,
     missing0.3, missing0.5, missing0.7, 
     columnwise_missing0.3, columnwise_missing0.5, response_continuous_missing0.1,
     response_continuous_missing0.3,
     file = "~/BayesianPMFWithGit/validation_results/validation_sim_5112022.rda")

# -----------------------------------------------------------------------------
# Creating the results table
# -----------------------------------------------------------------------------

load("~/BayesianPMFWithGit/validation_results/validation_sim_162022.rda", verbose = TRUE)

# Create the dataframe
validation_results <- data.frame(Condition = character(), Metric = character(), 
                                 Joint_Obs = numeric(), Indiv_Obs = numeric(), 
                                 Joint_Mis = numeric(), Indiv_Obs = numeric(),
                                 EY_Obs = numeric(), EY_Mis = numeric(), tau2 = numeric())

# No response, no missing
validation_results <- rbind.data.frame(validation_results, create_validation_table(no_response_no_missing, condition = "No Response, No Missing"))

# Continuous response, no missing
validation_results <- rbind.data.frame(validation_results, create_validation_table(response_continuous, condition = "Continuous Response"))

# Binary response, no missing
validation_results <- rbind.data.frame(validation_results, create_validation_table(response_binary, condition = "Binary Response"))

# Entrywise missingness in data
validation_results <- rbind.data.frame(validation_results, create_validation_table(missing0.3, condition = "30% Entrywise Missing"))
validation_results <- rbind.data.frame(validation_results, create_validation_table(missing0.5, condition = "50% Entrywise Missing"))
validation_results <- rbind.data.frame(validation_results, create_validation_table(missing0.7, condition = "70% Entrywise Missing"))

# Missing continuous response
validation_results <- rbind.data.frame(validation_results, create_validation_table(response_continuous_missing0.3, condition = "30% Continuous Response Missing"))
validation_results <- rbind.data.frame(validation_results, create_validation_table(response_continuous_missing0.5, condition = "50% Continuous Response Missing"))
validation_results <- rbind.data.frame(validation_results, create_validation_table(response_continuous_missing0.7, condition = "70% Continuous Response Missing"))

# Missing binary response
validation_results <- rbind.data.frame(validation_results, create_validation_table(response_binary_missing0.3, condition = "30% Binary Response Missing"))
validation_results <- rbind.data.frame(validation_results, create_validation_table(response_binary_missing0.5, condition = "50% Binary Response Missing"))
validation_results <- rbind.data.frame(validation_results, create_validation_table(response_binary_missing0.7, condition = "70% Binary Response Missing"))


# Sparsity
validation_results <- rbind.data.frame(validation_results, create_validation_table(continuous_response_sparsity, condition = "Continuous Response with Sparsity"))
validation_results <- rbind.data.frame(validation_results, create_validation_table(binary_response_sparsity, condition = "Binary Response with Sparsity"))

# Columnwise missingness
validation_results <- rbind.data.frame(validation_results, create_validation_table(columnwise_missing0.3, condition = "30% Columnwise Missingness"))
validation_results <- rbind.data.frame(validation_results, create_validation_table(columnwise_missing0.5, condition = "50% Columnwise Missingness"))


# Exporting the table to LaTeX
library(xtable)
library(magrittr)
library(tidyverse)

# Basic
print(xtable(validation_results, digits = 4), include.rownames = FALSE)


# -----------------------------------------------------------------------------
# Coverage simulations when ranks are 0
# -----------------------------------------------------------------------------

r <- 0
r.vec <- c(1, 0)
ranks <- c(r, r.vec)

# No response, no missingness
no_response_no_missing <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2nX = NULL, s2nY = NULL, center = FALSE, nninit = FALSE, ranks = ranks)

# Continuous response
response_continuous_joint_rank_only <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2nX = NULL, s2nY = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous")

# Binary response
response_binary_joint_rank_only <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2nX = NULL, s2nY = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary")

# -----------------------------------------------------------------------------
# Coverage simulations with 3 data sources
# -----------------------------------------------------------------------------

# Setting up the data
n <- 50
p.vec <- c(150, 100, 100)
r <- 1
r.vec <- c(1, 1, 1)
ranks <- c(r, r.vec)
q <- 3

# Setting up the model parameters
model_params <- true_params <- list(error_vars = c(1,1,1),
                                    joint_var = 1,
                                    indiv_vars = c(1,1,1),
                                    beta_vars = c(1, 1, rep(1, q)), # Use the same variance for all the effects from each source
                                    response_vars = c(shape = 1, rate = 1))

# No response, no missingness
no_response_no_missing <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2nX = NULL, s2nY = NULL, center = FALSE, nninit = FALSE, ranks = ranks)

# Continuous response
response_continuous <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2nX = NULL, s2nY = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous")

# Binary response
response_binary <- validation_simulation(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2nX = NULL, s2nY = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary")
