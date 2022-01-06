# -----------------------------------------------------------------------------
# Conducting model validation simulations to ensure coverage, MSE, and CI
# width look appropriate. 
#
# This script contains finalized results for validation studies. 
# -----------------------------------------------------------------------------

# Load in the helper functions
source("bpmf.R")

# Setting up the data
n <- 50
p.vec <- c(150, 100)
r <- 1
r.vec <- c(1, 1)
ranks <- c(r, r.vec)
q <- 2

# Setting up the model parameters
model_params <- true_params <- list(error_vars = c(1,1),
                                    joint_var = 1,
                                    indiv_vars = c(1,1),
                                    beta_vars = c(1, 1, rep(1, q)), # Use the same variance for all the effects from each source
                                    response_vars = c(shape = 1, rate = 1))

# -----------------------------------------------------------------------------
# Coverage simulations
# -----------------------------------------------------------------------------

# No response, no missingness
no_response_no_missing <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks)

# Continuous response
response_continuous <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous")

# Binary response
response_binary <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary")

# Missing continuous response
response_continuous_missing0.3 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100,  center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.3)
response_continuous_missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.5)
response_continuous_missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.7)

# Missing binary response
response_binary_missing0.3 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.3)
response_binary_missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.5)
response_binary_missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.7)

# Entrywise missing data
missing0.3 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.3)
missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.5)
missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.7)

# Columnwise missing data
columnwise_missing0.1 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = FALSE, prop_missing = 0.1)
columnwise_missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = FALSE, prop_missing = 0.5)
columnwise_missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = FALSE, prop_missing = 0.7)

# Sparsity
continuous_response_sparsity <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks, response = "continuous", sparsity = TRUE)
binary_response_sparsity <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks, response = "binary", sparsity = TRUE)

save(no_response_no_missing, response_continuous, response_binary, 
     response_continuous_missing0.3, response_continuous_missing0.5, response_continuous_missing0.7,
     missing0.3, missing0.5, missing0.7, continuous_response_sparsity, binary_response_sparsity, 
     response_binary_missing0.3, response_binary_missing0.5, response_binary_missing0.7,
     columnwise_missing0.3, columnwise_missing0.5,
     file = "~/BayesianPMFWithGit/validation_results/validation_sim_162022.rda")

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




