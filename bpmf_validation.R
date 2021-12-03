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
                                    beta_vars = c(10, 1, rep(1, q)), # Use the same variance for all the effects from each source
                                    response_vars = c(shape = 1,rate = 1))

# Generating example data
data <- bpmf_data(p.vec, n, ranks, true_params, response = "continuous", missingness = NULL, entrywise = NULL, prop_missing = NULL)
X <- data$data
y <- data$Y


# -----------------------------------------------------------------------------
# Coverage simulations
# -----------------------------------------------------------------------------

# No response, no missingness
no_response_no_missing <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks)

# Continuous response
response_continuous <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous")

# Binary response
response_binary <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary")

# Missing continuous response
response_continuous_missing0.1 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.1)
response_continuous_missing0.3 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.3)
response_continuous_missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.5)
response_continuous_missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.7)
response_continuous_missing0.9 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.9)

# Missing data
missing0.1 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.1)
missing0.3 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.3)
missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.5)
missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.7)
missing0.9 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.9)

# Save results
save(test, test_with_response_continuous, test_with_missing0.1, test_with_missing0.3, test_with_missing0.5,
     test_with_missing0.7, test_with_missing0.9, test_with_response_binary, file = "~/BayesianPMFWithGit/validation_testing/validation_results_1292021.rda")


# -----------------------------------------------------------------------------
# Creating the results table
# -----------------------------------------------------------------------------

# Create the dataframe
validation_results <- data.frame(Condition = character(), Metric = character(1), Joint = numeric(), 
                                 Indiv = numeric(), EY = numeric(), tau2 = numeric(), Xm = numeric(), 
                                 Ym = numeric())

# No response, no missing
validation_results <- rbind.data.frame(validation_results,
                                       create_validation_table(test, condition = "No Response, No Missing"))


# Continuous response, no missing
validation_results <- rbind.data.frame(validation_results,
                                       create_validation_table(response_continuous, condition = "Continuous Response"))
try <- validation_results %>% gather(key = Type, value = Measurement, Joint:Ym)

# Exporting the table to LaTeX
library(xtable)
library(magrittr)
library(tidyverse)

# Basic
xtable(validation_results, digits = 4)

validation_tbl <- ftable(Type ~ Condition + Metric, data =  try)
validation_ht <- as_hux(validation_tbl)
  



