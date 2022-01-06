# -----------------------------------------------------------------------------
# Comparing the performance of BPMF to other low-rank factorization 
# models. Performance is assessed based on (1) recovery of underlying
# structure and (2) predictive performance. 
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set up
# -----------------------------------------------------------------------------

# Load in helper functions
source("~/BayesianPMFWithGit/bpmf.R")

# Number of data sources
q <- 2

# Joint and individual ranks (ranks[1] = joint, ranks[2:r_total] = individual)
ranks <- c(2,2,2)

# Number of predictors for each source
p.vec = c(10, 15)

# Sample size 
n = 5

# Parameters for data generation
true_params <- model_params <- list(error_vars = c(1,1), # Error variance for each source
                                    joint_var = 1, # Variance for joint structure
                                    indiv_vars = c(1,1), # Variance for each individual structure
                                    beta_vars = c(10, 1, rep(1, q)), # Variance of intercept effect and each joint effect 
                                    response_vars = c(shape = 1, rate = 1)) # Hyperparameters for variance of response (if continuous)


# -----------------------------------------------------------------------------
# sJIVE
# -----------------------------------------------------------------------------

run_each_mod(mod = "sJIVE", p.vec, n, ranks, response = "continuous", true_params, 
             s2nX = 1, s2nY = 1, nsim = 100, nsample = 1000, n_clust = 10)

