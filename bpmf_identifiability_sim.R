# -----------------------------------------------------------------------------
# Simulation study to assess the effectiveness of our approach to undoing
# the permutation invariance of our Bayesian model
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set up
# -----------------------------------------------------------------------------

# Load in helper functions
source("~/BayesianPMFWithGit/bpmf.R")

# Number of data sources
q <- 2

# Joint and individual ranks (ranks[1] = joint, ranks[2:r_total] = individual)
ranks <- c(5,5,5)

# Number of predictors for each source
p.vec = c(100, 100)

# Sample size 
n = 300

# Parameters for data generation
true_params <- list(error_vars = c(1,1), # Error variance for each source
                    joint_var = 1, # Variance for joint structure
                    indiv_vars = c(1,1), # Variance for each individual structure
                    beta_vars = c(10, 1, rep(1, q)), # Variance of intercept effect and each joint effect 
                    response_vars = c(shape = 1, rate = 1)) # Hyperparameters for variance of response (if continuous)

# Setting the model variances
model_params <- list(error_vars = c(1,1), # Error variance for each source
                     joint_var = 1/(sqrt(n) + sqrt(sum(p.vec))), # Variance for joint structure
                     indiv_vars = c(1/(sqrt(n) + sqrt(p.vec[1])), 1/(sqrt(n) + sqrt(p.vec[2]))), # Variance for each individual structure
                     beta_vars = c(10, 1, rep(1, q)), # Variance of intercept effect and each joint effect 
                     response_vars = c(shape = 1, rate = 1)) # Hyperparameters for variance of response (if continuous)

# Parameters for the simulation
nsample <- 2000
nsim <- 100

# -----------------------------------------------------------------------------
# Run the simulation
# -----------------------------------------------------------------------------

ident_sim_results <- identifiability_sim(p.vec, n, ranks, response = "continuous", true_params,
                                         model_params, sparsity = TRUE, nsim = 100, nsample = 2000)
