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
n = 100

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

# Saving the different potential s2n's for X and Y
s2nX.list <- s2nY.list <- c(0.99/0.01, 0.5/0.5, 0.01/0.99)

# -----------------------------------------------------------------------------
# Run the simulation 
# -----------------------------------------------------------------------------

# Initialize a list to contain the results
ident_res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1 # Index the results

# Iterate through the different conditions
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    
    # -----------------------------------------------------------------------------
    # Run the simulation with no s2n adjustment in Y, initializing at the truth
    # -----------------------------------------------------------------------------
    
    res <- identifiability_sim(p.vec, n, ranks, response = "continuous", true_params,
                               model_params, sparsity = TRUE, s2nX = s2nX, s2nY = NULL,
                               init_at_truth = TRUE, nsim = nsim, nsample = nsample)
    ident_res[[ind]] <- res
    
    # Save the results
    settings <- list(s2nX = s2nX, s2nY = NULL, init_at_truth = TRUE, num_in_spike = c(2,2,2))
    save(res, settings,
         file = paste0("~/BayesianPMF/03Simulations/Identifiability/id_sim_s2nX_", s2nX, "_s2nY_NULL_spike222_initattruth_TRUE.rda"))
    
    # Update counter
    ind <-  ind + 1
    
    # Remove res to save on memory
    rm(res)
    gc()
  }
}

for (s2nX in s2nX.list) {
#  for (s2nY in s2nY.list) {
    
    # -----------------------------------------------------------------------------
    # Run the simulation with no s2n adjustment in Y, initializing at the truth
    # -----------------------------------------------------------------------------
    
    res <- identifiability_sim(p.vec, n, ranks, response = "continuous", true_params,
                               model_params, sparsity = TRUE, s2nX = s2nX, s2nY = NULL,
                               init_at_truth = FALSE, nsim = nsim, nsample = nsample)
    ident_res[[ind]] <- res
    
    # Save the results
    settings <- list(s2nX = s2nX, s2nY = NULL, init_at_truth = TRUE, num_in_spike = c(2,2,2))
    save(res, settings,
         file = paste0("~/BayesianPMF/03Simulations/Identifiability/id_sim_s2nX_", s2nX, "_s2nY_NULL_spike222_initattruth_FALSE.rda"))
    
    # Update counter
    ind <-  ind + 1
    
    # Remove res to save on memory
    rm(res)
    gc()
#  }
}

# -----------------------------------------------------------------------------
# Run the simulation with no s2n adjustment in Y, initializing at the truth
# -----------------------------------------------------------------------------

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    
    # -----------------------------------------------------------------------------
    # Run the simulation with no s2n adjustment in Y, initializing at the truth
    # -----------------------------------------------------------------------------
    
    res <- identifiability_sim(p.vec, n, ranks, response = "continuous", true_params,
                               model_params, sparsity = TRUE, s2nX = s2nX, s2nY = s2nY,
                               init_at_truth = TRUE, nsim = nsim, nsample = nsample)
    
    ident_res[[ind]] <- res
    
    # Save the results
    settings <- list(s2nX = s2nX, s2nY = NULL, init_at_truth = FALSE, num_in_spike = c(2,2,2))
    save(res, settings,
         file = paste0("~/BayesianPMF/03Simulations/Identifiability/id_sim_s2nX_", s2nX, "_s2nY_NULL_spike222_initattruth_TRUE.rda"))
    
    
    # Update counter
    ind <-  ind + 1
    
    # Remove res to save on memory
    rm(res)
    gc()
  }
}



# -----------------------------------------------------------------------------
# Run the simulation with s2n adjustment in Y, initializing at the truth
# -----------------------------------------------------------------------------

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    
    # -----------------------------------------------------------------------------
    # Run the simulation with no s2n adjustment in Y, without initializing at the truth
    # -----------------------------------------------------------------------------
    
    res <- identifiability_sim(p.vec, n, ranks, response = "continuous", true_params,
                               model_params, sparsity = TRUE, s2nX = s2nX, s2nY = s2nY,
                               init_at_truth = FALSE, nsim = nsim, nsample = nsample)
    
    ident_res[[ind]] <- res
    
    # Save the results
    settings <- list(s2nX = s2nX, s2nY = NULL, init_at_truth = FALSE, num_in_spike = c(2,2,2))
    save(res, settings,
         file = paste0("~/BayesianPMF/03Simulations/Identifiability/id_sim_s2nX_", s2nX, "_s2nY_NULL_spike222_initattruth_FALSE.rda"))
    
    
    # Update counter
    ind <-  ind + 1
    
    # Remove res to save on memory
    rm(res)
    gc()
  }
}

# -----------------------------------------------------------------------------
# Run the simulation with s2n adjustment in Y, without initializing at the truth
# -----------------------------------------------------------------------------

res <- identifiability_sim(p.vec, n, ranks, response = "continuous", true_params,
                           model_params, sparsity = TRUE, s2nX = s2nX, s2nY = s2nY,
                           init_at_truth = FALSE, nsim = nsim, nsample = nsample)

ident_res[[ind]] <- res

# Save the results
settings <- list(s2nX = s2nX, s2nY = s2nY, init_at_truth = FALSE, num_in_spike = c(2,2,2))
save(res, settings,
     file = paste0("~/BayesianPMF/03Simulations/Identifiability/id_sim_s2nX_", s2nX, "_s2nY_", s2nY, "_spike222_initattruth_FALSE.rda"))

# Update counter
ind <-  ind + 1
