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
# s2nX.list <- s2nY.list <- c(0.99/0.01, 0.5/0.5, 0.01/0.99)
s2nX.list <- s2nY.list <- c(1000, 75, 50, 25, 10, 5)

# -----------------------------------------------------------------------------
# Run the simulation with no s2n adjustment in Y
# -----------------------------------------------------------------------------

# Initialize a list to contain the results
ident_res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1 # Index the results

# Initialize at the truth
for (s2nX in s2nX.list) {
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

# Do not initialize at the truth
for (s2nX in s2nX.list) {
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
}

# -----------------------------------------------------------------------------
# Examining the results
# -----------------------------------------------------------------------------

# Initialize a table to contain results
identifiability_results <- data.frame(s2nX = rep(s2nX.list, 2),
                                      init_at_truth = rep(c(TRUE, FALSE), each = 3),
                                      CorrectedSSD = numeric(6),
                                      NonCorrectedSSD = numeric(6),
                                      Diff = numeric(6))

# Iterate through the conditions
for (init_at_truth in c(TRUE, FALSE)) { # Did we initialize at at the true value?
  for (s2nX in s2nX.list) { # What was the s2n for X?
    # Load in the results
    load(paste0("~/BayesianPMF/03Simulations/Identifiability/id_sim_s2nX_", s2nX, "_s2nY_NULL_spike222_initattruth_", init_at_truth,".rda"), verbose = TRUE)
    
    # Compute the SSDs
    corrected_ssd <- mean(res[,1])
    noncorrected_ssd <- mean(res[,2])
    
    # Save the SSDs in the table
    identifiability_results[identifiability_results$s2nX == s2nX & 
                              identifiability_results$init_at_truth == init_at_truth,]$CorrectedSSD <- corrected_ssd
    
    identifiability_results[identifiability_results$s2nX == s2nX & 
                              identifiability_results$init_at_truth == init_at_truth,]$NonCorrectedSSD <- noncorrected_ssd
  }
}

# Calculate the difference
identifiability_results$Diff <- identifiability_results$CorrectedSSD - identifiability_results$NonCorrectedSSD

# Display results
library(xtable)
print(xtable(identifiability_results, digits = 4))
