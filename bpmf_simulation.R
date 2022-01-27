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
ranks <- c(1,1,1)

# Number of predictors for each source
p.vec = c(100, 100)

# Sample size 
n = 100

# Parameters for data generation
true_params <- model_params <- list(error_vars = c(1,1), # Error variance for each source
                                    joint_var = 1, # Variance for joint structure
                                    indiv_vars = c(1,1), # Variance for each individual structure
                                    beta_vars = c(10, 1, rep(1, q)), # Variance of intercept effect and each joint effect 
                                    response_vars = c(shape = 1, rate = 1)) # Hyperparameters for variance of response (if continuous)

# Signal-to-noise ratios to consider
s2nX.list <- s2nY.list <- c(0.99/0.01, 0.9/0.1, 0.75/0.25, 0.5/0.5, 0.25/0.75, 0.1/0.9, 0.01/0.99)

# -----------------------------------------------------------------------------
# sJIVE
# -----------------------------------------------------------------------------

start <- Sys.time()

sJIVE.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    sJIVE.res[[ind]] <- run_each_mod(mod = "sJIVE", p.vec, n, ranks, response = "continuous", true_params, 
                                     s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = 100, nsample = 1000, n_clust = 10)
    ind <- ind + 1
  }
}

end <- Sys.time()
end-start

# Sparse condition
sJIVE.res[[ind]] <- run_each_mod(mod = "sJIVE", p.vec, n, ranks, response = "continuous", true_params, 
                                 s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = 100, nsample = 1000, n_clust = 10)

# -----------------------------------------------------------------------------
# BIDIFAC+
# -----------------------------------------------------------------------------

BIDIFAC.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BIDIFAC.res[[ind]] <- run_each_mod(mod = "BIDIFAC+", p.vec, n, ranks, response = "continuous", true_params, 
                                     s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = 100, nsample = 1000, n_clust = 10)
    ind <- ind + 1
  }
}

# Sparse condition
BIDIFAC.res[[ind]] <- run_each_mod(mod = "BIDIFAC+", p.vec, n, ranks, response = "continuous", true_params, 
                                 s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = 100, nsample = 1000, n_clust = 10)

# -----------------------------------------------------------------------------
# JIVE
# -----------------------------------------------------------------------------

JIVE.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    JIVE.res[[ind]] <- run_each_mod(mod = "JIVE", p.vec, n, ranks, response = "continuous", true_params, 
                                       s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = 100, nsample = 1000, n_clust = 10)
    ind <- ind + 1
  }
}

# Sparse condition
JIVE.res[[ind]] <- run_each_mod(mod = "BIDIFAC+", p.vec, n, ranks, response = "continuous", true_params, 
                                s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = 100, nsample = 1000, n_clust = 10)

# -----------------------------------------------------------------------------
# MOFA
# -----------------------------------------------------------------------------

MOFA.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    MOFA.res[[ind]] <- run_each_mod(mod = "MOFA", p.vec, n, ranks, response = "continuous", true_params, 
                                    s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = 100, nsample = 1000, n_clust = 10)
    ind <- ind + 1
  }
}

# Sparse condition
MOFA.res[[ind]] <- run_each_mod(mod = "MOFA", p.vec, n, ranks, response = "continuous", true_params, 
                                s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = 100, nsample = 1000, n_clust = 10)

# -----------------------------------------------------------------------------
# BPMF
# -----------------------------------------------------------------------------

BPMF.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BPMF.res[[ind]] <- run_each_mod(mod = "BPMF", p.vec, n, ranks, response = "continuous", true_params, 
                                    s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = 100, nsample = 1000, n_clust = 10)
    ind <- ind + 1
  }
}

# Sparse condition
BPMF.res[[ind]] <- run_each_mod(mod = "BPMF", p.vec, n, ranks, response = "continuous", true_params, 
                                s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = 100, nsample = 1000, n_clust = 10)