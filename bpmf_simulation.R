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

# Signal-to-noise ratios to consider
# s2nX.list <- s2nY.list <- c(0.99/0.01, 0.9/0.1, 0.75/0.25, 0.5/0.5, 0.25/0.75, 0.1/0.9, 0.01/0.99)
s2nX.list <- s2nY.list <- c(0.99/0.01, 0.5/0.5, 0.01/0.99)


# -----------------------------------------------------------------------------
# test
# -----------------------------------------------------------------------------

s2nX <- s2nY <- 99
test.res <- run_each_mod(mod = "test", p.vec, n, ranks, response = "continuous", true_params, model_params,
                         s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)

# -----------------------------------------------------------------------------
# sJIVE
# -----------------------------------------------------------------------------

start <- Sys.time()

sJIVE.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list[2:7]) {
  for (s2nY in s2nY.list) {
    if (s2nX == 9) {
      if (s2nY != 3) {
        sJIVE.res[[ind]] <- run_each_mod(mod = "sJIVE", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                         s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
        ind <- ind + 1
      }
    }
    if (s2nX != 9) {
      sJIVE.res[[ind]] <- run_each_mod(mod = "sJIVE", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                       s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
      ind <- ind + 1
    }
  }
}

end <- Sys.time()
end-start

# Sparse condition
sJIVE.res[[ind]] <- run_each_mod(mod = "sJIVE", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                 s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = nsim, nsample = nsample, n_clust = 10)

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/sJIVE")
al_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(al_files_split, function(file) (file[5] == s2nX) & (file[7] == paste0(s2nY, ".rda")))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# BIDIFAC+
# -----------------------------------------------------------------------------

BIDIFAC.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BIDIFAC.res[[ind]] <- run_each_mod(mod = "BIDIFAC+", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                     s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# Sparse condition
BIDIFAC.res[[ind]] <- run_each_mod(mod = "BIDIFAC+", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                 s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = nsim, nsample = nsample, n_clust = 10)

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/BIDIFAC+")
al_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(al_files_split, function(file) (file[5] == s2nX) & (file[7] == paste0(s2nY, ".rda")))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# JIVE
# -----------------------------------------------------------------------------

JIVE.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    JIVE.res[[ind]] <- run_each_mod(mod = "JIVE", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                       s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# Sparse condition
JIVE.res[[ind]] <- run_each_mod(mod = "JIVE", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = nsim, nsample = nsample, n_clust = 10)

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/JIVE")
al_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(al_files_split, function(file) (file[5] == s2nX) & (file[7] == paste0(s2nY, ".rda")))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# MOFA
# -----------------------------------------------------------------------------

MOFA.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    MOFA.res[[ind]] <- run_each_mod(mod = "MOFA", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                    s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# Sparse condition
MOFA.res[[ind]] <- run_each_mod(mod = "MOFA", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = nsim, nsample = nsample, n_clust = 10)

# -----------------------------------------------------------------------------
# BPMF
# -----------------------------------------------------------------------------

BPMF.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BPMF.res[[ind]] <- run_each_mod(mod = "BPMF", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                    s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# Sparse condition
BPMF.res[[ind]] <- run_each_mod(mod = "BPMF", p.vec, n, ranks, response = "continuous", true_params, model_params,
                                s2nX = NULL, s2nY = NULL, sparsity = TRUE, nsim = nsim, nsample = nsample, n_clust = 10)

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/BPMF")
al_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(al_files_split, function(file) (file[5] == s2nX) & (file[7] == paste0(s2nY, ".rda")))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos


# -----------------------------------------------------------------------------
# Results
# -----------------------------------------------------------------------------

# Set up the table 
simulation_results <- data.frame(s2nX = numeric(),
                                 s2nY = numeric(),
                                 Metric = character(),
                                 BPMF = numeric(),
                                 `BIDIFAC+` = numeric(),
                                 check.names = FALSE)

# Iterate through the s2ns and append to the table
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Current results
    res <- create_simulation_table(mod.list = c("BPMF", "BIDIFAC+"),
                                   path.list = list(BPMF = "~/BayesianPMF/03Simulations/BPMF/",
                                                    `BIDIFAC+` = "~/BayesianPMF/03Simulations/BIDIFAC+/"),
                                   s2nX = s2nX, 
                                   s2nY = s2nY)
    # Save
    simulation_results <- rbind.data.frame(simulation_results, res)
  }
}

