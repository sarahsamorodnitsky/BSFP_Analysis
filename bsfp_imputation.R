# -----------------------------------------------------------------------------
# Comparing the performance of BPMF to other approaches to imputation of 
# missing values in omics data
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set up
# -----------------------------------------------------------------------------

# Load in helper functions
source("~/BSFP/bpmf.R")

# Number of data sources
q <- 2

# Joint and individual ranks (ranks[1] = joint, ranks[2:r_total] = individual)
ranks <- c(1,1,1)
names(ranks) <- c("joint", paste("indiv", 1:q))

# Number of predictors for each source
p.vec = c(100, 100)

# Sample size 
n = 100

# Parameters for data generation
true_params <- list(error_vars = c(X1 = 1, X2 = 1, Y = 1), # Error variance for each source (all have error 1)
                    joint_var = 1, # Prior variance for joint structure
                    indiv_vars = c(1,1), #  Prior variance for each individual structure
                    beta_vars = c(intercept = 10, joint_effect = 1, indiv_effects = rep(1, q)) # Variance on the effect of each factor on y including the intercept
)

# Parameters for the simulation
nsample <- 2000
nsim <- 100

# Signal-to-noise ratios to consider
s2nX.list <- s2nY.list <- c(0.9/0.1, 0.75/0.25, 0.5/0.5, 0.25/0.75)

# The parameters for the Bayesian model
joint_var_no_y <- 1/(sqrt(n) + sqrt(sum(p.vec))) # Number of samples and number of features across all sources and Y
indiv_vars_no_y <- sapply(1:q, function(s) 1/(sqrt(n) + sqrt(p.vec[s]))) # Number of samples and number of features in each source plus Y

# Setting the model variances
model_params <- list(error_vars = c(X1 = 1, X2 = 1), # Error variance for each source
                    joint_var = joint_var_no_y, # Variance for joint structure
                    indiv_vars = indiv_vars_no_y, # Variance for each individual structure
                    beta_vars = c(intercept = 10, joint  = 1, indiv = rep(1,q)), # Variance of intercept effect and each joint effect 
                    response_vars = c(shape = 1, rate = 1) # Hyperparameters for prior on tau2
)   

# The missingness proportion
prop_missing <- 0.1

# Model simulation parameters 
n_clust <- 10

# -----------------------------------------------------------------------------
# BPMF (Data Mode) 
# -----------------------------------------------------------------------------

# Entrywise missingness
BPMF.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BPMF.entrywise[[ind]] <- imputation_simulation(mod = "BSFP", p.vec = p.vec, n = n, ranks = ranks, response = NULL, true_params, model_params = model_params, s2nX = s2nX, s2nY = s2nY, nsim = nsim, nsample = 2000, n_clust = n_clust, missing_data_type = "entrywise", prop_missing = prop_missing)  
    ind <- ind + 1
  }
}


# -----------------------------------------------------------------------------
# BIDIFAC
# -----------------------------------------------------------------------------

# Entrywise missingness
BIDIFAC.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BIDIFAC.entrywise[[ind]] <- model_comparison(mod = "UNIFAC", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_data,
                                           s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/BIDIFAC")
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
