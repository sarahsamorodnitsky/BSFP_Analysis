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
ranks <- c(5,5,5)
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
s2nX.list <- c(0.9/0.1, 0.75/0.25, 0.5/0.5, 0.25/0.75)

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

# Initialize lists to store the results
BPMF.entrywise <- UNIFAC.entrywise <- mean_imputation.entrywise <- svd_combined.entrywise <-
  svd_separate.entrywise <- knn_combined.entrywise <- knn_separate.entrywise <- 
  rf_combined.entrywise <- rf_separate.entrywise <-
  lapply(1:(length(s2nX.list)), function(rep) list())

BPMF.columnwise <- UNIFAC.columnwise <- mean_imputation.columnwise <- svd_combined.columnwise <-
  svd_separate.columnwise <- knn_combined.columnwise <- knn_separate.columnwise <- 
  rf_combined.columnwise <- rf_separate.columnwise <-
  lapply(1:(length(s2nX.list)), function(rep) list())

BPMF.MNAR <- UNIFAC.MNAR <- mean_imputation.MNAR <- svd_combined.MNAR <-
  svd_separate.MNAR <- knn_combined.MNAR <- knn_separate.MNAR <- 
  rf_combined.MNAR <- rf_separate.MNAR <-
  lapply(1:(length(s2nX.list)), function(rep) list())

# -----------------------------------------------------------------------------
# BSFP (Data Mode) 
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  BPMF.entrywise[[ind]] <- imputation_simulation(mod = "BSFP", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/BSFP", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/BSFP", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list[2]) {
  BPMF.columnwise[[ind]] <- imputation_simulation(mod = "BSFP", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/BSFP", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/BSFP", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  BPMF.MNAR[[ind]] <- imputation_simulation(mod = "BSFP", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/BSFP", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/BSFP", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))

# -----------------------------------------------------------------------------
# UNIFAC
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  UNIFAC.entrywise[[ind]] <- imputation_simulation(mod = "UNIFAC", p.vec = p.vec, n = n, ranks = ranks,
                                                   response = NULL, true_params, model_params = model_params,
                                                   s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                   missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/UNIFAC", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/UNIFAC", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  UNIFAC.columnwise[[ind]] <- imputation_simulation(mod = "UNIFAC", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/UNIFAC", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/UNIFAC", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  UNIFAC.MNAR[[ind]] <- imputation_simulation(mod = "UNIFAC", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/UNIFAC", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/UNIFAC", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))

# -----------------------------------------------------------------------------
# Mean Imputation
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  mean_imputation.entrywise[[ind]] <- imputation_simulation(mod = "Mean_Imputation", p.vec = p.vec, n = n, ranks = ranks,
                                                            response = NULL, true_params, model_params = model_params,
                                                            s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                            missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/Mean_Imputation", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/Mean_Imputation", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  mean_imputation.columnwise[[ind]] <- imputation_simulation(mod = "Mean_Imputation", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/Mean_Imputation", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/Mean_Imputation", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  mean_imputation.MNAR[[ind]] <- imputation_simulation(mod = "Mean_Imputation", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/Mean_Imputation", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/Mean_Imputation", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))

# -----------------------------------------------------------------------------
# SVD Combined Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  svd_combined.entrywise[[ind]] <- imputation_simulation(mod = "SVD_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                            response = NULL, true_params, model_params = model_params,
                                                            s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                            missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/SVD_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/SVD_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  BPMF.columnwise[[ind]] <- imputation_simulation(mod = "SVD_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/SVD_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/SVD_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  svd_combined.MNAR[[ind]] <- imputation_simulation(mod = "SVD_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/SVD_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/SVD_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))


# -----------------------------------------------------------------------------
# SVD Separate Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  svd_separate.entrywise[[ind]] <- imputation_simulation(mod = "SVD_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                         response = NULL, true_params, model_params = model_params,
                                                         s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                         missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/SVD_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/SVD_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness (doesn't work for this model!)
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  svd_separate.columnwise[[ind]] <- imputation_simulation(mod = "SVD_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
# all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/SVD_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
# all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/SVD_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  svd_separate.MNAR[[ind]] <- imputation_simulation(mod = "SVD_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/SVD_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/SVD_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))

# -----------------------------------------------------------------------------
# KNN Combined Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  knn_combined.entrywise[[ind]] <- imputation_simulation(mod = "KNN_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                         response = NULL, true_params, model_params = model_params,
                                                         s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                         missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/KNN_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/KNN_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  knn_combined.columnwise[[ind]] <- imputation_simulation(mod = "KNN_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/KNN_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/KNN_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  knn_combined.MNAR[[ind]] <- imputation_simulation(mod = "KNN_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/KNN_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/KNN_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))

# -----------------------------------------------------------------------------
# KNN Separate Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  knn_separate.entrywise[[ind]] <- imputation_simulation(mod = "KNN_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                         response = NULL, true_params, model_params = model_params,
                                                         s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                         missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/KNN_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/KNN_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list[-1]) {
  knn_separate.columnwise[[ind]] <- imputation_simulation(mod = "KNN_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/KNN_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/KNN_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  knn_separate.MNAR[[ind]] <- imputation_simulation(mod = "KNN_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/KNN_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/KNN_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))

# -----------------------------------------------------------------------------
# RF Combined Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  rf_combined.entrywise[[ind]] <- imputation_simulation(mod = "RF_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                        response = NULL, true_params, model_params = model_params,
                                                        s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                        missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/RF_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/RF_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  rf_combined.columnwise[[ind]] <- imputation_simulation(mod = "RF_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/RF_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/RF_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  rf_combined.MNAR[[ind]] <- imputation_simulation(mod = "RF_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/RF_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/RF_Combined_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))

# -----------------------------------------------------------------------------
# RF Separate Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  rf_separate.entrywise[[ind]] <- imputation_simulation(mod = "RF_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                         response = NULL, true_params, model_params = model_params,
                                                         s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                         missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/RF_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/RF_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "entrywise")))

# -------------------------------------
# Columnwise missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  rf_separate.columnwise[[ind]] <- imputation_simulation(mod = "RF_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks,
                                                 response = NULL, true_params, model_params = model_params,
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust,
                                                 missing_data_type = "columnwise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/RF_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/RF_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "columnwise")))

# -------------------------------------
# MNAR missingness
# -------------------------------------

ind <- 1

for (s2nX in s2nX.list) {
  rf_separate.MNAR[[ind]] <- imputation_simulation(mod = "RF_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "MNAR", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank3/RF_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))
all(sapply(s2nX.list, function(s2nX) check_all_sims(path = "~/BayesianPMF/03Simulations/Imputation_Rank15/RF_Separate_Sources", s2nX = s2nX, nsim = nsim, missing_data_type = "MNAR")))






