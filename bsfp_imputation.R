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

# -----------------------------------------------------------------------------
# BPMF (Data Mode) 
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

BPMF.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  BPMF.entrywise[[ind]] <- imputation_simulation(mod = "BSFP", p.vec = p.vec, n = n, ranks = ranks, 
                                                 response = NULL, true_params, model_params = model_params, 
                                                 s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                 missing_data_type = "entrywise", prop_missing = prop_missing)  
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/BSFP")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  # Select all the files corresponding to current s2nX and s2nY 
  files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) file[6] == s2nX)]
  all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
  combos[ind] <- paste(s2nX, "&", s2nY)
  ind <- ind + 1
}
names(all_s2n) <- combos
all(all_s2n)

# -----------------------------------------------------------------------------
# UNIFAC
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

UNIFAC.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  UNIFAC.entrywise[[ind]] <- imputation_simulation(mod = "UNIFAC", p.vec = p.vec, n = n, ranks = ranks, 
                                                   response = NULL, true_params, model_params = model_params, 
                                                   s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                   missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/UNIFAC")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) file[6] == s2nX))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos
all(all_s2n)

# -----------------------------------------------------------------------------
# Mean Imputation
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

mean_imputation.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  mean_imputation.entrywise[[ind]] <- imputation_simulation(mod = "Mean_Imputation", p.vec = p.vec, n = n, ranks = ranks, 
                                                            response = NULL, true_params, model_params = model_params, 
                                                            s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                            missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/Mean_Imputation")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) file[6] == s2nX))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos
all(all_s2n)

# -----------------------------------------------------------------------------
# SVD Combined Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

svd_combined.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  svd_combined.entrywise[[ind]] <- imputation_simulation(mod = "SVD_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                            response = NULL, true_params, model_params = model_params, 
                                                            s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                            missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/SVD_Combined_Sources")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) file[6] == s2nX))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos
all(all_s2n)

# -----------------------------------------------------------------------------
# SVD Separate Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

svd_combined.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  svd_combined.entrywise[[ind]] <- imputation_simulation(mod = "SVD_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                         response = NULL, true_params, model_params = model_params, 
                                                         s2nX = s2nX, s2nY = s2nY, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                         missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/SVD_Separate_Sources")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) file[6] == s2nX))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos
all(all_s2n)


# -----------------------------------------------------------------------------
# KNN Combined Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

knn_combined.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  knn_combined.entrywise[[ind]] <- imputation_simulation(mod = "KNN_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                         response = NULL, true_params, model_params = model_params, 
                                                         s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                         missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/KNN_Combined_Sources")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) file[6] == s2nX))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos
all(all_s2n)

# -----------------------------------------------------------------------------
# KNN Separate Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

knn_separate.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  knn_separate.entrywise[[ind]] <- imputation_simulation(mod = "KNN_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                         response = NULL, true_params, model_params = model_params, 
                                                         s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                         missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/KNN_Combined_Sources")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) file[6] == s2nX))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos
all(all_s2n)

# -----------------------------------------------------------------------------
# RF Combined Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

rf_combined.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  rf_combined.entrywise[[ind]] <- imputation_simulation(mod = "RF_Combined_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                        response = NULL, true_params, model_params = model_params, 
                                                        s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                        missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/RF_Combined_Sources")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) function(file) file[6] == s2nX))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos
all(all_s2n)

# -----------------------------------------------------------------------------
# RF Separate Sources
# -----------------------------------------------------------------------------

# -------------------------------------
# Entrywise missingness
# -------------------------------------

rf_separate.entrywise <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  rf_separate.entrywise[[ind]] <- imputation_simulation(mod = "RF_Separate_Sources", p.vec = p.vec, n = n, ranks = ranks, 
                                                         response = NULL, true_params, model_params = model_params, 
                                                         s2nX = s2nX, s2nY = NULL, nsim = nsim, nsample = 2000, n_clust = n_clust, 
                                                         missing_data_type = "entrywise", prop_missing = prop_missing)
  ind <- ind + 1
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/Imputation/RF_Separate_Sources")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) function(file) file[6] == s2nX))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos
all(all_s2n)







