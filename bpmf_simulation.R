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
# s2nX.list <- s2nY.list <- c(0.99/0.01, 0.5/0.5, 0.01/0.99)

# The parameters for the Bayesian model
joint_var_no_y <- 1/(sqrt(n) + sqrt(sum(p.vec))) # Number of samples and number of features across all sources and Y
indiv_vars_no_y <- sapply(1:q, function(s) 1/(sqrt(n) + sqrt(p.vec[s]))) # Number of samples and number of features in each source plus Y

# Setting the model variances
model_params_bpmf_data <- list(error_vars = c(X1 = 1, X2 = 1), # Error variance for each source
                               joint_var = joint_var_no_y, # Variance for joint structure
                               indiv_vars = indiv_vars_no_y, # Variance for each individual structure
                               beta_vars = c(intercept = 10, joint  = 1, indiv = rep(1,q)), # Variance of intercept effect and each joint effect 
                               response_vars = c(shape = 1, rate = 1) # Hyperparameters for prior on tau2
)   

# When we initialize with y
joint_var_with_y <- 1/(sqrt(n) + sqrt(sum(p.vec) + 1)) # Number of samples and number of features across all sources and Y
indiv_vars_with_y <- sapply(1:q, function(s) 1/(sqrt(n) + sqrt(p.vec[s] + 1))) # Number of samples and number of features in each source plus Y

# Setting the model variances
model_params_bpmf_test <- list(error_vars = c(X1 = 1, X2 = 1), # Error variance for each source
                               joint_var = joint_var_with_y, # Variance for joint structure
                               indiv_vars = indiv_vars_with_y, # Variance for each individual structure
                               beta_vars = c(intercept = 10, joint = joint_var_with_y, indiv = indiv_vars_with_y), # Variance of intercept effect and each joint effect 
                               response_vars = c(shape = 1, rate = 1) # These will not be used
)   

# -----------------------------------------------------------------------------
# test (troubleshooting what happens when I use the true scores in the linear model)
# -----------------------------------------------------------------------------

s2nX <- s2nY <- 99
test.res <- model_comparison(mod = "test", p.vec, n, ranks, response = "continuous", true_params, model_params,
                         s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)

# -----------------------------------------------------------------------------
# sJIVE (Estimated Ranks)
# -----------------------------------------------------------------------------

sJIVE.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    sJIVE.res[[ind]] <- model_comparison(mod = "sJIVE", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_data,
                                     s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

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
# sJIVE (Fixed Ranks)
# -----------------------------------------------------------------------------

sJIVE.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    sJIVE.res[[ind]] <- model_comparison(mod = "sJIVE", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_data,
                                         s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10, estim_ranks = FALSE)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/sJIVE_Fixed_Ranks")
al_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(al_files_split, function(file) (file[5] == s2nX) & (file[7] == s2nY))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# BIDIFAC+ (not sure I will keep this condition)
# -----------------------------------------------------------------------------

BIDIFAC.plus.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BIDIFAC.plus.res[[ind]] <- model_comparison(mod = "BIDIFAC+", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_test,
                                                s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

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
# BIDIFAC
# -----------------------------------------------------------------------------

BIDIFAC.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BIDIFAC.res[[ind]] <- model_comparison(mod = "BIDIFAC", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_data,
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

# -----------------------------------------------------------------------------
# JIVE (Estimated Ranks)
# -----------------------------------------------------------------------------

JIVE.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    JIVE.res[[ind]] <- model_comparison(mod = "JIVE", p.vec = p.vec, n = n, ranks = ranks, response = "continuous", true_params = true_params, model_params = model_params_bpmf_data,
                                       s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

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
# JIVE (Fixed Ranks)
# -----------------------------------------------------------------------------

JIVE.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    JIVE.res[[ind]] <- model_comparison(mod = "JIVE", p.vec = p.vec, n = n, ranks = ranks, response = "continuous", true_params = true_params, model_params = model_params_bpmf_data,
                                        s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10, estim_ranks = FALSE)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/JIVE_Fixed_Ranks")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[5] == s2nX) & (file[7] == s2nY))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos


# -----------------------------------------------------------------------------
# MOFA (Estimated Ranks)
# -----------------------------------------------------------------------------

MOFA.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    MOFA.res[[ind]] <- model_comparison(mod = "MOFA", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_data,
                                    s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/MOFA")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[5] == s2nX) & (file[7] == paste0(s2nY, ".rda")))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# MOFA (Fixed Ranks)
# -----------------------------------------------------------------------------

MOFA.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    MOFA.res[[ind]] <- model_comparison(mod = "MOFA", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_data,
                                        s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10, estim_ranks = FALSE)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/MOFA_Fixed_Ranks")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[5] == s2nX) & (file[7] == s2nY))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# BPMF (Data Mode)
# -----------------------------------------------------------------------------

BPMF.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BPMF.res[[ind]] <- model_comparison(mod = "BPMF_Data_Mode", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_data,
                                        s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/BPMF_Data_Mode")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[7] == s2nX) & (file[9] == paste0(s2nY, ".rda")))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# BPMF (Data Mode) (Fixed Ranks)
# -----------------------------------------------------------------------------

BPMF.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BPMF.res[[ind]] <- model_comparison(mod = "BPMF_Data_Mode", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_data,
                                        s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10, estim_ranks = FALSE)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/BPMF_Data_Mode_Fixed_Ranks")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[7] == s2nX) & (file[9] == s2nY))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# BPMF (TEST - initializing with BIDIFAC+ and Y, no scaling, fixing tau2 at 1)
# -----------------------------------------------------------------------------

# Setting the model variances
# model_params_bpmf_test <- list(error_vars = c(X1 = 1, X2 = 1), # Error variance for each source
#                                joint_var = joint_var_with_y, # Variance for joint structure
#                                indiv_vars = indiv_vars_with_y, # Variance for each individual structure
#                                beta_vars = c(intercept = 10, joint  = 1, indiv = rep(1,q)), # Variance of intercept effect and each joint effect 
#                                response_vars = c(shape = 1, rate = 1) # Hyperparameters for prior on tau2
# )   

BPMF.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BPMF.res[[ind]] <- model_comparison(mod = "BPMF_test", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_test,
                                        s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# -----------------------------------------------------------------------------
# BPMF (TEST - initializing with BIDIFAC+ and Y, scaling, fixing tau2 at 1)
# -----------------------------------------------------------------------------

BPMF.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BPMF.res[[ind]] <- model_comparison(mod = "BPMF_test_scale", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_test,
                                        s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# -----------------------------------------------------------------------------
# BIP (Estimated Ranks)
# -----------------------------------------------------------------------------

BIP.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BIP.res[[ind]] <- model_comparison(mod = "BIP", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_test,
                                        s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/BIP")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[5] == s2nX) & (file[7] == paste0(s2nY, ".rda")))]
    all_s2n[ind] <- length(files_for_s2nX_s2nY) == 100
    combos[ind] <- paste(s2nX, "&", s2nY)
    ind <- ind + 1
  }
}
names(all_s2n) <- combos

# -----------------------------------------------------------------------------
# BIP (Fixed Ranks)
# -----------------------------------------------------------------------------

BIP.res <- lapply(1:(length(s2nX.list) * length(s2nY.list)), function(rep) list())
ind <- 1

for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    BIP.res[[ind]] <- model_comparison(mod = "BIP", p.vec, n, ranks, response = "continuous", true_params, model_params_bpmf_test,
                                       s2nX = s2nX, s2nY = s2nY, sparsity = FALSE, nsim = nsim, nsample = nsample, n_clust = 10, estim_ranks = FALSE)
    ind <- ind + 1
  }
}

# Check that all conditions ran
all_s2n <- c()
combos <- c()
all_files <- list.files("~/BayesianPMF/03Simulations/BIP_Fixed_Ranks")
all_files_split <- strsplit(all_files, split = "_")
ind <- 1
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    # Select all the files corresponding to current s2nX and s2nY 
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[5] == s2nX) & (file[7] == s2nY))]
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
                                 BPMF_Full_Mode = numeric(),
                                 BPMF_Full_Mode_No_Scaling = numeric(),
                                 BPMF_Data_Mode = numeric(),
                                 BPMF_test = numeric(),
                                 BPMF_test_scale = numeric(),
                                 check.names = FALSE)

# Iterate through the s2ns and append to the table
for (s2nX in s2nX.list) {
  for (s2nY in s2nY.list) {
    temp <- data.frame(s2nX = numeric(12),
                       s2nY = numeric(12),
                       Metric = character(12),
                       BPMF_Full_Mode = numeric(12),
                       BPMF_Full_Mode_No_Scaling = numeric(12),
                       BPMF_Data_Mode = numeric(12),
                       BPMF_test = numeric(12),
                       BPMF_test_scale = numeric(12),
                       check.names = FALSE)
    
    # Current results
    res <- create_simulation_table(simulation_results = temp,
                                   mod.list = c("BPMF_Full_Mode", "BPMF_Full_Mode_No_Scaling", "BPMF_Data_Mode", "BPMF_test", "BPMF_test_scale"),
                                   path.list = list(BPMF_Full_Mode = "~/BayesianPMF/03Simulations/BPMF_Full_Mode/",
                                                    BPMF_Full_Mode_No_Scaling = "~/BayesianPMF/03Simulations/BPMF_Full_Mode_No_Scaling/",
                                                    BPMF_Data_Mode = "~/BayesianPMF/03Simulations/BPMF_Data_Mode/",
                                                    BPMF_test = "~/BayesianPMF/03Simulations/BPMF_test/",
                                                    BPMF_test_scale = "~/BayesianPMF/03Simulations/BPMF_test_scale/"),
                                   s2nX = s2nX, 
                                   s2nY = s2nY)
    # Save
    simulation_results <- rbind.data.frame(simulation_results, res)
    
    # Update
    print(paste0("s2nX:", s2nX, ", s2nY:", s2nY))
  }
}

