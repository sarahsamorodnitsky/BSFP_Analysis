# -----------------------------------------------------------------------------
# Applying the BPMF model to the HIV-COPD BALF Biocrates and Somascan data. 
# Comparing the model application when there is no sparsity vs. when
# there is sparsity. 
# * Note on V2: This script runs the data application again with centered and 
#   data NOT scaled to have overall variance 1. 
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set-up
# -----------------------------------------------------------------------------

# Packages
library(doParallel)
library(foreach)

# Working directories
data_wd <- "~/BSFP_Analysis/data/"
model_wd <- "~/BSFP_Analysis/"
results_wd <- "~/BayesianPMF/04DataApplication/"

# Loading in the model functions
source(paste0(model_wd, "bpmf.R"))

# -----------------------------------------------------------------------------
# Preparing the data
# -----------------------------------------------------------------------------

# Loading data in (must load in separately on 22 server)
load(paste0(data_wd, "BiocratesLavageProcessed.rda"))
load(paste0(data_wd, "HIV_COPD_SomaScan_Normalized_Clean.rda"))

# BALF Biocrates
lavage_processed_no_info <- lavage_processed[, -c(1:2)] # removing meta
bioc_meta_data <- lavage_processed[, 1:2] # saving meta elsewhere
lavage_processed_no_info_log_scale <- t(apply(lavage_processed_no_info, 1, function(row) {
  scale(log(row + 1), center = TRUE, scale = FALSE)
})) # scaling the rows
rownames(lavage_processed_no_info_log_scale) <- bioc_meta_data$Metabolite
colnames(lavage_processed_no_info_log_scale) <- colnames(lavage_processed_no_info)

# BALF Somascan
somascan_normalized_clean_no_info <- somascan_normalized_clean[, -c(1:29)] # removing meta
soma_meta_data <- somascan_normalized_clean[, 1:29] # saving meta elsewhere
soma_pids <- somascan_normalized_clean$PID
somascan_normalized_clean_no_info_transpose <- t(somascan_normalized_clean_no_info)
somascan_normalized_clean_no_info_transpose_scale <- t(apply(somascan_normalized_clean_no_info_transpose, 1, function(row) {
  scale(log(row + 1), center = TRUE, scale = FALSE)
})) # scaling the rows
colnames(somascan_normalized_clean_no_info_transpose_scale) <- soma_pids

# Reorder the Somascan cols to match the Biocrates cols
somascan_normalized_clean_no_info_transpose_reorder <- 
  somascan_normalized_clean_no_info_transpose_scale[, match(colnames(lavage_processed_no_info_log_scale),
                                                            colnames(somascan_normalized_clean_no_info_transpose_scale))]

# Check
all(colnames(somascan_normalized_clean_no_info_transpose_reorder) == colnames(lavage_processed_no_info_log_scale))
dim(somascan_normalized_clean_no_info_transpose_reorder)
dim(lavage_processed_no_info)

# Combining the sources together
hiv_copd_data <- matrix(list(), nrow = 2, ncol = 1)
hiv_copd_data[[1,1]] <- lavage_processed_no_info_log_scale
hiv_copd_data[[2,1]] <- somascan_normalized_clean_no_info_transpose_reorder

# Creating a list of the two sources
hiv_copd_data_list <- list(Metabolomics = lavage_processed_no_info_log_scale, Proteomics = somascan_normalized_clean_no_info_transpose_reorder)

# Preparing the response variable
fev1pp <- matrix(list(), nrow = 1, ncol = 1)
fev1pp[[1,1]] <- matrix(subject_processed$FEV1_percent_predicted, ncol = 1)

# -----------------------------------------------------------------------------
# Model parameters
# -----------------------------------------------------------------------------

# Data dimensions
p1 <- nrow(lavage_processed_no_info_log_scale)
p2 <- nrow(somascan_normalized_clean_no_info_transpose_reorder)
p <- p1 + p2
p.vec <- c(p1, p2)
n <- ncol(lavage_processed_no_info_log_scale)
q <- nrow(hiv_copd_data)

# Error variances
sigma21 <- 1 # error variance for Biocrates
sigma22 <- 1 # error variance for Somascan

# Variance of joint structure
lambda_joint <- sqrt(nrow(lavage_processed_no_info_log_scale) + nrow(somascan_normalized_clean_no_info_transpose_reorder)) + 
  sqrt(ncol(lavage_processed_no_info_log_scale))
sigma2_joint <- 1/(lambda_joint)

# Variance of individual structure for Biocrates
lambda_indiv1 <- sqrt(nrow(lavage_processed_no_info_log_scale)) + sqrt(ncol(lavage_processed_no_info_log_scale))
sigma2_indiv1 <- 1/(lambda_indiv1) 

# Variance of individual structure for Somascan
lambda_indiv2 <- sqrt(nrow(somascan_normalized_clean_no_info_transpose_reorder)) + sqrt(ncol(lavage_processed_no_info_log_scale))
sigma2_indiv2 <- 1/(lambda_indiv2)

# For the regression coefficients, beta
lambda2_intercept <- 1e6
lambda2_joint <- 1
lambda2_indiv1 <- 1
lambda2_indiv2 <- 1
beta_vars <- c(lambda2_intercept, lambda2_joint, lambda2_indiv1, lambda2_indiv2)

# For the response vector
shape <- 0.01
rate <- 0.01

# Putting the model parameters together
model_params <- list(error_vars = c(sigma21, sigma22),
                     joint_var = sigma2_joint,
                     indiv_vars = c(sigma2_indiv1, sigma2_indiv2),
                     beta_vars = beta_vars,
                     response_vars = c(shape = shape, rate = rate))

# Gibbs sampler parameters
nsample <- 5000
burnin <- 1000
nsample_after_burnin <- nsample-burnin
iters_burnin <- seq(burnin+1, nsample)
thinned_iters <- seq(1, nsample, by = 10)
thinned_iters_burnin <- seq(burnin, nsample, by = 10)

# -----------------------------------------------------------------------------
# Training Data Model Fit (BPMF, BIDIFAC, sJIVE, MOFA, BIP)
# -----------------------------------------------------------------------------

# Fitting BPMF
fev1pp_training_fit_nonsparse_V2 <- bpmf_data_mode(
  data = hiv_copd_data,
  Y = fev1pp,
  nninit = TRUE,
  model_params = model_params,
  sparsity = FALSE,
  nsample = nsample,
  progress = TRUE
)

# Save the results
save(fev1pp_training_fit_nonsparse_V2, file = paste0(results_wd, "BPMF/Training_Fit/training_data_fit_V2.rda"))

# Fitting BIDIFAC
BIDIFAC_training_fit <- BIDIFAC(hiv_copd_data, rmt = TRUE, pbar = FALSE, scale_back = TRUE)
save(BIDIFAC_training_fit, file = paste0(results_wd, "/BIDIFAC/BIDIFAC_training_data_fit.rda"))
  
# Fitting sJIVE
eta <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
sJIVE_training_fit <- sJIVE(X = hiv_copd_data_list, Y = c(fev1pp[[1,1]]), eta = eta, rankA = NULL, rankJ = NULL, method = "permute", threshold = 0.001, center.scale = FALSE, reduce.dim = TRUE)
save(sJIVE_training_fit, file = paste0(results_wd, "/sJIVE/sJIVE_training_data_fit.rda"))

# Fitting JIVE
JIVE_training_fit <- jive(hiv_copd_data_list, center = FALSE, scale = FALSE, method = "perm")
save(JIVE_training_fit, file = paste0(results_wd, "/JIVE/JIVE_training_data_fit.rda"))

# Fitting MOFA
mofa_pre_train <- create_mofa(hiv_copd_data_list)

# Set the data options so that the data is not additionally centered and fix the ranks if desired
data_opts <- get_default_data_options(mofa_pre_train)
model_opts <- get_default_model_options(mofa_pre_train)

# No centering
data_opts$center_groups <- FALSE

# If using fixed ranks
if (!estim_ranks) {
  model_opts$num_factors <- sum(ranks)
}

# Create the MOFA object
MOFAobject <- prepare_mofa(
  object = mofa_pre_train,
  model_options = model_opts
)

# Train the MOFA model
MOFA_training_fit <- run_mofa(MOFAobject)
save(MOFA_training_fit, file = paste0(results_wd, "/MOFA/MOFA_training_data_fit.rda"))

# Fitting BIP
library(BIPnet)

# Transposing the data matrices
hiv_copd_data_list_bip <- lapply(hiv_copd_data_list, function(data) t(data))

# Scaling the data to have sd 1
hiv_copd_data_list_bip <- lapply(1:q, function(s) scale(hiv_copd_data_list_bip[[s]], center = FALSE, scale = TRUE))

# Adding the response vector
hiv_copd_data_list_bip[[3]] <- fev1pp[[1,1]]

# Fitting the model
BIP_training_fit <- BIP(dataList = hiv_copd_data_list_bip, IndicVar = c(0,0,1), Method = "BIP", sample = burnin, burnin = burnin-1, nbrcomp = 35)

# -----------------------------------------------------------------------------
# Assessing convergence of BPMF
# -----------------------------------------------------------------------------

# Assessing convergence for both model fits
load(paste0(results_wd, "BPMF/Training_Fit/training_data_fit_V2.rda"), verbose = TRUE)

fev1pp_training_nonsparse_conv_V2 <- sapply(1:nsample, function(sim_iter) {
  # Calculate the log-joint density at each thinned iterations
  log_joint_density(data = hiv_copd_data, 
                    U.iter = fev1pp_training_fit_nonsparse_V2$U.draw[[sim_iter]], 
                    V.iter = fev1pp_training_fit_nonsparse_V2$V.draw[[sim_iter]], 
                    W.iter = fev1pp_training_fit_nonsparse_V2$W.draw[[sim_iter]], 
                    Vs.iter = fev1pp_training_fit_nonsparse_V2$Vs.draw[[sim_iter]],
                    model_params = model_params,
                    ranks = fev1pp_training_fit_nonsparse_V2$ranks,
                    Y = fev1pp,
                    beta.iter = fev1pp_training_fit_nonsparse_V2$beta.draw[[sim_iter]],
                    tau2.iter = fev1pp_training_fit_nonsparse_V2$tau2.draw[[sim_iter]])
})

# Save the log-joint density 
# save(fev1pp_training_nonsparse_conv_V2, file = "~/BayesianPMF/04DataApplication/BPMF/Training_Fit/Log_Joint_Density_NonSparse_V2.rda")

# Plotting the log-joint densities
load("~/BayesianPMF/04DataApplication/BPMF/Training_Fit/Log_Joint_Density_NonSparse.rda", verbose = TRUE)
load("~/BayesianPMF/04DataApplication/BPMF/Training_Fit/Log_Joint_Density_NonSparse_V2.rda", verbose = TRUE)

# Calculate the ranges for the convergences
ymin <- min(range(fev1pp_training_nonsparse_conv), range(fev1pp_training_nonsparse_conv_V2))
ymax <- max(range(fev1pp_training_nonsparse_conv), range(fev1pp_training_nonsparse_conv_V2))

png("~/BayesianPMF/04DataApplication/BPMF/Training_Fit/log_joint_density_plot_all_iters_V2.png")
plot(fev1pp_training_nonsparse_conv, ylab = "Log Joint Density", main = "Log Joint Density Across Gibbs Sampling Iterations",
     ylim = c(ymin, ymax))
points(fev1pp_training_nonsparse_conv_V2, col = 2)
dev.off()

# Calculate the ranges for the burnin
ymin.burnin <- min(range(fev1pp_training_nonsparse_conv[burnin:nsample]), range(fev1pp_training_nonsparse_conv_V2[burnin:nsample]))
ymax.burnin <- max(range(fev1pp_training_nonsparse_conv[burnin:nsample]), range(fev1pp_training_nonsparse_conv_V2[burnin:nsample]))

png("~/BayesianPMF/04DataApplication/BPMF/Training_Fit/log_joint_density_plot_burnin_V2.png")
plot(fev1pp_training_nonsparse_conv[burnin:nsample], ylab = "Log Joint Density", main = "Log Joint Density After Burn-in",
     ylim = c(ymin.burnin, ymax.burnin))
points(fev1pp_training_nonsparse_conv_V2[burnin:nsample], col = 2)
dev.off()

# -----------------------------------------------------------------------------
# Plotting a random selection of loadings and scores to check for mode-switching
# * Includes a trace plot of the underlying structure
# -----------------------------------------------------------------------------

set.seed(1)

# Randomly sample loadings
ranks <- fev1pp_training_fit_nonsparse_V2$ranks
rand.factor <- lapply(1:(q+1), function(s) sample(1:ranks[[s]], size = 1))
rand.load <- lapply(1:q, function(s) sample(1:p.vec[s], size = 3, replace = FALSE))

# Randomply sample scores
rand.score <- sample(1:n, size = 2, replace = FALSE)

# Plot
pdf("~/BayesianPMF/04DataApplication/BPMF/Figures/Unaligned_Trace_Plots_Structure_Loadings_Scores.pdf")

# Plot the structure
par(mfrow = c(3,2))
plot(sapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$J.draw[[iter]][[1,1]][rand.load[[1]][1], rand.score[[1]]]), ylab = "", main = "Joint Structure")
plot(sapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$J.draw[[iter]][[2,1]][rand.load[[2]][1], rand.score[[2]]]), ylab = "", main = "Joint Structure")
plot(sapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$A.draw[[iter]][[1,1]][rand.load[[1]][2], rand.score[[1]]]), ylab = "", main = "Metabolite Structure")
plot(sapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$A.draw[[iter]][[1,1]][rand.load[[1]][3], rand.score[[2]]]), ylab = "", main = "Metabolite Structure")
plot(sapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$A.draw[[iter]][[2,1]][rand.load[[2]][2], rand.score[[1]]]), ylab = "", main = "Protein Structure")
plot(sapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$A.draw[[iter]][[2,1]][rand.load[[2]][3], rand.score[[2]]]), ylab = "", main = "Protein Structure")
par(mfrow = c(1,1))

# Plot the Loadings
par(mfrow = c(3,2))
for (s in 1:q) {
  for (load in rand.load[[s]]) {
    plot(sapply(iters_burnin, function(iter) {
      fev1pp_training_fit_nonsparse_V2$U.draw[[iter]][[s,1]][load, rand.factor[[1]]]
    }), ylab = "", main = "Joint Loading")
  }
}
par(mfrow = c(1,1))
par(mfrow = c(3,1))
for (load in rand.load[[1]]) {
  plot(sapply(iters_burnin, function(iter) {
    fev1pp_training_fit_nonsparse_V2$W.draw[[iter]][[1,1]][load, rand.factor[[2]]]
  }), ylab = "", main = "Metabolite Loading")
}
par(mfrow = c(1,1))
par(mfrow = c(3,1))
for (load in rand.load[[2]]) {
  plot(sapply(iters_burnin, function(iter) {
    fev1pp_training_fit_nonsparse_V2$W.draw[[iter]][[2,2]][load, rand.factor[[3]]]
  }), ylab = "", main = "Protein Loading")
}
par(mfrow = c(1,1))

# Plot the scores
par(mfrow = c(3,2))
plot(sapply(iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_V2$V.draw[[iter]][[1,1]][rand.score[1], rand.factor[[1]]]
}), ylab = "", main = "Joint Score")
plot(sapply(iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_V2$V.draw[[iter]][[1,1]][rand.score[2], rand.factor[[1]]]
}), ylab = "", main = "Joint Score")
plot(sapply(iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_V2$Vs.draw[[iter]][[1,1]][rand.score[1], rand.factor[[2]]]
}), ylab = "", main = "Metabolite Score")
plot(sapply(iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_V2$Vs.draw[[iter]][[1,1]][rand.score[2], rand.factor[[2]]]
}), ylab = "", main = "Metabolite Score")
plot(sapply(iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_V2$Vs.draw[[iter]][[1,2]][rand.score[1], rand.factor[[3]]]
}), ylab = "", main = "Protein Score")
plot(sapply(iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_V2$Vs.draw[[iter]][[1,2]][rand.score[2], rand.factor[[3]]]
}), ylab = "", main = "Protein Score")
par(mfrow = c(1,1))

dev.off()

# -----------------------------------------------------------------------------
# Adjusting for factor switching
# * Using the Poworoznek et al. (2021) solution 
# -----------------------------------------------------------------------------

# Loading in the full training data fit
load(paste0(results_wd, "BPMF/Training_Fit/training_data_fit_V2.rda"), verbose = TRUE)

# Save the ranks
ranks <- fev1pp_training_fit_nonsparse_V2$ranks
joint.rank <- fev1pp_training_fit_nonsparse_V2$ranks[1]
indiv.ranks <- fev1pp_training_fit_nonsparse_V2$ranks[-1]

# Applying the Match Align algorithm to undo rotational invariance in the results
fev1pp_training_fit_nonsparse_aligned_V2 <- match_align_bpmf(fev1pp_training_fit_nonsparse_V2, y = fev1pp,
                                                             model_params = model_params, p.vec = p.vec, 
                                                             iters_burnin = iters_burnin)
# Save
joint.scores.final <- fev1pp_training_fit_nonsparse_aligned_V2$joint.scores.final
joint.loadings.final <- fev1pp_training_fit_nonsparse_aligned_V2$joint.loadings.final
joint.betas.final <- fev1pp_training_fit_nonsparse_aligned_V2$joint.betas.final

individual.scores.final <- fev1pp_training_fit_nonsparse_aligned_V2$individual.scores.final
individual.loadings.final <- fev1pp_training_fit_nonsparse_aligned_V2$individual.loadings.final
individual.betas.final <- fev1pp_training_fit_nonsparse_aligned_V2$individual.betas.final

# -------------------------------------
# Check the joint structure after the algorithm matches the original structure 
# -------------------------------------

# Create a vector for the indices of each rank
beta.ind <- lapply(1:length(ranks), function(i) {
  if (i == 1) { 
    (1:ranks[i]) + 1
  } else {
    ((sum(ranks[1:(i-1)])+1):sum(ranks[1:i])) + 1
  }
})

# Save the joint loadings
joint.loadings <- lapply(iters_burnin, function(iter) {
  rbind(do.call(rbind, fev1pp_training_fit_nonsparse_V2$U.draw[[iter]]), # Joint loadings  
        t(fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][beta.ind[[1]],])) # Joint regression coefficients
})

# Save joint scores
joint.scores <- lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$V.draw[[iter]][[1,1]])

# Save the original joint structure
joint.structure.original <- lapply(1:burnin, function(iter) {
  joint.loadings[[iter]] %*% t(joint.scores[[iter]])
})

# Calculate the new structure
joint.structure.new <- lapply(1:burnin, function(iter) {
  rbind(joint.loadings.final[[iter]], t(joint.betas.final[[iter]])) %*%
    t(joint.scores.final[[iter]])
})

# Check
all.equal(joint.structure.original, joint.structure.new) # TRUE!

# Check the individual structure after the algorithm matches the original structure --

# Combine the individual loadings and betas
individual.loadings <- lapply(1:q, function(s) lapply(iters_burnin, function(iter) {
  rbind(fev1pp_training_fit_nonsparse_V2$W.draw[[iter]][[s,s]], 
        fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][beta.ind[[s+1]],])
}))
individual.scores <- lapply(1:q, function(s) lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$Vs.draw[[iter]][[1,s]]))

# Calculate the original structure
individual.structure.original <- lapply(1:q, function(s) {
  lapply(1:burnin, function(iter) {
    individual.loadings[[s]][[iter]] %*% t(individual.scores[[s]][[iter]])
  })
})

# Calculate the new structure
individual.structure.new <- lapply(1:q, function(s) {
  lapply(1:burnin, function(iter) {
    rbind(individual.loadings.final[[s]][[iter]], t(individual.betas.final[[s]][[iter]])) %*%
      t(individual.scores.final[[s]][[iter]])
  })
})

# Check
lapply(1:q, function(s) all.equal(individual.structure.original[[s]], individual.structure.new[[s]])) # TRUE TRUE!

# Save the results
save(joint.scores.final, joint.loadings.final, joint.betas.final, individual.scores.final, individual.loadings.final, individual.betas.final,
     file = "/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Joint_Individual_Structures_Factor_Switching_Ordered_V2_NewPivot_Burnin.rda")

# -------------------------------------
# Run alignment again with different pivot as a test of sensitivity
# -------------------------------------

fev1pp_training_fit_nonsparse_aligned_V2_sensitivity <- lapply(c(-1,1), function(index) {
  match_align_bpmf(fev1pp_training_fit_nonsparse_V2, y = fev1pp,
                   model_params = model_params, p.vec = p.vec, 
                   iters_burnin = iters_burnin, index = index)
})

# Save the results
save(fev1pp_training_fit_nonsparse_aligned_V2_sensitivity,
     file = "/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Factor_Switching_Sensitivity_V2_NewPivot_Burnin.rda")


# -------------------------------------
# Calculate the difference between the posterior mean covariance and estimated covariance
# using aligned loadings. 
# -------------------------------------

p.ind <- lapply(1:q, function(s) {
  if (s == 1) {
    1:p.vec[s]
  } else {
    (p.vec[s-1]+1):(sum(p.vec))
  }
})

# Combining the loadings together into one big matrix
combined.loadings.final <- lapply(1:burnin, function(iter) {
  load <- matrix(list(), nrow = 3, ncol = 3)
  
  # Joint loadings
  load[[1,1]] <- joint.loadings.final[[iter]][p.ind[[1]],]
  load[[2,1]] <- joint.loadings.final[[iter]][p.ind[[2]],]
  load[[3,1]] <- t(joint.betas.final[[iter]])
  
  # Individual loadings
  load[[1,2]] <- individual.loadings.final[[1]][[iter]]
  load[[2,2]] <- matrix(0, nrow = nrow(load[[2,1]]), ncol = ncol(load[[1,2]]))
  load[[3,2]] <- t(individual.betas.final[[1]][[iter]])
  
  load[[2,3]] <- individual.loadings.final[[2]][[iter]] 
  load[[1,3]] <- matrix(0, nrow = nrow(load[[1,2]]), ncol = ncol(load[[2,3]]))
  load[[3,3]] <- t(individual.betas.final[[2]][[iter]])
  
  data.rearrange(load)$out
})

combined.loadings.unaligned <- lapply(iters_burnin, function(iter) {
  load <- matrix(list(), nrow = 3, ncol = 3)
  
  # Joint
  load[[1,1]] <- fev1pp_training_fit_nonsparse_V2$U.draw[[iter]][[1,1]]
  load[[2,1]] <- fev1pp_training_fit_nonsparse_V2$U.draw[[iter]][[2,1]]
  load[[3,1]] <- t(fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][beta.ind[[1]],,drop=FALSE])
  
  # Individual
  load[[1,2]] <- fev1pp_training_fit_nonsparse_V2$W.draw[[iter]][[1,1]]
  load[[2,2]] <- fev1pp_training_fit_nonsparse_V2$W.draw[[iter]][[2,1]]
  load[[3,2]] <- t(fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][beta.ind[[2]],,drop=FALSE])
  
  load[[1,3]] <- fev1pp_training_fit_nonsparse_V2$W.draw[[iter]][[1,2]]
  load[[2,3]] <- fev1pp_training_fit_nonsparse_V2$W.draw[[iter]][[2,2]]
  load[[3,3]] <- t(fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][beta.ind[[3]],,drop=FALSE])
  
  data.rearrange(load)$out
})

# Calculate the posterior mean of the posterior mean of the covariance
post.mean.covariance <- lmean(lapply(1:burnin, function(iter) combined.loadings.unaligned[[iter]] %*% t(combined.loadings.unaligned[[iter]])))

# Calculate the estimated covariance using the posterior mean of the aligned loadings
combined.loadings.final.mean <- lmean(combined.loadings.final)
est.covariance.aligned.load <- 

# -----------------------------------------------------------------------------
# Create PCA-like plots using ALIGNED factors
# -----------------------------------------------------------------------------

library(viridis)

# Load in the clustering information
load("~/HIV-COPD/Data/UnivariatePValuesUncombinedClusterComparison.rda", verbose = TRUE)

# Creating colors for the plots
ccstat_colors <- viridis(2)[clinical_data_soma$ccstat]
cluster_colors <- viridis(3)[1:2][clusters.soma]

# Creating shapes
cluster_shapes <- c(3,4)[clusters.soma]

# Load in the aligned factors
load("/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Joint_Individual_Structures_Factor_Switching_Ordered_V2_NewPivot_Burnin.rda", verbose = TRUE)

# -------------------------------------
# Joint structure
# -------------------------------------

# Calculate the posterior mean of the joint scores
joint.scores.final.mean <- Reduce("+", lapply(1:burnin, function(iter) joint.scores.final[[iter]]))/(nsample/2)

# Plot subjects against every combination of PCs
joint_rank <- fev1pp_training_fit_nonsparse_V2$ranks[1]
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Joint_Structure_Aligned_Ordered_PCA_Plot_V2_NewPivot_Burnin.pdf"), width = 15)
for (rs1 in 1:joint_rank) {
 for (rs2 in 1:joint_rank) {
   if (rs2 > rs1) {
     par(mfrow = c(1,2), mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
     # Plot with case-control label as colors
     plot(joint.scores.final.mean[,rs1], joint.scores.final.mean[,rs2], pch = 16, xlab = paste0("Joint Factor ", rs1),
          ylab = paste0("Joint Factor ", rs2), main = paste0("Joint Factors ", rs1, " vs ", rs2),
          col = ccstat_colors, lwd = 3)
     
     # Plot with stand-out subjects as colors
     plot(joint.scores.final.mean[,rs1], joint.scores.final.mean[,rs2], xlab = paste0("Joint Factor ", rs1),
          ylab = paste0("Joint Factor ", rs2), main = paste0("Joint Factors ", rs1, " vs ", rs2),
          col = cluster_colors, pch = cluster_shapes, lwd = 3)
     legend("topleft", col = c("#440154FF", "#21908CFF"), pch = c(3,4),
            legend = c("Cluster 1", "Cluster 2"))
     par(mfrow = c(1,1))
   }
 } 
}
dev.off()

# -------------------------------------
# Individual structure
# -------------------------------------

# Saving the reordered individual factors, ordered by variance explained in each source
individual.scores.final.mean <- lapply(1:q, function(s) Reduce("+", lapply(1:burnin, function(iter) individual.scores.final[[s]][[iter]]))/(nsample/2))

# Plot subjects against every combination of PCs for Biocrates
indiv_rank1 <- fev1pp_training_fit_nonsparse_V2$ranks[2]
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Individual_Structure_Biocrates_Aligned_Ordered_PCA_Plot_V2_NewPivot_Burnin.pdf"), width = 15)
for (rs1 in 1:indiv_rank1) {
  for (rs2 in 1:indiv_rank1) {
    if (rs2 > rs1) {
      par(mfrow = c(1,2))
      # Plot with case-control label as colors
      plot(individual.scores.final.mean[[1]][,rs1], individual.scores.final.mean[[1]][,rs2], pch = 16, xlab = paste0("Individual Factor ", rs1),
           ylab = paste0("Individual Factor ", rs2), main = paste0("Individual Factors Metabolomic ", rs1, " vs ", rs2),
           col = ccstat_colors)
      
      # Plot with stand-out subjects as colors
      plot(individual.scores.final.mean[[1]][,rs1], individual.scores.final.mean[[1]][,rs2], xlab = paste0("Individual Factor ", rs1),
           ylab = paste0("Individual Factor ", rs2), main = paste0("Individual Factors Metabolomic ", rs1, " vs ", rs2),
           col = cluster_colors, pch = cluster_shapes, lwd = 3)
      legend("bottomleft", col = c("#440154FF", "#21908CFF"), pch = c(3,4),
             legend = c("Cluster 1", "Cluster 2"))
      par(mfrow = c(1,1))
    }
  } 
}
dev.off()

# Plot subjects against every combination of PCs for Somascan
indiv_rank2 <- fev1pp_training_fit_nonsparse_V2$ranks[3]
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Individual_Structure_Somascan_Aligned_Ordered_PCA_Plot_V2_NewPivot_Burnin.pdf"), width = 15)
for (rs1 in 1:indiv_rank2) {
  for (rs2 in 1:indiv_rank2) {
    if (rs2 > rs1) {
      par(mfrow = c(1,2))
      # Plot with case-control label as colors
      plot(individual.scores.final.mean[[2]][,rs1], individual.scores.final.mean[[2]][,rs2], pch = 16, xlab = paste0("Individual Factor ", rs1),
           ylab = paste0("Individual Factor ", rs2), main = paste0("Individual Factors Proteomic ", rs1, " vs ", rs2),
           col = ccstat_colors)
      
      # Plot with stand-out subjects as colors
      plot(individual.scores.final.mean[[2]][,rs1], individual.scores.final.mean[[2]][,rs2], xlab = paste0("Individual Factor ", rs1),
           ylab = paste0("Individual Factor ", rs2), main = paste0("Individual Factors Proteomic ", rs1, " vs ", rs2),
           col = cluster_colors, pch = cluster_shapes, lwd = 3)
      par(mfrow = c(1,1))
    }
  } 
}
dev.off()

# -----------------------------------------------------------------------------
# Calculate the posterior mean of the scores, loadings, and betas for the
# aligned and unaligned versions of the models
# * Also, calculate the posterior mean of the joint and individual structures
# -----------------------------------------------------------------------------

# Load in the aligned results --
load("/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Joint_Individual_Structures_Factor_Switching_Ordered_V2_NewPivot_Burnin.rda", verbose = TRUE)

# Load in the original results --
load(paste0(results_wd, "BPMF/Training_Fit/training_data_fit_V2.rda"), verbose = TRUE)

# Load in infinitefactor
library(infinitefactor)

# Calculate the means
joint.scores.final.mean <- lmean(joint.scores.final)
joint.loadings.final.mean <- lmean(joint.loadings.final)
joint.betas.final.mean <- lmean(joint.betas.final)

individual.scores.final.mean <- lapply(1:q, function(s) lmean(individual.scores.final[[s]]))
individual.loadings.final.mean <- lapply(1:q, function(s) lmean(individual.loadings.final[[s]]))
individual.betas.final.mean <- lapply(1:q, function(s) lmean(individual.betas.final[[s]]))

joint.structure.mean <- lapply(1:q, function(s) {
  lmean(lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$J.draw[[iter]][[s,1]]))
})
indiv.structure.mean <- lapply(1:q, function(s) {
  lmean(lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$A.draw[[iter]][[s,1]]))
})

# Save
save(joint.scores.final.mean, joint.loadings.final.mean, joint.betas.final.mean,
     individual.scores.final.mean, individual.loadings.final.mean, individual.betas.final.mean,
     joint.structure.mean, indiv.structure.mean,
     file = "~/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Aligned_Posterior_Mean_NewPivot_Burnin.rda")

# Load in the unaligned results --
load(paste0(results_wd, "BPMF/Training_Fit/training_data_fit_V2.rda"), verbose = TRUE)

# Calculate the means
joint.scores.unaligned.mean <- lmean(lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$V.draw[[iter]][[1,1]]))
joint.loadings.unaligned.mean <- lmean(lapply(iters_burnin, function(iter) do.call(rbind, fev1pp_training_fit_nonsparse_V2$U.draw[[iter]])))
joint.betas.unaligned.mean <- lmean(lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][beta.ind[[1]],,drop=FALSE]))

individual.scores.unaligned.mean <- lapply(1:q, function(s) lmean(lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$Vs.draw[[iter]][[1,s]])))
individual.loadings.unaligned.mean <- lapply(1:q, function(s) lmean(lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$W.draw[[iter]][[s,s]])))
individual.betas.unaligned.mean <- lapply(1:q, function(s) lmean(lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][beta.ind[[s+1]],,drop=FALSE])))

# Save
save(joint.scores.unaligned.mean, joint.loadings.unaligned.mean, joint.betas.unaligned.mean,
     individual.scores.unaligned.mean, individual.loadings.unaligned.mean, individual.betas.unaligned.mean,
     file = "~/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Unaligned_Posterior_Mean.rda")

# -----------------------------------------------------------------------------
# Heatmaps 
# * Initially using JIVE heatmaps
# * Later, switched to sJIVE heatmaps to include response
# -----------------------------------------------------------------------------

# -------------------------------------
# Heatmaps using sup.r.jive
# -------------------------------------

# Load in the plotting library
source("~/BayesianPMF/04DataApplication/BPMF/Figures/plotHeatmap.R")

# Save the posterior mean of the joint scores
S_J <- t(Reduce("+", joint.scores.final)/burnin)

# Posterior mean of individual scores
S_I <- lapply(1:q, function(s) t(Reduce("+", individual.scores.final[[s]])/burnin))

# Posterior mean of joint loadings
U_I_full <- Reduce("+", joint.loadings.final)/ burnin
U_I <- list(U_I_full[1:p.vec[1],,drop=FALSE],
            U_I_full[(p.vec[1]+1):p,,drop=FALSE])

# Posterior mean of the individual loadings
W_I <- lapply(1:q, function(s) Reduce("+", individual.loadings.final[[s]])/burnin)

# Posterior mean of joint betas
theta1 <- t(Reduce("+", joint.betas.final)/burnin)

# Posterior mean of individual betas
theta2 <- lapply(1:q, function(s) t(Reduce("+", individual.betas.final[[s]]))/burnin)

# Combine into list
hiv_copd_data_list_fev1pp <- list(X = list(hiv_copd_data_list[[1]], hiv_copd_data_list[[2]]), 
                                  Y = unlist(fev1pp[[1,1]]))
heatmap_sup.r.jive <- list(data = hiv_copd_data_list_fev1pp,
                           S_J = S_J, S_I = S_I,
                           U_I = U_I, W_I = W_I,
                           theta1 = theta1, theta2 = theta2)

# Ordered by FEV1pp
png(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Heatmap_Ordered_by_FEV1pp_with_Outcome_V2_NewPivot_Burnin.png"), width = 800)
plotHeatmap(heatmap_sup.r.jive, ylab = "FEV1pp", xlab = c("Metabolomics", "Proteomics"))
dev.off()

# Ordered by Joint
png(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Heatmap_Ordered_by_Joint_with_Outcome_V2_NewPivot_Burnin.png"), width = 800)
plotHeatmap(heatmap_sup.r.jive, ylab = "FEV1pp", xlab = c("Metabolomics", "Proteomics"), order_by = 0)
dev.off()

# Ordered by Metabolite
png(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Heatmap_Ordered_by_Metabolite_with_Outcome_V2_NewPivot_Burnin.png"), width = 800)
plotHeatmap(heatmap_sup.r.jive, ylab = "FEV1pp", xlab = c("Metabolomics", "Proteomics"), order_by = 1)
dev.off()

# Ordered by Protein
png(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Heatmap_Ordered_by_Protein_with_Outcome_V2_NewPivot_Burnin.png"), width = 800)
plotHeatmap(heatmap_sup.r.jive, ylab = "FEV1pp", xlab = c("Metabolomics", "Proteomics"), order_by = 2)
dev.off()


# -----------------------------------------------------------------------------
# Credible Intervals for Loadings, Scores, and Betas
# * I create credible interval plots locally, but summarize the results here
# -----------------------------------------------------------------------------

# -------------------------------------
# Creating summaries for ALIGNED factors
# -------------------------------------

# Load in the aligned factors
load("/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Joint_Individual_Structures_Factor_Switching_Ordered_V2_NewPivot_Burnin.rda", verbose = TRUE)

library(dplyr)

# Save the working directory where the summaries will go
exploring_factors_wd <- "~/BayesianPMF/04DataApplication/BPMF/Exploring_Factors/V2/Aligned/NewPivot_Burnin"

# Save the joint and individual ranks
ranks <- fev1pp_training_fit_nonsparse_V2$ranks

# Calculate summaries of the factors 
factor_summaries(joint.loadings.final, joint.scores.final, joint.betas.final,
                 individual.loadings.final, individual.scores.final, individual.betas.final,
                 ranks = ranks, exploring_factors_wd = exploring_factors_wd,
                 nsample_after_burnin = nsample_after_burnin)


# -------------------------------------
# Creating summaries of the UNALIGNED factors
# -------------------------------------

# Adding a burnin
U.draw.burnin <- fev1pp_training_fit_nonsparse_V2$U.draw[(burnin+1):nsample]
W.draw.burnin <- fev1pp_training_fit_nonsparse_V2$W.draw[(burnin+1):nsample]
V.draw.burnin <- fev1pp_training_fit_nonsparse_V2$V.draw[(burnin+1):nsample]
Vs.draw.burnin <- fev1pp_training_fit_nonsparse_V2$Vs.draw[(burnin+1):nsample]

# For each factor,
for (i in 1:sum(ranks)) {
  
  # Save the results for the given factor
  if (i %in% rank.inds[[1]]) {
    rank_index <- i
    factor.final.by.source <- lapply(1:q, function(s) do.call(cbind, lapply(1:burnin, function(iter) U.draw.burnin[[iter]][[s,1]][,rank_index,drop=FALSE])))
    scores.final <- do.call(cbind, lapply(1:burnin, function(iter) V.draw.burnin[[iter]][[1,1]][,rank_index,drop=FALSE]))
    
    # Calculating means and CIs for each loading
    factor.final.summary.by.source <- lapply(1:q, function(s) t(apply(factor.final.by.source[[s]], 1, function(load) {
      c(mean = mean(load), lower.ci = quantile(load, 0.025), upper.ci = quantile(load, 0.975))
    })))
    
    # Name the rows by the respective biomarker
    rownames(factor.final.summary.by.source[[1]]) <- rownames(lavage_processed_no_info_log_scale)
    rownames(factor.final.summary.by.source[[2]]) <- rownames(somascan_normalized_clean_no_info_transpose_scale)
    
    # Separate the loadings by metabolites and proteins
    factor.final.summary <- NULL
    factor.final.summary.metabolites <- factor.final.summary.by.source[[1]]
    factor.final.summary.proteins <- factor.final.summary.by.source[[2]]
    
    # Order the loadings within the metabolites and proteins by the posterior mean
    factor.final.summary.metabolites.order <- 
      factor.final.summary.metabolites[order(factor.final.summary.metabolites[,1], decreasing = FALSE),]
    
    factor.final.summary.proteins.order <- 
      factor.final.summary.proteins[order(factor.final.summary.proteins[,1], decreasing = FALSE),]
    
    # Create the file name
    file_name <- paste0(exploring_factors_wd, "Unaligned/Joint_Factor", rank_index, "_Unaligned_Ordered_Summary.rda")
  }
  
  if (i %in% rank.inds[[2]]) {
    rank_index <- i - ranks[1]
    factor.final <- do.call(cbind, lapply(1:burnin, function(iter) W.draw.burnin[[iter]][[1,1]][,rank_index,drop=FALSE]))
    
    # Calculating means and CIs for each loading
    factor.final.summary <- t(apply(factor.final, 1, function(load) {
      c(mean = mean(load), lower.ci = quantile(load, 0.025), upper.ci = quantile(load, 0.975))
    }))
    
    # Name the rows by the respective biomarker
    rownames(factor.final.summary) <- c(rownames(lavage_processed_no_info_log_scale))
    
    # Separate the loadings by metabolites and proteins (only metabolites here)
    factor.final.summary.metabolites <- factor.final.summary
    factor.final.summary.proteins <- NULL
    
    # Order the loadings within the metabolites and proteins by the posterior mean
    factor.final.summary.metabolites.order <- 
      factor.final.summary.metabolites[order(factor.final.summary.metabolites[,1], decreasing = FALSE),]
    
    # Create a file name
    file_name <- paste0(exploring_factors_wd, "Unaligned/Metabolite_Indiv_Factor", rank_index, "_Unaligned_Ordered_Summary.rda")
    
  }
  
  if (i %in% rank.inds[[3]]) {
    rank_index <- i - cumsum(ranks[1:2])[2]
    factor.final <- do.call(cbind, lapply(1:burnin, function(iter) W.draw.burnin[[iter]][[2,2]][,rank_index,drop=FALSE]))
    
    # Calculating means and CIs for each loading
    factor.final.summary <- t(apply(factor.final, 1, function(load) {
      c(mean = mean(load), lower.ci = quantile(load, 0.025), upper.ci = quantile(load, 0.975))
    }))
    
    # Name the rows by the respective biomarker
    rownames(factor.final.summary) <- rownames(somascan_normalized_clean_no_info_transpose_scale)
    
    # Separate the loadings by metabolites and proteins (only proteins here)
    factor.final.summary.metabolites.order <- NULL
    factor.final.summary.proteins <- factor.final.summary
    
    # Order the loadings within the metabolites and proteins by the posterior mean
    factor.final.summary.proteins.order <- 
      factor.final.summary.proteins[order(factor.final.summary.proteins[,1], decreasing = FALSE),]
    
    # Create the file name
    file_name <- paste0(exploring_factors_wd, "Unaligned/Protein_Indiv_Factor", rank_index, "_Unaligned_Ordered_Summary.rda")
    
  }
  
  # Save the results
  save(factor.final.summary, factor.final.summary.metabolites.order, factor.final.summary.proteins.order,
       file = file_name)
}

# -----------------------------------------------------------------------------
# Contributions of each structure (joint, individual) to explaining variation in 
# the omics sources and in Y
# -----------------------------------------------------------------------------

# Calculating a summary of the variance explained
var_exp_summary <- var_explained(BPMF.fit = fev1pp_training_fit_nonsparse_V2,
                                 hiv_copd_data = hiv_copd_data, iters_burnin = iters_burnin)

# Joint structure 
var_exp_summary$Joint

# Individual structure
var_exp_summary$Individual

# -----------------------------------------------------------------------------
# Contributions of each factor (joint, individual) to predicting Y
# -----------------------------------------------------------------------------

# Calculating the contribution of the joint factors to predicting FEV1pp
joint_contribution_to_fev1pp <- sapply(1:nsample_after_burnin, function(iter) {
  var(joint.scores.final[[iter]] %*% joint.betas.final[[iter]])/var(fev1pp[[1,1]])
})

# Calculating the contribution of the individual factors to predicting FEV1pp
individual_contribution_to_fev1pp <- lapply(1:q, function(s) {
  sapply(1:nsample_after_burnin, function(iter) {
    var(individual.scores.final[[s]][[iter]] %*% individual.betas.final[[s]][[iter]])/var(fev1pp[[1,1]])
  })
})

# Overall of Joint 
mean(joint_contribution_to_fev1pp)
c(quantile(joint_contribution_to_fev1pp, 0.025), quantile(joint_contribution_to_fev1pp, 0.975))

# Overall contribution of Indiviudal
sapply(individual_contribution_to_fev1pp, mean)
sapply(individual_contribution_to_fev1pp, function(indiv) {
  c(quantile(indiv, 0.025), quantile(indiv, 0.975))
})

# Calculate the contribution of the joint structure to predicting FEV1pp for each sample
joint_contribution_to_fev1pp_by_sample <- sapply(1:nsample_after_burnin, function(iter) {
  joint.scores.final[[iter]] %*% joint.betas.final[[iter]]
})

# Calculate the contribution of the individual structures to predicting FEV1pp for each sample
individual_contribution_to_fev1pp_by_sample <- lapply(1:q, function(s) {
  sapply(1:nsample_after_burnin, function(iter) {
    individual.scores.final[[s]][[iter]] %*% individual.betas.final[[s]][[iter]]
  })
})

# Calculate the overall contribution of the factors in explaining FEV1pp
overall_contribution_to_fev1pp <- sapply(1:nsample_after_burnin, function(iter) {
    fev1pp_training_fit_nonsparse_V2$beta.draw[iters_burnin][[iter]][[1,1]][1,] +  # Intercept
    joint.scores.final[[iter]] %*% joint.betas.final[[iter]] +                     # Joint
    individual.scores.final[[1]][[iter]] %*% individual.betas.final[[1]][[iter]] + # Metabolome
    individual.scores.final[[2]][[iter]] %*% individual.betas.final[[2]][[iter]]   # Proteome
})

# Save the results
save(joint_contribution_to_fev1pp_by_sample, individual_contribution_to_fev1pp_by_sample, overall_contribution_to_fev1pp,
     file = "~/BayesianPMF/04DataApplication/BPMF/Exploring_Factors/V2/Aligned/NewPivot_Burnin/Structure_Contribution_to_FEV1pp_NewPivot_Burnin.rda")


# -----------------------------------------------------------------------------
# K-Means on the Posterior Sampling Iterations
# -----------------------------------------------------------------------------

# Loading in the joint structure with burnin
joint.structure.burnin <- lapply(iters_burnin, function(iter) {
  
  # Add in E(Y) for the joint structure
  J.draw.EY.combined <- matrix(list(), nrow = q+1, ncol = 1)
  J.draw.EY.combined[[1,1]] <- fev1pp_training_fit_nonsparse_V2$J.draw[[iter]][[1,1]]
  J.draw.EY.combined[[2,1]] <- fev1pp_training_fit_nonsparse_V2$J.draw[[iter]][[2,1]]
  J.draw.EY.combined[[3,1]] <- t(fev1pp_training_fit_nonsparse_V2$V.draw[[iter]][[1,1]] %*% (fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][-1,,drop=FALSE][rank.inds[[1]],,drop=FALSE]))
  
  data.rearrange(J.draw.EY.combined)$out
})

# For each posterior sampling iteration, apply the K-means algorithm
joint.structure.kmeans <- sapply(1:length(iters_burnin), function(iter) {
  # For reproducibility
  set.seed(iter)
  
  # Run K-means
  res <- kmeans(x = t(joint.structure.burnin[[iter]]),
                centers = 2,
                iter.max = 50,
                nstart = 10)
  
  # Return cluster
  res$cluster
})

# Transpose
joint.structure.kmeans <- t(joint.structure.kmeans)

# Add subject IDs
colnames(joint.structure.kmeans) <- colnames(somascan_normalized_clean_no_info_transpose_scale)

# Check cluster sizes
cluster.sizes <- t(apply(joint.structure.kmeans, 1, table))

# Rename the smaller cluster to 1 and larger to 2
joint.structure.kmeans.rename <- joint.structure.kmeans

for (iter in 1:length(iters_burnin)) {
  # Save the current clustering scheme
  current_cluster <- joint.structure.kmeans[iter,]
  
  # If the larger cluster is a 1, rename to 2
  if (which.min(table(current_cluster)) == 2) {
    joint.structure.kmeans.rename[iter,][current_cluster == 2] <- 1
    joint.structure.kmeans.rename[iter,][current_cluster == 1] <- 2
  }
}

# Save the cluster sizes
cluster.sizes.rename <- t(apply(joint.structure.kmeans.rename, 1, table))

# Plot the sizes the smaller and larger clusters over time
plot(cluster.sizes.rename[,2], col = 1, pch = 16, ylim = c(0, 52), ylab = "Cluster Size", main = "Cluster Sizes Across Sampling Iterations")
abline(a = mean(cluster.sizes.rename[,2]), b = 0, col = 1, lty = 2, lwd = 2)
points(cluster.sizes.rename[,1], col = 5, pch = 16)
abline(a = mean(cluster.sizes.rename[,1]), b = 0, col = 5, lty = 2, lwd = 2)

# Calculating the proportion of time each individual spent in Cluster 1
prop.cluster1 <- apply(joint.structure.kmeans.rename, 2, function(sample) {
  sum(sample == 1)/length(iters_burnin)
})
sort(prop.cluster1, decreasing = TRUE)

# -----------------------------------------------------------------------------
# Cross-Validated Model Fit
# -----------------------------------------------------------------------------

# Saving the index for each pair 
ind_of_pairs <- seq(1, n, by = 2)
n_pair <- length(ind_of_pairs)

# For running in parallel
funcs <- c("bpmf_data", "center_data", "bpmf_data_mode", "get_results", "BIDIFAC",
           "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
           "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum")
packs <- c("MASS", "truncnorm", "EnvStats", "svMisc", "Matrix")

# -------------------------------------
# Running BPMF with cross validation
# -------------------------------------

cl <- makeCluster(2)
registerDoParallel(cl)
fev1pp_cv <- foreach(pair = ind_of_pairs, .packages = packs, .export = funcs, .verbose = TRUE) %dopar% {
  # Create a new vector of the outcome with the current pair set to NA
  fev1pp_cv <- fev1pp
  fev1pp_cv[[1,1]][pair:(pair+1),] <- NA
  
  # Run the model with the above pair's continuous outcome missing
  fev1pp_cv_fit_nonsparse <- bpmf_data_mode(
    data = hiv_copd_data,
    Y = fev1pp_cv,
    nninit = TRUE,
    model_params = model_params,
    sparsity = FALSE,
    nsample = nsample,
    progress = TRUE
  )
  
  # -------------------------------------
  # Calculate the log joint density for held-out pair
  # -------------------------------------
  
  # Subset the data to just this pair
  hiv_copd_data_pair <- matrix(list(), nrow = q, ncol = 1)
  for (s in 1:q) {
    hiv_copd_data_pair[[s,1]] <- hiv_copd_data[[s,1]][, pair:(pair+1)]
  }
  
  # Subset the posterior samples for the scores for just these subjects
  V.draw_pair <- lapply(fev1pp_cv_fit_nonsparse$V.draw, function(iter) {
    # Init the subsetted matrix
    V <- matrix(list(), nrow = 1, ncol = 1)
    
    # Fill in 
    V[[1,1]] <- iter[[1,1]][pair:(pair+1),]
    
    # Return
    V
  })
  
  Vs.draw_pair <- lapply(fev1pp_cv_fit_nonsparse$Vs.draw, function(iter) {
    # Init the subsetted matrix
    Vs <- matrix(list(), nrow = 1, ncol = q)
    
    # Fill in 
    for (s in 1:q) {
      Vs[[1,s]] <- iter[[1,s]][pair:(pair+1),]
    }
    
    # Return
    Vs
  })
  
  # Save the ranks
  ranks <- fev1pp_cv_fit_nonsparse$ranks
  
  # Saving the true outcomes for the missing subjects
  Y_pair <- matrix(list(), nrow = 1, ncol = 1)
  Y_pair[[1,1]] <- fev1pp_cv[[1,1]][pair:(pair+1),,drop = FALSE]
  
  # Save the imputed outcomes
  Ym.draw_pair <- fev1pp_cv_fit_nonsparse$Ym.draw
  
  # Calculating the log-joint density after burn-in
  convergence <- sapply(thinned_iters_burnin, function(sim_iter) {
    # Calculate the log-joint density at each thinned iterations
    log_joint_density(data = hiv_copd_data_pair, 
                      U.iter = fev1pp_cv_fit_nonsparse$U.draw[[sim_iter]], 
                      V.iter = V.draw_pair[[sim_iter]], 
                      W.iter = fev1pp_cv_fit_nonsparse$W.draw[[sim_iter]], 
                      Vs.iter = Vs.draw_pair[[sim_iter]],
                      model_params = model_params,
                      ranks = fev1pp_cv_fit_nonsparse$ranks,
                      Y = Y_pair,
                      Ym.iter = Ym.draw_pair[[sim_iter]],
                      beta.iter = fev1pp_cv_fit_nonsparse$beta.draw[[sim_iter]],
                      tau2.iter = fev1pp_cv_fit_nonsparse$tau2.draw[[sim_iter]])
  })
  
  # Save just the relevant output
  save(Ym.draw_pair, ranks, convergence, file = paste0(results_wd, "/BPMF/Cross_Validation/FEV1pp_CV_NonSparse_Pair", pair, ".rda"))
  
  # Remove large objects
  rm(hiv_copd_data_pair, fev1pp_cv_fit_nonsparse, V.draw_pair, Vs.draw_pair, Ym.draw_pair)
  
  # Garbage collection
  gc()
}
stopCluster(cl)

# -------------------------------------
# Running BIDIFAC with cross validation 
# -------------------------------------

run_model_with_cv(mod = "BIDIFAC", hiv_copd_data = hiv_copd_data, outcome = fev1pp,
                  outcome_name = "FEV1pp", ind_of_pairs = ind_of_pairs, 
                  model_params = model_params, nsample = nsample)


# -------------------------------------
# Running JIVE with cross validation
# -------------------------------------

run_model_with_cv(mod = "JIVE", hiv_copd_data = hiv_copd_data, outcome = fev1pp,
                  outcome_name = "FEV1pp", ind_of_pairs = ind_of_pairs, 
                  model_params = model_params, nsample = nsample)

# -------------------------------------
# Running MOFA with cross validation
# -------------------------------------

run_model_with_cv(mod = "MOFA", hiv_copd_data = hiv_copd_data, outcome = fev1pp,
                  outcome_name = "FEV1pp", ind_of_pairs = ind_of_pairs,
                  model_params = model_params, nsample = nsample)

# -------------------------------------
# Running sJIVE with cross validation
# -------------------------------------

run_model_with_cv(mod = "sJIVE", hiv_copd_data = hiv_copd_data, outcome = fev1pp,
                  outcome_name = "FEV1pp", ind_of_pairs = ind_of_pairs, 
                  model_params = model_params, nsample = nsample)

# -------------------------------------
# Running LASSO with cross validation on combined data
# -------------------------------------

run_model_with_cv(mod = "LASSO_Combined_Sources", hiv_copd_data = hiv_copd_data, outcome = fev1pp,
                  outcome_name = "FEV1pp", ind_of_pairs = ind_of_pairs, 
                  model_params = model_params, nsample = nsample)

# -------------------------------------
# Running LASSO with cross validation on metabolites only
# -------------------------------------

run_model_with_cv(mod = "LASSO_Metabolite_Only", hiv_copd_data = hiv_copd_data, outcome = fev1pp,
                  outcome_name = "FEV1pp", ind_of_pairs = ind_of_pairs, 
                  model_params = model_params, nsample = nsample)

# -------------------------------------
# Running LASSO with cross validation on proteins only
# -------------------------------------

run_model_with_cv(mod = "LASSO_Protein_Only", hiv_copd_data = hiv_copd_data, outcome = fev1pp,
                  outcome_name = "FEV1pp", ind_of_pairs = ind_of_pairs, 
                  model_params = model_params, nsample = nsample)

# -------------------------------------
# Cross-Validated Model Fit Results
# -------------------------------------

# Create a vector with model names
models <- c("BPMF", "BIDIFAC", "JIVE", "MOFA", "sJIVE", "LASSO_Combined_Sources", "LASSO_Metabolite_Only", "LASSO_Protein_Only")

# Create a vector of cross-validated FEV1pp results for each model
fev1pp_cv <- lapply(models, function(mod) c())
fev1pp_cv_ci <- lapply(models, function(mod) matrix(nrow = n, ncol = 2))
names(fev1pp_cv) <- names(fev1pp_cv_ci) <- models

# Iterate through the case-control pairs
for (mod in models) {
  for (pair in ind_of_pairs) {
    
    # Load in the results
    if (mod == "BPMF") {
      load(paste0(results_wd, "BPMF/Cross_Validation/FEV1pp_CV_NonSparse_Pair", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "BIDIFAC") {
      load(paste0(results_wd, "BIDIFAC/Cross_Validation/FEV1pp_CV_BIDIFAC_Pair_", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "JIVE") {
      load(paste0(results_wd, "JIVE/Cross_Validation/FEV1pp_CV_JIVE_Pair_", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "MOFA") {
      load(paste0(results_wd, "MOFA/Cross_Validation/FEV1pp_CV_MOFA_Pair_", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "sJIVE") {
      load(paste0(results_wd, "sJIVE/Cross_Validation/FEV1pp_CV_sJIVE_Pair_", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "LASSO_Combined_Sources") {
      load(paste0(results_wd, "LASSO_Combined_Sources/Cross_Validation/FEV1pp_CV_LASSO_Combined_Sources_Pair_", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "LASSO_Metabolite_Only") {
      load(paste0(results_wd, "LASSO_Metabolite_Only/Cross_Validation/FEV1pp_CV_LASSO_Metabolite_Only_Pair_", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "LASSO_Protein_Only") {
      load(paste0(results_wd, "LASSO_Protein_Only/Cross_Validation/FEV1pp_CV_LASSO_Protein_Only_Pair_", pair, ".rda"), verbose = TRUE)
    }
    
    # Combine the samples
    if (mod != "sJIVE" & !grepl("LASSO", mod)) {
      samps <- do.call(cbind, do.call(cbind, Ym.draw_pair))
      
      # Take a burn-in
      samps_burnin <- samps[,thinned_iters_burnin]
      
      # Save in the vector
      fev1pp_cv[[mod]][pair:(pair+1)] <- rowMeans(samps_burnin)
      
      # Calculate the 95% credible interval for each held-out individual
      fev1pp_cv_ci[[mod]][pair,] <- c(quantile(samps_burnin[1,], 0.025), quantile(samps_burnin[1,], 0.975))
      fev1pp_cv_ci[[mod]][pair+1,] <- c(quantile(samps_burnin[2,], 0.025), quantile(samps_burnin[2,], 0.975))
    }
    
    if (mod == "sJIVE" | grepl("LASSO", mod)) {
      # Save in the vector
      fev1pp_cv[[mod]][pair:(pair+1)] <- Ym.draw_pair
    }
  }
}

# Plotting the results from each model against the truth
plot(fev1pp_cv$BPMF, c(fev1pp[[1,1]]), xlab = "Predicted FEV1pp", ylab = "Observed FEV1pp", main = "Cross-Validated FEV1pp vs. True FEV1pp", pch = 16)
points(fev1pp_cv$BIDIFAC, c(fev1pp[[1,1]]), col = 2, pch = 16) # BIDIFAC
points(fev1pp_cv$JIVE, c(fev1pp[[1,1]]), col = 3, pch = 16) # JIVE
points(fev1pp_cv$MOFA, c(fev1pp[[1,1]]), col = 4, pch = 16) # MOFA
points(fev1pp_cv$sJIVE, c(fev1pp[[1,1]]), col = 4, pch = 16) # sJIVE
abline(a=0, b=1, lwd = 2)
legend("bottomright", legend = c("BPMF", "BIDIFAC", "JIVE", "MOFA"), col = c(1, 2, 3, 4), pch = rep(16, 4), cex = 0.5)

# Calculating a correlation test for each cross validated outcome vs. true FEV1pp
cor.test(fev1pp_cv$BPMF, c(fev1pp[[1,1]])) # BPMF
cor.test(fev1pp_cv$BIDIFAC, c(fev1pp[[1,1]])) # BIDIFAC
cor.test(fev1pp_cv$sJIVE, c(fev1pp[[1,1]])) # sJIVE
cor.test(fev1pp_cv$JIVE, c(fev1pp[[1,1]])) # JIVE
cor.test(fev1pp_cv$MOFA, c(fev1pp[[1,1]])) # MOFA
cor.test(fev1pp_cv$LASSO_Combined_Sources, c(fev1pp[[1,1]])) # LASSO Combined Sources
cor.test(fev1pp_cv$LASSO_Metabolite_Only, c(fev1pp[[1,1]])) # LASSO Combined Sources
cor.test(fev1pp_cv$LASSO_Protein_Only, c(fev1pp[[1,1]])) # LASSO Combined Sources

# Compare the coverage rates for FEV1pp
mean(sapply(1:n, function(i) fev1pp_cv_ci$BPMF[i,1] <= fev1pp[[1,1]][i] & fev1pp[[1,1]][i] <= fev1pp_cv_ci$BPMF[i,2])) # BPMF
mean(sapply(1:n, function(i) fev1pp_cv_ci$BIDIFAC[i,1] <= fev1pp[[1,1]][i] & fev1pp[[1,1]][i] <= fev1pp_cv_ci$BIDIFAC[i,2])) # BIDIFAC
mean(sapply(1:n, function(i) fev1pp_cv_ci$JIVE[i,1] <= fev1pp[[1,1]][i] & fev1pp[[1,1]][i] <= fev1pp_cv_ci$JIVE[i,2])) # JIVE
mean(sapply(1:n, function(i) fev1pp_cv_ci$MOFA[i,1] <= fev1pp[[1,1]][i] & fev1pp[[1,1]][i] <= fev1pp_cv_ci$MOFA[i,2])) # MOFA

# Save the imputed values and credible intervals to plot locally
save(fev1pp_cv, fev1pp_cv_ci, file = "~/BayesianPMF/04DataApplication/FEV1pp_EY_CV.rda")

# Create table
library(xtable)

cv_results <- data.frame(Model = character(), Correlation = numeric(), PValue = numeric())

for (mod in models) {
  res <- cor.test(fev1pp_cv[[mod]], c(fev1pp[[1,1]]))
  cv_results[nrow(cv_results) + 1,] <- data.frame(mod, res$estimate, res$p.value)
}

cv_results[order(cv_results$Correlation, decreasing = TRUE),] %>%
  xtable(digits = 4) %>%
    print.xtable(include.rownames=FALSE)

# -----------------------------------------------------------------------------
# Missing data imputation comparing several methods
# -----------------------------------------------------------------------------

# -------------------------------------
# Missing data imputation using BPMF
# -------------------------------------

# Setting the proportion of missing values (columns) in each dataset
prop_missing <- 0.1

# Save the number of replications
nsim <- 20

# Running BPMF with entrywise missing data
hiv_copd_bpmf_imputation <- model_imputation(mod = "BPMF", hiv_copd_data = hiv_copd_data,
                                             outcome = NULL, outcome_name = "fev1pp", nsample = 5000,
                                             nsim = nsim, prop_missing = prop_missing, entrywise = TRUE,
                                             nclust = 5)

# Running BPMF with columnwise missing data
hiv_copd_bpmf_imputation <- model_imputation(mod = "BPMF", hiv_copd_data = hiv_copd_data,
                                             outcome = NULL, outcome_name = "fev1pp", nsample = 5000,
                                             nsim = nsim, prop_missing = prop_missing, entrywise = FALSE,
                                             nclust = 5)

# -------------------------------------
# Missing data imputation using BIDIFAC
# -------------------------------------

# Running BIDIFAC with entrywise missing data
hiv_copd_bidifac_imputation <- model_imputation(mod = "BIDIFAC", hiv_copd_data = hiv_copd_data,
                                                 outcome = fev1pp, outcome_name = "fev1pp", p.vec = p.vec,
                                                 nsim = nsim, prop_missing = prop_missing, entrywise = TRUE)

# Running BIDIFAC with columnwise missing data
hiv_copd_bidifac_imputation <- model_imputation(mod = "BIDIFAC", hiv_copd_data = hiv_copd_data,
                                                 outcome = fev1pp, outcome_name = "fev1pp", p.vec = p.vec,
                                                 nsim = nsim, prop_missing = prop_missing, entrywise = FALSE)

# -------------------------------------
# Missing data imputation using SVD
# -------------------------------------

# Running SVDmiss with entrywise missing data, combined sources
hiv_copd_svd_imputation <- model_imputation(mod = "SVD_Combined_Sources", hiv_copd_data = hiv_copd_data,
                                              outcome = fev1pp, outcome_name = "fev1pp", p.vec = p.vec,
                                              nsim = nsim, prop_missing = prop_missing, entrywise = TRUE)

# Running SVDmiss with entrywise missing data, separate sources
hiv_copd_svdmiss_imputation <- model_imputation(mod = "SVD_Separate_Sources", hiv_copd_data = hiv_copd_data,
                                              outcome = fev1pp, outcome_name = "fev1pp", p.vec = p.vec,
                                              nsim = nsim, prop_missing = prop_missing, entrywise = TRUE)

# Running SVDmiss with columnwise missing data, combined sources
hiv_copd_svd_imputation <- model_imputation(mod = "SVD_Combined_Sources", hiv_copd_data = hiv_copd_data,
                                             outcome = fev1pp, outcome_name = "fev1pp", p.vec = p.vec,
                                             nsim = nsim, prop_missing = prop_missing, entrywise = FALSE)

# Running SVDmiss with columnwise missing data, separate sources - Unable to complete matrix, too much missing data
hiv_copd_svdmiss_imputation <- model_imputation(mod = "SVD_Separate_Sources", hiv_copd_data = hiv_copd_data,
                                                 outcome = fev1pp, outcome_name = "fev1pp", p.vec = p.vec,
                                                 nsim = nsim, prop_missing = prop_missing, entrywise = FALSE)

# -------------------------------------
# Missing data imputation using mean imputation
# -------------------------------------

# Entrywise
hiv_copd_mean_imputation <- model_imputation(mod = "Mean_Imputation", hiv_copd_data = hiv_copd_data,
                                             outcome = fev1pp, outcome_name = "fev1pp", p.vec = p.vec,
                                             nsim = nsim, prop_missing = prop_missing, entrywise = TRUE)

# Columnwise
hiv_copd_mean_imputation <- model_imputation(mod = "Mean_Imputation", hiv_copd_data = hiv_copd_data,
                                             outcome = fev1pp, outcome_name = "fev1pp", p.vec = p.vec,
                                             nsim = nsim, prop_missing = prop_missing, entrywise = FALSE)

# -------------------------------------
# Missing data imputation using k-Nearest Neighbors
# -------------------------------------

# Entrywise
hiv_copd_knn_combined <- model_imputation(mod = "KNN_Combined_Sources", hiv_copd_data = hiv_copd_data,
                                          outcome = fev1pp, outcome_name = "fev1pp", 
                                          nsim = nsim, prop_missing = prop_missing, entrywise = TRUE,
                                          nclust = 10)

# Columnwise
hiv_copd_knn_combined <- model_imputation(mod = "KNN_Combined_Sources", hiv_copd_data = hiv_copd_data,
                                          outcome = fev1pp, outcome_name = "fev1pp",
                                          nsim = nsim, prop_missing = prop_missing, entrywise = FALSE,
                                          nclust = 10)

# Entrywise
hiv_copd_knn_combined <- model_imputation(mod = "KNN_Separate_Sources", hiv_copd_data = hiv_copd_data,
                                          outcome = fev1pp, outcome_name = "fev1pp",
                                          nsim = nsim, prop_missing = prop_missing, entrywise = TRUE,
                                          nclust = 10)

# Columnwise
hiv_copd_knn_combined <- model_imputation(mod = "KNN_Separate_Sources", hiv_copd_data = hiv_copd_data,
                                          outcome = fev1pp, outcome_name = "fev1pp", 
                                          nsim = nsim, prop_missing = prop_missing, entrywise = FALSE,
                                          nclust = 10)

# -------------------------------------
# Missing data imputation using random forest
# -------------------------------------

# Entrywise
hiv_copd_rf_combined <- model_imputation(mod = "RF_Combined_Sources", hiv_copd_data = hiv_copd_data,
                                          outcome = fev1pp, outcome_name = "fev1pp", 
                                          nsim = nsim, prop_missing = prop_missing, entrywise = TRUE,
                                          nclust = 10)

# Columnwise
hiv_copd_rf_combined <- model_imputation(mod = "RF_Combined_Sources", hiv_copd_data = hiv_copd_data,
                                          outcome = fev1pp, outcome_name = "fev1pp",
                                          nsim = nsim, prop_missing = prop_missing, entrywise = FALSE,
                                          nclust = 10)

# Entrywise
hiv_copd_rf_separate <- model_imputation(mod = "RF_Separate_Sources", hiv_copd_data = hiv_copd_data,
                                          outcome = fev1pp, outcome_name = "fev1pp",
                                          nsim = nsim, prop_missing = prop_missing, entrywise = TRUE,
                                          nclust = 10)

# Columnwise
hiv_copd_rf_separate <- model_imputation(mod = "RF_Separate_Sources", hiv_copd_data = hiv_copd_data,
                                          outcome = fev1pp, outcome_name = "fev1pp", 
                                          nsim = nsim, prop_missing = prop_missing, entrywise = FALSE,
                                          nclust = 10)


# -------------------------------------
# Tabulating the results
# -------------------------------------

# Initialize table with results
hiv_copd_imputation_results <- data.frame(Model = character(), 
                                          Missingness = character(),
                                          Metabolome_MSE = character(),
                                          Proteome_MSE = character(),
                                          Metabolome_Coverage = character(),
                                          Proteome_Coverage = character(),
                                          Metabolome_CI_Width = character(),
                                          Proteome_CI_Width = character())

# For each model
mods <- c("BPMF", "BIDIFAC", "Mean_Imputation", "SVD_Combined_Sources", "SVD_Separate_Sources",
          "KNN_Combined_Sources", "KNN_Separate_Sources", "RF_Combined_Sources", "RF_Separate_Sources")

for (mod in mods) {
  
  # List the results for the ENTRYWISE imputations
  entrywise_files <- list.files(paste0("~/BayesianPMF/04DataApplication/", mod,"/Imputation/Entrywise/"))
  columnwise_files <- list.files(paste0("~/BayesianPMF/04DataApplication/", mod,"/Imputation/Columnwise/"))
  nsim <- length(entrywise_files)
  
  # Initialize a data.frame for the interim results
  current_mod_entrywise_results <- current_mod_columnwise_results <-
    data.frame(Metabolome_MSE = numeric(),
               Proteome_MSE = numeric(),
               Metabolome_Coverage = numeric(),
               Proteome_Coverage = numeric(),
               Metabolome_CI_Width = numeric(),
               Proteome_CI_Width = numeric())
  
  # Iterate through the results and load in
  for (i in 1:nsim) {
    
    # Load in ith entrywise result
    entrywise_out <- load(paste0("~/BayesianPMF/04DataApplication/", mod,"/Imputation/Entrywise/", entrywise_files[i]))
    
    # Save entrywise result
    current_mod_entrywise_results[i,] <- c(metabolome_mse, proteome_mse, metabolome_coverage,
                                           proteome_coverage, metabolome_ci_width, proteome_ci_width)
    
    # Avoid confusion
    rm(entrywise_out)
    
    # If there are columnwise results
    if (length(columnwise_files) > 0) {
      # Load in ith columnwise result
      columnwise_out <- load(paste0("~/BayesianPMF/04DataApplication/", mod,"/Imputation/Columnwise/", columnwise_files[i]))
      
      # Save columnwise result
      current_mod_columnwise_results[i,] <- c(metabolome_mse, proteome_mse, metabolome_coverage,
                                              proteome_coverage, metabolome_ci_width, proteome_ci_width)
      
      # Avoid confusion
      rm(columnwise_out)
    } 
  }
  
  # Calculate the mean for entrywise and columnwise results
  current_mod_entrywise_mean <- colMeans(current_mod_entrywise_results)
  current_mod_columnwise_mean <- colMeans(current_mod_columnwise_results)
  
  current_mod_entrywise_sd <- apply(current_mod_entrywise_results, 2, sd)
  current_mod_columnwise_sd <- apply(current_mod_columnwise_results, 2, sd)
  
  # Construct vector to save results
  entrywise_vec <- c(mod, "Entrywise", paste0(round(current_mod_entrywise_mean, 3), " (", round(current_mod_entrywise_sd, 3), ")"))
  columnwise_vec <- c(mod, "Columnwise", paste0(round(current_mod_columnwise_mean, 3), " (", round(current_mod_columnwise_sd, 3), ")"))
  
  # Save
  hiv_copd_imputation_results <- rbind.data.frame(hiv_copd_imputation_results, entrywise_vec)
  hiv_copd_imputation_results <- rbind.data.frame(hiv_copd_imputation_results, columnwise_vec)
  
}

# Adjust the column names
colnames(hiv_copd_imputation_results) <- c("Model", "Missingness", "Metabolome_MSE",
                                           "Proteome_MSE", "Metabolome_Coverage", "Proteome_Coverage",
                                           "Metabolome_CI_Width", "Proteome_CI_Width")

# Select just the first four columns
library(xtable)

hiv_copd_imputation_results %>%
  select(Model, Missingness, Metabolome_MSE, Proteome_MSE) %>%
  arrange(Missingness) %>%
  xtable() %>%
  print(include.rownames = FALSE)
