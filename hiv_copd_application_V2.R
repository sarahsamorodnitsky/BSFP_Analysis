# -----------------------------------------------------------------------------
# Applying the BPMF model to the HIV-COPD BALF Biocrates and Somascan data. 
# Comparing the model application when there is no sparsity vs. when
# there is sparsity. 
# * Note on V2: This is called V2 because I wanted to start a clean version of 
#   HIV-COPD application script. 
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set-up
# -----------------------------------------------------------------------------

# Packages
library(doParallel)
library(foreach)

# Working directories
data_wd <- "~/BayesianPMFWithGit/data/"
model_wd <- "~/BayesianPMFWithGit/"
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
burnin <- nsample/2
iters_burnin <- seq(burnin+1, nsample)
thinned_iters <- seq(1, nsample, by = 10)
thinned_iters_burnin <- seq(burnin, nsample, by = 10)

# -----------------------------------------------------------------------------
# Training Data Model Fit (BPMF, BIDIFAC, sJIVE, MOFA, BIP)
# -----------------------------------------------------------------------------

# Fitting BPMF
fev1pp_training_fit_nonsparse <- bpmf(
  data = hiv_copd_data,
  Y = fev1pp,
  nninit = TRUE,
  model_params = model_params,
  sparsity = FALSE,
  nsample = nsample,
  progress = TRUE
)

# Save the results
# save(fev1pp_training_fit_nonsparse, fev1pp_training_fit_sparse, file = paste0(results_wd, "training_data_fit.rda"))

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
BIP_training_fit <- BIP(dataList = hiv_copd_data_list_bip, IndicVar = c(0,0,1), Method = "BIP", sample = burnin, burnin = burnin-1, nbrcomp = 50)

# -----------------------------------------------------------------------------
# Assessing convergence of BPMF
# -----------------------------------------------------------------------------

# Assessing convergence for both model fits
load(paste0(results_wd, "BPMF/Training_Fit/training_data_fit.rda"), verbose = TRUE)

fev1pp_training_nonsparse_conv <- sapply(1:nsample, function(sim_iter) {
  # Calculate the log-joint density at each thinned iterations
  log_joint_density(data = hiv_copd_data, 
                    U.iter = fev1pp_training_fit_nonsparse$U.draw[[sim_iter]], 
                    V.iter = fev1pp_training_fit_nonsparse$V.draw[[sim_iter]], 
                    W.iter = fev1pp_training_fit_nonsparse$W.draw[[sim_iter]], 
                    Vs.iter = fev1pp_training_fit_nonsparse$Vs.draw[[sim_iter]],
                    model_params = model_params,
                    ranks = fev1pp_training_fit_nonsparse$ranks,
                    Y = fev1pp,
                    beta.iter = fev1pp_training_fit_nonsparse$beta.draw[[sim_iter]],
                    tau2.iter = fev1pp_training_fit_nonsparse$tau2.draw[[sim_iter]])
})

# Save the log-joint density 
# save(fev1pp_training_nonsparse_conv, file = "~/BayesianPMF/04DataApplication/BPMF/Training_Fit/Log_Joint_Density_NonSparse.rda")

# Plotting the log-joint densities
load("~/BayesianPMF/04DataApplication/BPMF/Training_Fit/Log_Joint_Density_NonSparse.rda", verbose = TRUE)

png("~/BayesianPMF/04DataApplication/BPMF/Training_Fit/log_joint_density_plot_all_iters.png")
plot(fev1pp_training_nonsparse_conv, ylab = "Log Joint Density", main = "Log Joint Density Across Gibbs Sampling Iterations")
dev.off()

png("~/BayesianPMF/04DataApplication/BPMF/Training_Fit/log_joint_density_plot_burnin.png")
plot(fev1pp_training_nonsparse_conv[burnin:nsample], ylab = "Log Joint Density", main = "Log Joint Density After Burn-in")
dev.off()

# -----------------------------------------------------------------------------
# Adjusting for factor switching
# * Using the Poworoznek et al. (2021) solution 
# -----------------------------------------------------------------------------

# Loading in the full training data fit
load(paste0(results_wd, "BPMF/Training_Fit/training_data_fit.rda"), verbose = TRUE)

# Save the ranks
ranks <- fev1pp_training_fit_nonsparse$ranks
joint.rank <- fev1pp_training_fit_nonsparse$ranks[1]
indiv.ranks <- fev1pp_training_fit_nonsparse$ranks[-1]

# Applying the Match Align algorithm to undo rotational invariance in the results
fev1pp_training_fit_nonsparse_aligned <- match_align_bpmf(fev1pp_training_fit_nonsparse, y = fev1pp,
                                                          model_params = model_params, p.vec = p.vec)

# Save the results from the label switching algorithm
joint.scores.aligned <- fev1pp_training_fit_nonsparse_aligned$joint.scores.final
joint.loadings.aligned <- fev1pp_training_fit_nonsparse_aligned$joint.loadings.final
joint.betas.aligned <- fev1pp_training_fit_nonsparse_aligned$joint.betas.final

individual.scores.aligned <- fev1pp_training_fit_nonsparse_aligned$individual.scores.final
individual.loadings.aligned <- fev1pp_training_fit_nonsparse_aligned$individual.loadings.final
individual.betas.aligned <- fev1pp_training_fit_nonsparse_aligned$individual.betas.final

# Order the components by squared Frobenius-norm of the corresponding rank-1 structure after burn-in --

# Save the aligned components after burn-in
joint.scores.aligned.burnin <- joint.scores.aligned[iters_burnin]
joint.loadings.aligned.burnin <- joint.loadings.aligned[iters_burnin]
joint.betas.aligned.burnin <- joint.betas.aligned[iters_burnin]

# For each component, calculate the posterior mean of the corresponding rank-1 structure
joint.rank1.structure <- lapply(1:joint.rank, function(r) {
  # Calculate the rank-1 structure for the given component at each iteration
  factor_r_structure <- lapply(1:burnin, function(iter) {
    rbind(joint.loadings.aligned.burnin[[iter]][,r,drop=FALSE], t(joint.betas.aligned.burnin[[iter]][r,,drop=FALSE])) %*% 
      t(joint.scores.aligned.burnin[[iter]][,r,drop=FALSE])
  })
  
  # Calculate the posterior mean
  Reduce("+", factor_r_structure)/length(factor_r_structure)
})

# Calculate the norm of each rank-1 structure
joint.structure.norm <- sapply(joint.rank1.structure, function(str) frob(str))

# Order the factors
joint.factor.order <- order(joint.structure.norm, decreasing = TRUE)

# Reorder the joint scores, loadings, and betas after burn-in
joint.scores.aligned.burnin.order <- lapply(1:burnin, function(iter) {
  joint.scores.aligned.burnin[[iter]][,joint.factor.order,drop=FALSE]
})
joint.loadings.aligned.burnin.order <- lapply(1:burnin, function(iter) {
  joint.loadings.aligned.burnin[[iter]][,joint.factor.order,drop=FALSE]
})
joint.betas.aligned.burnin.order <- lapply(1:burnin, function(iter) {
  joint.betas.aligned.burnin[[iter]][joint.factor.order,,drop=FALSE]
})

# Order the factors in each individual structure by the Frobenius norm of their corresponding rank-1 structure --

# Save the aligned components after burn-in
individual.scores.aligned.burnin <- lapply(1:q, function(s) individual.scores.aligned[[s]][iters_burnin])
individual.loadings.aligned.burnin <- lapply(1:q, function(s) individual.loadings.aligned[[s]][iters_burnin])
individual.betas.aligned.burnin <- lapply(1:q, function(s) individual.betas.aligned[[s]][iters_burnin])

# For each source and each component, calculate the posterior mean of the corresponding rank-1 structure
individual.rank1.structure <- lapply(1:q, function(s) {
  lapply(1:indiv.ranks[s], function(r) {
    # Calculate the rank-1 structure for the given component at each iteration
    factor_r_structure <- lapply(1:burnin, function(iter) {
      rbind(individual.loadings.aligned.burnin[[s]][[iter]][,r,drop=FALSE], t(individual.betas.aligned.burnin[[s]][[iter]][r,,drop=FALSE])) %*% 
        t(individual.scores.aligned.burnin[[s]][[iter]][,r,drop=FALSE])
    })
    
    # Calculate the posterior mean
    Reduce("+", factor_r_structure)/length(factor_r_structure)
  })
})

# Calculate the norm of each rank-1 structure for each source
individual.structure.norm <- lapply(1:q, function(s) {
  sapply(individual.rank1.structure[[s]], function(str) frob(str))
})

# Order the factors
individual.factor.order <- lapply(1:q, function(s) {
  order(individual.structure.norm[[s]], decreasing = TRUE)
})

# Reorder the individual scores, loadings, and betas after burn-in
individual.scores.aligned.burnin.order <- lapply(1:q, function(s) {
  lapply(1:burnin, function(iter) {
    individual.scores.aligned.burnin[[s]][[iter]][,individual.factor.order[[s]],drop=FALSE]
  })
})
individual.loadings.aligned.burnin.order <- lapply(1:q, function(s) {
  lapply(1:burnin, function(iter) {
    individual.loadings.aligned.burnin[[s]][[iter]][,individual.factor.order[[s]],drop=FALSE]
  })
})
individual.betas.aligned.burnin.order <- lapply(1:q, function(s) {
  lapply(1:burnin, function(iter) {
    individual.betas.aligned.burnin[[s]][[iter]][individual.factor.order[[s]],,drop=FALSE]
  })
})

# Check the joint structure after the algorithm matches the original structure --

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
  rbind(do.call(rbind, fev1pp_training_fit_nonsparse$U.draw[[iter]]), # Joint loadings  
        t(fev1pp_training_fit_nonsparse$beta.draw[[iter]][[1,1]][beta.ind[[1]],])) # Joint regression coefficients
})

# Save joint scores
joint.scores <- lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse$V.draw[[iter]][[1,1]])

# Save the original joint structure
joint.structure.original <- lapply(1:burnin, function(iter) {
  joint.loadings[[iter]] %*% t(joint.scores[[iter]])
})

# Calculate the new structure
joint.structure.new <- lapply(1:burnin, function(iter) {
  rbind(joint.loadings.aligned.burnin.order[[iter]], t(joint.betas.aligned.burnin.order[[iter]])) %*%
    t(joint.scores.aligned.burnin.order[[iter]])
})

# Check
all.equal(joint.structure.original, joint.structure.new) # TRUE!

# Check the individual structure after the algorithm matches the original structure --

# Combine the individual loadings and betas
individual.loadings <- lapply(1:q, function(s) lapply(iters_burnin, function(iter) {
  rbind(fev1pp_training_fit_nonsparse$W.draw[[iter]][[s,s]], 
        fev1pp_training_fit_nonsparse$beta.draw[[iter]][[1,1]][beta.ind[[s+1]],])
}))
individual.scores <- lapply(1:q, function(s) lapply(iters_burnin, function(iter) fev1pp_training_fit_nonsparse$Vs.draw[[iter]][[1,s]]))

# Calculate the original structure
individual.structure.original <- lapply(1:q, function(s) {
  lapply(1:burnin, function(iter) {
    individual.loadings[[s]][[iter]] %*% t(individual.scores[[s]][[iter]])
  })
})

# Calculate the new structure
individual.structure.new <- lapply(1:q, function(s) {
  lapply(1:burnin, function(iter) {
    rbind(individual.loadings.aligned.burnin.order[[s]][[iter]], t(individual.betas.aligned.burnin.order[[s]][[iter]])) %*%
      t(individual.scores.aligned.burnin.order[[s]][[iter]])
  })
})

# Check
lapply(1:q, function(s) all.equal(individual.structure.original[[s]], individual.structure.new[[s]])) # TRUE TRUE!

# Rename
joint.scores.final <- joint.scores.aligned.burnin.order
joint.loadings.final <- joint.loadings.aligned.burnin.order
joint.betas.final <- joint.betas.aligned.burnin.order

individual.scores.final <- individual.scores.aligned.burnin.order
individual.loadings.final <- individual.loadings.aligned.burnin.order
individual.betas.final <- individual.betas.aligned.burnin.order

# Save the results
save(joint.scores.final, joint.loadings.final, joint.betas.final, individual.scores.final, individual.loadings.final, individual.betas.final,
     file = "/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/FEV1pp_Joint_Individual_Structures_Factor_Switching_Ordered.rda")

# -----------------------------------------------------------------------------
# Create PCA-like plots using ALIGNED factors from non-sparse model.
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
load("/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/FEV1pp_Joint_Individual_Structures_Factor_Switching_Ordered.rda", verbose = TRUE)

# -------------------------------------
# Joint structure
# -------------------------------------

# Calculate the posterior mean of the joint scores
joint.scores.final.mean <- Reduce("+", lapply(1:burnin, function(iter) joint.scores.final[[iter]]))/(nsample/2)

# Plot subjects against every combination of PCs
joint_rank <- fev1pp_training_fit_nonsparse$ranks[1]
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Joint_Structure_Aligned_Ordered_PCA_Plot.pdf"), width = 15)
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
indiv_rank1 <- fev1pp_training_fit_nonsparse$ranks[2]
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Individual_Structure_Biocrates_Aligned_Ordered_PCA_Plot.pdf"), width = 15)
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
indiv_rank2 <- fev1pp_training_fit_nonsparse$ranks[3]
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Individual_Structure_Somascan_Aligned_Ordered_PCA_Plot.pdf"), width = 15)
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
# Heatmaps 
# * Initially using JIVE heatmaps
# * Later, switched to sJIVE heatmaps to include response
# -----------------------------------------------------------------------------

# -------------------------------------
# Heatmaps using r.jive package
# -------------------------------------

# Save the ranks
ranks <- fev1pp_training_fit_nonsparse$ranks
joint_rank <- ranks[1]
indiv_rank1 <- ranks[2]
indiv_rank2 <- ranks[3]

# Calculating the posterior mean of the joint structure (Biocrates)
# joint_biocrates_mean <- Reduce("+", lapply(iters_burnin, function(iter) {
#   fev1pp_training_fit_nonsparse$U.draw[[iter]][[1,1]] %*% t(fev1pp_training_fit_nonsparse$V.draw[[iter]][[1,1]])
# }))/length(iters_burnin)
# 
# # Calculating the posterior mean of the joint structure (Somascan)
# joint_somascan_mean <- Reduce("+", lapply(iters_burnin, function(iter) {
#   fev1pp_training_fit_nonsparse$U.draw[[iter]][[2,1]] %*% t(fev1pp_training_fit_nonsparse$V.draw[[iter]][[1,1]])
# }))/length(iters_burnin)
# 
# # Calculating the posterior mean of the individual structure (Biocrates)
# indiv_biocrates_mean <- Reduce("+", lapply(iters_burnin, function(iter) {
#   fev1pp_training_fit_nonsparse$W.draw[[iter]][[1,1]] %*% t(fev1pp_training_fit_nonsparse$Vs.draw[[iter]][[1,1]])
# }))/length(iters_burnin)
# 
# # Calculating the posterior mean of the individual structure (Somascan)
# indiv_somascan_mean <- Reduce("+", lapply(iters_burnin, function(iter) {
#   fev1pp_training_fit_nonsparse$W.draw[[iter]][[2,2]] %*% t(fev1pp_training_fit_nonsparse$Vs.draw[[iter]][[1,2]])
# }))/length(iters_burnin)

# Calculating the joint structure mean
joint_structure_mean <- Reduce("+", lapply(1:burnin, function(iter) {
  joint.loadings.final[[iter]] %*% t(joint.scores.final[[iter]])
}))/burnin

metabolite_structure_mean <-  Reduce("+", lapply(1:burnin, function(iter) {
  individual.loadings.final[[1]][[iter]] %*% t(individual.scores.final[[1]][[iter]])
}))/burnin

protein_structure_mean <-  Reduce("+", lapply(1:burnin, function(iter) {
  individual.loadings.final[[2]][[iter]] %*% t(individual.scores.final[[2]][[iter]])
}))/burnin

# Calculating E(Y) for joint contribution
EY_joint_mean <- Reduce("+", lapply(1:burnin, function(iter) {
  joint.scores.final[[iter]] %*% joint.betas.final[[iter]]
}))/burnin

EY_indiv_mean <-lapply(1:q, function(s) {
  Reduce("+", lapply(1:burnin, function(iter) {
    individual.scores.final[[s]][[iter]] %*% individual.betas.final[[s]][[iter]]
  }))/burnin
})

# Reorganizing the results to fit the JIVE structure
heatmap_nonsparse_jive <- list(data = list(hiv_copd_data[[1,1]], hiv_copd_data[[2,1]]),
                              joint = list(joint_structure_mean[1:p.vec[1],,drop=FALSE], 
                                           joint_structure_mean[(p.vec[1]+1):p,,drop=FALSE]),
                              individual = list(metabolite_structure_mean, protein_structure_mean),
                              rankJ = joint_rank,
                              rankA = c(indiv_rank1, indiv_rank2))

# Applying the showHeatmaps function (ordered by joint structure)
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Heatmap_Ordered_by_Joint.pdf"), width = 15)
showHeatmaps(heatmap_nonsparse_jive, order_by = 0)
dev.off()

# Applying the showHeatmaps function (ordered by Biocrates structure)
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Heatmap_Ordered_by_Biocrates.pdf"), width = 15)
showHeatmaps(heatmap_nonsparse_jive, order_by = 1)
dev.off()

# Applying the showHeatmaps function (ordered by Somascan structure)
pdf(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Heatmap_Ordered_by_Somascan.pdf"), width = 15)
showHeatmaps(heatmap_nonsparse_jive, order_by = 2)
dev.off()

# -------------------------------------
# Heatmaps using sup.r.jive
# -------------------------------------

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

png(paste0("~/BayesianPMF/04DataApplication/BPMF/Figures/Heatmap_Ordered_by_FEV1pp_with_Outcome.png"), width = 800)
plotHeatmap(heatmap_sup.r.jive, ylab = "FEV1pp", xlab = c("Metabolomics", "Proteomics"))
dev.off()

# -----------------------------------------------------------------------------
# Credible Intervals for Coefficients Using Non-Sparse Model
# -----------------------------------------------------------------------------

# -------------------------------------
# Creating summaries for ALIGNED factors
# -------------------------------------

# Load in the aligned factors
load("/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Joint_Individual_Structures_Factor_Switching_Ordered.rda", verbose = TRUE)

library(dplyr)

# Save the working directory where the summaries will go
exploring_factors_wd <- "~/BayesianPMF/04DataApplication/BPMF/Exploring_Factors/"

# Save the joint and individual ranks
ranks <- fev1pp_training_fit_nonsparse$ranks

# Save the rank indices
rank.inds <- lapply(1:(q+1), function(s) {
  if (s == 1) {
    1:ranks[s]
  } else {
    (cumsum(ranks[1:(s-1)])[s-1] + 1):cumsum(ranks[1:s])[s]
  }
})

# For each factor, 
for (i in 1:sum(ranks)) {
  
  # Save the results for the given factor
  if (i %in% rank.inds[[1]]) {
    rank_index <- i
    factor.final <- do.call(cbind, lapply(1:burnin, function(iter) joint.loadings.final[[iter]][,rank_index,drop=FALSE]))
    
    # Calculating means and CIs for each loading
    factor.final.summary <- t(apply(factor.final, 1, function(load) {
      c(mean = mean(load), lower.ci = quantile(load, 0.025), upper.ci = quantile(load, 0.975))
    }))
    
    # Name the rows by the respective biomarker
    rownames(factor.final.summary) <- c(rownames(lavage_processed_no_info_log_scale), rownames(somascan_normalized_clean_no_info_transpose_scale))
    
    # Separate the loadings by metabolites and proteins
    factor.final.summary.metabolites <- factor.final.summary[1:p.vec[1],,drop=FALSE]
    factor.final.summary.proteins <- factor.final.summary[(p.vec[1]+1):p,,drop=FALSE]
    
    # Order the loadings within the metabolites and proteins by the posterior mean
    factor.final.summary.metabolites.order <- 
      factor.final.summary.metabolites[order(factor.final.summary.metabolites[,1], decreasing = FALSE),]
    
    factor.final.summary.proteins.order <- 
      factor.final.summary.proteins[order(factor.final.summary.proteins[,1], decreasing = FALSE),]
    
    # Create the file name
    file_name <- paste0(exploring_factors_wd, "Aligned/Joint_Factor", rank_index, "_Aligned_Ordered_Summary.rda")
    
    # Save the results
    save(factor.final.summary, factor.final.summary.metabolites.order, factor.final.summary.proteins.order, file = file_name)
  }
  
  if (i %in% rank.inds[[2]]) {
    rank_index <- i - ranks[1]
    factor.final <- do.call(cbind, lapply(1:burnin, function(iter) individual.loadings.final[[1]][[iter]][,rank_index,drop=FALSE]))
    
    # Calculating means and CIs for each loading
    factor.final.summary <- t(apply(factor.final, 1, function(load) {
      c(mean = mean(load), lower.ci = quantile(load, 0.025), upper.ci = quantile(load, 0.975))
    }))
    
    # Name the rows by the respective biomarker
    rownames(factor.final.summary) <- c(rownames(lavage_processed_no_info_log_scale))
    
    # Separate the loadings by metabolites and proteins (only metabolites here)
    factor.final.summary.metabolites <- factor.final.summary
    factor.final.summary.proteins <- NA
    
    # Order the loadings within the metabolites and proteins by the posterior mean
    factor.final.summary.metabolites.order <- 
      factor.final.summary.metabolites[order(factor.final.summary.metabolites[,1], decreasing = FALSE),]
    
    # Create a file name
    file_name <- paste0(exploring_factors_wd, "Aligned/Metabolite_Indiv_Factor", rank_index, "_Aligned_Ordered_Summary.rda")
    
    # Save the results
    save(factor.final.summary.metabolites.order, file = file_name)
  }
  
  if (i %in% rank.inds[[3]]) {
    rank_index <- i - cumsum(ranks[1:2])[2]
    factor.final <- do.call(cbind, lapply(1:burnin, function(iter) individual.loadings.final[[2]][[iter]][,rank_index,drop=FALSE]))
    
    # Calculating means and CIs for each loading
    factor.final.summary <- t(apply(factor.final, 1, function(load) {
      c(mean = mean(load), lower.ci = quantile(load, 0.025), upper.ci = quantile(load, 0.975))
    }))
    
    # Name the rows by the respective biomarker
    rownames(factor.final.summary) <- rownames(somascan_normalized_clean_no_info_transpose_scale)
    
    # Separate the loadings by metabolites and proteins (only proteins here)
    factor.final.summary.metabolites <- NA
    factor.final.summary.proteins <- factor.final.summary
    
    # Order the loadings within the metabolites and proteins by the posterior mean
    factor.final.summary.proteins.order <- 
      factor.final.summary.proteins[order(factor.final.summary.proteins[,1], decreasing = FALSE),]
    
    # Create the file name
    file_name <- paste0(exploring_factors_wd, "Aligned/Protein_Indiv_Factor", rank_index, "_Aligned_Ordered_Summary.rda")
    
    # Save the results
    save(factor.final.summary.proteins.order, file = file_name)
  }
  
}

# -------------------------------------
# Creating summaries of the UNALIGNED factors
# -------------------------------------

# Adding a burnin
U.draw.burnin <- fev1pp_training_fit_nonsparse$U.draw[(burnin+1):nsample]
W.draw.burnin <- fev1pp_training_fit_nonsparse$W.draw[(burnin+1):nsample]

# For each factor,
for (i in 1:sum(ranks)) {
  
  # Save the results for the given factor
  if (i %in% rank.inds[[1]]) {
    rank_index <- i
    factor.final.by.source <- lapply(1:q, function(s) do.call(cbind, lapply(1:burnin, function(iter) U.draw.burnin[[iter]][[s,1]][,rank_index,drop=FALSE])))
    
    # Calculating means and CIs for each loading
    factor.final.summary.by.source <- lapply(1:q, function(s) t(apply(factor.final.by.source[[s]], 1, function(load) {
      c(mean = mean(load), lower.ci = quantile(load, 0.025), upper.ci = quantile(load, 0.975))
    })))
    
    # Name the rows by the respective biomarker
    rownames(factor.final.summary.by.source[[1]]) <- rownames(lavage_processed_no_info_log_scale)
    rownames(factor.final.summary.by.source[[2]]) <- rownames(somascan_normalized_clean_no_info_transpose_scale)
    
    # Separate the loadings by metabolites and proteins
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
# Cross-Validated Model Fit
# -----------------------------------------------------------------------------

# Saving the index for each pair 
ind_of_pairs <- seq(1, n, by = 2)
n_pair <- length(ind_of_pairs)

# For running in parallel
funcs <- c("bpmf_data", "center_data", "bpmf", "get_results", "BIDIFAC",
           "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
           "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum")
packs <- c("MASS", "truncnorm", "EnvStats", "svMisc", "Matrix")

# -------------------------------------
# Running BPMF with cross validation
# -------------------------------------

cl <- makeCluster(3)
registerDoParallel(cl)
fev1pp_cv <- foreach(pair = ind_of_pairs, .packages = packs, .export = funcs, .verbose = TRUE) %dopar% {
  # Create a new vector of the outcome with the current pair set to NA
  fev1pp_cv <- fev1pp
  fev1pp_cv[[1,1]][pair:(pair+1),] <- NA
  
  # Run the model with the above pair's continuous outcome missing
  fev1pp_cv_fit_nonsparse <- bpmf(
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
  save(Ym.draw_pair, ranks, convergence, file = paste0(results_wd, "FEV1pp_CV_NonSparse_Pair", pair, ".rda"))
  
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

# -----------------------------------------------------------------------------
# Cross-Validated Model Fit Results
# -----------------------------------------------------------------------------

# Create a vector with model names
models <- c("BPMF", "BIDIFAC", "JIVE", "MOFA")

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
      load(paste0(results_wd, "BIDIFAC/FEV1pp_CV_BIDIFAC_Pair_", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "JIVE") {
      load(paste0(results_wd, "JIVE/FEV1pp_CV_JIVE_Pair_", pair, ".rda"), verbose = TRUE)
    }
    if (mod == "MOFA") {
      load(paste0(results_wd, "MOFA/FEV1pp_CV_MOFA_Pair_", pair, ".rda"), verbose = TRUE)
    }
    
    # Combine the samples
    samps <- do.call(cbind, do.call(cbind, Ym.draw_pair))
    
    # Take a burn-in
    samps_burnin <- samps[,thinned_iters_burnin]
    
    # Save in the vector
    fev1pp_cv[[mod]][pair:(pair+1)] <- rowMeans(samps_burnin)
    
    # Calculate the 95% credible interval for each held-out individual
    fev1pp_cv_ci[[mod]][pair,] <- c(quantile(samps_burnin[1,], 0.025), quantile(samps_burnin[1,], 0.975))
    fev1pp_cv_ci[[mod]][pair+1,] <- c(quantile(samps_burnin[2,], 0.025), quantile(samps_burnin[2,], 0.975))
  }
}

# Plotting the results from each model against the truth
plot(fev1pp_cv$BPMF, c(fev1pp[[1,1]]), xlab = "Predicted FEV1pp", ylab = "Observed FEV1pp", main = "Cross-Validated FEV1pp vs. True FEV1pp", pch = 16)
points(fev1pp_cv$BIDIFAC, c(fev1pp[[1,1]]), col = 2, pch = 16) # BIDIFAC
points(fev1pp_cv$JIVE, c(fev1pp[[1,1]]), col = 3, pch = 16) # JIVE
points(fev1pp_cv$MOFA, c(fev1pp[[1,1]]), col = 4, pch = 16) # MOFA
abline(a=0, b=1, lwd = 2)
legend("bottomright", legend = c("BPMF", "BIDIFAC", "JIVE", "MOFA"), col = c(1, 2, 3, 4), pch = rep(16, 4), cex = 0.5)

# Calculating a correlation test for each cross validated outcome vs. true FEV1pp
cor.test(fev1pp_cv$BPMF, c(fev1pp[[1,1]])) # BPMF
cor.test(fev1pp_cv$BIDIFAC, c(fev1pp[[1,1]])) # BIDIFAC
cor.test(fev1pp_cv$JIVE, c(fev1pp[[1,1]])) # JIVE
cor.test(fev1pp_cv$MOFA, c(fev1pp[[1,1]])) # MOFA

# Compare the coverage rates for FEV1pp
mean(sapply(1:n, function(i) fev1pp_cv_ci$BPMF[i,1] <= fev1pp[[1,1]][i] & fev1pp[[1,1]][i] <= fev1pp_cv_ci$BPMF[i,2])) # BPMF
mean(sapply(1:n, function(i) fev1pp_cv_ci$BIDIFAC[i,1] <= fev1pp[[1,1]][i] & fev1pp[[1,1]][i] <= fev1pp_cv_ci$BIDIFAC[i,2])) # BIDIFAC
mean(sapply(1:n, function(i) fev1pp_cv_ci$JIVE[i,1] <= fev1pp[[1,1]][i] & fev1pp[[1,1]][i] <= fev1pp_cv_ci$JIVE[i,2])) # JIVE
mean(sapply(1:n, function(i) fev1pp_cv_ci$MOFA[i,1] <= fev1pp[[1,1]][i] & fev1pp[[1,1]][i] <= fev1pp_cv_ci$MOFA[i,2])) # MOFA


# -----------------------------------------------------------------------------
# Missing data imputation
# -----------------------------------------------------------------------------

# Setting the proportion of missing values (columns) in each dataset
prop_missing <- 0.1

# Save the number of replications
nsim <- 100

# For 100 replications of entrywise missingness:
cl <- makeCluster(2)
registerDoParallel(cl)
# For running in parallel
funcs <- c("bpmf_data", "center_data", "bpmf_data_mode", "get_results", "BIDIFAC",
           "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
           "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum")
packs <- c("MASS", "truncnorm", "EnvStats", "svMisc", "Matrix")
hiv_copd_imputation <- foreach(sim_iter = 13:20, .packages = packs, .export = funcs, .verbose = TRUE) %dopar% { 
  
  # Set a seed
  set.seed(sim_iter)
  
  # Randomly sample samples to remove from each source
  metabolome_missing <- sort(sample(1:length(hiv_copd_data[[1,1]]), size = prop_missing * length(hiv_copd_data[[1,1]]), replace = FALSE))
  proteome_missing <- sort(sample(1:length(hiv_copd_data[[2,1]]), size = prop_missing * length(hiv_copd_data[[2,1]]), replace = FALSE))
  
  # Setting these samples to missing 
  hiv_copd_data_missing <- hiv_copd_data
  hiv_copd_data_missing[[1,1]][metabolome_missing] <- NA
  hiv_copd_data_missing[[2,1]][proteome_missing] <- NA
  
  # Save the missing samples
  metabolome_missing_samples <- hiv_copd_data[[1,1]][metabolome_missing]
  proteome_missing_samples <- hiv_copd_data[[2,1]][proteome_missing]
  
  # Running the model
  bpmf_impute <- bpmf_data_mode(
      data = hiv_copd_data_missing,
      Y = fev1pp,
      nninit = TRUE,
      model_params = model_params,
      sparsity = FALSE,
      nsample = nsample,
      progress = TRUE
  )
  
  # Saving the imputed values for each source with burn-in
  metabolome_impute_burnin <- lapply(burnin:nsample, function(iter) {
    matrix(bpmf_impute$Xm.draw[[iter]][[1,1]], nrow = 1) 
  })
  
  proteome_impute_burnin <- lapply(burnin:nsample, function(iter) {
    matrix(bpmf_impute$Xm.draw[[iter]][[2,1]], nrow = 1)
  })
  
  # Create a matrix of imputed values
  metabolome_impute_burnin <- do.call(rbind, metabolome_impute_burnin)
  proteome_impute_burnin <- do.call(rbind, proteome_impute_burnin)
  
  # Calculate credible intervals for each imputed values
  metabolome_cis <- apply(metabolome_impute_burnin, 2, function(row) {
    c(quantile(row, 0.025), quantile(row, 0.975))
  })
  
  proteome_cis <- apply(proteome_impute_burnin, 2, function(row) {
    c(quantile(row, 0.025), quantile(row, 0.975))
  })
  
  # Calculate coverage
  metabolome_coverage <- mean(sapply(1:length(metabolome_missing), function(i) {
    metabolome_cis[1,i] <= metabolome_missing_samples[i] & metabolome_missing_samples[i] <= metabolome_cis[2,i]
  }))
  
  proteome_coverage <- mean(sapply(1:length(proteome_missing), function(i) {
    proteome_cis[1,i] <= proteome_missing_samples[i] & proteome_missing_samples[i] <= proteome_cis[2,i]
  }))
  
  # Calculate the posterior mean of the imputed values
  metabolome_impute_mean <- colMeans(metabolome_impute_burnin)
  proteome_impute_mean <- colMeans(proteome_impute_burnin)
  
  # Calculate MSE with posterior mean and true observed values
  metabolome_mse <- frob(metabolome_missing_samples - metabolome_impute_mean)/frob(metabolome_missing_samples)
  proteome_mse <- frob(proteome_missing_samples - proteome_impute_mean)/frob(proteome_missing_samples)
  
  # Calculate the CI width
  metabolome_ci_width <- mean(abs(metabolome_cis[2,] - metabolome_cis[1,]))
  proteome_ci_width <- mean(abs(proteome_cis[2,] - proteome_cis[1,]))
  
  # Save the results
  save(metabolome_missing, proteome_missing,
       metabolome_impute_mean, proteome_impute_mean,
       metabolome_coverage, proteome_coverage, 
       metabolome_ci_width, proteome_ci_width, 
       metabolome_mse, proteome_mse,
       file = paste0("~/BayesianPMF/04DataApplication/BPMF/Imputation/Entrywise/", "FEV1pp_", prop_missing, "_Entrywise_Impute_", sim_iter, ".rda"))
  
  # Remove large matrices and clean
  rm(bpmf_impute, metabolome_cis, proteome_cis)
  gc()
}

# Load in the results so far
entrywise_imputation_files <- list.files("~/BayesianPMF/04DataApplication/BPMF/Imputation/Entrywise")

# Create a dataframe with the results
entrywise_imputation_results <- matrix(nrow = 20, ncol = 6)

# Iteratively load in each replication and add to table
for (i in 1:length(entrywise_imputation_files)) {
  # Load in the file
  file <- entrywise_imputation_files[i]
  load(file)
  
  # Add results to table
  entrywise_imputation_results[i,] <- 
    c(metabolome_coverage, metabolome_mse, metabolome_ci_width,
      proteome_coverage, proteome_mse, proteome_ci_width)
}

# Add the column names
colnames(entrywise_imputation_results) <-
  c("Metabolome Coverage", "Metabolome MSE", "Metabolome CI Width",
    "Proteome Coverage", "Proteome MSE", "Proteome CI Width")

# Average the results
colMeans(entrywise_imputation_results)

# For 100 replications of COLUMNWISE missingness:
cl <- makeCluster(2)
registerDoParallel(cl)
# For running in parallel
funcs <- c("bpmf_data", "center_data", "bpmf_data_mode", "get_results", "BIDIFAC",
           "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
           "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum")
packs <- c("MASS", "truncnorm", "EnvStats", "svMisc", "Matrix")
hiv_copd_imputation <- foreach(sim_iter = 1:20, .packages = packs, .export = funcs, .verbose = TRUE) %dopar% { 
  
  # Set a seed
  set.seed(sim_iter)
  
  # Randomly sample samples to remove from each source
  metabolome_obs_missing <- sort(sample(1:n, size = prop_missing * n, replace = FALSE))
  
  # To prevent the same sample missing from both sources
  avail_obs <- c(1:n)[!(c(1:n) %in% metabolome_obs_missing)]
  
  proteome_obs_missing <- sort(sample(avail_obs, size = prop_missing * n, replace = FALSE))
  
  # Setting these samples to missing 
  hiv_copd_data_missing <- hiv_copd_data
  hiv_copd_data_missing[[1,1]][,metabolome_obs_missing] <- NA
  hiv_copd_data_missing[[2,1]][,proteome_obs_missing] <- NA
  
  # Save the entries that are NA
  metabolome_missing <- which(is.na(hiv_copd_data_missing[[1,1]]))
  proteome_missing <- which(is.na(hiv_copd_data_missing[[2,1]]))
  
  # Save the missing samples
  metabolome_missing_samples <- hiv_copd_data[[1,1]][metabolome_missing]
  proteome_missing_samples <- hiv_copd_data[[2,1]][proteome_missing]
  
  # Running the model
  bpmf_impute <- bpmf_data_mode(
    data = hiv_copd_data_missing,
    Y = fev1pp,
    nninit = TRUE,
    model_params = model_params,
    sparsity = FALSE,
    nsample = nsample,
    progress = TRUE
  )
  
  # Saving the imputed values for each source with burn-in
  metabolome_impute_burnin <- lapply(burnin:nsample, function(iter) {
    matrix(bpmf_impute$Xm.draw[[iter]][[1,1]], nrow = 1) 
  })
  
  proteome_impute_burnin <- lapply(burnin:nsample, function(iter) {
    matrix(bpmf_impute$Xm.draw[[iter]][[2,1]], nrow = 1)
  })
  
  # Create a matrix of imputed values
  metabolome_impute_burnin <- do.call(rbind, metabolome_impute_burnin)
  proteome_impute_burnin <- do.call(rbind, proteome_impute_burnin)
  
  # Calculate credible intervals for each imputed values
  metabolome_cis <- apply(metabolome_impute_burnin, 2, function(row) {
    c(quantile(row, 0.025), quantile(row, 0.975))
  })
  
  proteome_cis <- apply(proteome_impute_burnin, 2, function(row) {
    c(quantile(row, 0.025), quantile(row, 0.975))
  })
  
  # Calculate coverage
  metabolome_coverage <- mean(sapply(1:length(metabolome_missing), function(i) {
    metabolome_cis[1,i] <= metabolome_missing_samples[i] & metabolome_missing_samples[i] <= metabolome_cis[2,i]
  }))
  
  proteome_coverage <- mean(sapply(1:length(proteome_missing), function(i) {
    proteome_cis[1,i] <= proteome_missing_samples[i] & proteome_missing_samples[i] <= proteome_cis[2,i]
  }))
  
  # Calculate the posterior mean of the imputed values
  metabolome_impute_mean <- colMeans(metabolome_impute_burnin)
  proteome_impute_mean <- colMeans(proteome_impute_burnin)
  
  # Calculate MSE with posterior mean and true observed values
  metabolome_mse <- frob(metabolome_missing_samples - metabolome_impute_mean)/frob(metabolome_missing_samples)
  proteome_mse <- frob(proteome_missing_samples - proteome_impute_mean)/frob(proteome_missing_samples)
  
  # Calculate the CI width
  metabolome_ci_width <- mean(abs(metabolome_cis[2,] - metabolome_cis[1,]))
  proteome_ci_width <- mean(abs(proteome_cis[2,] - proteome_cis[1,]))
  
  # Save the results
  save(metabolome_missing, proteome_missing,
       metabolome_impute_mean, proteome_impute_mean,
       metabolome_coverage, proteome_coverage, 
       metabolome_ci_width, proteome_ci_width, 
       metabolome_mse, proteome_mse,
       file = paste0("~/BayesianPMF/04DataApplication/BPMF/Imputation/Columnwise/", "FEV1pp_", prop_missing, "_Columnwise_Impute_", sim_iter, ".rda"))
  
  # Remove large matrices and clean
  rm(bpmf_impute, metabolome_cis, proteome_cis)
  gc()
}

# Load in the results so far
columnwise_imputation_files <- list.files("~/BayesianPMF/04DataApplication/BPMF/Imputation/Columnwise")

# Create a dataframe with the results
columnwise_imputation_results <- matrix(nrow = 20, ncol = 6)

# Iteratively load in each replication and add to table
for (i in 1:length(columnwise_imputation_files)) {
  # Load in the file
  file <- columnwise_imputation_files[i]
  load(file)
  
  # Add results to table
  columnwise_imputation_results[i,] <- 
    c(metabolome_coverage, metabolome_mse, metabolome_ci_width,
      proteome_coverage, proteome_mse, proteome_ci_width)
}

# Add the column names
colnames(columnwise_imputation_results) <-
  c("Metabolome Coverage", "Metabolome MSE", "Metabolome CI Width",
    "Proteome Coverage", "Proteome MSE", "Proteome CI Width")

# Average the results
colMeans(columnwise_imputation_results)