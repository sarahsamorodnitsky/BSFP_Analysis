# -----------------------------------------------------------------------------
# Applying the BPMF model to the HIV-COPD BALF Biocrates and Somascan data. 
# Comparing the model application when there is no sparsity vs. when
# there is sparsity. 
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

# Loading data in 
load(paste0(data_wd, "BiocratesLavageProcessed.rda"))
load(paste0(data_wd, "HIV_COPD_SomaScan_Normalized_Clean.rda"))

# BALF Biocrates
lavage_processed_no_info <- lavage_processed[, -c(1:2)] # removing meta
bioc_meta_data <- lavage_processed[, 1:2] # saving meta elsewhere
lavage_processed_no_info_log_scale <- t(apply(lavage_processed_no_info, 1, function(row) {
  scale(log(row + 1), center = TRUE, scale = TRUE)
})) # scaling the rows
rownames(lavage_processed_no_info_log_scale) <- bioc_meta_data$Metabolite
colnames(lavage_processed_no_info_log_scale) <- colnames(lavage_processed_no_info)

# BALF Somascan
somascan_normalized_clean_no_info <- somascan_normalized_clean[, -c(1:29)] # removing meta
soma_meta_data <- somascan_normalized_clean[, 1:29] # saving meta elsewhere
soma_pids <- somascan_normalized_clean$PID
somascan_normalized_clean_no_info_transpose <- t(somascan_normalized_clean_no_info)
somascan_normalized_clean_no_info_transpose_scale <- t(apply(somascan_normalized_clean_no_info_transpose, 1, function(row) {
  scale(log(row + 1), center = TRUE, scale = TRUE)
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

# Preparing the response variable
fev1pp <- matrix(list(), nrow = 1, ncol = 1)
fev1pp[[1,1]] <- matrix(subject_processed$FEV1_percent_predicted, ncol = 1)

# -----------------------------------------------------------------------------
# Model parameters
# -----------------------------------------------------------------------------

# Data dimensions
p1 <- nrow(lavage_processed_no_info_log_scale)
p2 <- nrow(somascan_normalized_clean_no_info_transpose_reorder)
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
thinned_iters <- seq(1, nsample, by = 10)
thinned_iters_burnin <- seq(burnin, nsample, by = 10)

# -----------------------------------------------------------------------------
# Training Data Model Fit
# -----------------------------------------------------------------------------

# Fitting the model without sparsity first
fev1pp_training_fit_nonsparse <- bpmf(
  data = hiv_copd_data,
  Y = fev1pp,
  nninit = TRUE,
  model_params = model_params,
  sparsity = FALSE,
  nsample = nsample,
  progress = TRUE
)

# Fitting the model with sparsity
fev1pp_training_fit_sparse <- bpmf(
  data = hiv_copd_data,
  Y = fev1pp,
  nninit = TRUE,
  model_params = model_params,
  sparsity = TRUE,
  nsample = nsample,
  progress = TRUE
)

# Save the results
save(fev1pp_training_fit_nonsparse, fev1pp_training_fit_sparse, file = paste0(results_wd, "training_data_fit.rda"))

# -----------------------------------------------------------------------------
# Investigating training data fit results
# -----------------------------------------------------------------------------

# -------------------------------------
# Convergence
# -------------------------------------

# Assessing convergence for both model fits
load(paste0(results_wd, "training_data_fit.rda"))

fev1pp_training_nonsparse_conv <- sapply(thinned_iters, function(sim_iter) {
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

fev1pp_training_sparse_conv <- sapply(thinned_iters, function(sim_iter) {
  # Calculate the log-joint density at each thinned iterations
  log_joint_density(data = hiv_copd_data, 
                    U.iter = fev1pp_training_fit_sparse$U.draw[[sim_iter]], 
                    V.iter = fev1pp_training_fit_sparse$V.draw[[sim_iter]], 
                    W.iter = fev1pp_training_fit_sparse$W.draw[[sim_iter]], 
                    Vs.iter = fev1pp_training_fit_sparse$Vs.draw[[sim_iter]],
                    model_params = model_params,
                    ranks = fev1pp_training_fit_sparse$ranks,
                    Y = fev1pp,
                    beta.iter = fev1pp_training_fit_sparse$beta.draw[[sim_iter]],
                    tau2.iter = fev1pp_training_fit_sparse$tau2.draw[[sim_iter]])
})

# Plotting the log-joint densities
plot(fev1pp_training_nonsparse_conv, ylab = "Log Joint Density", main = "Log Joint Density for Non-Sparse Model")
plot(fev1pp_training_sparse_conv, ylab = "Log Joint Density", main = "Log Joint Density for Sparse Model")

# Applying the label switching algorithm to the training data fits
fev1pp_training_fit_sparse_ls <- label_switching(U.draw = fev1pp_training_fit_sparse$U.draw,
                                                 V.draw = fev1pp_training_fit_sparse$V.draw,
                                                 W.draw = fev1pp_training_fit_sparse$W.draw,
                                                 Vs.draw = fev1pp_training_fit_sparse$Vs.draw,
                                                 betas = fev1pp_training_fit_sparse$beta.draw,
                                                 gammas = fev1pp_training_fit_sparse$gamma.draw,
                                                 r = fev1pp_training_fit_sparse$ranks[1],
                                                 r.vec = fev1pp_training_fit_sparse$ranks[-1],
                                                 nsample = nsample,
                                                 thinned_iters_burnin = thinned_iters_burnin)

fev1pp_training_fit_nonsparse_ls <- label_switching(U.draw = fev1pp_training_fit_nonsparse$U.draw,
                                                    V.draw = fev1pp_training_fit_nonsparse$V.draw,
                                                    W.draw = fev1pp_training_fit_nonsparse$W.draw,
                                                    Vs.draw = fev1pp_training_fit_nonsparse$Vs.draw,
                                                    betas = fev1pp_training_fit_nonsparse$beta.draw,
                                                    gammas = NULL,
                                                    r = fev1pp_training_fit_nonsparse$ranks[1],
                                                    r.vec = fev1pp_training_fit_nonsparse$ranks[-1],
                                                    nsample = nsample,
                                                    thinned_iters_burnin = thinned_iters_burnin)

# Check that label switching correctly swapped the columns
correct_swaps_structure <- correct_swaps_lm <- c()
for (iter in 1:nsample) {
  correct_swaps_structure_s <- correct_swaps_lm_s <- c()
  for (s in 1:q) {
    # Checking the structure
    before_ls <- fev1pp_training_fit_sparse$U.draw[[iter]][[s,1]] %*% t(fev1pp_training_fit_sparse$V.draw[[iter]][[1,1]]) +
      fev1pp_training_fit_sparse$W.draw[[iter]][[s,s]] %*% t(fev1pp_training_fit_sparse$Vs.draw[[iter]][[1,s]])
    
    after_ls <- fev1pp_training_fit_sparse_ls$swapped_U.draw[[iter]][[s,1]] %*% t(fev1pp_training_fit_sparse_ls$swapped_V.draw[[iter]][[1,1]]) +
      fev1pp_training_fit_sparse_ls$swapped_W.draw[[iter]][[s,s]] %*% t(fev1pp_training_fit_sparse_ls$swapped_Vs.draw[[iter]][[1,s]])
    
   # Checking the gammas * betas
   before_ls_lm <- t(fev1pp_training_fit_sparse$gamma.draw[[iter]][[1,1]]) %*% abs(fev1pp_training_fit_sparse$beta.draw[[iter]][[1,1]])
   
   after_ls_lm <- t(fev1pp_training_fit_sparse_ls$swapped_gammas[[iter]][[1,1]]) %*% abs(fev1pp_training_fit_sparse_ls$swapped_betas[[iter]][[1,1]])
   
   # Checking that everything matched
   correct_swaps_structure_s[s] <- all.equal(before_ls, after_ls)
   correct_swaps_lm_s[s] <- all.equal(before_ls_lm, after_ls_lm)
  }
  # Add the results from both sources
  correct_swaps_structure[iter] <- all(correct_swaps_structure_s)
  correct_swaps_lm[iter] <- all(correct_swaps_lm_s)
}
all(correct_swaps_structure)
all(correct_swaps_lm)

correct_swaps_structure <- correct_swaps_lm <- c()
for (iter in 1:nsample) {
  correct_swaps_structure_s <- correct_swaps_lm_s <- c()
  for (s in 1:q) {
    # Checking the structure
    before_ls <- fev1pp_training_fit_nonsparse$U.draw[[iter]][[s,1]] %*% t(fev1pp_training_fit_nonsparse$V.draw[[iter]][[1,1]]) +
      fev1pp_training_fit_nonsparse$W.draw[[iter]][[s,s]] %*% t(fev1pp_training_fit_nonsparse$Vs.draw[[iter]][[1,s]])
    
    after_ls <- fev1pp_training_fit_nonsparse_ls$swapped_U.draw[[iter]][[s,1]] %*% t(fev1pp_training_fit_nonsparse_ls$swapped_V.draw[[iter]][[1,1]]) +
      fev1pp_training_fit_nonsparse_ls$swapped_W.draw[[iter]][[s,s]] %*% t(fev1pp_training_fit_nonsparse_ls$swapped_Vs.draw[[iter]][[1,s]])
    
    # Checking the VStar * betas
    before.VStar.iter <- cbind(1, fev1pp_training_fit_nonsparse$V.draw[[iter]][[1,1]], do.call(cbind, fev1pp_training_fit_nonsparse$Vs.draw[[iter]]))
    before_ls_lm <- before.VStar.iter %*% (fev1pp_training_fit_nonsparse$beta.draw[[iter]][[1,1]])
    
    after.VStar.iter <- cbind(1, fev1pp_training_fit_nonsparse_ls$swapped_V.draw[[iter]][[1,1]], do.call(cbind, fev1pp_training_fit_nonsparse_ls$swapped_Vs.draw[[iter]]))
    after_ls_lm <- after.VStar.iter %*% (fev1pp_training_fit_nonsparse_ls$swapped_betas[[iter]][[1,1]])
    
    # Checking that everything matched
    correct_swaps_structure_s[s] <- all.equal(before_ls, after_ls)
    correct_swaps_lm_s[s] <- all.equal(before_ls_lm, after_ls_lm)
  }
  # Add the results from both sources
  correct_swaps_structure[iter] <- all(correct_swaps_structure_s)
  correct_swaps_lm[iter] <- all(correct_swaps_lm_s)
}
all(correct_swaps_structure)
all(correct_swaps_lm)


# -------------------------------------
# Selection of factors in sparse model
# -------------------------------------

# Saving posterior probabilities ordered by variance explained
gammas_ls <- fev1pp_training_fit_sparse_ls$gammas.mean.reorder 
order_joint_ls <- fev1pp_training_fit_sparse_ls$order_joint
order_indiv_ls <- fev1pp_training_fit_sparse_ls$order_indiv
var_exp_joint_ls <- fev1pp_training_fit_sparse_ls$var_exp_joint
var_exp_indiv_ls <- fev1pp_training_fit_sparse_ls$var_exp_indiv

# Creating table
library(xtable)
selection_df <- data.frame(Factor = character(length(gammas_ls) - 1),
                           Post.Prob = numeric(length(gammas_ls) - 1),
                           Var.Exp = numeric(length(gammas_ls) - 1))

# Saving the ranks from the sparse model
ranks_ls <- fev1pp_training_fit_sparse$ranks

# Filling in the table
selection_df$Post.Prob <- gammas_ls[-1,]

# Ordering the names of the factors
joint_factor_names <- paste0("Joint ", 1:ranks_ls[1])
joint_factor_names <- joint_factor_names[order_joint_ls]

indiv_factor_names <- list(paste0("Metabolomics ", 1:ranks_ls[2]),
                           paste0("Proteomics ", 1:ranks_ls[3]))
indiv_factor_names <- lapply(1:q, function(s) {
  indiv_factor_names[[s]][order_indiv_ls[[s,1]]]
})

# Combining them together
factor_names <- c(joint_factor_names, unlist(indiv_factor_names))

# Adding them to the table
selection_df$Factor <- factor_names

# Adding variance explained
var_exp_joint_ls <- var_exp_joint_ls[order_joint_ls]
var_exp_indiv_ls <- lapply(1:q, function(s) {
  var_exp_indiv_ls[[s,1]][order_indiv_ls[[s]]]
})
var_exp <- c(var_exp_joint_ls, unlist(var_exp_indiv_ls))

# Adding them to the table
selection_df$Var.Exp <- var_exp

# xtable
print(xtable(selection_df, digits = 3), include.rownames = FALSE)


# -------------------------------------
# Comparing the training data fits
# -------------------------------------

latent_structure_sparse <- lapply(1:nsample, function(iter) {
  lapply(1:q, function(s) {
    fev1pp_training_fit_sparse$U.draw[[iter]][[s,1]] %*% t(fev1pp_training_fit_sparse$V.draw[[iter]][[1,1]]) + fev1pp_training_fit_sparse$W.draw[[iter]][[s,s]] %*% t(fev1pp_training_fit_sparse$Vs.draw[[iter]][[1,s]])
  })
})
latent_structure_sparse <- list(Reduce("+", lapply(latent_structure_sparse, function(iter) iter[[1]]))/nsample,
                                Reduce("+", lapply(latent_structure_sparse, function(iter) iter[[2]]))/nsample)

latent_structure_nonsparse <- lapply(1:nsample, function(iter) {
  lapply(1:q, function(s) {
    fev1pp_training_fit_nonsparse$U.draw[[iter]][[s,1]] %*% t(fev1pp_training_fit_nonsparse$V.draw[[iter]][[1,1]]) + fev1pp_training_fit_nonsparse$W.draw[[iter]][[s,s]] %*% t(fev1pp_training_fit_nonsparse$Vs.draw[[iter]][[1,s]])
  })
})
latent_structure_nonsparse <- list(Reduce("+", lapply(latent_structure_nonsparse, function(iter) iter[[1]]))/nsample,
                                   Reduce("+", lapply(latent_structure_nonsparse, function(iter) iter[[2]]))/nsample)

frob(hiv_copd_data[[1,1]] - latent_structure_sparse[[1]]) # Sparse model fit for Biocrates
frob(hiv_copd_data[[2,1]] - latent_structure_sparse[[2]]) # Sparse model fit for Somascan

frob(hiv_copd_data[[1,1]] - latent_structure_nonsparse[[1]]) # Sparse model fit for Biocrates
frob(hiv_copd_data[[2,1]] - latent_structure_nonsparse[[2]]) # Sparse model fit for Somascan

# Joint Scores
V_nonsparse <- Reduce("+", lapply(fev1pp_training_fit_nonsparse_ls$swapped_V.draw, function(iter) iter[[1,1]]))/nsample
V_sparse <- Reduce("+", lapply(fev1pp_training_fit_sparse_ls$swapped_V.draw, function(iter) iter[[1,1]]))/nsample

# Individual scores
Vs_nonsparse <- Reduce("+", lapply(fev1pp_training_fit_nonsparse_ls$swapped_Vs.draw, function(iter) do.call(cbind, iter)))/nsample
Vs_sparse <- Reduce("+", lapply(fev1pp_training_fit_sparse_ls$swapped_Vs.draw, function(iter) do.call(cbind, iter)))/nsample

# Betas
betas_nonsparse <- Reduce("+", lapply(fev1pp_training_fit_nonsparse_ls$swapped_betas, function(iter) iter[[1,1]]))/nsample
betas_sparse <- Reduce("+", lapply(fev1pp_training_fit_sparse_ls$swapped_betas, function(iter) iter[[1,1]]))/nsample

# E(Y)
EY_nonsparse <- cbind(1, V_nonsparse, Vs_nonsparse) %*% betas_nonsparse
EY_sparse <- cbind(1, V_sparse, Vs_sparse) %*% betas_sparse

# Compare to FEV1pp
cor.test(EY_nonsparse, fev1pp[[1,1]])
cor.test(EY_sparse, fev1pp[[1,1]])

# Plotting E(Y) against FEV1pp
png("~/BayesianPMF/04DataApplication/Figures/EY_vs_FEV1pp_Nonsparse.png")
par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), xpd = FALSE)
plot(EY_nonsparse, fev1pp[[1,1]], xlab = "Estimated FEV1pp (%)", ylab = "Observed FEV1pp (%)",
     main = "Predicted vs. Observed FEV1pp", pch = 16)
abline(a=0, b=1, lwd = 2, col = viridis(1))
dev.off()

# -----------------------------------------------------------------------------
# PCA-like plots using non-sparse model
# -----------------------------------------------------------------------------

library(viridis)

# Creating colors for the plots
ccstat_colors <- viridis(2)[clinical_data_soma$ccstat]
cluster_colors <- viridis(2)[clusters.soma]

# -------------------------------------
# Joint structure
# -------------------------------------

# Saving the posterior mean of the joint structure, ordered by variance explained 
V_ls_nonsparse_mean_ordered <- fev1pp_training_fit_nonsparse_ls$V.mean.reorder[[1,1]]

# Load in the clustering information
load("~/HIV-COPD/Data/UnivariatePValuesUncombinedClusterComparison.rda", verbose = TRUE)

# Plot subjects against every combination of PCs
joint_rank <- fev1pp_training_fit_nonsparse$ranks[1]
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Joint_Structure_PCA_Plot.pdf"), width = 15)
for (rs1 in 1:joint_rank) {
 for (rs2 in 1:joint_rank) {
   if (rs1 != rs2) {
     par(mfrow = c(1,2), mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
     # Plot with case-control label as colors
     plot(V_ls_nonsparse_mean_ordered[,rs1], V_ls_nonsparse_mean_ordered[,rs2], pch = 16, xlab = paste0("Joint Score ", rs1),
          ylab = paste0("Joint Score ", rs2), main = paste0("Joint Factors ", rs1, " vs ", rs2),
          col = ccstat_colors)
     
     # Plot with stand-out subjects as colors
     plot(V_ls_nonsparse_mean[,rs1], V_ls_nonsparse_mean[,rs2], pch = 16, xlab = paste0("Joint Score ", rs1),
          ylab = paste0("Joint Score ", rs2), main = paste0("Joint Factors ", rs1, " vs ", rs2),
          col = cluster_colors)
     par(mfrow = c(1,1))
   }
 } 
}
dev.off()

# -------------------------------------
# Individual structure
# -------------------------------------

# Saving the reordered individual factors, ordered by variance explained in each source
Vs_ls_nonsparse_mean_ordered <- fev1pp_training_fit_nonsparse_ls$Vs.mean.reorder

# Plot subjects against every combination of PCs for Biocrates
indiv_rank1 <- fev1pp_training_fit_nonsparse$ranks[2]
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Individual_Structure_Biocrates_PCA_Plot.pdf"), width = 15)
for (rs1 in 1:indiv_rank1) {
  for (rs2 in 1:indiv_rank1) {
    if (rs1 != rs2) {
      par(mfrow = c(1,2))
      # Plot with case-control label as colors
      plot(Vs_ls_nonsparse_mean_ordered[[1,1]][,rs1], Vs_ls_nonsparse_mean_ordered[[1,1]][,rs2], pch = 16, xlab = paste0("Individual Score ", rs1),
           ylab = paste0("Individual Score ", rs2), main = paste0("Individual Factors Biocrates ", rs1, " vs ", rs2),
           col = ccstat_colors)
      
      # Plot with stand-out subjects as colors
      plot(Vs_ls_nonsparse_mean_ordered[[1,1]][,rs1], Vs_ls_nonsparse_mean_ordered[[1,1]][,rs2], pch = 16, xlab = paste0("Individual Score ", rs1),
           ylab = paste0("Individual Score ", rs2), main = paste0("Individual Factors Biocrates ", rs1, " vs ", rs2),
           col = cluster_colors)
      par(mfrow = c(1,1))
    }
  } 
}
dev.off()

# Plot subjects against every combination of PCs for Somascan
indiv_rank2 <- fev1pp_training_fit_nonsparse$ranks[3]
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Individual_Structure_Somascan_PCA_Plot.pdf"), width = 15)
for (rs1 in 1:indiv_rank2) {
  for (rs2 in 1:indiv_rank2) {
    if (rs1 != rs2) {
      par(mfrow = c(1,2))
      # Plot with case-control label as colors
      plot(Vs_ls_nonsparse_mean_ordered[[1,2]][,rs1], Vs_ls_nonsparse_mean_ordered[[1,2]][,rs2], pch = 16, xlab = paste0("Individual Score ", rs1),
           ylab = paste0("Individual Score ", rs2), main = paste0("Individual Factors Somascan ", rs1, " vs ", rs2),
           col = ccstat_colors)
      
      # Plot with stand-out subjects as colors
      plot(Vs_ls_nonsparse_mean_ordered[[1,2]][,rs1], Vs_ls_nonsparse_mean_ordered[[1,2]][,rs2], pch = 16, xlab = paste0("Individual Score ", rs1),
           ylab = paste0("Individual Score ", rs2), main = paste0("Individual Factors Somascan ", rs1, " vs ", rs2),
           col = cluster_colors)
      par(mfrow = c(1,1))
    }
  } 
}
dev.off()

# -----------------------------------------------------------------------------
# Heatmaps on non-sparse model
# -----------------------------------------------------------------------------

# Calculating the posterior mean of the joint structure (Biocrates)
joint_biocrates_ls_nonsparse_mean <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_ls$swapped_U.draw[[iter]][[1,1]] %*% t(fev1pp_training_fit_nonsparse_ls$swapped_V.draw[[iter]][[1,1]])
}))/length(thinned_iters_burnin)

# Calculating the posterior mean of the joint structure (Somascan)
joint_somascan_ls_nonsparse_mean <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_ls$swapped_U.draw[[iter]][[2,1]] %*% t(fev1pp_training_fit_nonsparse_ls$swapped_V.draw[[iter]][[1,1]])
}))/length(thinned_iters_burnin)

# Calculating the posterior mean of the individual structure (Biocrates)
indiv_biocrates_ls_nonsparse_mean <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_ls$swapped_W.draw[[iter]][[1,1]] %*% t(fev1pp_training_fit_nonsparse_ls$swapped_Vs.draw[[iter]][[1,1]])
}))/length(thinned_iters_burnin)

# Calculating the posterior mean of the individual structure (Somascan)
indiv_somascan_ls_nonsparse_mean <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_ls$swapped_W.draw[[iter]][[2,2]] %*% t(fev1pp_training_fit_nonsparse_ls$swapped_Vs.draw[[iter]][[1,2]])
}))/length(thinned_iters_burnin)


# Reorganizing the results to fit the JIVE structure
heatmap_ls_nonsparse_jive <- list(data = list(hiv_copd_data[[1,1]], hiv_copd_data[[2,1]]),
                                  joint = list(joint_biocrates_ls_nonsparse_mean, joint_somascan_ls_nonsparse_mean),
                                  individual = list(indiv_biocrates_ls_nonsparse_mean, indiv_somascan_ls_nonsparse_mean),
                                  rankJ = joint_rank,
                                  rankA = c(indiv_rank1, indiv_rank2))

# Applying the showHeatmaps function (ordered by joint structure)
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Heatmap_Ordered_by_Joint.pdf"), width = 15)
showHeatmaps(heatmap_ls_nonsparse_jive, order_by = 0)
dev.off()

# Applying the showHeatmaps function (ordered by Biocrates structure)
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Heatmap_Ordered_by_Biocrates.pdf"), width = 15)
showHeatmaps(heatmap_ls_nonsparse_jive, order_by = 1)
dev.off()

# Applying the showHeatmaps function (ordered by Somascan structure)
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Heatmap_Ordered_by_Somascan.pdf"), width = 15)
showHeatmaps(heatmap_ls_nonsparse_jive, order_by = 2)
dev.off()

# -----------------------------------------------------------------------------
# Calculating meta-loadings for each feature
# -----------------------------------------------------------------------------

# Saving the joint loadings after taking the posterior mean and ordering by variance explained
U.mean.reorder <- fev1pp_training_fit_nonsparse_ls$U.mean.reorder

# Saving the individual loadings after taking the posterior mean and ordering by variance explained
W.mean.reorder <- fev1pp_training_fit_nonsparse_ls$W.mean.reorder

# Saving the betas vector ordered by variance explained
betas.mean.reorder <- fev1pp_training_fit_nonsparse_ls$betas.mean.reorder

# Iterating through the Biocrates metabolite factors and calculating the linear combination between
# each joint and individual loading with the corresponding beta

# Save the meta loadings vector
n_metab <- nrow(lavage_processed_no_info_log_scale)
n_protein <- nrow(somascan_normalized_clean_no_info_transpose_scale)
meta_loadings <- matrix(list(), nrow = q, ncol = 1)

# Save the indiv ranks
indiv_ranks <- fev1pp_training_fit_nonsparse$ranks[-1]

# Iterate through each metabolite, take linear combination with betas
for (s in 1:q) {
  if (s == 1) {
    for (i in 1:n_metab) {
      meta_loadings[[s,1]][i] <- 
        U.mean.reorder[[s,1]][i,,drop = FALSE] %*% betas.mean.reorder[2:(joint_rank+1),,drop = FALSE] + 
        W.mean.reorder[[s,s]][i,,drop = FALSE] %*% betas.mean.reorder[(joint_rank+2):(indiv_rank1+joint_rank+1),,drop = FALSE] 
    }
  }
  
  if (s == 2) {
    for (i in 1:n_protein) {
      meta_loadings[[s,1]][i] <- 
        U.mean.reorder[[s,1]][i,,drop = FALSE] %*% betas.mean.reorder[2:(joint_rank+1),,drop = FALSE] + 
        W.mean.reorder[[s,s]][i,,drop = FALSE] %*% betas.mean.reorder[(indiv_rank1+joint_rank+2):(indiv_rank1+indiv_rank2+joint_rank+1),,drop = FALSE] 
    }
  }
}

# -----------------------------------------------------------------------------
# Credible Intervals for Coefficients Using Non-Sparse Model
# -----------------------------------------------------------------------------

# Creating a matrix with the posterior samples for the betas after burnin and thinning
betas_ls_nonsparse_mat <- do.call(cbind, lapply(thinned_iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_ls$swapped_betas[[iter]][[1,1]]
}))

# Calculate the mean and CIs
betas_ls_nonsparse_mat_results <- t(apply(betas_ls_nonsparse_mat, 1, function(beta) {
  c(mean(beta), quantile(beta, 0.025), quantile(beta, 0.975))
}))

# Add factor labels
betas_ls_nonsparse_mat_results <- as.data.frame(betas_ls_nonsparse_mat_results)
betas_ls_nonsparse_mat_results <- cbind.data.frame(0:(nrow(betas_ls_nonsparse_mat_results)-1),
                                                   betas_ls_nonsparse_mat_results)
colnames(betas_ls_nonsparse_mat_results) <- c("Factor", "Posterior Mean", "Lower (2.5%)", "Upper (97.5%)")
betas_ls_nonsparse_mat_results <- betas_ls_nonsparse_mat_results[-1,] # Remove intercept

# Create plot
library(plotrix)
plotCI(x = betas_ls_nonsparse_mat_results$Factor,               # plotrix plot with confidence intervals
       y = betas_ls_nonsparse_mat_results$`Posterior Mean`,
       li = betas_ls_nonsparse_mat_results$`Lower (2.5%)`,
       ui = betas_ls_nonsparse_mat_results$`Upper (97.5%)`,
       main = "95% Credible Intervals for Factor Effects \n on FEV1pp (Non-Sparse)",
       xlab = "Factor", ylab = "Posterior Mean",
       scol = c(rep(1, joint_rank), rep(2, indiv_rank1), rep(3, indiv_rank2)))


# -----------------------------------------------------------------------------
# Credible Intervals for Coefficients Using Sparse Model
# -----------------------------------------------------------------------------

# Creating a matrix with the posterior samples for the betas after burnin and thinning
betas_ls_sparse_mat <- do.call(cbind, lapply(thinned_iters_burnin, function(iter) {
  fev1pp_training_fit_sparse_ls$swapped_betas[[iter]][[1,1]]
}))

# Calculate the mean and CIs
betas_ls_sparse_mat_results <- t(apply(betas_ls_sparse_mat, 1, function(beta) {
  c(mean(beta), quantile(beta, 0.025), quantile(beta, 0.975))
}))

# Add factor labels
betas_ls_sparse_mat_results <- as.data.frame(betas_ls_sparse_mat_results)
betas_ls_sparse_mat_results <- cbind.data.frame(0:(nrow(betas_ls_sparse_mat_results)-1),
                                                   betas_ls_sparse_mat_results)
colnames(betas_ls_sparse_mat_results) <- c("Factor", "Posterior Mean", "Lower (2.5%)", "Upper (97.5%)")
betas_ls_sparse_mat_results <- betas_ls_sparse_mat_results[-1,] # Remove intercept

# Create plot
library(plotrix)
plotCI(x = betas_ls_sparse_mat_results$Factor,               # plotrix plot with confidence intervals
       y = betas_ls_sparse_mat_results$`Posterior Mean`,
       li = betas_ls_sparse_mat_results$`Lower (2.5%)`,
       ui = betas_ls_sparse_mat_results$`Upper (97.5%)`,
       main = "95% Credible Intervals for Factor Effects \n on FEV1pp (Sparse)",
       xlab = "Factor", ylab = "Posterior Mean",
       scol = c(rep(1, joint_rank), rep(2, indiv_rank1), rep(3, indiv_rank2)))


# -----------------------------------------------------------------------------
# Cross Validated Model Fit
# -----------------------------------------------------------------------------

# Saving the index for each pair 
ind_of_pairs <- seq(1, n, by = 2)
n_pair <- length(ind_of_pairs)

# For running in parallel
funcs <- c("bpmf_data", "center_data", "bpmf", "get_results", "BIDIFAC",
           "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
           "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum")
packs <- c("MASS", "truncnorm", "EnvStats", "svMisc", "Matrix")

# Running each training and test run in parallel
cl <- makeCluster(3)
registerDoParallel(cl)
fev1pp_cv <- foreach(pair = ind_of_pairs, .packages = packs, .export = funcs, .verbose = TRUE) %dopar% {
  # Create a new vector of the outcome with the current pair set to NA
  fev1pp_cv <- fev1pp
  fev1pp_cv[[1,1]][pair:(pair+1),] <- NA
  
  # Run the model with the above pair's continuous outcome missing
  fev1pp_cv_fit_sparse <- bpmf(
    data = hiv_copd_data,
    Y = fev1pp_cv,
    nninit = TRUE,
    model_params = model_params,
    sparsity = TRUE,
    nsample = nsample,
    progress = TRUE
  )
  
  # Saving the predicted outcomes for the missing subjects
  Ym.draw <- fev1pp_cv_fit_sparse$Ym.draw
  ranks <- fev1pp_cv_fit_sparse$ranks
  
  # Calculating the log-joint density after burn-in
  convergence <- sapply(thinned_iters_burnin, function(sim_iter) {
    # Calculate the log-joint density at each thinned iterations
    log_joint_density(data = hiv_copd_data, 
                      U.iter = fev1pp_cv_fit_sparse$U.draw[[sim_iter]], 
                      V.iter = fev1pp_cv_fit_sparse$V.draw[[sim_iter]], 
                      W.iter = fev1pp_cv_fit_sparse$W.draw[[sim_iter]], 
                      Vs.iter = fev1pp_cv_fit_sparse$Vs.draw[[sim_iter]],
                      model_params = model_params,
                      ranks = fev1pp_cv_fit_sparse$ranks,
                      Y = fev1pp,
                      Ym.iter = fev1pp_cv_fit_sparse$Ym.draw[[sim_iter]],
                      beta.iter = fev1pp_cv_fit_sparse$beta.draw[[sim_iter]],
                      tau2.iter = fev1pp_cv_fit_sparse$tau2.draw[[sim_iter]],
                      gamma.iter = fev1pp_cv_fit_sparse$gamma.draw[[sim_iter]],
                      p.iter = fev1pp_cv_fit_sparse$p.draw[[sim_iter]])
  })

  # Save just the relevant output
  save(Ym.draw, ranks, convergence, file = paste0(results_wd, "FEV1pp_CV_Sparse_Pair", pair, ".rda"))
}
stopCluster(cl)

# Running the cross validation algorithm again, this time without sparsity
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
  
  # Saving the predicted outcomes for the missing subjects
  Ym.draw <- fev1pp_cv_fit_nonsparse$Ym.draw
  ranks <- fev1pp_cv_fit_nonsparse$ranks
  
  # Calculating the log-joint density after burn-in
  convergence <- sapply(thinned_iters_burnin, function(sim_iter) {
    # Calculate the log-joint density at each thinned iterations
    log_joint_density(data = hiv_copd_data, 
                      U.iter = fev1pp_cv_fit_nonsparse$U.draw[[sim_iter]], 
                      V.iter = fev1pp_cv_fit_nonsparse$V.draw[[sim_iter]], 
                      W.iter = fev1pp_cv_fit_nonsparse$W.draw[[sim_iter]], 
                      Vs.iter = fev1pp_cv_fit_nonsparse$Vs.draw[[sim_iter]],
                      model_params = model_params,
                      ranks = fev1pp_cv_fit_nonsparse$ranks,
                      Y = fev1pp,
                      Ym.iter = fev1pp_cv_fit_nonsparse$Ym.draw[[sim_iter]],
                      beta.iter = fev1pp_cv_fit_nonsparse$beta.draw[[sim_iter]],
                      tau2.iter = fev1pp_cv_fit_nonsparse$tau2.draw[[sim_iter]])
  })
  
  # Save just the relevant output
  save(Ym.draw, ranks, convergence, file = paste0(results_wd, "FEV1pp_CV_NonSparse_Pair", pair, ".rda"))
}
stopCluster(cl)

# Loading in the predicted outcome iteratively for each pair and computing
# the posterior mean predicted outcome. Comparing to the true FEV1pp 

fev1pp_cv_sparsity <- c()
convergence_cv_sparsity <- c()
for (pair in ind_of_pairs) {
  # Load in the results
  load(paste0(results_wd, "FEV1pp_CV_Sparse_Pair", pair, ".rda"))
  
  # Combine the samples
  samps <- do.call(cbind, do.call(cbind, Ym.draw))
  
  # Take a burn-in
  samps_burnin <- samps[,thinned_iters_burnin]
  
  # Save in the vector
  fev1pp_cv_sparsity[pair:(pair+1)] <- rowMeans(samps_burnin)
  
  # Save the vector of log-joint densities after burn-in
  convergence_cv_sparsity <- c(convergence_cv_sparsity, mean(convergence))
}

# Plotting the results
plot(fev1pp_cv_sparsity, c(fev1pp[[1,1]]), xlab = "Predicted FEV1pp", ylab = "Observed FEV1pp", main = "Cross Validated FEV1pp from Sparse Model")
abline(a=0, b=1)
cor.test(fev1pp_cv_sparsity, c(fev1pp[[1,1]]))

# Loading in the predicted outcome iteratively for each pair and computing
# the posterior mean predicted outcome. Comparing to the true FEV1pp 
# Using the non-sparse model

fev1pp_cv <- c()
convergence_cv_nonsparse <- c()
for (pair in ind_of_pairs) {
  # Load in the results
  load(paste0(results_wd, "FEV1pp_CV_NonSparse_Pair", pair, ".rda"))
  
  # Combine the samples
  samps <- do.call(cbind, do.call(cbind, Ym.draw))
  
  # Take a burn-in
  samps_burnin <- samps[,thinned_iters_burnin]
  
  # Save in the vector
  fev1pp_cv[pair:(pair+1)] <- rowMeans(samps_burnin)
  
  # Save the vector of log-joint densities after burn-in
  convergence_cv_nonsparse <- c(convergence_cv_nonsparse, mean(convergence))
}

# Plotting the results
plot(fev1pp_cv, c(fev1pp[[1,1]]), xlab = "Predicted FEV1pp", ylab = "Observed FEV1pp", main = "Cross Validated FEV1pp from Model Without Sparsity")
abline(a=0, b=1)
cor.test(fev1pp_cv, c(fev1pp[[1,1]]))



