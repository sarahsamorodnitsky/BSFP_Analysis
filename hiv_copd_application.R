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
lavage_processed_no_info_log_scale <- scale(log(lavage_processed_no_info + 1), # log/scale
                                            center = TRUE, scale = TRUE)
rownames(lavage_processed_no_info_log_scale) <- bioc_meta_data$Metabolite

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
                                                 nsample = nsample)

fev1pp_training_fit_nonsparse_ls <- label_switching(U.draw = fev1pp_training_fit_nonsparse$U.draw,
                                                   V.draw = fev1pp_training_fit_nonsparse$V.draw,
                                                   W.draw = fev1pp_training_fit_nonsparse$W.draw,
                                                   Vs.draw = fev1pp_training_fit_nonsparse$Vs.draw,
                                                   betas = fev1pp_training_fit_nonsparse$beta.draw,
                                                   gammas = NULL,
                                                   r = fev1pp_training_fit_nonsparse$ranks[1],
                                                   r.vec = fev1pp_training_fit_nonsparse$ranks[-1],
                                                   nsample = nsample)

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

# In sparse model, which factors were selected?
mat_of_gammas_ls <- do.call(cbind, lapply(fev1pp_training_fit_sparse_ls$swapped_gammas, function(iter) iter[[1,1]]))
mat_of_gammas <- do.call(cbind, lapply(fev1pp_training_fit_sparse$gamma.draw, function(iter) iter[[1,1]]))

post_prob_ls <- rowMeans(mat_of_gammas_ls)
post_prob <- rowMeans(mat_of_gammas)

# Comparing the training data fits
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

# -----------------------------------------------------------------------------
# PCA-like plots using non-sparse model
# -----------------------------------------------------------------------------

# Calculate posterior mean for joint structure -- 
V_ls_nonsparse_mean <- Reduce("+", lapply(thinned_iters_burnin, function(iter) fev1pp_training_fit_nonsparse_ls$swapped_V.draw[[iter]][[1,1]]))/length(thinned_iters_burnin)

# Load in the clustering information
load("~/HIV-COPD/Data/UnivariatePValuesUncombinedClusterComparison.rda", verbose = TRUE)

# Plot subjects against every combination of PCs
joint_rank <- fev1pp_training_fit_nonsparse$ranks[1]
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Joint_Structure_PCA_Plot.pdf"), width = 15)
for (rs1 in 1:joint_rank) {
 for (rs2 in 1:joint_rank) {
   if (rs1 != rs2) {
     par(mfrow = c(1,2))
     # Plot with case-control label as colors
     plot(V_ls_nonsparse_mean[,rs1], V_ls_nonsparse_mean[,rs2], pch = 16, xlab = paste0("Joint Score ", rs1),
          ylab = paste0("Joint Score ", rs2), main = paste0("Joint Factors ", rs1, " vs ", rs2),
          col = clinical_data_soma$ccstat)
     # Plot with stand-out subjects as colors
     plot(V_ls_nonsparse_mean[,rs1], V_ls_nonsparse_mean[,rs2], pch = 16, xlab = paste0("Joint Score ", rs1),
          ylab = paste0("Joint Score ", rs2), main = paste0("Joint Factors ", rs1, " vs ", rs2),
          col = clusters.soma)
     par(mfrow = c(1,1))
   }
 } 
}
dev.off()

# Calculate the posterior mean for the individual structure --
V1_ls_nonsparse_mean <- Reduce("+", lapply(thinned_iters_burnin, function(iter) fev1pp_training_fit_nonsparse_ls$swapped_Vs.draw[[iter]][[1,1]]))/length(thinned_iters_burnin)
V2_ls_nonsparse_mean <- Reduce("+", lapply(thinned_iters_burnin, function(iter) fev1pp_training_fit_nonsparse_ls$swapped_Vs.draw[[iter]][[1,2]]))/length(thinned_iters_burnin)

# Plot subjects against every combination of PCs for Biocrates
indiv_rank1 <- fev1pp_training_fit_nonsparse$ranks[2]
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Individual_Structure_Biocrates_PCA_Plot.pdf"), width = 15)
for (rs1 in 1:indiv_rank1) {
  for (rs2 in 1:indiv_rank1) {
    if (rs1 != rs2) {
      par(mfrow = c(1,2))
      # Plot with case-control label as colors
      plot(V1_ls_nonsparse_mean[,rs1], V1_ls_nonsparse_mean[,rs2], pch = 16, xlab = paste0("Individual Score ", rs1),
           ylab = paste0("Individual Score ", rs2), main = paste0("Individual Factors Biocrates ", rs1, " vs ", rs2),
           col = clinical_data_soma$ccstat)
      # Plot with stand-out subjects as colors
      plot(V1_ls_nonsparse_mean[,rs1], V1_ls_nonsparse_mean[,rs2], pch = 16, xlab = paste0("Individual Score ", rs1),
           ylab = paste0("Individual Score ", rs2), main = paste0("Individual Factors Biocrates ", rs1, " vs ", rs2),
           col = clusters.soma)
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
      plot(V2_ls_nonsparse_mean[,rs1], V2_ls_nonsparse_mean[,rs2], pch = 16, xlab = paste0("Individual Score ", rs1),
           ylab = paste0("Individual Score ", rs2), main = paste0("Individual Factors Somascan ", rs1, " vs ", rs2),
           col = clinical_data_soma$ccstat)
      # Plot with stand-out subjects as colors
      plot(V2_ls_nonsparse_mean[,rs1], V2_ls_nonsparse_mean[,rs2], pch = 16, xlab = paste0("Individual Score ", rs1),
           ylab = paste0("Individual Score ", rs2), main = paste0("Individual Factors Somascan ", rs1, " vs ", rs2),
           col = clusters.soma)
      par(mfrow = c(1,1))
    }
  } 
}
dev.off()

# -----------------------------------------------------------------------------
# Investigating factor 2
# -----------------------------------------------------------------------------

# Selecting the loadings for factor 2 from Biocrates --
U_Biocrates_Factor2_nonsparse_ls <- do.call(cbind, lapply(thinned_iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_ls$swapped_U.draw[[iter]][[1,1]][,2]
}))
rownames(U_Biocrates_Factor2_nonsparse_ls) <- rownames(hiv_copd_data[[1,1]])

# Calculate the posterior mean loading for each protein
U_Biocrates_Factor2_nonsparse_ls_mean <- rowMeans(U_Biocrates_Factor2_nonsparse_ls)
names(U_Biocrates_Factor2_nonsparse_ls_mean) <- rownames(hiv_copd_data[[1,1]])
U_Biocrates_Factor2_nonsparse_ls_mean_order <- U_Biocrates_Factor2_nonsparse_ls_mean[order(U_Biocrates_Factor2_nonsparse_ls_mean, decreasing = TRUE)]

# Waterfall plot
barplot(U_Biocrates_Factor2_nonsparse_ls_mean_order, col="blue", border="blue", space=0.5,
        main = "Posterior Mean of Biocrates Factor 2 Loadings", ylab = "Mean Loading",
        cex.axis=1.2, cex.lab=1.4, names.arg = NULL)



# Constructing credible intervals
U_Biocrates_Factor2_nonsparse_ls_mat <- t(apply(U_Biocrates_Factor2_nonsparse_ls, 1, function(protein) {
  c(mean(protein), quantile(protein, 0.025), quantile(protein, 0.975))
}))
U_Biocrates_Factor2_nonsparse_ls_mat <- cbind.data.frame(rownames(U_Biocrates_Factor2_nonsparse_ls_mat), U_Biocrates_Factor2_nonsparse_ls_mat)
colnames(U_Biocrates_Factor2_nonsparse_ls_mat) <- c("Protein", "Mean", "Lower", "Upper")
U_Biocrates_Factor2_nonsparse_ls_mat$Protein <- 1:nrow(U_Biocrates_Factor2_nonsparse_ls_mat)

# Order by the mean
U_Biocrates_Factor2_nonsparse_ls_mat <- U_Biocrates_Factor2_nonsparse_ls_mat[order(U_Biocrates_Factor2_nonsparse_ls_mat$Mean, decreasing = TRUE),]

# Creating the same plot with error bars
error.bar <- function(x, y, upper, lower, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

pdf("~/BayesianPMF/04DataApplication/Figures/Biocrates_Factor2_Loadings.pdf")
waterfall_factor2 <- barplot(U_Biocrates_Factor2_nonsparse_ls_mat$Mean, col="blue", border="blue", space=0.5,
                             main = "Posterior Mean of Biocrates Factor 2 Loadings", ylab = "Mean Loading",
                             cex.axis=1.2, cex.lab=1.4, names.arg = NULL, ylim = c(-1,1))
error.bar(waterfall_factor2, U_Biocrates_Factor2_nonsparse_ls_mat$Mean, U_Biocrates_Factor2_nonsparse_ls_mat$Upper-U_Biocrates_Factor2_nonsparse_ls_mat$Mean, U_Biocrates_Factor2_nonsparse_ls_mat$Mean-U_Biocrates_Factor2_nonsparse_ls_mat$Lower)
dev.off()


# Selecting the loadings for factor 2 from Somascan --
U_Somascan_Factor2_nonsparse_ls <- do.call(cbind, lapply(thinned_iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_ls$swapped_U.draw[[iter]][[2,1]][,2]
}))

# Calculate the posterior mean loading for each protein
U_Somascan_Factor2_nonsparse_ls_mean <- rowMeans(U_Somascan_Factor2_nonsparse_ls)
U_Somascan_Factor2_nonsparse_ls_mean_order <- U_Somascan_Factor2_nonsparse_ls_mean[order(U_Somascan_Factor2_nonsparse_ls_mean, decreasing = TRUE)]

# Waterfall plot
barplot(U_Somascan_Factor2_nonsparse_ls_mean_order, col="blue", border="blue",
        main = "Posterior Mean of Somascan Factor 2 Loadings", ylab = "Mean Loading",
        cex.axis=1.2, cex.lab=1.4)

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
cl <- makeCluster(5)
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

  # Save just the relevant output
  save(Ym.draw, ranks, file = paste0(results_wd, "FEV1pp_CV_Sparse_Pair", pair, ".rda"))
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
  
  # Save just the relevant output
  save(Ym.draw, ranks, file = paste0(results_wd, "FEV1pp_CV_NonSparse_Pair", pair, ".rda"))
}
stopCluster(cl)

# Loading in the predicted outcome iteratively for each pair and computing
# the posterior mean predicted outcome. Comparing to the true FEV1pp 

fev1pp_cv_sparsity <- c()
for (pair in ind_of_pairs) {
  # Load in the results
  load(paste0(results_wd, "FEV1pp_CV_Sparse_Pair", pair, ".rda"))
  
  # Combine the samples
  samps <- do.call(cbind, do.call(cbind, Ym.draw))
  
  # Take a burn-in
  samps_burnin <- samps[,thinned_iters_burnin]
  
  # Save in the vector
  fev1pp_cv_sparsity[pair:(pair+1)] <- rowMeans(samps_burnin)
}

# Plotting the results
plot(fev1pp_cv_sparsity, c(fev1pp[[1,1]]), xlab = "Predicted FEV1pp", ylab = "Observed FEV1pp", main = "Cross Validated FEV1pp from Sparse Model")
abline(a=0, b=1)
cor.test(fev1pp_cv_sparsity, c(fev1pp[[1,1]]))

# Loading in the predicted outcome iteratively for each pair and computing
# the posterior mean predicted outcome. Comparing to the true FEV1pp 
# Using the non-sparse model

fev1pp_cv <- c()
for (pair in ind_of_pairs) {
  # Load in the results
  load(paste0(results_wd, "FEV1pp_CV_NonSparse_Pair", pair, ".rda"))
  
  # Combine the samples
  samps <- do.call(cbind, do.call(cbind, Ym.draw))
  
  # Take a burn-in
  samps_burnin <- samps[,thinned_iters_burnin]
  
  # Save in the vector
  fev1pp_cv[pair:(pair+1)] <- rowMeans(samps_burnin)
}

# Plotting the results
plot(fev1pp_cv, c(fev1pp[[1,1]]), xlab = "Predicted FEV1pp", ylab = "Observed FEV1pp", main = "Cross Validated FEV1pp from Model Without Sparsity")
abline(a=0, b=1)
cor.test(fev1pp_cv, c(fev1pp[[1,1]]))



