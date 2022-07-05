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

# -----------------------------------------------------------------------------
# Investigating training data fit results
# -----------------------------------------------------------------------------

# -------------------------------------
# Convergence
# -------------------------------------

# Assessing convergence for both model fits
load(paste0(results_wd, "training_data_fit.rda"), verbose = TRUE)

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

# -------------------------------------
# Applying our label switching algorithm
# to the results
# -------------------------------------

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
                                                    thinned_iters_burnin = thinned_iters_burnin,
                                                    nninit = TRUE)

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
# Calculation of variance explained in data and response
# -------------------------------------

# Save the parameters each ordered by variance explained
ranks <- fev1pp_training_fit_nonsparse$ranks
betas_ls <- fev1pp_training_fit_nonsparse_ls$betas.mean.reorder
V_ls <- fev1pp_training_fit_nonsparse_ls$V.mean.reorder
U_ls <- fev1pp_training_fit_nonsparse_ls$U.mean.reorder
Vs_ls <- fev1pp_training_fit_nonsparse_ls$Vs.mean.reorder
W_ls <- fev1pp_training_fit_nonsparse_ls$W.mean.reorder

# -----------------
# Variance exp in X
# -----------------

# Calculate the Frobenius norm of each data source
ss_biocrates <- frob(hiv_copd_data[[1,1]])
ss_somascan <- frob(hiv_copd_data[[2,1]])

# Calculate the contribution of joint structure for Biocrates
ss_biocrates_joint <- frob(U_ls[[1,1]] %*% t(V_ls[[1,1]]))

# Calculate the contribution of joint structure for Somascan
ss_somascan_joint <- frob(U_ls[[2,1]] %*% t(V_ls[[1,1]]))

# Calculate the contribution of individual structure for Biocrates
ss_biocrates_indiv <- frob(W_ls[[1,1]] %*% t(Vs_ls[[1,1]]))

# Calculate the contribution of individual structure for Somascan
ss_somascan_indiv <- frob(W_ls[[2,2]] %*% t(Vs_ls[[1,2]]))

# Calculate the variance explained
prop_var_exp_biocrates_joint <- ss_biocrates_joint/ss_biocrates # Joint Biocrates
prop_var_exp_somascan_joint <- ss_somascan_joint/ss_somascan # Joint Somascan
prop_var_exp_biocrates_indiv <- ss_biocrates_indiv/ss_biocrates
prop_var_exp_somascan_indiv <- ss_somascan_indiv/ss_somascan

# -----------------
# Variance exp in Y
# -----------------

# Calculate the sum of squared Ys
ssY <- sum(fev1pp[[1,1]])

# Calculate the contribution from joint prediction
ss_joint_Y <- sum((V_ls[[1,1]] %*% betas_ls[2:(1+ranks[1]),,drop = FALSE])^2)
ss_biocrates_Y <- sum((Vs_ls[[1,1]] %*% betas_ls[(2+ranks[1]):(1+ranks[1]+ranks[2]),,drop = FALSE])^2)
ss_somascan_Y <- sum((Vs_ls[[1,2]] %*% betas_ls[(2+ranks[1]+ranks[2]):(1+sum(ranks)),,drop = FALSE])^2)

# Calculate proportion variance explained
var_exp_Y_joint <- ss_joint_Y/ssY
var_exp_Y_biocrates <- ss_biocrates_Y/ssY
var_exp_Y_somascan <- ss_somascan_Y/ssY

# Creating table
library(xtable)
var_exp_df <- data.frame(Source = c("Metabolite", "Protein", "FEV1pp"),
                         Joint = numeric(3),
                         Indiv.Metab = numeric(3),
                         Indiv.Protein = numeric(3))
var_exp_df[1,2:4] <- c(prop_var_exp_biocrates_joint, prop_var_exp_biocrates_indiv, 0) # 0 because Somascan does not explain var in Biocrates
var_exp_df[2,2:4] <- c(prop_var_exp_somascan_joint, 0, prop_var_exp_somascan_indiv) # 0 because Biocrates does not explain var in Somascan
var_exp_df[3,2:4] <- c(var_exp_Y_joint, var_exp_Y_biocrates, var_exp_Y_somascan)

# Print table
print(xtable(var_exp_df, digits = 4), include.rownames = FALSE)

# -----------------------------------------------------------------------------
# PCA-like plots using non-sparse model
# -----------------------------------------------------------------------------

library(viridis)

# Load in the clustering information
load("~/HIV-COPD/Data/UnivariatePValuesUncombinedClusterComparison.rda", verbose = TRUE)

# Creating colors for the plots
ccstat_colors <- viridis(2)[clinical_data_soma$ccstat]
cluster_colors <- viridis(3)[1:2][clusters.soma]

# Creating shapes
cluster_shapes <- c(3,4)[clusters.soma]

# -------------------------------------
# Joint structure
# -------------------------------------

# Saving the posterior mean of the joint structure, ordered by variance explained 
V_ls_nonsparse_mean_ordered <- fev1pp_training_fit_nonsparse_ls$V.mean.reorder[[1,1]]

# Plot subjects against every combination of PCs
joint_rank <- fev1pp_training_fit_nonsparse$ranks[1]
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Joint_Structure_PCA_Plot.pdf"), width = 15)
for (rs1 in 1:joint_rank) {
 for (rs2 in 1:joint_rank) {
   if (rs1 != rs2) {
     par(mfrow = c(1,2), mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
     # Plot with case-control label as colors
     plot(V_ls_nonsparse_mean_ordered[,rs1], V_ls_nonsparse_mean_ordered[,rs2], pch = 16, xlab = paste0("Joint Factor ", rs1),
          ylab = paste0("Joint Factor ", rs2), main = paste0("Joint Factors ", rs1, " vs ", rs2),
          col = ccstat_colors)
     
     # Plot with stand-out subjects as colors
     plot(V_ls_nonsparse_mean_ordered[,rs1], V_ls_nonsparse_mean_ordered[,rs2], xlab = paste0("Joint Factor ", rs1),
          ylab = paste0("Joint Factor ", rs2), main = paste0("Joint Factors ", rs1, " vs ", rs2),
          col = cluster_colors, pch = cluster_shapes)
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
Vs_ls_nonsparse_mean_ordered <- fev1pp_training_fit_nonsparse_ls$Vs.mean.reorder

# Plot subjects against every combination of PCs for Biocrates
indiv_rank1 <- fev1pp_training_fit_nonsparse$ranks[2]
pdf(paste0("~/BayesianPMF/04DataApplication/Figures/Individual_Structure_Biocrates_PCA_Plot.pdf"), width = 15)
for (rs1 in 1:indiv_rank1) {
  for (rs2 in 1:indiv_rank1) {
    if (rs1 != rs2) {
      par(mfrow = c(1,2))
      # Plot with case-control label as colors
      plot(Vs_ls_nonsparse_mean_ordered[[1,1]][,rs1], Vs_ls_nonsparse_mean_ordered[[1,1]][,rs2], pch = 16, xlab = paste0("Individual Factor ", rs1),
           ylab = paste0("Individual Factor ", rs2), main = paste0("Individual Factors Metabolomic ", rs1, " vs ", rs2),
           col = ccstat_colors)
      
      # Plot with stand-out subjects as colors
      plot(Vs_ls_nonsparse_mean_ordered[[1,1]][,rs1], Vs_ls_nonsparse_mean_ordered[[1,1]][,rs2], xlab = paste0("Individual Factor ", rs1),
           ylab = paste0("Individual Factor ", rs2), main = paste0("Individual Factors Metabolomic ", rs1, " vs ", rs2),
           col = cluster_colors, pch = cluster_shapes)
      legend("bottomleft", col = c("#440154FF", "#21908CFF"), pch = c(3,4),
             legend = c("Cluster 1", "Cluster 2"))
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
      plot(Vs_ls_nonsparse_mean_ordered[[1,2]][,rs1], Vs_ls_nonsparse_mean_ordered[[1,2]][,rs2], pch = 16, xlab = paste0("Individual Factor ", rs1),
           ylab = paste0("Individual Factor ", rs2), main = paste0("Individual Factors Proteomic ", rs1, " vs ", rs2),
           col = ccstat_colors)
      
      # Plot with stand-out subjects as colors
      plot(Vs_ls_nonsparse_mean_ordered[[1,2]][,rs1], Vs_ls_nonsparse_mean_ordered[[1,2]][,rs2], xlab = paste0("Individual Factor ", rs1),
           ylab = paste0("Individual Factor ", rs2), main = paste0("Individual Factors Proteomic ", rs1, " vs ", rs2),
           col = cluster_colors, pch = cluster_shapes)
      par(mfrow = c(1,1))
    }
  } 
}
dev.off()

# -----------------------------------------------------------------------------
# Kernel density plots for particular factors
# -----------------------------------------------------------------------------

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
      load(paste0(results_wd, "BPMF/FEV1pp_CV_NonSparse_Pair", pair, ".rda"), verbose = TRUE)
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
