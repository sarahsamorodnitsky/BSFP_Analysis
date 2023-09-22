# -----------------------------------------------------------------------------
# Applying the BSFP model to the BRCA data from TCGA. 
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set-up
# -----------------------------------------------------------------------------

# Packages
library(doParallel)
library(foreach)
library(r.jive)

# Working directories
data_wd <- "~/BayesianPMF/data/"
model_wd <- "~/BSFP_Analysis/"
results_wd <- "~/BayesianPMF/04DataApplication/TCGA_BRCA"

# Loading in the model functions
source(paste0(model_wd, "bsfp.R"))

# -----------------------------------------------------------------------------
# Preparing the data
# -----------------------------------------------------------------------------

# Load in the omics data for TCGA BRCA from r.jive
data(BRCA_data)

# Load in the clinical data
clin.data <- read.csv("~/BayesianPMF/04DataApplication/data/BRCA_clin_data.csv")
colnames(clin.data) <- clin.data[1,]
clin.data <- clin.data[-1,]

# Loads in an object called `Data`. Rename to BRCA. 
# Includes Expression, Methylation, miRNA
# 348 samples

# Save the sample IDs from the r.jive package
sample.ids <- substr(colnames(Data$Expression), 1, 12)
sample.ids <- gsub('\\.', '-', sample.ids) # Convert . to -

# Match omics data to clinical data
Match <- match(sample.ids, clin.data$`Complete TCGA ID`)
clin.data.matched <- clin.data[Match,]

# CHECK
all(clin.data.matched$`Complete TCGA ID` == sample.ids) # TRUE!

# Use estrogen receptor status as the outcome
# clin.data.matched$ER.Status #use this as outcome?

# Data processing from BCC:
# Impute missing values in gene expression using KNN.
# Remove miRNAs with over 50% missing values
# Select only most variable genes in gene expression 
# Take square root of methylation data
# Log(1+x) transformation of miRNA data
# Center the features but do not scale. 

# Center the data, do not scale the data
Data$Expression <- t(scale(t(Data$Expression), center = TRUE, scale = FALSE))
Data$Methylation <- t(scale(t(Data$Methylation), center = TRUE, scale = FALSE))
Data$miRNA <- t(scale(t(Data$miRNA), center = TRUE, scale = FALSE))

# Combining the sources together
tcga_brca_data <- matrix(list(), nrow = 3, ncol = 1)
tcga_brca_data[[1,1]] <- Data$Expression
tcga_brca_data[[2,1]] <- Data$Methylation
tcga_brca_data[[3,1]] <- Data$miRNA

# Preparing the response variable
er_status <- matrix(list(), nrow = 1, ncol = 1)
er_status[[1,1]] <- matrix(clin.data.matched$`ER Status`, ncol = 1)
er_status[[1,1]][!(er_status[[1,1]] %in% c("Positive", "Negative"))] <- NA
er_status[[1,1]][er_status[[1,1]] == "Positive"] <- 1
er_status[[1,1]][er_status[[1,1]] == "Negative"] <- 0
er_status[[1,1]] <- matrix(as.numeric(er_status[[1,1]]), ncol = 1)
rownames(er_status[[1,1]]) <- clin.data.matched$`Complete TCGA ID`

# -----------------------------------------------------------------------------
# Model parameters
# -----------------------------------------------------------------------------

# Data dimensions
p1 <- nrow(Data$Expression)
p2 <- nrow(Data$Methylation)
p3 <- nrow(Data$miRNA)
p <- p1 + p2 + p3
p.vec <- c(p1, p2, p3)
n <- ncol(Data$Expression)
q <- nrow(tcga_brca_data)

# Error variances
sigma21 <- 1 # error variance for Expression
sigma22 <- 1 # error variance for Methylation
sigma23 <- 1 # error variance for miRNA

# Variance of joint structure
lambda_joint <- sqrt(p1 + p2 + p3) + sqrt(n)
sigma2_joint <- 1/(lambda_joint)

# Variance of individual structure for expression
lambda_indiv1 <- sqrt(p1) + sqrt(n)
sigma2_indiv1 <- 1/(lambda_indiv1) 

# Variance of individual structure for methylation
lambda_indiv2 <- sqrt(p2) + sqrt(n)
sigma2_indiv2 <- 1/(lambda_indiv2)

# Variance of individual structure for methylation
lambda_indiv3 <- sqrt(p3) + sqrt(n)
sigma2_indiv3 <- 1/(lambda_indiv3)

# For the regression coefficients, beta
lambda2_intercept <- 4
lambda2_joint <- 1
lambda2_indiv1 <- 1
lambda2_indiv2 <- 1
lambda2_indiv3 <- 1
beta_vars <- c(lambda2_intercept, lambda2_joint, lambda2_indiv1, lambda2_indiv2, lambda2_indiv3)

# For the response vector
shape <- 0.01
rate <- 0.01

# Putting the model parameters together
model_params <- list(error_vars = c(sigma21, sigma22, sigma23),
                     joint_var = sigma2_joint,
                     indiv_vars = c(sigma2_indiv1, sigma2_indiv2, sigma2_indiv3),
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
# Training Data Model Fit 
# -----------------------------------------------------------------------------

# Fitting BSFP
tcga_brca_er_training <- BSFP::bsfp(
  data = tcga_brca_data,
  Y = er_status,
  nninit = TRUE,
  previous_init = NULL,
  model_params = model_params,
  save_structures = TRUE,
  save_loadings_scores = TRUE,
  save_predictive_model = TRUE,
  save_imputations = TRUE,
  nsample = nsample,
  thinning = 10,
  progress = TRUE,
  save_last_sample = TRUE
)

# Save
save(tcga_brca_er_training,
     file = "~/BayesianPMF/04DataApplication/TCGA_BRCA/TCGA_BRCA_Training_V2_thinned.rda")

# -----------------------------------------------------------------------------
# Assessing convergence 
# -----------------------------------------------------------------------------

# Use trace plots
plot(sapply(seq(length(tcga_brca_er_training$J.draw)), function(i) 
  tcga_brca_er_training$J.draw[[i]][[1,1]][1,1]))
plot(sapply(seq(length(tcga_brca_er_training$J.draw)), function(i) 
  tcga_brca_er_training$J.draw[[i]][[2,1]][50,101]))
plot(sapply(seq(length(tcga_brca_er_training$J.draw)), function(i) 
  tcga_brca_er_training$J.draw[[i]][[3,1]][10,101]))

plot(sapply(seq(length(tcga_brca_er_training$A.draw)), function(i) 
  tcga_brca_er_training$A.draw[[i]][[1,1]][1,1]))
plot(sapply(seq(length(tcga_brca_er_training$A.draw)), function(i) 
  tcga_brca_er_training$A.draw[[i]][[2,1]][50,10]))
plot(sapply(seq(length(tcga_brca_er_training$A.draw)), function(i) 
  tcga_brca_er_training$A.draw[[i]][[3,1]][10,20]))

plot(sapply(seq(length(tcga_brca_er_training$Z.draw))[-1], function(i) 
  tcga_brca_er_training$Z.draw[[i]][1,1]))
plot(sapply(seq(length(tcga_brca_er_training$Z.draw))[-1], function(i) 
  tcga_brca_er_training$Z.draw[[i]][2,1]))
plot(sapply(seq(length(tcga_brca_er_training$Z.draw))[-1], function(i) 
  tcga_brca_er_training$Z.draw[[i]][10,1]))

# Check prediction accuracy
EY.draw <- do.call(cbind, lapply(tcga_brca_er_training$EY.draw, function(i) i[[1,1]]))

# Compare the number of 0s and 1s in truth vs. predicted (about the same)
table(er_status)
table(ifelse(rowMeans(EY.draw)<0.5, 0, 1))

# Plot as comparison
plot(rowMeans(EY.draw), er_status[[1,1]])

# Using log-joint density
fev1pp_training_nonsparse_conv_V2 <- sapply(1:nsample, function(sim_iter) {
  # Calculate the log-joint density at each thinned iterations
  log_joint_density(data = tcga_brca_data, 
                    U.iter = tcga_brca_er_training_thinned$U.draw[[sim_iter]], 
                    V.iter = tcga_brca_er_training_thinned$V.draw[[sim_iter]], 
                    W.iter = tcga_brca_er_training_thinned$W.draw[[sim_iter]], 
                    Vs.iter = tcga_brca_er_training_thinned$Vs.draw[[sim_iter]],
                    model_params = model_params,
                    ranks = tcga_brca_er_training_thinned$ranks,
                    Y = er_status,
                    beta.iter = tcga_brca_er_training_thinned$beta.draw[[sim_iter]])
})

# -----------------------------------------------------------------------------
# Alignment
# -----------------------------------------------------------------------------

# Save the ranks
ranks <- tcga_brca_er_training$ranks
joint.rank <- tcga_brca_er_training$ranks[1]
indiv.ranks <- tcga_brca_er_training$ranks[-1]

# Applying the Match Align algorithm to undo rotational invariance in the results
er_status_training_fit_aligned <- match_align_bpmf(tcga_brca_er_training, 
                                                   y = er_status,
                                                   model_params = model_params,
                                                   p.vec = p.vec, 
                                                   iters_burnin = seq(length(tcga_brca_er_training$J.draw)))

# Save the results
save(er_status_training_fit_aligned, 
     file = "~/BayesianPMF/04DataApplication/TCGA_BRCA/TCGA_BRCA_Training_Thinned_Aligned_V2.rda")

# Exact the scores
joint.scores.final <- er_status_training_fit_aligned$joint.scores.final
joint.loadings.final <- er_status_training_fit_aligned$joint.loadings.final
joint.betas.final <- er_status_training_fit_aligned$joint.betas.final

individual.scores.final <- er_status_training_fit_aligned$individual.scores.final
individual.loadings.final <- er_status_training_fit_aligned$individual.loadings.final
individual.betas.final <- er_status_training_fit_aligned$individual.betas.final

# -----------------------------------------------------------------------------
# Heatmaps using sup.r.jive
# -----------------------------------------------------------------------------

# Load in the plotting library
source("~/BayesianPMF/04DataApplication/BPMF/Figures/plotHeatmap.R")

# Save the posterior mean of the joint scores
S_J <- t(Reduce("+", joint.scores.final)/burnin)

# Posterior mean of individual scores
S_I <- lapply(1:q, function(s) t(Reduce("+", individual.scores.final[[s]])/burnin))

# Posterior mean of joint loadings
inds <- cumsum(p.vec)
U_I_full <- Reduce("+", joint.loadings.final)/ burnin
U_I <- list(U_I_full[1:inds[1],,drop=FALSE],
            U_I_full[(inds[1]+1):inds[2],,drop=FALSE],
            U_I_full[(inds[2]+1):inds[3],,drop=FALSE])

# Posterior mean of the individual loadings
W_I <- lapply(1:q, function(s) Reduce("+", individual.loadings.final[[s]])/burnin)

# Posterior mean of joint betas
theta1 <- t(Reduce("+", joint.betas.final)/burnin)

# Posterior mean of individual betas
theta2 <- lapply(1:q, function(s) t(Reduce("+", individual.betas.final[[s]]))/burnin)

# Combine into list
tcga_brca_data_list_er_status <- list(X = list(tcga_brca_data[[1]], tcga_brca_data[[2]], tcga_brca_data[[3]]), 
                                      Y = unlist(er_status[[1,1]]))

heatmap_sup.r.jive <- list(data = tcga_brca_data_list_er_status,
                           S_J = S_J, S_I = S_I,
                           U_I = U_I, W_I = W_I,
                           theta1 = theta1, theta2 = theta2)

plotHeatmap(heatmap_sup.r.jive)

# -----------------------------------------------------------------------------
# Plot the credible intervals for the factors
# -----------------------------------------------------------------------------

library(dplyr)
library(BSFP)

# Save the working directory where the summaries will go
exploring_factors_wd <- "~/BayesianPMF/04DataApplication/TCGA_BRCA/"

# Calculate the summaries
summaries <- summarize_factors(data = tcga_brca_data, Y = er_status,
                               iters_burnin = 1:length(tcga_brca_er_training$J.draw), 
                               aligned_results = er_status_training_fit_aligned,
                               ranks = tcga_brca_er_training$ranks,
                               Ym.draw = tcga_brca_er_training$Ym.draw,
                               tau2.draw = NULL)

# Save results
save(summaries, er_status,
     file = "~/BayesianPMF/04DataApplication/TCGA_BRCA/Exploring_Factors/TCGA_BRCA_Factor_Summaries_Burnin_Aligned_V2.rda")

# -----------------------------------------------------------------------------
# Calculating variance expained 
# -----------------------------------------------------------------------------

# Calculating a summary of the variance explained
var_exp_summary <- var_explained(BPMF.fit = tcga_brca_er_training,
                                 iters_burnin = 1:length(tcga_brca_er_training$J.draw))

# Joint structure 
var_exp_summary$Joint

# Individual structure
var_exp_summary$Individual

# -----------------------------------------------------------------------------
# Contributions of each factor (joint, individual) to predicting Y
# -----------------------------------------------------------------------------

# Calculating the contribution of the joint factors to predicting FEV1pp
joint_contribution_to_er <- sapply(1:length(tcga_brca_er_training$J.draw), function(iter) {
  var(joint.scores.final[[iter]] %*% joint.betas.final[[iter]])/var(er_status[[1,1]], na.rm = TRUE)
})

# Calculating the contribution of the individual factors to predicting FEV1pp
individual_contribution_to_er <- lapply(1:q, function(s) {
  sapply(1:length(tcga_brca_er_training$J.draw), function(iter) {
    var(individual.scores.final[[s]][[iter]] %*% individual.betas.final[[s]][[iter]])/var(er_status[[1,1]], na.rm = TRUE)
  })
})

# Overall of Joint 
mean(joint_contribution_to_er)
c(quantile(joint_contribution_to_er, 0.025), quantile(joint_contribution_to_er, 0.975))

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


# -----------------------------------------------------------------------------
# 5-fold Cross Validation with BSFP
# -----------------------------------------------------------------------------

library(caret)

# Which outcomes have non-missing values?
outcomes_to_keep <- which(!is.na(er_status[[1,1]]))

# Keep the observations with non-missing values
er_status_no_missing <- matrix(list(), nrow = 1, ncol = 1)
er_status_no_missing[[1,1]] <- er_status[[1,1]][outcomes_to_keep,,drop=FALSE]

# Retain these observations in the data
tcga_brca_data_no_missing <- tcga_brca_data
tcga_brca_data_no_missing[[1,1]] <- tcga_brca_data[[1,1]][,outcomes_to_keep]
tcga_brca_data_no_missing[[2,1]] <- tcga_brca_data[[2,1]][,outcomes_to_keep]
tcga_brca_data_no_missing[[3,1]] <- tcga_brca_data[[3,1]][,outcomes_to_keep]

# Check
sapply(tcga_brca_data_no_missing, dim)

sample.ids.no.missing <- substr(colnames(tcga_brca_data_no_missing[[1,1]]), 1, 12)
sample.ids.no.missing <- gsub('\\.', '-', sample.ids.no.missing) # Convert . to -

all(sample.ids.no.missing == rownames(er_status_no_missing[[1,1]]))
all(sample.ids.no.missing == rownames(er_status_no_missing[[1,1]]))
all(sample.ids.no.missing == rownames(er_status_no_missing[[1,1]]))

# Fit BIDIFAC to the subsetted data
# rank_init <- BIDIFAC(tcga_brca_data_no_missing, rmt = TRUE, pbar = TRUE, scale_back = FALSE)
# save(rank_init, file = "~/BayesianPMF/04DataApplication/TCGA_BRCA/Cross_Validation/TCGA_BRCA_BIDIFAC_Init_No_Missing_CV.rda")

# Set up the folds
# set.seed(1)
# test_folds <- createFolds(er_status_no_missing[[1,1]], k = 5, list = TRUE, returnTrain = FALSE)
# 
# # Save the test folds
# save(test_folds, file = paste0(results_wd, "/Cross_Validation/ER_Status_Test_Folds_09072023.rda"))

# Load in folds
load(paste0(results_wd, "/Cross_Validation/ER_Status_Test_Folds_09072023.rda"), verbose = TRUE)

# Run cross validation
for (test_fold in 1:length(test_folds)) {

  # Save the indices for the test fold
  inds <- test_folds[[test_fold]]

  # Create a new vector of the outcome with the current pair set to NA
  er_status_cv <- er_status_no_missing
  er_status_cv[[1,1]][inds,] <- NA

  # Run the model with the above pair's continuous outcome missing
  er_status_cv_fit <- BSFP::bsfp(
    data = tcga_brca_data_no_missing,
    Y = er_status_cv,
    nninit = TRUE,
    previous_init = "~/BayesianPMF/04DataApplication/TCGA_BRCA/Cross_Validation/TCGA_BRCA_BIDIFAC_Init_No_Missing_CV.rda",
    model_params = model_params,
    save_structures = FALSE,
    save_loadings_scores = FALSE,
    save_predictive_model = FALSE,
    save_imputations = TRUE,
    nsample = nsample,
    thinning = 10,
    progress = TRUE,
    save_last_sample = TRUE
  )

  # Save the imputed outcomes
  Ym.draw_pair <- er_status_cv_fit$Ym.draw
  
  # Save the last iteration
  last.iter <- er_status_cv_fit$last.iter

  # Save just the relevant output
  save(Ym.draw_pair, last.iter, file = paste0(results_wd, "/Cross_Validation/ER_Status_CV_Fold", test_fold, "_V2.rda"))

  # Remove large objects
  rm(er_status_cv_fit)

  # Garbage collection
  gc()
}

# -------------------------------------
# Results
# -------------------------------------

# Create a matrix with the estimated outcome for each sample
pred.out <- rep(NA, length(unlist(test_folds)))
names(pred.out) <- rownames(er_status_no_missing[[1,1]])

# Load in the results for each fold and calculate the posterior mean after burn-in
for (test_fold in 1:length(test_folds)) {
 
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/Cross_Validation/ER_Status_CV_Fold", test_fold, "_V2.rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # Create a list of vectors from the results
  Ym.draw_pair_burnin_thinned_list <- lapply(101:500, function(iter) Ym.draw_pair[[iter]])
  
  # For testing
  test <- er_status_no_missing[[1,1]][rownames(er_status_no_missing[[1,1]]) %in% current_fold,]
  
  # Aggregate results into a matrix
  Ym.draw_pair_burnin_thinned_matrix <- t(do.call(cbind, Ym.draw_pair_burnin_thinned_list))
  colnames(Ym.draw_pair_burnin_thinned_matrix) <- current_fold
  
  # Calculate the mean response
  Ym.draw_pair_burnin_thinned_mean <- colMeans(Ym.draw_pair_burnin_thinned_matrix)
  
  # Check
  all(colnames(pred.out[names(pred.out) %in% current_fold]) == colnames(Ym.draw_pair_burnin_thinned_matrix)) # TRUE!
  
  # Save these in the matrix
  pred.out[names(pred.out) %in% current_fold] <- Ym.draw_pair_burnin_thinned_mean
}

# Check
all(rownames(er_status_no_missing[[1,1]][unlist(test_folds[1:3])]) == names(pred.out)) # TRUE!

# Compare to the truth
frob(er_status_no_missing[[1,1]]- pred.out)/frob(er_status_no_missing[[1,1]])

# -----------------------------------------------------------------------------
# 5-fold Cross Validation with other models
# -----------------------------------------------------------------------------

# -------------------------------------
# esJIVE
# -------------------------------------

# -------------------------------------
# JIVE
# -------------------------------------

# Fitting JIVE
# tcga_brca_data_no_missing_list <- lapply(1:q, function(s) tcga_brca_data_no_missing[[s,1]])
# JIVE_training_fit <- jive(tcga_brca_data_no_missing_list, center = FALSE, scale = TRUE, method = "perm")
# save(JIVE_training_fit, file = paste0(results_wd, "/JIVE/JIVE_training_data_fit_no_missing.rda"))

# Running with cross validation
run_model_with_cv_tcga(mod = "JIVE", data = tcga_brca_data_no_missing,
                       outcome = er_status_no_missing, outcome_name = "ER_Status",
                       test_folds = test_folds, model_params = model_params,
                       nsample = 5000)

# Load in the results for each fold and calculate the posterior mean after burn-in
pred.out.jive <- rep(NA, nrow(er_status_no_missing[[1,1]]))
names(pred.out.jive) <- rownames(er_status_no_missing[[1,1]])

for (test_fold in 1:length(test_folds)) {
  
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/JIVE/ER_Status_CV_JIVE_Fold_", test_fold, ".rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # CHECK
  all(names(pred.out.jive[names(pred.out.jive) %in% current_fold]) == current_fold) # TRUE!
  
  # Create a list of vectors from the results
  Ym.draw_pair_burnin_thinned_list <- lapply(thinned_iters_burnin, function(iter) Ym.draw_pair[[iter]][[1,1]])
  
  # Aggregate results into a matrix
  Ym.draw_pair_burnin_thinned_matrix <- t(do.call(cbind, Ym.draw_pair_burnin_thinned_list))
  colnames(Ym.draw_pair_burnin_thinned_matrix) <- current_fold
  
  # Calculate the mean response
  Ym.draw_pair_burnin_thinned_mean <- colMeans(Ym.draw_pair_burnin_thinned_matrix)
  
  # Save these in the matrix
  pred.out.jive[names(pred.out.jive) %in% current_fold] <- Ym.draw_pair_burnin_thinned_mean
}

# Compare to the truth
frob(er_status_no_missing[[1,1]] - pred.out.jive)/frob(er_status_no_missing[[1,1]])

# -------------------------------------
# MOFA -- NEED TO CHECK HOW I'VE DEFINED JOINT AND INDIVIDUAL FACTORS
# -------------------------------------

library(MOFA2)

# # Fitting MOFA
# colnames(tcga_brca_data_no_missing_list[[1]]) <-
#   colnames(tcga_brca_data_no_missing_list[[2]]) <-
#   colnames(tcga_brca_data_no_missing_list[[3]]) <- colnames(tcga_brca_data_no_missing_list[[1]])
# 
# mofa_pre_train <- create_mofa(tcga_brca_data_no_missing_list)
# 
# # Set the data options so that the data is not additionally centered and fix the ranks if desired
# data_opts <- get_default_data_options(mofa_pre_train)
# model_opts <- get_default_model_options(mofa_pre_train)
# 
# # No centering
# data_opts$center_groups <- FALSE
# 
# # Create the MOFA object
# MOFAobject <- prepare_mofa(
#   object = mofa_pre_train,
#   model_options = model_opts
# )
# 
# # Train the MOFA model
# MOFA_training_fit <- run_mofa(MOFAobject)
# save(MOFA_training_fit, file = paste0(results_wd, "/MOFA/MOFA_training_data_fit_no_missing.rda"))

# Run with cross validation
run_model_with_cv_tcga(mod = "MOFA", data = tcga_brca_data_no_missing,
                       outcome = er_status_no_missing, outcome_name = "ER_Status",
                       test_folds = test_folds, model_params = model_params,
                       nsample = 5000)

# Load in the results for each fold and calculate the posterior mean after burn-in
pred.out.mofa <- rep(NA, nrow(er_status_no_missing[[1,1]]))
names(pred.out.mofa) <- rownames(er_status_no_missing[[1,1]])

for (test_fold in 1:length(test_folds)) {
  
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/MOFA/ER_Status_CV_MOFA_Fold_", test_fold, ".rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # CHECK
  all(names(pred.out.mofa[names(pred.out.mofa) %in% current_fold]) == current_fold) # TRUE!
  
  # Create a list of vectors from the results
  Ym.draw_pair_burnin_thinned_list <- lapply(thinned_iters_burnin, function(iter) Ym.draw_pair[[iter]][[1,1]])
  
  # Aggregate results into a matrix
  Ym.draw_pair_burnin_thinned_matrix <- t(do.call(cbind, Ym.draw_pair_burnin_thinned_list))
  colnames(Ym.draw_pair_burnin_thinned_matrix) <- current_fold
  
  # Calculate the mean response
  Ym.draw_pair_burnin_thinned_mean <- colMeans(Ym.draw_pair_burnin_thinned_matrix)
  
  # Save these in the matrix
  pred.out.mofa[names(pred.out.mofa) %in% current_fold] <- Ym.draw_pair_burnin_thinned_mean
}

# Compare to the truth
frob(er_status_no_missing[[1,1]] - pred.out.mofa)/frob(er_status_no_missing[[1,1]])

# -------------------------------------
# BIDIFAC
# -------------------------------------

# Fitting BIDIFAC
# BIDIFAC_training_fit <- BIDIFAC(tcga_brca_data_no_missing, rmt = TRUE, pbar = FALSE, scale_back = TRUE)
# save(BIDIFAC_training_fit, file = paste0(results_wd, "/BIDIFAC/BIDIFAC_training_data_fit_no_missing.rda"))

# Run with cross validation
run_model_with_cv_tcga(mod = "BIDIFAC", data = tcga_brca_data_no_missing,
                       outcome = er_status_no_missing, outcome_name = "ER_Status",
                       test_folds = test_folds, model_params = model_params,
                       nsample = 5000)

# Load in the results for each fold and calculate the posterior mean after burn-in
pred.out.bidifac <- rep(NA, nrow(er_status_no_missing[[1,1]]))
names(pred.out.bidifac) <- rownames(er_status_no_missing[[1,1]])

for (test_fold in 1:length(test_folds)) {
  
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/BIDIFAC/ER_Status_CV_BIDIFAC_Fold_", test_fold, ".rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # CHECK
  all(names(pred.out.bidifac[names(pred.out.bidifac) %in% current_fold]) == current_fold) # TRUE!
  
  # Create a list of vectors from the results
  Ym.draw_pair_burnin_thinned_list <- lapply(thinned_iters_burnin, function(iter) Ym.draw_pair[[iter]][[1,1]])
  
  # Aggregate results into a matrix
  Ym.draw_pair_burnin_thinned_matrix <- t(do.call(cbind, Ym.draw_pair_burnin_thinned_list))
  colnames(Ym.draw_pair_burnin_thinned_matrix) <- current_fold
  
  # Calculate the mean response
  Ym.draw_pair_burnin_thinned_mean <- colMeans(Ym.draw_pair_burnin_thinned_matrix)
  
  # Save these in the matrix
  pred.out.bidifac[names(pred.out.bidifac) %in% current_fold] <- Ym.draw_pair_burnin_thinned_mean
}

# Compare to the truth
frob(er_status_no_missing[[1,1]] - pred.out.bidifac)/frob(er_status_no_missing[[1,1]])

# -------------------------------------
# LASSO (combined)
# -------------------------------------

run_model_with_cv_tcga(mod = "LASSO_Combined_Sources", data = tcga_brca_data_no_missing,
                       outcome = er_status_no_missing, outcome_name = "ER_Status",
                       test_folds = test_folds, model_params = model_params,
                       nsample = 5000)

# Load in the results for each fold and calculate the posterior mean after burn-in
pred.out.lasso.combined <- rep(NA, nrow(er_status_no_missing[[1,1]]))
names(pred.out.lasso.combined) <- rownames(er_status_no_missing[[1,1]])

for (test_fold in 1:length(test_folds)) {
  
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/Cross_Validation/LASSO_Combined_Sources/ER_Status_CV_LASSO_Combined_Sources_Fold_", test_fold, ".rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # CHECK
  all(names(pred.out.lasso.combined[names(pred.out.lasso.combined) %in% current_fold]) == current_fold) # TRUE!

  # Save these in the matrix
  pred.out.lasso.combined[names(pred.out.lasso.combined) %in% current_fold] <- Ym.draw_pair
}

# Compare to the truth
frob(er_status_no_missing[[1,1]] - pred.out.lasso.combined)/frob(er_status_no_missing[[1,1]])

# -------------------------------------
# LASSO (expression only)
# -------------------------------------

run_model_with_cv_tcga(mod = "LASSO_Expression_Only", data = tcga_brca_data_no_missing,
                       outcome = er_status_no_missing, outcome_name = "ER_Status",
                       test_folds = test_folds, model_params = model_params,
                       nsample = 5000)

# Load in the results for each fold and calculate the posterior mean after burn-in
pred.out.lasso.expression <- rep(NA, nrow(er_status_no_missing[[1,1]]))
names(pred.out.lasso.expression) <- rownames(er_status_no_missing[[1,1]])

for (test_fold in 1:length(test_folds)) {
  
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/Cross_Validation/LASSO_Expression_Only/ER_Status_CV_LASSO_Expression_Only_Fold_", test_fold, ".rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # CHECK
  all(names(pred.out.lasso.expression[names(pred.out.lasso.expression) %in% current_fold]) == current_fold) # TRUE!
  
  # Save these in the matrix
  pred.out.lasso.expression[names(pred.out.lasso.expression) %in% current_fold] <- Ym.draw_pair
}

# Compare to the truth
frob(er_status_no_missing[[1,1]] - pred.out.lasso.expression)/frob(er_status_no_missing[[1,1]])

# -------------------------------------
# LASSO (methylation only)
# -------------------------------------

run_model_with_cv_tcga(mod = "LASSO_Methylation_Only", data = tcga_brca_data_no_missing,
                       outcome = er_status_no_missing, outcome_name = "ER_Status",
                       test_folds = test_folds, model_params = model_params,
                       nsample = 5000)

# Load in the results for each fold and calculate the posterior mean after burn-in
pred.out.lasso.methylation <- rep(NA, nrow(er_status_no_missing[[1,1]]))
names(pred.out.lasso.methylation) <- rownames(er_status_no_missing[[1,1]])

for (test_fold in 1:length(test_folds)) {
  
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/Cross_Validation/LASSO_Methylation_Only/ER_Status_CV_LASSO_Methylation_Only_Fold_", test_fold, ".rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # CHECK
  all(names(pred.out.lasso.methylation[names(pred.out.lasso.methylation) %in% current_fold]) == current_fold) # TRUE!
  
  # Save these in the matrix
  pred.out.lasso.methylation[names(pred.out.lasso.methylation) %in% current_fold] <- Ym.draw_pair
}

# Compare to the truth
frob(er_status_no_missing[[1,1]] - pred.out.lasso.methylation)/frob(er_status_no_missing[[1,1]])

# -------------------------------------
# LASSO (miRNA only)
# -------------------------------------

run_model_with_cv_tcga(mod = "LASSO_miRNA_Only", data = tcga_brca_data_no_missing,
                       outcome = er_status_no_missing, outcome_name = "ER_Status",
                       test_folds = test_folds, model_params = model_params,
                       nsample = 5000)

# Load in the results for each fold and calculate the posterior mean after burn-in
pred.out.lasso.miRNA <- rep(NA, nrow(er_status_no_missing[[1,1]]))
names(pred.out.lasso.miRNA) <- rownames(er_status_no_missing[[1,1]])

for (test_fold in 1:length(test_folds)) {
  
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/Cross_Validation/LASSO_miRNA_Only/ER_Status_CV_LASSO_miRNA_Only_Fold_", test_fold, ".rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # CHECK
  all(names(pred.out.lasso.miRNA[names(pred.out.lasso.miRNA) %in% current_fold]) == current_fold) # TRUE!
  
  # Save these in the matrix
  pred.out.lasso.miRNA[names(pred.out.lasso.miRNA) %in% current_fold] <- Ym.draw_pair
}

# Compare to the truth
frob(er_status_no_missing[[1,1]] - pred.out.lasso.miRNA)/frob(er_status_no_missing[[1,1]])

# -------------------------------------
# multiview
# -------------------------------------

run_model_with_cv_tcga(mod = "multiview", data = tcga_brca_data_no_missing,
                       outcome = er_status_no_missing, outcome_name = "ER_Status",
                       test_folds = test_folds, model_params = model_params,
                       nsample = 5000)

# Load in the results for each fold and calculate the posterior mean after burn-in
pred.out.multiview <- rep(NA, nrow(er_status_no_missing[[1,1]]))
names(pred.out.multiview) <- rownames(er_status_no_missing[[1,1]])

for (test_fold in 1:length(test_folds)) {
  
  # Load in the results
  load(paste0("~/BayesianPMF/04DataApplication/TCGA_BRCA/multiview/ER_Status_CV_multiview_Fold_", test_fold, ".rda"), verbose = TRUE)
  
  # Save the patient IDs for this test fold
  current_fold <- rownames(er_status_no_missing[[1,1]])[test_folds[[test_fold]]]
  
  # CHECK
  all(names(pred.out.multiview[names(pred.out.multiview) %in% current_fold]) == current_fold) # TRUE!
  
  # Save these in the matrix
  pred.out.multiview[names(pred.out.multiview) %in% current_fold] <- Ym.draw_pair
}

# Compare to the truth
frob(er_status_no_missing[[1,1]] - pred.out.multiview)/frob(er_status_no_missing[[1,1]])
