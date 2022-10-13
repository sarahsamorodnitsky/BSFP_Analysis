# -----------------------------------------------------------------------------
# Applying the BPMF model to the HIV-COPD BALF Biocrates and Somascan data. 
# Studying the sensitivity of the results to the alignment by studying the
# factor loadings and scores summaries after aligning to one posterior 
# sampling iteration before and after the pivot used in the analysis. 
# * This analysis is based on the V2 HIV-COPD results, where the data was 
#   centered but not scaled to overall variance 1. 
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set-up
# -----------------------------------------------------------------------------

# Working directories
data_wd <- "~/BSFP/data/"
model_wd <- "~/BSFP/"
results_wd <- "~/BayesianPMF/04DataApplication/"

# Loading in the model functions
source(paste0(model_wd, "bpmf.R"))

# Load in the original, unaligned results 
load(paste0("~/BayesianPMF/04DataApplication/BPMF/Training_Fit/training_data_fit_V2.rda"), verbose = TRUE)

# Load in the sensitivity analysis results
load("/home/samorodnitsky/BayesianPMF/04DataApplication/BPMF/Training_Fit/FEV1pp_Factor_Switching_Sensitivity_V2_NewPivot_Burnin.rda", verbose = TRUE)

# Loading data in
data_wd <- "~/BSFP/data/"
load(paste0(data_wd, "BiocratesLavageProcessed.rda"), verbose = TRUE)
load(paste0(data_wd, "HIV_COPD_SomaScan_Normalized_Clean.rda"), verbose = TRUE)

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

# Data dimensions
p1 <- nrow(lavage_processed_no_info_log_scale)
p2 <- nrow(somascan_normalized_clean_no_info_transpose_reorder)
p <- p1 + p2
p.vec <- c(p1, p2)
n <- ncol(lavage_processed_no_info_log_scale)
q <- length(p.vec)

# Gibbs sampler parameters
nsample <- 5000
burnin <- 1000
iters_burnin <- seq(burnin+1, nsample)
nsample_after_burnin <- nsample-burnin

# Set up the rank indices
ranks <- fev1pp_training_fit_nonsparse_V2$ranks
rank.inds <- lapply(1:(q+1), function(s) {
  if (s == 1) {
    1:ranks[s]
  } else {
    (cumsum(ranks[1:(s-1)])[s-1] + 1):cumsum(ranks[1:s])[s]
  }
})

# -----------------------------------------------------------------------------
# Constructing summaries of loadings and scores from the different alignments 
# -----------------------------------------------------------------------------

# Summaries for pivot -1
factor_summaries(joint.loadings.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$joint.loadings.final,
                 joint.scores.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$joint.scores.final,
                 joint.betas.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$joint.betas.final,
                 individual.loadings.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$individual.loadings.final,
                 individual.scores.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$individual.scores.final,
                 individual.betas.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$individual.betas.final,
                 ranks = fev1pp_training_fit_nonsparse_V2$ranks,
                 nsample_after_burnin = nsample_after_burnin,
                 exploring_factors_wd = "~/BayesianPMF/04DataApplication/BPMF/Exploring_Factors/V2/Aligned/NewPivot_Burnin_Sensitivity/Minus_1",
                 lavage_processed_no_info_log_scale = lavage_processed_no_info_log_scale,
                 somascan_normalized_clean_no_info_transpose_scale = somascan_normalized_clean_no_info_transpose_scale)

# Summaries for pivot +1
factor_summaries(joint.loadings.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[2]]$joint.loadings.final,
                 joint.scores.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[2]]$joint.scores.final,
                 joint.betas.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[2]]$joint.betas.final,
                 individual.loadings.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[2]]$individual.loadings.final,
                 individual.scores.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[2]]$individual.scores.final,
                 individual.betas.final = fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[2]]$individual.betas.final,
                 ranks = fev1pp_training_fit_nonsparse_V2$ranks,
                 nsample_after_burnin = nsample_after_burnin,
                 exploring_factors_wd = "~/BayesianPMF/04DataApplication/BPMF/Exploring_Factors/V2/Aligned/NewPivot_Burnin_Sensitivity/Plus_1",
                 lavage_processed_no_info_log_scale = lavage_processed_no_info_log_scale,
                 somascan_normalized_clean_no_info_transpose_scale = somascan_normalized_clean_no_info_transpose_scale)


# -----------------------------------------------------------------------------
# Calculating the decomposition of the covariance matrix as suggested in 
# Poworoznek et al. 8
# -----------------------------------------------------------------------------

# For alignment to pivot -1 iteration
U.draw <- fev1pp_training_fit_nonsparse_V2$U.draw[iters_burnin]
W.draw <- fev1pp_training_fit_nonsparse_V2$W.draw[iters_burnin]
joint.betas.draw <- lapply(iters_burnin, function(iter) {
  fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][rank.inds[[1]],,drop=FALSE]
})
individual.betas.draw <- lapply(1:q, function(s) {
  lapply(iters_burnin, function(iter) {
    fev1pp_training_fit_nonsparse_V2$beta.draw[[iter]][[1,1]][rank.inds[[s+1]],,drop=FALSE]
  })
})

U.draw.aligned <- fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$joint.loadings.final
W.draw.aligned <- fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$individual.loadings.final
joint.betas.draw.aligned <-fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$joint.betas.final
individual.betas.draw.aligned <- fev1pp_training_fit_nonsparse_aligned_V2_sensitivity[[1]]$individual.betas.final
