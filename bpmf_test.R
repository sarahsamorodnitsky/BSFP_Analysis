# -----------------------------------------------------------------------------
# Model testing for various scenarios
# * Assessing the coverage of the various parameters informally 
# * Assessing if the log joint density function works
# -----------------------------------------------------------------------------

# Load in the helper functions
source("bpmf.R")

# Setting up the data
n <- 50
p.vec <- c(150, 100)
r <- 1
r.vec <- c(1, 1)
ranks <- c(r, r.vec)
q <- 2

# Setting up the model parameters
model_params <- true_params <- list(error_vars = c(1,1),
                                    joint_var = 1,
                                    indiv_vars = c(1,1),
                                    beta_vars = c(10, 1, rep(1, q)), # Use the same variance for all the effects from each source
                                    response_vars = c(shape = 1,rate = 1))

# Generating example data
data <- bpmf_data(p.vec, n, ranks, true_params, response = "continuous", missingness = NULL, entrywise = NULL, prop_missing = NULL)
X <- data$data
y <- data$Y

# -----------------------------------------------------------------------------
# Testing the coverage
# -----------------------------------------------------------------------------

# No response
test <- bpmf_sim(nsample = 5, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks)

# Missing data
test_with_missing0.1 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.1)
test_with_missing0.3 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.3)
test_with_missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.5)
test_with_missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.7)
test_with_missing0.9 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.9)

# Continuous response
test_with_response_continuous <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous")

# Binary response
test_with_response_binary <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary")

# Missing continuous response
test_with_response_continuous_missing0.1 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.1)
test_with_response_continuous_missing0.3 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.3)
test_with_response_continuous_missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.5)
test_with_response_continuous_missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.7)
test_with_response_continuous_missing0.9 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous", missingness = "missingness_in_response", entrywise = NULL, prop_missing = 0.9)


# Save results
save(test, test_with_response_continuous, test_with_missing0.1, test_with_missing0.3, test_with_missing0.5,
     test_with_missing0.7, test_with_missing0.9, file = "~/BayesianPMFWithGit/validation_testing/validation_results_1292021.rda")


# -----------------------------------------------------------------------------
# Testing the convergence function
# -----------------------------------------------------------------------------

# Running the Gibbs sampler
bpmf_test <- bpmf(data, Y = NULL, nninit = TRUE, model_params = model_params, nsample = 2000)

# Saving the parameters
U.draw <- bpmf_test$U.draw
V.draw <- bpmf_test$V.draw
W.draw <- bpmf_test$W.draw
Vs.draw <- bpmf_test$Vs.draw
ranks <- bpmf_test$ranks

# For each iteration, calculate the log joint density
log_joint_density_over_iters <- c()
for (iter in 1:2000) {
  log_joint_density_over_iters[iter] <- log_joint_density(data = data, U.iter = U.draw[[iter]], V.iter = V.draw[[iter]], W.iter = W.draw[[iter]], Vs.iter = Vs.draw[[iter]], model_params = model_params, ranks = ranks)
}

# I think it looks good!
plot(log_joint_density_over_iters)

# How much difference is there between each point at the tail end?
diff(log_joint_density_over_iters)

# -----------------------------------------------------------------------------
# Testing obtaining structure from a different method and then inputting it into
# the model function
# -----------------------------------------------------------------------------

# Running the BIDIFAC model
bidifac.test <- BIDIFAC(X, rmt = TRUE, pbar = FALSE)

# Saving the column structure (the joint structure)
bidifac.joint <- lapply(1:q, function(source) {
  bidifac.test$C[[source,1]]
})
bidifac.joint.rank <- rankMatrix(bidifac.test$C[[1,1]])[1]

# Saving the individual structure 
bidifac.individual <- lapply(1:q, function(source) {
  bidifac.test$I[[source,1]]
})
bidifac.indiv.rank <- sapply(bidifac.test$I, function(source) {
  rankMatrix(source)[1]
}) # Not sure this is the most generalizable approach

# Obtaining the joint scores
bidifac.joint.scores <- svd(bidifac.joint[[1]])$v[,1:bidifac.joint.rank]

# Obtaining the individual scores
bidifac.indiv.scores <- do.call(cbind, lapply(1:q, function(source) {
  svd.source <- svd(bidifac.individual[[source]])
  svd.source$v[,1:bidifac.indiv.rank[source]]
}))

# Combining the scores together
bidifac.all.scores <- cbind.data.frame(y[[1,1]], bidifac.joint.scores, bidifac.indiv.scores)
colnames(bidifac.all.scores) <- c("y", "joint1", "indiv1", "indiv2")

# Fitting the Bayesian linear model
bidifac.bayes <- bpmf(data = X, Y = y, nninit = FALSE, model_params = model_params, 
                      ranks = ranks, scores = as.matrix(bidifac.all.scores[,-1]), nsample = 1000)

