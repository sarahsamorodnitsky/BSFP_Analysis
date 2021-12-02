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
data <- matrix(list(), ncol = 1, nrow = q)

V <- matrix(list(), nrow = 1, ncol = 1)
V[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(10)), nrow = n, ncol = r)

U <- matrix(list(), nrow = q, ncol = 1)
Vs <- matrix(list(), nrow = 1, ncol = q)
W <- matrix(list(), nrow = q, ncol = q)

for (s in 1:q) {
  U[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(10)), nrow = p.vec[s], ncol = r)

  Vs[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(5)), nrow = n, ncol = r.vec[s])

  W[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(5)), nrow = p.vec[s], ncol = r.vec[s])

  for (ss in 1:q) {
    if (ss != s) {
      W[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
    }
  }
}

for (s in 1:q) {
  data[s,1][[1]] <- U[[s,1]] %*% t(V[[1,1]]) + W[[s,s]] %*% t(Vs[[1,s]]) + matrix(rnorm(p.vec[s]*n), nrow = p.vec[s], ncol = n)
}

# -----------------------------------------------------------------------------
# Testing the coverage
# -----------------------------------------------------------------------------

# No response
test <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks)

# Continuous response
test_with_response_continuous <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "continuous")

# Binary response
test_with_response_binary <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, response = "binary")

# Missing data
test_with_missing0.1 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.1)
test_with_missing0.3 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.3)
test_with_missing0.5 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.5)
test_with_missing0.7 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.7)
test_with_missing0.9 <- bpmf_sim(nsample = 2000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = NULL, center = FALSE, nninit = FALSE, ranks = ranks, missingness = "missingness_in_data", entrywise = TRUE, prop_missing = 0.9)

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




