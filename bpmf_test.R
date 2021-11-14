# -----------------------------------------------------------------------------
# Model testing for various scenarios
# -----------------------------------------------------------------------------

# Load in the helper functions
source("bpmf.R")

# Setting up the data
n <- 20
p.vec <- c(100, 150)
r <- 5
r.vec <- c(2, 5)
ranks <- c(r, r.vec)
q <- 2
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

# Setting up the model parameters
model_params <- list(error_vars = c(1,1),
                     joint_var = 1,
                     indiv_vars = c(1,1),
                     beta_vars = c(10, 1, rep(1, q)), # Use the same variance for all the effects from each source
                     response_vars = c(1,1))

true_params <- model_params

nsample <- 10
n_clust <- 10
nsim <- 10
s2n <- 1
nninit <- FALSE
response <- NULL
missingness <- NULL
prop_missing <- 0.1
entrywise = TRUE

# -----------------------------------------------------------------------------
# No outcome
# -----------------------------------------------------------------------------

Y <- NULL
no.out <- bpmf(data, Y = Y, nninit = TRUE, model_params, ranks = NULL, nsample = 10, progress = TRUE)

# -----------------------------------------------------------------------------
# Continuous outcome
# -----------------------------------------------------------------------------

Y <- matrix(rnorm(n), nrow = n)
continuous.out <- bpmf(data, Y = Y, nninit = TRUE, model_params, ranks = NULL, nsample = 10, progress = TRUE)

# -----------------------------------------------------------------------------
# Binary outcome
# -----------------------------------------------------------------------------

Y <- matrix(rbinom(n, size = 1, p = 0.5), nrow = n)
bin.out <- bpmf(data, Y = Y, nninit = TRUE, model_params, ranks = NULL, nsample = 10, progress = TRUE)

# -----------------------------------------------------------------------------
# Missing data 
# -----------------------------------------------------------------------------

data.missing <- data
for (s in 1:q) {
  n_s <- length(data[[s,1]])
  random_missing <- sample(1:n_s, size = 0.1 * n_s, replace = FALSE)
  data.missing[[s,1]][random_missing] <- NA
}

Y <- NULL
missing.data.out <- bpmf(data.missing, Y = Y, nninit = TRUE, model_params, ranks = NULL, nsample = 10, progress = TRUE)

# -----------------------------------------------------------------------------
# Missing continuous response
# -----------------------------------------------------------------------------

Y.missing <- matrix(rnorm(n), nrow = n)
Y.missing[sample(1:n, size = 0.1 * n, replace = FALSE)] <- NA
continuous.out <- bpmf(data, Y = Y.missing, nninit = TRUE, model_params, ranks = NULL, nsample = 10, progress = TRUE)

# -----------------------------------------------------------------------------
# Missing binary response
# -----------------------------------------------------------------------------

Y.missing <- matrix(rbinom(n, size = 1, p = 0.5), nrow = n)
Y.missing[sample(1:n, size = 0.1 * n, replace = FALSE)] <- NA
continuous.out <- bpmf(data, Y = Y.missing, nninit = TRUE, model_params, ranks = NULL, nsample = 10, progress = TRUE)

# -----------------------------------------------------------------------------
# Both
# -----------------------------------------------------------------------------

# Missing data and missing continuous
Y.missing <- matrix(rnorm(n), nrow = n)
Y.missing[sample(1:n, size = 0.1 * n, replace = FALSE)] <- NA
continous.out.missing <- bpmf(data.missing, Y = Y.missing, nninit = TRUE, model_params, ranks = NULL, nsample = 10, progress = TRUE)

# Missing data and missing binary
Y.missing <- matrix(rbinom(n, size = 1, p = 0.5), nrow = n)
Y.missing[sample(1:n, size = 0.1 * n, replace = FALSE)] <- NA
continous.out.missing <- bpmf(data.missing, Y = Y.missing, nninit = TRUE, model_params, ranks = NULL, nsample = 10, progress = TRUE)


# -----------------------------------------------------------------------------
# Testing the coverage
# -----------------------------------------------------------------------------

test <- bpmf_sim(nsample = 1000, n_clust = 10, p.vec, n, true_params, model_params, nsim = 100, s2n = s2n, nninit = FALSE, ranks)

test_with_resp <- bpmf_sim(nsample = 10, n_clust = 10, p.vec, n, true_params, model_params, nsim = 10, s2n = s2n, nninit = FALSE, ranks, response = "continuous")

test_with_resp_bin <- bpmf_sim(nsample = 10, n_clust = 10, p.vec, n, true_params, model_params, nsim = 10, s2n = s2n, nninit = FALSE, ranks, response = "binary")

test_with_resp_missing <- bpmf_sim(nsample = 10, n_clust = 10, p.vec, n, true_params, model_params, nsim = 10, s2n = s2n, nninit = FALSE, ranks, 
                                   missingness = "missingness_in_data", prop_missing = 0.1, entrywise = TRUE)

