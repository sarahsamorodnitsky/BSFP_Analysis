# -----------------------------------------------------------------------------
# Example application of BSPF with alignment and posterior summaries 
# using simulated data. 
# -----------------------------------------------------------------------------

# Load in the functions
source("~/BSFP_Analysis/BSFP_functions.R")

# -----------------------------------------------------------------------------
# Defining parameters for testing the functions
# -----------------------------------------------------------------------------

# Setting up the data
n <- 50
p.vec <- c(75, 100)
q <- 2

# Setting up the model parameters
true_params <- list(error_vars = c(1,1),
                    joint_var = 1,
                    indiv_vars = c(1,1),
                    beta_vars = c(1, 1, rep(1, q)), 
                    response_vars = c(shape = 1, rate = 1))

# Choose ranks
r <- 3
r.vec <- c(3, 3)
ranks <- c(r, r.vec)

# Number of posterior sampling iterations
nsample <- 1000
burnin <- nsample/2
iters_burnin <- (burnin+1):nsample

# -----------------------------------------------------------------------------
# Example 1: No missingness
# -----------------------------------------------------------------------------

# Generate data
data.c1 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE)

# Run BSFP for 1000 iterations
bsfp.c1 <- bsfp(data = data.c1$data, Y = data.c1$Y, nsample = nsample)

# Check convergence
log_joint_density_by_iter.c1 <- sapply(1:nsample, function(iter) {
  log_joint_density(data = data.c1$data, Y = data.c1$Y, 
                    U.iter = bsfp.c1$U.draw[[iter]],
                    V.iter = bsfp.c1$V.draw[[iter]],
                    W.iter = bsfp.c1$W.draw[[iter]],
                    Vs.iter = bsfp.c1$Vs.draw[[iter]],
                    model_params = bsfp.c1$model_params,
                    ranks = bsfp.c1$ranks,
                    beta.iter = bsfp.c1$beta.draw[[iter]],
                    tau2.iter = bsfp.c1$tau2.draw[[iter]])
})

# Expect to see increase and stabilize over convergence
plot(log_joint_density_by_iter.c1)

# Check some trace plots of structures 

# Source 1 joint structure, (1,1) entry
plot(sapply(1:nsample, function(iter) {
  bsfp.c1$J.draw[[iter]][[1,1]][1,1]
}))

# Source 2 individual structure, (5,10) entry
plot(sapply(1:nsample, function(iter) {
  bsfp.c1$A.draw[[iter]][[2,1]][5,10]
}))

# Run the alignment algorithm
alignment.c1 <- match_align_bsfp(BSFP.fit = bsfp.c1, y = data.c1$Y, 
                                 model_params = bsfp.c1$model_params,
                                 p.vec = p.vec, iters_burnin = iters_burnin)

# Summarize aligned factors
summary.aligned.c1 <- summarize_factors(data = data.c1$data, Y = data.c1$Y,
                                        iters_burnin = iters_burnin,
                                        aligned_results = alignment.c1,
                                        ranks = bsfp.c1$ranks, tau2.draw = bsfp.c1$tau2.draw)

# -----------------------------------------------------------------------------
# Example 2: Entrywise missingness
# -----------------------------------------------------------------------------

# Generate data
data.c2 <- bsfp_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = "continuous", sparsity = FALSE,
                     missingness = "missingness_in_data", missing_data_type = "entrywise", prop_missing = 0.1)

# Run BSFP for 1000 iterations
bsfp.c2 <- bsfp(data = data.c2$missing_data, Y = data.c2$Y, nsample = nsample)

# Check convergence
log_joint_density_by_iter.c2 <- sapply(1:nsample, function(iter) {
  log_joint_density(data = data.c2$data, Y = data.c2$Y, 
                    U.iter = bsfp.c2$U.draw[[iter]],
                    V.iter = bsfp.c2$V.draw[[iter]],
                    W.iter = bsfp.c2$W.draw[[iter]],
                    Vs.iter = bsfp.c2$Vs.draw[[iter]],
                    model_params = bsfp.c2$model_params,
                    ranks = bsfp.c2$ranks,
                    beta.iter = bsfp.c2$beta.draw[[iter]],
                    tau2.iter = bsfp.c2$tau2.draw[[iter]],
                    Xm.iter = bsfp.c2$Xm.draw[[iter]])
})

# Expect to see increase and stabilize over convergence
plot(log_joint_density_by_iter.c2)

# Check some trace plots of structures 

# Source 1 joint structure, (1,1) entry
plot(sapply(1:nsample, function(iter) {
  bsfp.c2$J.draw[[iter]][[1,1]][1,1]
}))

# Source 2 individual structure, (5,10) entry
plot(sapply(1:nsample, function(iter) {
  bsfp.c2$A.draw[[iter]][[2,1]][5,10]
}))

# Run the alignment algorithm
alignment.c2 <- match_align_bsfp(BSFP.fit = bsfp.c2, y = data.c2$Y, 
                                 model_params = bsfp.c2$model_params,
                                 p.vec = p.vec, iters_burnin = iters_burnin)

# Summarize aligned factors
summary.aligned.c2 <- summarize_factors(data = data.c2$missing_data, Y = data.c2$Y,
                                        iters_burnin = iters_burnin,
                                        aligned_results = alignment.c2,
                                        ranks = bsfp.c2$ranks, tau2.draw = bsfp.c2$tau2.draw,
                                        Xm.draw = bsfp.c2$Xm.draw)
