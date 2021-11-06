# -----------------------------------------------------------------------------
# Contains some fake parameters to use for testing the model
# -----------------------------------------------------------------------------

# Setting up the data
n <- 20
p.vec <- c(100, 150)
r <- 5
r.vec <- c(2, 5)
data <- matrix(list(), ncol = 1, nrow = 2)

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

Y <- matrix(rnorm(n), nrow = n)

# Setting up the model parameters
model_params <- list(error_vars = c(1,1),
                     joint_var = 1,
                     indiv_vars = c(1,1),
                     beta_vars = NULL,
                     response_vars = NULL)
