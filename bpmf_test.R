# -----------------------------------------------------------------------------
# Contains some fake parameters to use for testing the model
# -----------------------------------------------------------------------------

# Setting up the data
data <- matrix(list(), ncol = 1, nrow = 2)
dat1 <- matrix(rnorm(1000), nrow = 100, ncol = 10)
dat2 <- matrix(rnorm(1500), nrow = 150, ncol = 10)
data[1,1][[1]] <- dat1
data[2,1][[1]] <- dat2

Y <- NULL

# Setting up the model parameters
model_params <- list(error_vars = c(1,1),
                     joint_var = 1,
                     indiv_vars = c(1,1),
                     beta_vars = NULL,
                     response_vars = NULL)
