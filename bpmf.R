# -----------------------------------------------------------------------------
# Helper functions for Bayesian PMF
# -----------------------------------------------------------------------------

# Packages
library(doParallel)
library(foreach)
library(Matrix)
library(MASS)
library(truncnorm)
# library(natural)
# library(RSpectra)
# library(MCMCpack)

# -----------------------------------------------------------------------------
# Bayesian PMF functions
# -----------------------------------------------------------------------------

# This version of BPMF is initialized with BIDIFAC+ with Y as a source
bpmf_full_mode <- function(data, Y, nninit = TRUE, model_params, ranks = NULL, scores = NULL, nsample, progress = TRUE, starting_values = NULL, err.y.est = TRUE) {
  # Gibbs sampling algorithm for sampling the underlying structure and the 
  # regression coefficient vector for a response vector, Y, if given
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # 
  # data = matrix with list entries corresponding to each data source. assumed to be row-centered
  # Y = column vector with outcome or NULL. assumed to be centered
  # nninit = should the model be initialized with a nuclear norm penalized objective? if FALSE, provide ranks
  # model_params = (error_vars, joint_vars, indiv_vars, beta_vars = NULL, response_vars)
  # ranks = vec of specific ranks if not nninit. (1st entry = joint rank, kth for k>1 is the individual rank for the k-1'st source)
  # scores = if using structure from another method to fit a linear model, provide joint and individual scores here.
  #   in this case, ranks should be provided and nninit = FALSE, Y != NULL
  # nsample = number of Gibbs sampling iterations
  # progress = should the progress of the sampler be displayed?
  # starting_values = list of starting values for V, U, W, Vs. If NULL and nninit = TRUE, init with BIDIFAC+,
  #    if NULL and nninit = FALSE, init from prior. If not NULL, will init with provided starting values unless 
  #    nninit. 
  # err.y.est (logical): should the data be scaled to have error variance 1?
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  q <- nrow(data) # Number of sources
  p.vec <- apply(data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(data[[1,1]]) # Number of subjects
  
  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------
  
  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  # response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response (only needed if scores are provided by another method)
  
  # ---------------------------------------------------------------------------
  # Check for missingness in data
  # ---------------------------------------------------------------------------
  
  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
  # If there is missingness, initialize the missing values with 0s
  for (s in 1:q) {
    data[[s,1]][missing_obs[[s]]] <- 0
  }
  
  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------
  
  response_given <- !is.null(Y[[1,1]]) 
  
  # If so, what kind of response is it?
  if (response_given) {
    Y <- matrix(unlist(Y))
    
    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"
    
    # If there is a response, is there missingness in the outcome?
    missingness_in_response <- any(is.na(Y))
    
    # Which entries are missing?
    missing_obs_Y <- which(is.na(Y))
    
    # If there are missing entries, initialize them with 0s
    Y[missing_obs_Y,] <- 0
  }
  
  if (!response_given) {
    response_type <- missing_obs_Y <- NULL
    missingness_in_response <- FALSE
  }
  
  # ---------------------------------------------------------------------------
  # Scaling the data to have error variance 1
  # ---------------------------------------------------------------------------
  
  # Duplicate the data to have a scaled version
  scaled_data <- data
  scaled_Y <- Y
  
  # If initializing with BIDIFAC+, scale data to have error variance 1
  if (nninit & err.y.est) {
    
    # Create a matrix to store the estimated error sd
    if (!response_given) {
      sigma.mat <- matrix(nrow = q, ncol = 1)
    }
    
    if (response_given) {
      sigma.mat <- matrix(nrow = q+1, ncol = 1)
    }
    
    # Iterate through the sources and estimate the error variance
    for (s in 1:q) {
      
      # Save the estimated error variance
      sigma.mat[s,] <- sigma.rmt(data[[s,1]])
      
      # Scale the data
      scaled_data[[s,1]] <- data[[s,1]]/sigma.mat[s,]
    }
    
    # If a response vector is given, scale Y
    if (response_given) {
      
      # Save the estimated error variance using the natural lasso
      data_combined <- do.call(rbind, scaled_data)
      nl_cv <- nlasso_cv(x = t(data_combined), y = Y)
      sigma.mat[s+1,] <- nl_cv$sig_obj
      
      # Scale the response
      scaled_Y <- Y/sigma.mat[s+1,]
    }
  }
  
  # If not initializing with BIDIFAC+, return a sigma matrix with 1s
  if (!nninit | (nninit & !err.y.est)) {
    
    # Create a matrix of 1s
    if (!response_given) {
      sigma.mat <- matrix(1, nrow = q, ncol = 1)
    }
    
    if (response_given) {
      sigma.mat <- matrix(1, nrow = q+1, ncol = 1)
    }
  } 
  
  # ---------------------------------------------------------------------------
  # Obtaining the ranks 
  # ---------------------------------------------------------------------------
  
  # Initializing with BIDIFAC+
  if (nninit) {
    
    # Combine the data together 
    data_combined <- matrix(list(), nrow = (q+1), ncol = 1)
    
    for (s in 1:q) {
      # Append the scaled data
      data_combined[[s,1]] <- scaled_data[[s,1]]
      
      # Add back in missing data if any
      data_combined[[s,1]][missing_obs[[1]]] <- NA
    }
    
    # Initialize the indices of observations in each source 
    p.ind <- lapply(1:q, function(s) {
      if (s == 1) {
        1:p.vec[s]
      } else {
        (p.vec[s-1] + 1):cumsum(p.vec)[s]
      }
    })
    
    # Save indices of samples per source
    n.ind <- list(1:n)
    
    # If no response is given, p.ind.list should only include the sources
    if (!response_given) {
      
      # Save the indices for features for each identified structure
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- p.ind[[s]]
      }
    }
    
    # Include the response as a source if a response is given and add it to p.ind, p.ind.list
    if (response_given) {
      
      # Append the response to the sources
      data_combined[[q+1,1]] <- t(scaled_Y)
      
      # Add back in any missingness if exists
      data_combined[[q+1,1]][missing_obs_Y] <- NA
      
      # Add the response indices to p.ind
      p.ind[[q+1]] <- cumsum(p.vec)[q] + 1 # For the response
      
      # Save the indices for features for each identified structure
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- c(p.ind[[s]], p.ind[[q+1]])
      }
      
    }
    
    # Collapse the data into one matrix
    data_combined <- do.call(rbind, data_combined)
    
    # Save the indices for samples in each identified structure
    n.ind.list <- lapply(1:(q+1), function(s) c(1:n))
    
    # Run BIDIFAC+ --
    
    # If there is no missing data
    if (!missingness_in_data & !missingness_in_response) {
      rank_init <- bidifac.plus.given(data_combined, p.ind = p.ind, n.ind = n.ind,
                                      p.ind.list = p.ind.list, n.ind.list = n.ind.list)
    }
    
    # If there is missing data
    if (missingness_in_data | missingness_in_response) {
      # Save the indices of the missing values
      all.miss <- which(is.na(data_combined))
      
      # Initialize
      rank_init <- bidifac.plus.impute(data_combined, p.ind = p.ind, n.ind = n.ind,
                                      p.ind.list = p.ind.list, n.ind.list = n.ind.list,
                                      all.miss = all.miss)
    }
    
    # Print when finished
    print("Posterior mode obtained, ranks determined.")
    
    # Obtain the structures and their ranks
    S <- rank_init$S
    
    # Joint rank
    r <- Matrix::rankMatrix(S[[1]])[[1]]
    
    # Individual ranks if no response or pairwise-shared ranks of each source with Y if response given
    r.vec <- sapply(2:(q+1), function(s) Matrix::rankMatrix(S[[s]])[[1]])
    
  }
  
  if (!nninit) {
    r <- ranks[1]
    r.vec <- ranks[-1]
  }
  
  r_total <- n_beta <- r + sum(r.vec)
  
  # Warning if the rank is 0
  if (n_beta == 0) {
    warning("Rank of decomposition is 0. No structure will be generated.")
  }
  
  # If a response is given, set up the variance matrix for the prior of the betas using the ranks and 
  # set up rank indices for betas
  if (response_given) {
    
    # Setting up prior covariance matrix
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    beta_vars <- rep(beta_vars, c(r, r.vec))
    diag(Sigma_beta) <- beta_vars
    
    # Setting up rank indices
    r.ind <- lapply(1:(q+1), function(s) list())
    r.ind[[1]] <- 1:r # Indices for joint factors
    
    for (s in 1:q) {
      
      if (s == 1) {
        r.ind[[s+1]] <- r + (1:r.vec[s])
      }
      
      if (s > 1) {
        r.ind[[s+1]] <- r + cumsum(r.vec[s-1]) + (1:r.vec[s])
      }
      
    }
  }
  
  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  U.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  W.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = q))
  
  if (!response_given) {
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (!missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  if (response_given) {
    beta.draw <- Z.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
  }
  
  if (missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  # ---------------------------------------------------------------------------
  # Initialize V, U, V, W
  # ---------------------------------------------------------------------------
  
  # If initializing with nuclear norm, initialize sampling at posterior mode
  if (nninit) {
    
    # Initialize joint scores, V
    
    V0 <- matrix(list(), nrow = 1, ncol = 1)
    beta_joint0 <- matrix(list(), nrow = 1, ncol = 1)
    
    # If there is joint structure
    if (r > 0) {
      svd.joint <- svd(S[[1]])
      V0[[1,1]] <- (svd.joint$v[,1:r, drop = FALSE]) %*% diag(svd.joint$d[1:r], nrow = r)
      
      # Save beta_joint if response is given
      if (response_given) {
        beta_joint0[[1,1]] <- t(svd(S[[1]])$u[p.ind[[q+1]],1:r, drop = FALSE])
      }
    } 
    
    # If there is no joint structure
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      
      # Save beta joint if response is given
      if (response_given) {
        beta_joint0[[1,1]] <- matrix(0, nrow = r, ncol = 1)
      }
    }
    
    U0 <- matrix(list(), nrow = q, ncol = 1)
    Vs0 <- matrix(list(), nrow = 1, ncol = q)
    W0 <- matrix(list(), nrow = q, ncol = q)
    beta_indiv0 <- matrix(list(), nrow = q, ncol = 1)
    
    for (s in 1:q) {
      
      # Initialize joint loadings, U 
      if (r > 0) {
        U0[[s,1]] <- svd(S[[1]])$u[p.ind[[s]],1:r, drop = FALSE]
      } 
      if (r == 0) {
        U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
      }
      
      # Initialize individual loadings, W, and individual scores, Vs, and beta_indiv if response is given
      if (r.vec[s] > 0) {
        
        # Compute SVD
        svd.indiv.s <- svd(S[[s+1]])
        
        # Save scores and loadings
        Vs0[[1,s]] <- (svd.indiv.s$v[,1:r.vec[s], drop = FALSE]) %*% diag(svd.indiv.s$d[1:r.vec[s]], nrow = r.vec[s])
        W0[[s,s]] <- svd.indiv.s$u[p.ind[[s]],1:r.vec[s], drop = FALSE]
        
        # Save beta_indiv
        if (response_given) {
          beta_indiv0[[s,1]] <- t(svd.indiv.s$u[p.ind[[q+1]],1:r.vec[s], drop = FALSE])
        }
        
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      } 
      
      # If there is no individual/pairwise-shared structure, set to 0
      if (r.vec[s] == 0) {
        
        # Saved scores and laodings
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
        W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        
        # Save beta_indiv if response is given
        if (response_given) {
          beta_indiv0[[s,1]] <- matrix(0, nrow = r.vec[s], ncol = 1)
        }
        
        # Fill in off-diagonal W with 0s
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      }
      
    }
    
    # Combine the coefficients
    beta0 <- rbind(beta_joint0[[1,1]], do.call(rbind, beta_indiv0))
    
    # Combining the scores together 
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) V0.star[[1,1]] <- V0[[1,1]] else V0.star[[1,1]] <- matrix(nrow = n, ncol = r)
    
    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) Vs0.star[[1,s]] <- Vs0[[1,s]] else Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }
    VStar0 <- cbind(do.call(cbind, V0.star), do.call(cbind, Vs0.star))
    
    # Initialize the latent variable for a binary outcome
    if (response_given) {
      Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
    } 
    
  }
  
  # If ranks provided, initialize with prior or use given starting values
  if (!nninit) {
    
    # If no starting values were provided, initialize from priors
    if (is.null(starting_values)) {
      V0 <- matrix(list(), nrow = 1, ncol = 1)
      if (r > 0) {
        V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
      } 
      if (r == 0) {
        V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      U0 <- matrix(list(), nrow = q, ncol = 1)
      Vs0 <- matrix(list(), nrow = 1, ncol = q)
      W0 <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        
        # Initialize U
        if (r > 0) {
          U0[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
        } 
        if (r == 0) {
          U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
        
        # Initialize W and V
        if (r.vec[s] > 0) {
          Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
          W0[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        } 
        if (r.vec[s] == 0) {
          Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
          W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
      }
    }
    
    # If starting values were provided, use them as initial values
    if (!is.null(starting_values)) {
      V0 <- starting_values$V
      U0 <- starting_values$U
      W0 <- starting_values$W
      Vs0 <- starting_values$Vs
    }
    
    # Initialize the regression parameters
    
    # Combining the scores together
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) {
      V0.star[[1,1]] <- V0[[1,1]] 
    }
    
    if (r == 0) {
      V0.star[[1,1]] <- matrix(nrow = n, ncol = r)
    } 
    
    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) {
        Vs0.star[[1,s]] <- Vs0[[1,s]] 
      }
      
      if (r.vec[s] == 0) {
        Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
      } 
    }
    
    if (response_given) {
      VStar0 <- cbind(do.call(cbind, V0.star), do.call(cbind, Vs0.star))
      
      if (n_beta != 0) {
        beta0 <- matrix(mvrnorm(1, mu = c(rep(0, n_beta)), Sigma = Sigma_beta))
      }
      
      if (n_beta == 0) {
        beta0 <- matrix(nrow = 0, ncol = 1)
      }
      
      # If beta0 is empty, this is generated from a N(0,1)
      Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
  
    }
    
  }
  
  # If imputing missingness in Y
  if (response_given) {
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(error_vars[q+1])))[missing_obs_Y,, drop = FALSE]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }
  
  # If there is missingness in the data, generate starting values for the missing entries
  if (missingness_in_data) {
    
    # If not initializing with posterior mode, initialize missing data with 0s
    if (!nninit) {
      Xm0 <- matrix(list(), ncol = 1, nrow = q)
      for (s in 1:q) {
        Xm0[[s,1]] <- rep(0, length(missing_obs[[s]]))
      }
    }
    
    # If initializing with posterior mode, initialize missing data with BIDIFAC+
    if (nninit) {
      Xm0 <- matrix(list(), ncol = 1, nrow = q)
      for (s in 1:q) {
        # Generate random noise
        Es0 <- rnorm(length(missing_obs[[s]]), mean = 0, sd = sqrt(error_vars[s]))
        
        # Save the imputed values
        Xm0[[s,1]] <- (U0[[s,1]] %*% t(V0[[1,1]]) + W0[[s,s]] %*% t(Vs0[[1,s]]))[missing_obs[[s]]] + Es0
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Storing the initial values 
  # ---------------------------------------------------------------------------
  
  V.draw[[1]] <- V0
  U.draw[[1]] <- U0
  Vs.draw[[1]] <- Vs0
  W.draw[[1]] <- W0
  
  if (response_given) {
    beta.draw[[1]][[1,1]] <- beta0
    Z.draw[[1]][[1,1]] <- Z0
    VStar.draw[[1]][[1,1]] <- VStar0
    
    if (missingness_in_response) {
      Ym.draw[[1]][[1,1]] <- Ym0
    }
    
  }
  
  if (missingness_in_data) {
    Xm.draw[[1]] <- Xm0
  }
  
  # ---------------------------------------------------------------------------
  # Computing the inverses 
  # ---------------------------------------------------------------------------
  
  if (!response_given) {
    # Error variance for X. 
    SigmaVInv <- diag(rep(1/error_vars, p.vec))
  }
  
  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined precisions between data and Z
      SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))
      
      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
      }
    } 
    
    if (response_type == "continuous") {
      # For V - Combined error variances between X1, X2, and Y
      SigmaVInv <- diag(c(rep(1/error_vars, c(p.vec,1))))
      
      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/error_vars[q+1]))
      }
    }
    
    # For beta - Combined precisions between intercept and all betas
    if (n_beta != 0) {
      SigmaBetaInv <- solve(Sigma_beta)
    }
  }
  
  # ---------------------------------------------------------------------------
  # If structure from another method is given, save the scores as VStar
  # ---------------------------------------------------------------------------
  
  if (!is.null(scores)) {
    VStar <- scores
  }
  
  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------
  
  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))
    
    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------
    
    V.iter <- V.draw[[iter]]
    U.iter <- U.draw[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.draw[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.draw[[iter]][[1,1]] 
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        beta_joint.iter <- beta.iter[r.ind[[1]],, drop = FALSE]
      } 
      
      if (r == 0) {
        beta_joint.iter <- matrix(0)
      }
      
      # Individual effects
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) {
          beta_indiv.iter[[s, 1]] <- matrix(0)
        }
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          beta_indiv.iter[[s,1]] <- beta.iter[r.ind[[s+1]],,drop=FALSE]
        }
      }
      
      # If binary outcome, save latent response
      if (response_type == "binary") {
        Z.iter <- Z.draw[[iter]][[1,1]]
      }
      
      # If there is missingness, fill in missing values with latest imputed values
      if (missingness_in_response) {
        # Save the current imputations for the missing values
        Ym.iter <- Ym.draw[[iter]][[1,1]]
        
        # Creating the completed outcome vector
        Y_complete <- scaled_Y
        
        # Filling in the missing entries for R1 and R2. 
        Y_complete[missing_obs_Y,] <- Ym.iter
      }
      
      # If there is no missingness in the response, store in Y_complete
      if (!missingness_in_response) {
        Y_complete <- scaled_Y
      }
      
    }
    
    if (missingness_in_data) {
      # Creating the completed matrices. 
      X_complete <- scaled_data
      
      # Fill in the completed matrices with the imputed values
      for (s in 1:q) {
        X_complete[[s,1]][missing_obs[[s]]] <- Xm.draw[[iter]][[s,1]]
      }
    }
    
    if (!missingness_in_data) {
      X_complete <- scaled_data
    }
    
    # -------------------------------------------------------------------------
    # Computing the inverse that changes with tau2
    # -------------------------------------------------------------------------
    
    # if (response_given) {
    #   if (response_type == "continuous") {
    #     # For V - Combined error variances between X1, X2, and Y
    #     SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1/tau2.iter[[1,1]]))
    # 
    #     # For Vs
    #     SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
    # 
    #     for (s in 1:q) {
    #       SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/tau2.iter[[1,1]]))
    #     }
    #   }
    # }
    
    # If estimating the underlying structure
    if (is.null(scores)) {
      # -------------------------------------------------------------------------
      # Posterior sample for V
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        if (!response_given) {
          # Concatenating Ui's together
          U.iter.combined <- do.call(rbind, U.iter)
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          # The combined centered Xis with the latent response vector
          X.iter <- do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <-  tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
        if (response_given) {
          # Concatenating Ui's together
          U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv) 
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          if (response_type == "binary") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Z.iter - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          if (response_type == "continuous") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Y_complete - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <- tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
      }
      
      if (r == 0) {
        V.draw[[iter+1]][[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      # Updating the value of V
      V.iter <- V.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Us
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        for (s in 1:q) {
          Xs.iter <- X_complete[[s,1]] - W.iter[[s,s]] %*% t(Vs.iter[[1,s]])
          Bu <- solve((1/error_vars[s]) * t(V.iter[[1,1]]) %*% V.iter[[1,1]] + (1/sigma2_joint) * diag(r))
          U.draw[[iter+1]][[s,1]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bu <- (1/error_vars[s]) * t(V.iter[[1,1]]) %*% Xs.iter[j, ]
            
            U1j <- mvrnorm(1, mu = Bu %*% bu, Sigma = Bu)
            U1j
          }), nrow = r))
        }
      }
      
      if (r == 0) {
        for (s in 1:q) {
          U.draw[[iter+1]][[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
      }
      
      U.iter <- U.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Vs, s=1,...,q
      # -------------------------------------------------------------------------
      
      if (!response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
            Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      if (response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            # Combined Ws and beta
            W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))
            
            tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
            tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)
            
            Bvs <- solve(tW_Sigma_W + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            if (response_type == "binary") {
              # Combined centered Xs and Z
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Z.iter - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            if (response_type == "continuous") {
              # Combined centered Xs and Y
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Y_complete - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- tW_Sigma %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      # Update the current value of V
      Vs.iter <- Vs.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Combine the scores
      # -------------------------------------------------------------------------
      
      # Combine current values of V and V.
      V.iter.star.joint <- V.iter
      if (r == 0) {
        V.iter.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
      } 
      
      # Combine the scores together for response prediction
      Vs.iter.star <- Vs.iter
      for (s in 1:q) {
        if (r.vec[s] == 0) {
          Vs.iter.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
        } 
      }
      
      VStar.iter <- cbind(do.call(cbind, V.iter.star.joint), do.call(cbind, Vs.iter.star))
      
      # Save the current VStar
      VStar.draw[[iter+1]][[1,1]] <- VStar.iter
      
      # -------------------------------------------------------------------------
      # Posterior sample for W
      # -------------------------------------------------------------------------
      
      for (s in 1:q) {
        if (r.vec[s] > 0) {
          Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
          Bws <- solve((1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Vs.iter[[1,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
          
          W.draw[[iter+1]][[s,s]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bws <- (1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Xs.iter[j,] 
            
            Wsj <- mvrnorm(1, mu = Bws %*% bws, Sigma = Bws)
            Wsj
          }), nrow = r.vec[s]))
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
        if (r.vec[s] == 0) {
          W.draw[[iter+1]][[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
      }
      
      # Update the current value of W
      W.iter <- W.draw[[iter+1]]
      
    }
    
    # If structure from another method is provided
    if (!is.null(scores)) {
      VStar.iter <- VStar
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------
    
    if (response_given) {
      
      # If the rank of the whole factorization is > 0
      if (n_beta != 0) {
        
        if (response_type == "binary") {
          Bbeta <- solve(t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
          bbeta <- t(VStar.iter) %*% Z.iter
        }
        
        if (response_type == "continuous") {
          Bbeta <- solve((1/error_vars[q+1]) * t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
          bbeta <- (1/error_vars[q+1]) * t(VStar.iter) %*% Y_complete
        }
        
        beta.draw[[iter+1]][[1,1]] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      }
      
      # If the rank of the factorization is 0
      if (n_beta == 0) {
        beta.draw[[iter+1]][[1,1]] <- matrix(nrow = 0, ncol = 1)
      }
      
      # Update the current value of beta
      beta.iter <- beta.draw[[iter+1]][[1,1]]
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Joint effect
      if (r != 0) {
        beta_joint.iter <- beta.iter[r.ind[[1]],, drop = FALSE]
      }
      
      if (r == 0) {
        beta_joint.iter <- matrix(0)
      }
      
      # Individual effects
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          beta_indiv.iter[[s, 1]] <- beta.iter[r.ind[[s+1]],,drop = FALSE]
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for tau2
    # -------------------------------------------------------------------------
    
    # if (response_given & !is.null(scores)) {
    #   if (response_type == "continuous") {
    #     tau2.draw[[iter+1]][[1,1]] <- matrix(1/rgamma(1, shape = shape + (n/2), rate = rate + 0.5 * sum((Y_complete - VStar.iter %*% beta.iter)^2)))
    # 
    #     # Update the current value of tau2
    #     tau2.iter <- tau2.draw[[iter+1]][[1,1]]
    #   }
    # }
    
    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Z.draw[[iter+1]][[1,1]] <- matrix(sapply(1:n, function(i) {
          if (Y_complete[i,] == 1) {
            rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          } else {
            rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          }
        }), ncol = 1)
      }
    }
    
    # -------------------------------------------------------------------------
    # Impute missing data
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (missingness_in_response) {
        if (response_type == "continuous") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(error_vars[q+1])), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
        
        if (response_type == "binary") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
      }
    }
    
    if (missingness_in_data) {
      for (s in 1:q) {
        Es <-  matrix(rnorm(p.vec[s]*n, 0, sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
        Xm.draw[[iter+1]][[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the joint and individual/pairwise-shared structure, scaled to 
  # the original data. 
  # ---------------------------------------------------------------------------
  
  # Storing the joint and individual/pairwise-shared structure at each Gibbs sampling iteration
  J.draw <- A.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  # Storing the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  
  if (is.null(scores)) {
    for (iter in 1:nsample) {
      for (s in 1:q) {
        # Calculating the joint structure and scaling by sigma.mat
        J.draw[[iter]][[s,1]] <- (U.draw[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]])) * sigma.mat[s,1]
        
        # Calculating the individual structure and scaling by sigma.mat
        A.draw[[iter]][[s,1]] <- (W.draw[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) * sigma.mat[s,1]
      }
      
      # Calculate the structure for Y
      if (response_given) {
        EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]] * sigma.mat[q+1,1]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  
  list(scaled_data = scaled_data, # Returning the scaled version of the data
       scaled_Y = scaled_Y, # Return the scaled version of the response vector
       sigma.mat = sigma.mat, # Scaling factors
       J.draw = J.draw, A.draw = A.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
       VStar.draw = VStar.draw, # Components that predict Y
       Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
       scores = scores, # Scores if provided by another method 
       ranks = c(r, r.vec), # Ranks
       beta.draw = beta.draw # Regression parameters
  )
  
}

# This version initializes with BIDIFAC WITHOUT y as a source
bpmf_data_mode <- function(data, Y, nninit = TRUE, model_params, ranks = NULL, scores = NULL, sparsity = FALSE, nsample, progress = TRUE, starting_values = NULL) {
  # Gibbs sampling algorithm for sampling the underlying structure and the 
  # regression coefficient vector for a response vector. 
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # 
  # data = matrix with list entries corresponding to each data source
  # Y = column vector with outcome or NULL
  # nninit = should the model be initialized with a nuclear norm penalized objective? if FALSE, provide ranks
  # model_params = (error_vars, joint_vars, indiv_vars, beta_vars = NULL, response_vars)
  # ranks = joint and individual ranks (1st entry = joint rank, kth for k>1 is the individual rank for the k-1'st source)
  # scores = if using structure from another method to fit a linear model, provide joint and individual scores here.
  #   in this case, ranks should be provided and nninit = FALSE, Y != NULL
  # nsample = number of Gibbs sampling iterations
  # progress = should the progress of the sampler be displayed?
  # starting_values = list of starting values for V, U, W, Vs. If NULL and nninit = TRUE, init with BIDIFAC+,
  #    if NULL and nninit = FALSE, init from prior. If not NULL, will init with provided starting values unless 
  #    nninit. 
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  q <- nrow(data) # Number of sources
  p.vec <- apply(data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(data[[1,1]]) # Number of subjects
  
  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------
  
  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
  
  # ---------------------------------------------------------------------------
  # Check for missingness in data
  # ---------------------------------------------------------------------------
  
  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
  # If there is missingness, initialize the missing values with 0s
  for (s in 1:q) {
    data[[s,1]][missing_obs[[s]]] <- 0
  }
  
  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------

  response_given <- !is.null(Y[[1,1]]) 
  
  # If so, what kind of response is it?
  if (response_given) {
    Y <- matrix(unlist(Y))
    
    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"
    
    # If there is a response, is there missingness in the outcome?
    missingness_in_response <- any(is.na(Y))
    
    # Which entries are missing?
    missing_obs_Y <- which(is.na(Y))
  }
  
  # ---------------------------------------------------------------------------
  # Obtaining the ranks 
  # ---------------------------------------------------------------------------
  
  if (nninit) {
    rank_init <- BIDIFAC(data, rmt = TRUE, pbar = FALSE, scale_back = FALSE)
    
    # Print when finished
    print("Posterior mode obtained: joint and individual ranks determined.")
    
    # Saving the results
    sigma.mat <- rank_init$sigma.mat
    C <- rank_init$C
    r <- rankMatrix(C[[1,1]])[[1]] # Joint rank
    I <- rank_init$I
    r.vec <- sapply(1:q, function(s) rankMatrix(I[[s,1]])) # Individual ranks
    
    # Scaling the data
    for (s in 1:q) {
      data[[s,1]] <- data[[s,1]]/sigma.mat[s,]
    }
  }
  
  if (!nninit) {
    r <- ranks[1]
    r.vec <- ranks[-1]
    sigma.mat <- matrix(nrow = q, 1)
    for (s in 1:q) {
      sigma.mat[s,1] <- 1
    }
  }
  
  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total
  
  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    beta_vars <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    diag(Sigma_beta) <- beta_vars
  }
  
  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  U.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  W.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = q))
  
  if (!response_given) {
    beta.draw <- Z.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (!missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  if (!sparsity) {
    gamma.draw <- p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (response_given) {
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    if (sparsity) {
      gamma.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
      p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
    }
  }
  
  if (missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  # ---------------------------------------------------------------------------
  # Initialize V, U, V, W
  # ---------------------------------------------------------------------------
  
  # If initializing with nuclear norm, set initialize values to posterior mode
  if (nninit) {

    V0 <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) {
      svd.joint <- svd(rank_init$C[[1,1]])
      V0[[1,1]] <- (svd.joint$v[,1:r, drop = FALSE]) %*% diag(svd.joint$d[1:r], nrow = r)
    } 
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
    }
    
    U0 <- matrix(list(), nrow = q, ncol = 1)
    Vs0 <- matrix(list(), nrow = 1, ncol = q)
    W0 <- matrix(list(), nrow = q, ncol = q)
    
    for (s in 1:q) {
      
      # Initialize U
      if (r > 0) {
        U0[[s,1]] <- svd(rank_init$C[[s,1]])$u[,1:r, drop = FALSE]
      } 
      if (r == 0) {
        U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
      }
      
      # Initialize W and V
      if (r.vec[s] > 0) {
        svd.indiv.s <- svd(rank_init$I[[s,1]])
        Vs0[[1,s]] <- (svd.indiv.s$v[,1:r.vec[s], drop = FALSE]) %*% diag(svd.indiv.s$d[1:r.vec[s]], nrow = r.vec[s])
        W0[[s,s]] <- svd.indiv.s$u[,1:r.vec[s], drop = FALSE]
        
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      } 
      if (r.vec[s] == 0) {
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
        W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      }
      
    }
  }
  
  # If ranks provided, initialize with prior or use given starting values
  if (!nninit) {
    
    # If no starting values were provided, initialize from priors
    if (is.null(starting_values)) {
      V0 <- matrix(list(), nrow = 1, ncol = 1)
      if (r > 0) {
        V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
      } 
      if (r == 0) {
        V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      U0 <- matrix(list(), nrow = q, ncol = 1)
      Vs0 <- matrix(list(), nrow = 1, ncol = q)
      W0 <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        
        # Initialize U
        if (r > 0) {
          U0[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
        } 
        if (r == 0) {
          U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
        
        # Initialize W and V
        if (r.vec[s] > 0) {
          Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
          W0[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        } 
        if (r.vec[s] == 0) {
          Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
          W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
      }
    }
    
    # If starting values were provided, use them as initial values
    if (!is.null(starting_values)) {
      V0 <- starting_values$V
      U0 <- starting_values$U
      W0 <- starting_values$W
      Vs0 <- starting_values$Vs
    }

  }
  
  if (response_given) {
    # Combining the scores together
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) V0.star[[1,1]] <- V0[[1,1]] else V0.star[[1,1]] <- matrix(nrow = n, ncol = r)
    
    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) Vs0.star[[1,s]] <- Vs0[[1,s]] else Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }
    
    VStar0 <- cbind(1, do.call(cbind, V0.star), do.call(cbind, Vs0.star))
    
    beta0 <- matrix(mvrnorm(1, mu = c(rep(0, n_beta)), Sigma = Sigma_beta))
    Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
    tau20 <- matrix(1/rgamma(1, shape = shape, rate = rate))
    
    if (sparsity) {
      p0 <- 0.5
      gamma0 <- matrix(rbinom(n_beta, size = 1, prob = p0), ncol = 1)
    }
  }
  
  # If there is missingness in the data, generate starting values for the missing entries
  if (missingness_in_data) {
    Xm0 <- matrix(list(), ncol = 1, nrow = q)
    for (s in 1:q) {
      Xm0[[s,1]] <- rep(0, length(missing_obs[s]))
    }
  }
  
  # If there is missingness in Y, generate starting values for the missing entries
  if (response_given) {
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20)))[missing_obs_Y,, drop = FALSE]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Storing the initial values 
  # ---------------------------------------------------------------------------
  
  V.draw[[1]] <- V0
  U.draw[[1]] <- U0
  Vs.draw[[1]] <- Vs0
  W.draw[[1]] <- W0
  
  if (response_given) {
    beta.draw[[1]][[1,1]] <- beta0
    Z.draw[[1]][[1,1]] <- Z0
    tau2.draw[[1]][[1,1]] <- tau20
    VStar.draw[[1]][[1,1]] <- VStar0

    if (missingness_in_response) {
      Ym.draw[[1]][[1,1]] <- Ym0
    }
    
    if (sparsity) {
      gamma.draw[[1]][[1,1]] <- gamma0
      p.draw[[1]][[1,1]] <- p0
    }
  }
  
  if (missingness_in_data) {
    Xm.draw[[1]] <- Xm0
  }
  
  # ---------------------------------------------------------------------------
  # Computing the inverses that don't change from iteration to iteration
  # ---------------------------------------------------------------------------
  
  if (!response_given) {
    # Error variance for X. 
    SigmaVInv <- diag(rep(1/error_vars, p.vec))
  }
  
  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined precisions between data and Z
      SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))
      
      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
      }
    } 
    
    # For beta - Combined precisions between intercept and all betas
    SigmaBetaInv <- solve(Sigma_beta)
  }
  
  # ---------------------------------------------------------------------------
  # If structure from another method is given, save the scores as VStar
  # ---------------------------------------------------------------------------
  
  if (!is.null(scores)) {
    VStar <- cbind(1, scores)
  }
  
  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------
  
  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))
    
    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------
    
    V.iter <- V.draw[[iter]]
    U.iter <- U.draw[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.draw[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.draw[[iter]][[1,1]] 
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), nrow = q, ncol = 1)
      
      # Breaking beta down into the intercept,
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE] 
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
      
      if (response_type == "binary") {
        Z.iter <- Z.draw[[iter]][[1,1]]
      }
      
      if (response_type == "continuous") {
        tau2.iter <- tau2.draw[[iter]][[1,1]]
      }
      
      if (missingness_in_response) {
        # Save the current imputations for the missing values
        Ym.iter <- Ym.draw[[iter]][[1,1]]
        
        # Creating the completed outcome vector
        Y_complete <- Y
        
        # Filling in the missing entries for R1 and R2. 
        Y_complete[missing_obs_Y,] <- Ym.iter
      }
      
      if (!missingness_in_response) {
        Y_complete <- Y
      }
      
      if (sparsity) {
        gamma.iter <- gamma.draw[[iter]][[1,1]]
        p.iter <- p.draw[[iter]][[1,1]]
      }
    }
    
    if (missingness_in_data) {
      # Creating the completed matrices. 
      X_complete <- data
      
      # Fill in the completed matrices with the imputed values
      for (s in 1:q) {
        X_complete[[s,1]][missing_obs[[s]]] <- Xm.draw[[iter]][[s,1]]
      }
    }
    
    if (!missingness_in_data) {
      X_complete <- data
    }
    
    # -------------------------------------------------------------------------
    # Computing the inverse that changes with tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        # For V - Combined error variances between X1, X2, and Y
        SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1/tau2.iter[[1,1]]))
        
        # For Vs
        SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
        
        for (s in 1:q) {
          SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/tau2.iter[[1,1]]))
        }
      }
    }
    
    # If estimating the underlying structure
    if (is.null(scores)) {
      # -------------------------------------------------------------------------
      # Posterior sample for V
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        if (!response_given) {
          # Concatenating Ui's together
          U.iter.combined <- do.call(rbind, U.iter)
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          # The combined centered Xis with the latent response vector
          X.iter <- do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <-  tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
        if (response_given) {
          # Concatenating Ui's together
          U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv) 
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          if (response_type == "binary") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          if (response_type == "continuous") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <- tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
      }
      
      if (r == 0) {
        V.draw[[iter+1]][[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      # Updating the value of V
      V.iter <- V.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Us
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        for (s in 1:q) {
          Xs.iter <- X_complete[[s,1]] - W.iter[[s,s]] %*% t(Vs.iter[[1,s]])
          Bu <- solve((1/error_vars[s]) * t(V.iter[[1,1]]) %*% V.iter[[1,1]] + (1/sigma2_joint) * diag(r))
          U.draw[[iter+1]][[s,1]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bu <- (1/error_vars[s]) * t(V.iter[[1,1]]) %*% Xs.iter[j, ]
            
            U1j <- mvrnorm(1, mu = Bu %*% bu, Sigma = Bu)
            U1j
          }), nrow = r))
        }
      }
      
      if (r == 0) {
        for (s in 1:q) {
          U.draw[[iter+1]][[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
      }
  
      U.iter <- U.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Vs, s=1,...,q
      # -------------------------------------------------------------------------
      
      if (!response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
            Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      if (response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            # Combined Ws and beta
            W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))
            
            tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
            tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)
            
            Bvs <- solve(tW_Sigma_W + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            if (response_type == "binary") {
              # Combined centered Xs and Z
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            if (response_type == "continuous") {
              # Combined centered Xs and Y
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Y_complete - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- tW_Sigma %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      # Update the current value of V
      Vs.iter <- Vs.draw[[iter+1]]
      
      # Combine current values of V and V.
      V.iter.star.joint <- V.iter
      if (r == 0) {
        V.iter.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
      } 
      
      Vs.iter.star <- Vs.iter
      for (s in 1:q) {
        if (r.vec[s] == 0) {
          Vs.iter.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
        } 
      }
      
      VStar.iter <- cbind(1, do.call(cbind, V.iter.star.joint), do.call(cbind, Vs.iter.star))
      
      # Save the current VStar
      VStar.draw[[iter+1]][[1,1]] <- VStar.iter
      
      # -------------------------------------------------------------------------
      # Posterior sample for W
      # -------------------------------------------------------------------------
      
      for (s in 1:q) {
        if (r.vec[s] > 0) {
          Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
          Bws <- solve((1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Vs.iter[[1,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
          
          W.draw[[iter+1]][[s,s]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bws <- (1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Xs.iter[j,] 
            
            Wsj <- mvrnorm(1, mu = Bws %*% bws, Sigma = Bws)
            Wsj
          }), nrow = r.vec[s]))
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
        if (r.vec[s] == 0) {
          W.draw[[iter+1]][[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
      }
      
      # Update the current value of W
      W.iter <- W.draw[[iter+1]]
    
    }
    
    # If structure from another method is provided
    if (!is.null(scores)) {
      VStar.iter <- VStar
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------
    
    if (response_given) {
      
      if (sparsity) {
        # Change the precision for the intercept
        diag(SigmaBetaInv)[1] <- 1/beta_vars[1]
        
        # Change the precision of those betas under the slab, excluding the intercept
        diag(SigmaBetaInv)[-1][gamma.iter[-1] == 1] <- 1/(beta_vars[-1][gamma.iter[-1] == 1])
        
        # Change the precision of those betas under the spike
        diag(SigmaBetaInv)[-1][gamma.iter[-1] == 0] <- 1000 # Can change this precision
      }
      
      if (response_type == "binary") {
        Bbeta <- solve(t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- t(VStar.iter) %*% Z.iter
      }
      
      if (response_type == "continuous") {
        Bbeta <- solve((1/tau2.iter[[1,1]]) * t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- (1/tau2.iter[[1,1]]) * t(VStar.iter) %*% Y_complete
      }
      
      beta.draw[[iter+1]][[1,1]] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      
      # Update the current value of beta
      beta.iter <- beta.draw[[iter+1]][[1,1]]
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE] 
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for spike-and-slab indicators and probabilities
    # -------------------------------------------------------------------------
    
    if (response_given & sparsity) {
      # Calculate the likelihood for beta under the spike and slab
      like_spike <- dnorm(beta.iter, mean = 0, sd = sqrt(1/1000), log = TRUE)
      like_slab <- dnorm(beta.iter, mean = 0, sd = sqrt(beta_vars), log = TRUE)
      
      # Calculating the probability that each gamma equals 1
      prob = sapply(1:n_beta, function(rs) {
        x = log(p.iter) + like_slab[rs,]
        y = log(1 - p.iter) + like_spike[rs,]
        exp(x - logSum(c(x,y)))
      })
      prob[1] <- 1 # Always include the intercept
      
      # Generating the gammas
      gamma.draw[[iter+1]][[1,1]] <- matrix(rbinom(n_beta, size = 1, prob = prob), ncol = 1)
      
      # Saving the current value of gamma
      gamma.iter <- gamma.draw[[iter+1]][[1,1]]
      
      # Generating the prior probabilities (excluding the intercept, hence n_beta - 1 and sum(gamma[-1]))
      p.draw[[iter+1]][[1,1]] <- rbeta(1, 1 + sum(gamma.iter[-1,]), 1 + (n_beta-1) - sum(gamma.iter[-1,]))
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        tau2.draw[[iter+1]][[1,1]] <- matrix(1/rgamma(1, shape = shape + (n/2), rate = rate + 0.5 * sum((Y_complete - VStar.iter %*% beta.iter)^2)))
        
        # Update the current value of tau2
        tau2.iter <- tau2.draw[[iter+1]][[1,1]]
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Z.draw[[iter+1]][[1,1]] <- matrix(sapply(1:n, function(i) {
          if (Y_complete[i,] == 1) {
            rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          } else {
            rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          }
        }), ncol = 1)
      }
    }
    
    # -------------------------------------------------------------------------
    # Impute missing data
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (missingness_in_response) {
        if (response_type == "continuous") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter[[1,1]])), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
        
        if (response_type == "binary") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
      }
    }
    
    if (missingness_in_data) {
      for (s in 1:q) {
        Es <-  matrix(rnorm(p.vec[s]*n, 0, sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
        Xm.draw[[iter+1]][[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the joint and individual structure, scaled to the data
  # ---------------------------------------------------------------------------
  
  # Storing the structures at each Gibbs sampling iteration
  J.draw <- A.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  # Calculating the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  if (is.null(scores)) {
    for (iter in 1:nsample) {
      for (s in 1:q) {
        # Calculating the joint structure and scaling by sigma.mat
        J.draw[[iter]][[s,1]] <- (U.draw[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]])) * sigma.mat[s,1]
        
        # Calculating the individual structure and scaling by sigma.mat
        A.draw[[iter]][[s,1]] <- (W.draw[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) * sigma.mat[s,1]
      }
      
      # Calculate the structure for Y
      if (response_given) {
        EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]]
      }
    }
  }
  
  # Return
  list(data = data, # Returning the scaled version of the data
        Y = Y, # Return the response vector
        sigma.mat = sigma.mat, # Scaling factors
        J.draw = J.draw, A.draw = A.draw, EY.draw = EY.draw, # Underlying structure
        V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
        VStar.draw = VStar.draw, # Components that predict Y,
        Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
        scores = scores, # Scores if provided by another method 
        ranks = c(r, r.vec), # Ranks
        tau2.draw = tau2.draw, beta.draw = beta.draw, # Regression parameters
        gamma.draw = gamma.draw, p.draw = p.draw) # Sparsity parameters

}

# Troubleshooting the initialization with Y as a source
bpmf_test_v1 <- function(data, Y, nninit = TRUE, model_params, ranks = NULL, scores = NULL, sparsity = FALSE, nsample, progress = TRUE, starting_values = NULL) {
  # Gibbs sampling algorithm for sampling the underlying structure and the 
  # regression coefficient vector for a response vector. 
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # 
  # data = matrix with list entries corresponding to each data source
  # Y = column vector with outcome or NULL
  # nninit = should the model be initialized with a nuclear norm penalized objective? if FALSE, provide ranks
  # model_params = (error_vars, joint_vars, indiv_vars, beta_vars = NULL, response_vars)
  # ranks = joint and individual ranks (1st entry = joint rank, kth for k>1 is the individual rank for the k-1'st source)
  # scores = if using structure from another method to fit a linear model, provide joint and individual scores here.
  #   in this case, ranks should be provided and nninit = FALSE, Y != NULL
  # nsample = number of Gibbs sampling iterations
  # progress = should the progress of the sampler be displayed?
  # starting_values = list of starting values for V, U, W, Vs. If NULL and nninit = TRUE, init with BIDIFAC+,
  #    if NULL and nninit = FALSE, init from prior. If not NULL, will init with provided starting values unless 
  #    nninit. 
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  q <- nrow(data) # Number of sources
  p.vec <- apply(data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(data[[1,1]]) # Number of subjects
  
  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------
  
  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
  
  # ---------------------------------------------------------------------------
  # Check for missingness in data
  # ---------------------------------------------------------------------------
  
  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
  # If there is missingness, initialize the missing values with 0s
  for (s in 1:q) {
    data[[s,1]][missing_obs[[s]]] <- 0
  }
  
  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------
  
  response_given <- !is.null(Y[[1,1]]) 
  
  # If so, what kind of response is it?
  if (response_given) {
    Y <- matrix(unlist(Y))
    
    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"
    
    # If there is a response, is there missingness in the outcome?
    missingness_in_response <- any(is.na(Y))
    
    # Which entries are missing?
    missing_obs_Y <- which(is.na(Y))
  }
  
  # ---------------------------------------------------------------------------
  # Obtaining the ranks 
  # ---------------------------------------------------------------------------
  
  scaled_data <- data
  scaled_Y <- Y
  
  # Initializing with BIDIFAC+
  if (nninit) {
    
    # Combine the data together 
    data_combined <- matrix(list(), nrow = (q+1), ncol = 1)
    
    for (s in 1:q) {
      # Append the scaled data
      data_combined[[s,1]] <- scaled_data[[s,1]]
      
      # Add back in missing data if any
      data_combined[[s,1]][missing_obs[[1]]] <- NA
    }
    
    # Initialize the indices of observations in each source 
    p.ind <- lapply(1:q, function(s) {
      if (s == 1) {
        1:p.vec[s]
      } else {
        (p.vec[s-1] + 1):cumsum(p.vec)[s]
      }
    })
    
    # Save indices of samples per source
    n.ind <- list(1:n)
    
    # If no response is given, p.ind.list should only include the sources
    if (!response_given) {
      
      # Save the indices for features for each identified structure
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- p.ind[[s]]
      }
    }
    
    # Include the response as a source if a response is given and add it to p.ind, p.ind.list
    if (response_given) {
      
      # Append the response to the sources
      data_combined[[q+1,1]] <- t(scaled_Y)
      
      # Add back in any missingness if exists
      data_combined[[q+1,1]][missing_obs_Y] <- NA
      
      # Add the response indices to p.ind
      p.ind[[q+1]] <- cumsum(p.vec)[q] + 1 # For the response
      
      # Save the indices for features for each identified structure
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- c(p.ind[[s]], p.ind[[q+1]])
      }
      
    }
    
    # Collapse the data into one matrix
    data_combined <- do.call(rbind, data_combined)
    
    # Save the indices for samples in each identified structure
    n.ind.list <- lapply(1:(q+1), function(s) c(1:n))
    
    # Run BIDIFAC+ --
    
    # If there is no missing data
    if (!missingness_in_data & !missingness_in_response) {
      rank_init <- bidifac.plus.given(data_combined, p.ind = p.ind, n.ind = n.ind,
                                      p.ind.list = p.ind.list, n.ind.list = n.ind.list)
    }
    
    # If there is missing data
    if (missingness_in_data | missingness_in_response) {
      # Save the indices of the missing values
      all.miss <- which(is.na(data_combined))
      
      # Initialize
      rank_init <- bidifac.plus.impute(data_combined, p.ind = p.ind, n.ind = n.ind,
                                       p.ind.list = p.ind.list, n.ind.list = n.ind.list,
                                       all.miss = all.miss)
    }
    
    # Print when finished
    print("Posterior mode obtained, ranks determined.")
    
    # Obtain the structures and their ranks
    S <- rank_init$S
    
    # Joint rank
    r <- Matrix::rankMatrix(S[[1]])[[1]]
    
    # Individual ranks if no response or pairwise-shared ranks of each source with Y if response given
    r.vec <- sapply(2:(q+1), function(s) Matrix::rankMatrix(S[[s]])[[1]])
    
  }
  
  if (!nninit) {
    r <- ranks[1]
    r.vec <- ranks[-1]
    sigma.mat <- matrix(nrow = q, 1)
    for (s in 1:q) {
      sigma.mat[s,1] <- 1
    }
  }
  
  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total
  
  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    beta_vars <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    diag(Sigma_beta) <- beta_vars
  }
  
  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  U.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  W.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = q))
  
  if (!response_given) {
    beta.draw <- Z.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (!missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  if (!sparsity) {
    gamma.draw <- p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (response_given) {
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    if (sparsity) {
      gamma.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
      p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
    }
  }
  
  if (missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  # ---------------------------------------------------------------------------
  # Initialize V, U, V, W
  # ---------------------------------------------------------------------------
  
  # If initializing with nuclear norm, initialize sampling at posterior mode
  if (nninit) {
    
    # Initialize joint scores, V
    
    V0 <- matrix(list(), nrow = 1, ncol = 1)
    beta_joint0 <- matrix(list(), nrow = 1, ncol = 1)
    
    # If there is joint structure
    if (r > 0) {
      svd.joint <- svd(S[[1]])
      V0[[1,1]] <- (svd.joint$v[,1:r, drop = FALSE]) %*% diag(svd.joint$d[1:r], nrow = r)
      
      # Save beta_joint if response is given
      if (response_given) {
        beta_joint0[[1,1]] <- t(svd(S[[1]])$u[p.ind[[q+1]],1:r, drop = FALSE])
      }
    } 
    
    # If there is no joint structure
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      
      # Save beta joint if response is given
      if (response_given) {
        beta_joint0[[1,1]] <- matrix(0, nrow = r, ncol = 1)
      }
    }
    
    U0 <- matrix(list(), nrow = q, ncol = 1)
    Vs0 <- matrix(list(), nrow = 1, ncol = q)
    W0 <- matrix(list(), nrow = q, ncol = q)
    beta_indiv0 <- matrix(list(), nrow = q, ncol = 1)
    
    for (s in 1:q) {
      
      # Initialize joint loadings, U 
      if (r > 0) {
        U0[[s,1]] <- svd(S[[1]])$u[p.ind[[s]],1:r, drop = FALSE]
      } 
      if (r == 0) {
        U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
      }
      
      # Initialize individual loadings, W, and individual scores, Vs, and beta_indiv if response is given
      if (r.vec[s] > 0) {
        
        # Compute SVD
        svd.indiv.s <- svd(S[[s+1]])
        
        # Save scores and loadings
        Vs0[[1,s]] <- (svd.indiv.s$v[,1:r.vec[s], drop = FALSE]) %*% diag(svd.indiv.s$d[1:r.vec[s]], nrow = r.vec[s])
        W0[[s,s]] <- svd.indiv.s$u[p.ind[[s]],1:r.vec[s], drop = FALSE]
        
        # Save beta_indiv
        if (response_given) {
          beta_indiv0[[s,1]] <- t(svd.indiv.s$u[p.ind[[q+1]],1:r.vec[s], drop = FALSE])
        }
        
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      } 
      
      # If there is no individual/pairwise-shared structure, set to 0
      if (r.vec[s] == 0) {
        
        # Saved scores and laodings
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
        W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        
        # Save beta_indiv if response is given
        if (response_given) {
          beta_indiv0[[s,1]] <- matrix(0, nrow = r.vec[s], ncol = 1)
        }
        
        # Fill in off-diagonal W with 0s
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      }
      
    }
    
    # Combine the coefficients
    beta0 <- rbind(rnorm(1, 0, Sigma_beta[1,1]), beta_joint0[[1,1]], do.call(rbind, beta_indiv0))
    
    # Combining the scores together 
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) V0.star[[1,1]] <- V0[[1,1]] else V0.star[[1,1]] <- matrix(nrow = n, ncol = r)
    
    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) Vs0.star[[1,s]] <- Vs0[[1,s]] else Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }
    VStar0 <- cbind(1, do.call(cbind, V0.star), do.call(cbind, Vs0.star))
    
    # Initialize the latent variable for a binary outcome
    if (response_given) {
      Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
    } 
    
  }
  
  # If ranks provided, initialize with prior or use given starting values
  if (!nninit) {
    
    # If no starting values were provided, initialize from priors
    if (is.null(starting_values)) {
      V0 <- matrix(list(), nrow = 1, ncol = 1)
      if (r > 0) {
        V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
      } 
      if (r == 0) {
        V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      U0 <- matrix(list(), nrow = q, ncol = 1)
      Vs0 <- matrix(list(), nrow = 1, ncol = q)
      W0 <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        
        # Initialize U
        if (r > 0) {
          U0[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
        } 
        if (r == 0) {
          U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
        
        # Initialize W and V
        if (r.vec[s] > 0) {
          Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
          W0[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        } 
        if (r.vec[s] == 0) {
          Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
          W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
      }
    }
    
    # If starting values were provided, use them as initial values
    if (!is.null(starting_values)) {
      V0 <- starting_values$V
      U0 <- starting_values$U
      W0 <- starting_values$W
      Vs0 <- starting_values$Vs
    }
    
  }
  
  # If imputing missingness in Y
  if (response_given) {
    tau20 <- matrix(1/rgamma(1, shape = shape, rate = rate))  
  
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20[[1,1]])))[missing_obs_Y,, drop = FALSE]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }
  
  # If there is missingness in the data, generate starting values for the missing entries
  if (missingness_in_data) {
    Xm0 <- matrix(list(), ncol = 1, nrow = q)
    for (s in 1:q) {
      Xm0[[s,1]] <- rep(0, length(missing_obs[s]))
    }
  }
  
  # If there is missingness in Y, generate starting values for the missing entries
  if (response_given) {
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20)))[missing_obs_Y,, drop = FALSE]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Storing the initial values 
  # ---------------------------------------------------------------------------
  
  V.draw[[1]] <- V0
  U.draw[[1]] <- U0
  Vs.draw[[1]] <- Vs0
  W.draw[[1]] <- W0
  
  if (response_given) {
    beta.draw[[1]][[1,1]] <- beta0
    Z.draw[[1]][[1,1]] <- Z0
    tau2.draw[[1]][[1,1]] <- tau20
    VStar.draw[[1]][[1,1]] <- VStar0
    
    if (missingness_in_response) {
      Ym.draw[[1]][[1,1]] <- Ym0
    }
    
    if (sparsity) {
      gamma.draw[[1]][[1,1]] <- gamma0
      p.draw[[1]][[1,1]] <- p0
    }
  }
  
  if (missingness_in_data) {
    Xm.draw[[1]] <- Xm0
  }
  
  # ---------------------------------------------------------------------------
  # Computing the inverses that don't change from iteration to iteration
  # ---------------------------------------------------------------------------
  
  if (!response_given) {
    # Error variance for X. 
    SigmaVInv <- diag(rep(1/error_vars, p.vec))
  }
  
  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined precisions between data and Z
      SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))
      
      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
      }
    } 
    
    # For beta - Combined precisions between intercept and all betas
    SigmaBetaInv <- solve(Sigma_beta)
  }
  
  # ---------------------------------------------------------------------------
  # If structure from another method is given, save the scores as VStar
  # ---------------------------------------------------------------------------
  
  if (!is.null(scores)) {
    VStar <- cbind(1, scores)
  }
  
  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------
  
  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))
    
    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------
    
    V.iter <- V.draw[[iter]]
    U.iter <- U.draw[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.draw[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.draw[[iter]][[1,1]] 
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), nrow = q, ncol = 1)
      
      # Breaking beta down into the intercept,
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE] 
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
      
      if (response_type == "binary") {
        Z.iter <- Z.draw[[iter]][[1,1]]
      }
      
      if (response_type == "continuous") {
        tau2.iter <- tau2.draw[[iter]][[1,1]]
      }
      
      if (missingness_in_response) {
        # Save the current imputations for the missing values
        Ym.iter <- Ym.draw[[iter]][[1,1]]
        
        # Creating the completed outcome vector
        Y_complete <- Y
        
        # Filling in the missing entries for R1 and R2. 
        Y_complete[missing_obs_Y,] <- Ym.iter
      }
      
      if (!missingness_in_response) {
        Y_complete <- Y
      }
      
      if (sparsity) {
        gamma.iter <- gamma.draw[[iter]][[1,1]]
        p.iter <- p.draw[[iter]][[1,1]]
      }
    }
    
    if (missingness_in_data) {
      # Creating the completed matrices. 
      X_complete <- data
      
      # Fill in the completed matrices with the imputed values
      for (s in 1:q) {
        X_complete[[s,1]][missing_obs[[s]]] <- Xm.draw[[iter]][[s,1]]
      }
    }
    
    if (!missingness_in_data) {
      X_complete <- data
    }
    
    # -------------------------------------------------------------------------
    # Computing the inverse that changes with tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        # For V - Combined error variances between X1, X2, and Y
        SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1/tau2.iter[[1,1]]))
        
        # For Vs
        SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
        
        for (s in 1:q) {
          SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/tau2.iter[[1,1]]))
        }
      }
    }
    
    # If estimating the underlying structure
    if (is.null(scores)) {
      # -------------------------------------------------------------------------
      # Posterior sample for V
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        if (!response_given) {
          # Concatenating Ui's together
          U.iter.combined <- do.call(rbind, U.iter)
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          # The combined centered Xis with the latent response vector
          X.iter <- do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <-  tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
        if (response_given) {
          # Concatenating Ui's together
          U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv) 
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          if (response_type == "binary") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          if (response_type == "continuous") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <- tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
      }
      
      if (r == 0) {
        V.draw[[iter+1]][[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      # Updating the value of V
      V.iter <- V.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Us
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        for (s in 1:q) {
          Xs.iter <- X_complete[[s,1]] - W.iter[[s,s]] %*% t(Vs.iter[[1,s]])
          Bu <- solve((1/error_vars[s]) * t(V.iter[[1,1]]) %*% V.iter[[1,1]] + (1/sigma2_joint) * diag(r))
          U.draw[[iter+1]][[s,1]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bu <- (1/error_vars[s]) * t(V.iter[[1,1]]) %*% Xs.iter[j, ]
            
            U1j <- mvrnorm(1, mu = Bu %*% bu, Sigma = Bu)
            U1j
          }), nrow = r))
        }
      }
      
      if (r == 0) {
        for (s in 1:q) {
          U.draw[[iter+1]][[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
      }
      
      U.iter <- U.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Vs, s=1,...,q
      # -------------------------------------------------------------------------
      
      if (!response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
            Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      if (response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            # Combined Ws and beta
            W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))
            
            tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
            tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)
            
            Bvs <- solve(tW_Sigma_W + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            if (response_type == "binary") {
              # Combined centered Xs and Z
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            if (response_type == "continuous") {
              # Combined centered Xs and Y
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Y_complete - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- tW_Sigma %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      # Update the current value of V
      Vs.iter <- Vs.draw[[iter+1]]
      
      # Combine current values of V and V.
      V.iter.star.joint <- V.iter
      if (r == 0) {
        V.iter.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
      } 
      
      Vs.iter.star <- Vs.iter
      for (s in 1:q) {
        if (r.vec[s] == 0) {
          Vs.iter.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
        } 
      }
      
      VStar.iter <- cbind(1, do.call(cbind, V.iter.star.joint), do.call(cbind, Vs.iter.star))
      
      # Save the current VStar
      VStar.draw[[iter+1]][[1,1]] <- VStar.iter
      
      # -------------------------------------------------------------------------
      # Posterior sample for W
      # -------------------------------------------------------------------------
      
      for (s in 1:q) {
        if (r.vec[s] > 0) {
          Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
          Bws <- solve((1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Vs.iter[[1,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
          
          W.draw[[iter+1]][[s,s]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bws <- (1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Xs.iter[j,] 
            
            Wsj <- mvrnorm(1, mu = Bws %*% bws, Sigma = Bws)
            Wsj
          }), nrow = r.vec[s]))
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
        if (r.vec[s] == 0) {
          W.draw[[iter+1]][[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
      }
      
      # Update the current value of W
      W.iter <- W.draw[[iter+1]]
      
    }
    
    # If structure from another method is provided
    if (!is.null(scores)) {
      VStar.iter <- VStar
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------
    
    if (response_given) {
      
      if (sparsity) {
        # Change the precision for the intercept
        diag(SigmaBetaInv)[1] <- 1/beta_vars[1]
        
        # Change the precision of those betas under the slab, excluding the intercept
        diag(SigmaBetaInv)[-1][gamma.iter[-1] == 1] <- 1/(beta_vars[-1][gamma.iter[-1] == 1])
        
        # Change the precision of those betas under the spike
        diag(SigmaBetaInv)[-1][gamma.iter[-1] == 0] <- 1000 # Can change this precision
      }
      
      if (response_type == "binary") {
        Bbeta <- solve(t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- t(VStar.iter) %*% Z.iter
      }
      
      if (response_type == "continuous") {
        Bbeta <- solve((1/tau2.iter[[1,1]]) * t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- (1/tau2.iter[[1,1]]) * t(VStar.iter) %*% Y_complete
      }
      
      beta.draw[[iter+1]][[1,1]] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      
      # Update the current value of beta
      beta.iter <- beta.draw[[iter+1]][[1,1]]
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE] 
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for spike-and-slab indicators and probabilities
    # -------------------------------------------------------------------------
    
    if (response_given & sparsity) {
      # Calculate the likelihood for beta under the spike and slab
      like_spike <- dnorm(beta.iter, mean = 0, sd = sqrt(1/1000), log = TRUE)
      like_slab <- dnorm(beta.iter, mean = 0, sd = sqrt(beta_vars), log = TRUE)
      
      # Calculating the probability that each gamma equals 1
      prob = sapply(1:n_beta, function(rs) {
        x = log(p.iter) + like_slab[rs,]
        y = log(1 - p.iter) + like_spike[rs,]
        exp(x - logSum(c(x,y)))
      })
      prob[1] <- 1 # Always include the intercept
      
      # Generating the gammas
      gamma.draw[[iter+1]][[1,1]] <- matrix(rbinom(n_beta, size = 1, prob = prob), ncol = 1)
      
      # Saving the current value of gamma
      gamma.iter <- gamma.draw[[iter+1]][[1,1]]
      
      # Generating the prior probabilities (excluding the intercept, hence n_beta - 1 and sum(gamma[-1]))
      p.draw[[iter+1]][[1,1]] <- rbeta(1, 1 + sum(gamma.iter[-1,]), 1 + (n_beta-1) - sum(gamma.iter[-1,]))
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        tau2.draw[[iter+1]][[1,1]] <- matrix(1/rgamma(1, shape = shape + (n/2), rate = rate + 0.5 * sum((Y_complete - VStar.iter %*% beta.iter)^2)))
        
        # Update the current value of tau2
        tau2.iter <- tau2.draw[[iter+1]][[1,1]]
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Z.draw[[iter+1]][[1,1]] <- matrix(sapply(1:n, function(i) {
          if (Y_complete[i,] == 1) {
            rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          } else {
            rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          }
        }), ncol = 1)
      }
    }
    
    # -------------------------------------------------------------------------
    # Impute missing data
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (missingness_in_response) {
        if (response_type == "continuous") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter[[1,1]])), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
        
        if (response_type == "binary") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
      }
    }
    
    if (missingness_in_data) {
      for (s in 1:q) {
        Es <-  matrix(rnorm(p.vec[s]*n, 0, sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
        Xm.draw[[iter+1]][[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the joint and individual structure, scaled to the data
  # ---------------------------------------------------------------------------
  
  # Storing the structures at each Gibbs sampling iteration
  J.draw <- A.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  # Calculating the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  if (is.null(scores)) {
    for (iter in 1:nsample) {
      for (s in 1:q) {
        # Calculating the joint structure and scaling by sigma.mat
        J.draw[[iter]][[s,1]] <- (U.draw[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]]))
        
        # Calculating the individual structure and scaling by sigma.mat
        A.draw[[iter]][[s,1]] <- (W.draw[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) 
      }
      
      # Calculate the structure for Y
      if (response_given) {
        EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]]
      }
    }
  }
  
  # Return
  list(data = data, # Returning the scaled version of the data
       Y = Y, # Return the response vector
       sigma.mat = sigma.mat, # Scaling factors
       J.draw = J.draw, A.draw = A.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
       VStar.draw = VStar.draw, # Components that predict Y,
       Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
       scores = scores, # Scores if provided by another method 
       ranks = c(r, r.vec), # Ranks
       tau2.draw = tau2.draw, beta.draw = beta.draw, # Regression parameters
       gamma.draw = gamma.draw, p.draw = p.draw) # Sparsity parameters
  
}

# Initialize with Y as a source but do not scale the data or Y
bpmf_test <- function(data, Y, nninit = TRUE, model_params, ranks = NULL, scores = NULL, sparsity = FALSE, nsample, progress = TRUE, starting_values = NULL) {
  # Gibbs sampling algorithm for sampling the underlying structure and the 
  # regression coefficient vector for a response vector. 
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # 
  # data = matrix with list entries corresponding to each data source
  # Y = column vector with outcome or NULL
  # nninit = should the model be initialized with a nuclear norm penalized objective? if FALSE, provide ranks
  # model_params = (error_vars, joint_vars, indiv_vars, beta_vars = NULL, response_vars)
  # ranks = joint and individual ranks (1st entry = joint rank, kth for k>1 is the individual rank for the k-1'st source)
  # scores = if using structure from another method to fit a linear model, provide joint and individual scores here.
  #   in this case, ranks should be provided and nninit = FALSE, Y != NULL
  # nsample = number of Gibbs sampling iterations
  # progress = should the progress of the sampler be displayed?
  # starting_values = list of starting values for V, U, W, Vs. If NULL and nninit = TRUE, init with BIDIFAC+,
  #    if NULL and nninit = FALSE, init from prior. If not NULL, will init with provided starting values unless 
  #    nninit. 
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  q <- nrow(data) # Number of sources
  p.vec <- apply(data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(data[[1,1]]) # Number of subjects
  
  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------
  
  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
  
  # ---------------------------------------------------------------------------
  # Check for missingness in data
  # ---------------------------------------------------------------------------
  
  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
  # If there is missingness, initialize the missing values with 0s
  for (s in 1:q) {
    data[[s,1]][missing_obs[[s]]] <- 0
  }
  
  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------
  
  response_given <- !is.null(Y[[1,1]]) 
  
  # If so, what kind of response is it?
  if (response_given) {
    Y <- matrix(unlist(Y))
    
    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"
    
    # If there is a response, is there missingness in the outcome?
    missingness_in_response <- any(is.na(Y))
    
    # Which entries are missing?
    missing_obs_Y <- which(is.na(Y))
  }
  
  # ---------------------------------------------------------------------------
  # Obtaining the ranks 
  # ---------------------------------------------------------------------------
  
  scaled_data <- data
  scaled_Y <- Y
  
  # Initializing with BIDIFAC+
  if (nninit) {
    
    # Combine the data together 
    data_combined <- matrix(list(), nrow = (q+1), ncol = 1)
    
    for (s in 1:q) {
      # Append the scaled data
      data_combined[[s,1]] <- scaled_data[[s,1]]
      
      # Add back in missing data if any
      data_combined[[s,1]][missing_obs[[1]]] <- NA
    }
    
    # Initialize the indices of observations in each source 
    p.ind <- lapply(1:q, function(s) {
      if (s == 1) {
        1:p.vec[s]
      } else {
        (p.vec[s-1] + 1):cumsum(p.vec)[s]
      }
    })
    
    # Save indices of samples per source
    n.ind <- list(1:n)
    
    # If no response is given, p.ind.list should only include the sources
    if (!response_given) {
      
      # Save the indices for features for each identified structure
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- p.ind[[s]]
      }
    }
    
    # Include the response as a source if a response is given and add it to p.ind, p.ind.list
    if (response_given) {
      
      # Append the response to the sources
      data_combined[[q+1,1]] <- t(scaled_Y)
      
      # Add back in any missingness if exists
      data_combined[[q+1,1]][missing_obs_Y] <- NA
      
      # Add the response indices to p.ind
      p.ind[[q+1]] <- cumsum(p.vec)[q] + 1 # For the response
      
      # Save the indices for features for each identified structure
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- c(p.ind[[s]], p.ind[[q+1]])
      }
      
    }
    
    # Collapse the data into one matrix
    data_combined <- do.call(rbind, data_combined)
    
    # Save the indices for samples in each identified structure
    n.ind.list <- lapply(1:(q+1), function(s) c(1:n))
    
    # Run BIDIFAC+ --
    
    # If there is no missing data
    if (!missingness_in_data & !missingness_in_response) {
      rank_init <- bidifac.plus.given(data_combined, p.ind = p.ind, n.ind = n.ind,
                                      p.ind.list = p.ind.list, n.ind.list = n.ind.list)
    }
    
    # If there is missing data
    if (missingness_in_data | missingness_in_response) {
      # Save the indices of the missing values
      all.miss <- which(is.na(data_combined))
      
      # Initialize
      rank_init <- bidifac.plus.impute(data_combined, p.ind = p.ind, n.ind = n.ind,
                                       p.ind.list = p.ind.list, n.ind.list = n.ind.list,
                                       all.miss = all.miss)
    }
    
    # Print when finished
    print("Posterior mode obtained, ranks determined.")
    
    # Obtain the structures and their ranks
    S <- rank_init$S
    
    # Joint rank
    r <- Matrix::rankMatrix(S[[1]])[[1]]
    
    # Individual ranks if no response or pairwise-shared ranks of each source with Y if response given
    r.vec <- sapply(2:(q+1), function(s) Matrix::rankMatrix(S[[s]])[[1]])
    
    # Save the estimated error variances (which are not estimated and fixed to 1)
    sigma.mat <- matrix(1, nrow = q+1, ncol = 1)
  }
  
  if (!nninit) {
    r <- ranks[1]
    r.vec <- ranks[-1]
    sigma.mat <- matrix(nrow = q, 1)
    for (s in 1:q) {
      sigma.mat[s,1] <- 1
    }
  }
  
  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total
  
  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    beta_vars <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    diag(Sigma_beta) <- beta_vars
  }
  
  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  U.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  W.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = q))
  
  if (!response_given) {
    beta.draw <- Z.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (!missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  if (!sparsity) {
    gamma.draw <- p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (response_given) {
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    if (sparsity) {
      gamma.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
      p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
    }
  }
  
  if (missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  # ---------------------------------------------------------------------------
  # Initialize V, U, V, W
  # ---------------------------------------------------------------------------
  
  # If initializing with nuclear norm, initialize sampling at posterior mode
  if (nninit) {
    
    # Initialize joint scores, V
    
    V0 <- matrix(list(), nrow = 1, ncol = 1)
    beta_joint0 <- matrix(list(), nrow = 1, ncol = 1)
    
    # If there is joint structure
    if (r > 0) {
      svd.joint <- svd(S[[1]])
      V0[[1,1]] <- (svd.joint$v[,1:r, drop = FALSE]) %*% diag(svd.joint$d[1:r], nrow = r)
      
      # Save beta_joint if response is given
      if (response_given) {
        beta_joint0[[1,1]] <- t(svd(S[[1]])$u[p.ind[[q+1]],1:r, drop = FALSE])
      }
    } 
    
    # If there is no joint structure
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      
      # Save beta joint if response is given
      if (response_given) {
        beta_joint0[[1,1]] <- matrix(0, nrow = r, ncol = 1)
      }
    }
    
    U0 <- matrix(list(), nrow = q, ncol = 1)
    Vs0 <- matrix(list(), nrow = 1, ncol = q)
    W0 <- matrix(list(), nrow = q, ncol = q)
    beta_indiv0 <- matrix(list(), nrow = q, ncol = 1)
    
    for (s in 1:q) {
      
      # Initialize joint loadings, U 
      if (r > 0) {
        U0[[s,1]] <- svd(S[[1]])$u[p.ind[[s]],1:r, drop = FALSE]
      } 
      if (r == 0) {
        U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
      }
      
      # Initialize individual loadings, W, and individual scores, Vs, and beta_indiv if response is given
      if (r.vec[s] > 0) {
        
        # Compute SVD
        svd.indiv.s <- svd(S[[s+1]])
        
        # Save scores and loadings
        Vs0[[1,s]] <- (svd.indiv.s$v[,1:r.vec[s], drop = FALSE]) %*% diag(svd.indiv.s$d[1:r.vec[s]], nrow = r.vec[s])
        W0[[s,s]] <- svd.indiv.s$u[p.ind[[s]],1:r.vec[s], drop = FALSE]
        
        # Save beta_indiv
        if (response_given) {
          beta_indiv0[[s,1]] <- t(svd.indiv.s$u[p.ind[[q+1]],1:r.vec[s], drop = FALSE])
        }
        
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      } 
      
      # If there is no individual/pairwise-shared structure, set to 0
      if (r.vec[s] == 0) {
        
        # Saved scores and laodings
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
        W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        
        # Save beta_indiv if response is given
        if (response_given) {
          beta_indiv0[[s,1]] <- matrix(0, nrow = r.vec[s], ncol = 1)
        }
        
        # Fill in off-diagonal W with 0s
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      }
      
    }
    
    # Combine the coefficients
    beta0 <- rbind(rnorm(1, 0, Sigma_beta[1,1]), beta_joint0[[1,1]], do.call(rbind, beta_indiv0))
    
    # Combining the scores together 
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) V0.star[[1,1]] <- V0[[1,1]] else V0.star[[1,1]] <- matrix(nrow = n, ncol = r)
    
    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) Vs0.star[[1,s]] <- Vs0[[1,s]] else Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }
    VStar0 <- cbind(1, do.call(cbind, V0.star), do.call(cbind, Vs0.star))
    
    # Initialize the latent variable for a binary outcome
    if (response_given) {
      Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
    } 
    
  }
  
  # If ranks provided, initialize with prior or use given starting values
  if (!nninit) {
    
    # If no starting values were provided, initialize from priors
    if (is.null(starting_values)) {
      V0 <- matrix(list(), nrow = 1, ncol = 1)
      if (r > 0) {
        V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
      } 
      if (r == 0) {
        V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      U0 <- matrix(list(), nrow = q, ncol = 1)
      Vs0 <- matrix(list(), nrow = 1, ncol = q)
      W0 <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        
        # Initialize U
        if (r > 0) {
          U0[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
        } 
        if (r == 0) {
          U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
        
        # Initialize W and V
        if (r.vec[s] > 0) {
          Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
          W0[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        } 
        if (r.vec[s] == 0) {
          Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
          W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
      }
    }
    
    # If starting values were provided, use them as initial values
    if (!is.null(starting_values)) {
      V0 <- starting_values$V
      U0 <- starting_values$U
      W0 <- starting_values$W
      Vs0 <- starting_values$Vs
    }
    
  }
  
  # If imputing missingness in Y
  if (response_given) {
    tau20 <- matrix(1) # matrix(1/rgamma(1, shape = shape, rate = rate))  
    
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20[[1,1]])))[missing_obs_Y,, drop = FALSE]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }
  
  # If there is missingness in the data, generate starting values for the missing entries
  if (missingness_in_data) {
    Xm0 <- matrix(list(), ncol = 1, nrow = q)
    for (s in 1:q) {
      Xm0[[s,1]] <- rep(0, length(missing_obs[s]))
    }
  }
  
  # If there is missingness in Y, generate starting values for the missing entries
  if (response_given) {
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20)))[missing_obs_Y,, drop = FALSE]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Storing the initial values 
  # ---------------------------------------------------------------------------
  
  V.draw[[1]] <- V0
  U.draw[[1]] <- U0
  Vs.draw[[1]] <- Vs0
  W.draw[[1]] <- W0
  
  if (response_given) {
    beta.draw[[1]][[1,1]] <- beta0
    Z.draw[[1]][[1,1]] <- Z0
    tau2.draw[[1]][[1,1]] <- tau20
    VStar.draw[[1]][[1,1]] <- VStar0
    
    if (missingness_in_response) {
      Ym.draw[[1]][[1,1]] <- Ym0
    }
    
    if (sparsity) {
      gamma.draw[[1]][[1,1]] <- gamma0
      p.draw[[1]][[1,1]] <- p0
    }
  }
  
  if (missingness_in_data) {
    Xm.draw[[1]] <- Xm0
  }
  
  # ---------------------------------------------------------------------------
  # Computing the inverses that don't change from iteration to iteration
  # ---------------------------------------------------------------------------
  
  if (!response_given) {
    # Error variance for X. 
    SigmaVInv <- diag(rep(1/error_vars, p.vec))
  }
  
  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined precisions between data and Z
      SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))
      
      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
      }
    } 
    
    # For beta - Combined precisions between intercept and all betas
    SigmaBetaInv <- solve(Sigma_beta)
  }
  
  # ---------------------------------------------------------------------------
  # If structure from another method is given, save the scores as VStar
  # ---------------------------------------------------------------------------
  
  if (!is.null(scores)) {
    VStar <- cbind(1, scores)
  }
  
  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------
  
  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))
    
    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------
    
    V.iter <- V.draw[[iter]]
    U.iter <- U.draw[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.draw[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.draw[[iter]][[1,1]] 
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), nrow = q, ncol = 1)
      
      # Breaking beta down into the intercept,
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE] 
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
      
      if (response_type == "binary") {
        Z.iter <- Z.draw[[iter]][[1,1]]
      }
      
      if (response_type == "continuous") {
        tau2.iter <- tau2.draw[[iter]][[1,1]]
      }
      
      if (missingness_in_response) {
        # Save the current imputations for the missing values
        Ym.iter <- Ym.draw[[iter]][[1,1]]
        
        # Creating the completed outcome vector
        Y_complete <- Y
        
        # Filling in the missing entries for R1 and R2. 
        Y_complete[missing_obs_Y,] <- Ym.iter
      }
      
      if (!missingness_in_response) {
        Y_complete <- Y
      }
      
      if (sparsity) {
        gamma.iter <- gamma.draw[[iter]][[1,1]]
        p.iter <- p.draw[[iter]][[1,1]]
      }
    }
    
    if (missingness_in_data) {
      # Creating the completed matrices. 
      X_complete <- data
      
      # Fill in the completed matrices with the imputed values
      for (s in 1:q) {
        X_complete[[s,1]][missing_obs[[s]]] <- Xm.draw[[iter]][[s,1]]
      }
    }
    
    if (!missingness_in_data) {
      X_complete <- data
    }
    
    # -------------------------------------------------------------------------
    # Computing the inverse that changes with tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        # For V - Combined error variances between X1, X2, and Y
        SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1/tau2.iter[[1,1]]))
        
        # For Vs
        SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
        
        for (s in 1:q) {
          SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/tau2.iter[[1,1]]))
        }
      }
    }
    
    # If estimating the underlying structure
    if (is.null(scores)) {
      # -------------------------------------------------------------------------
      # Posterior sample for V
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        if (!response_given) {
          # Concatenating Ui's together
          U.iter.combined <- do.call(rbind, U.iter)
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          # The combined centered Xis with the latent response vector
          X.iter <- do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <-  tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
        if (response_given) {
          # Concatenating Ui's together
          U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv) 
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          if (response_type == "binary") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          if (response_type == "continuous") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <- tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
      }
      
      if (r == 0) {
        V.draw[[iter+1]][[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      # Updating the value of V
      V.iter <- V.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Us
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        for (s in 1:q) {
          Xs.iter <- X_complete[[s,1]] - W.iter[[s,s]] %*% t(Vs.iter[[1,s]])
          Bu <- solve((1/error_vars[s]) * t(V.iter[[1,1]]) %*% V.iter[[1,1]] + (1/sigma2_joint) * diag(r))
          U.draw[[iter+1]][[s,1]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bu <- (1/error_vars[s]) * t(V.iter[[1,1]]) %*% Xs.iter[j, ]
            
            U1j <- mvrnorm(1, mu = Bu %*% bu, Sigma = Bu)
            U1j
          }), nrow = r))
        }
      }
      
      if (r == 0) {
        for (s in 1:q) {
          U.draw[[iter+1]][[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
      }
      
      U.iter <- U.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Vs, s=1,...,q
      # -------------------------------------------------------------------------
      
      if (!response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
            Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      if (response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            # Combined Ws and beta
            W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))
            
            tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
            tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)
            
            Bvs <- solve(tW_Sigma_W + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            if (response_type == "binary") {
              # Combined centered Xs and Z
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            if (response_type == "continuous") {
              # Combined centered Xs and Y
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Y_complete - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- tW_Sigma %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      # Update the current value of V
      Vs.iter <- Vs.draw[[iter+1]]
      
      # Combine current values of V and V.
      V.iter.star.joint <- V.iter
      if (r == 0) {
        V.iter.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
      } 
      
      Vs.iter.star <- Vs.iter
      for (s in 1:q) {
        if (r.vec[s] == 0) {
          Vs.iter.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
        } 
      }
      
      VStar.iter <- cbind(1, do.call(cbind, V.iter.star.joint), do.call(cbind, Vs.iter.star))
      
      # Save the current VStar
      VStar.draw[[iter+1]][[1,1]] <- VStar.iter
      
      # -------------------------------------------------------------------------
      # Posterior sample for W
      # -------------------------------------------------------------------------
      
      for (s in 1:q) {
        if (r.vec[s] > 0) {
          Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
          Bws <- solve((1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Vs.iter[[1,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
          
          W.draw[[iter+1]][[s,s]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bws <- (1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Xs.iter[j,] 
            
            Wsj <- mvrnorm(1, mu = Bws %*% bws, Sigma = Bws)
            Wsj
          }), nrow = r.vec[s]))
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
        if (r.vec[s] == 0) {
          W.draw[[iter+1]][[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
      }
      
      # Update the current value of W
      W.iter <- W.draw[[iter+1]]
      
    }
    
    # If structure from another method is provided
    if (!is.null(scores)) {
      VStar.iter <- VStar
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------
    
    if (response_given) {
      
      if (sparsity) {
        # Change the precision for the intercept
        diag(SigmaBetaInv)[1] <- 1/beta_vars[1]
        
        # Change the precision of those betas under the slab, excluding the intercept
        diag(SigmaBetaInv)[-1][gamma.iter[-1] == 1] <- 1/(beta_vars[-1][gamma.iter[-1] == 1])
        
        # Change the precision of those betas under the spike
        diag(SigmaBetaInv)[-1][gamma.iter[-1] == 0] <- 1000 # Can change this precision
      }
      
      if (response_type == "binary") {
        Bbeta <- solve(t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- t(VStar.iter) %*% Z.iter
      }
      
      if (response_type == "continuous") {
        Bbeta <- solve((1/tau2.iter[[1,1]]) * t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- (1/tau2.iter[[1,1]]) * t(VStar.iter) %*% Y_complete
      }
      
      beta.draw[[iter+1]][[1,1]] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      
      # Update the current value of beta
      beta.iter <- beta.draw[[iter+1]][[1,1]]
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE] 
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for spike-and-slab indicators and probabilities
    # -------------------------------------------------------------------------
    
    if (response_given & sparsity) {
      # Calculate the likelihood for beta under the spike and slab
      like_spike <- dnorm(beta.iter, mean = 0, sd = sqrt(1/1000), log = TRUE)
      like_slab <- dnorm(beta.iter, mean = 0, sd = sqrt(beta_vars), log = TRUE)
      
      # Calculating the probability that each gamma equals 1
      prob = sapply(1:n_beta, function(rs) {
        x = log(p.iter) + like_slab[rs,]
        y = log(1 - p.iter) + like_spike[rs,]
        exp(x - logSum(c(x,y)))
      })
      prob[1] <- 1 # Always include the intercept
      
      # Generating the gammas
      gamma.draw[[iter+1]][[1,1]] <- matrix(rbinom(n_beta, size = 1, prob = prob), ncol = 1)
      
      # Saving the current value of gamma
      gamma.iter <- gamma.draw[[iter+1]][[1,1]]
      
      # Generating the prior probabilities (excluding the intercept, hence n_beta - 1 and sum(gamma[-1]))
      p.draw[[iter+1]][[1,1]] <- rbeta(1, 1 + sum(gamma.iter[-1,]), 1 + (n_beta-1) - sum(gamma.iter[-1,]))
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        tau2.draw[[iter+1]][[1,1]] <- matrix(1) # matrix(1/rgamma(1, shape = shape + (n/2), rate = rate + 0.5 * sum((Y_complete - VStar.iter %*% beta.iter)^2)))
        
        # Update the current value of tau2
        tau2.iter <- tau2.draw[[iter+1]][[1,1]]
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Z.draw[[iter+1]][[1,1]] <- matrix(sapply(1:n, function(i) {
          if (Y_complete[i,] == 1) {
            rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          } else {
            rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          }
        }), ncol = 1)
      }
    }
    
    # -------------------------------------------------------------------------
    # Impute missing data
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (missingness_in_response) {
        if (response_type == "continuous") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter[[1,1]])), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
        
        if (response_type == "binary") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
      }
    }
    
    if (missingness_in_data) {
      for (s in 1:q) {
        Es <-  matrix(rnorm(p.vec[s]*n, 0, sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
        Xm.draw[[iter+1]][[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the joint and individual structure, scaled to the data
  # ---------------------------------------------------------------------------
  
  # Storing the structures at each Gibbs sampling iteration
  J.draw <- A.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  # Calculating the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  if (is.null(scores)) {
    for (iter in 1:nsample) {
      for (s in 1:q) {
        # Calculating the joint structure and scaling by sigma.mat
        J.draw[[iter]][[s,1]] <- (U.draw[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]]))
        
        # Calculating the individual structure and scaling by sigma.mat
        A.draw[[iter]][[s,1]] <- (W.draw[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) 
      }
      
      # Calculate the structure for Y
      if (response_given) {
        EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]]
      }
    }
  }
  
  # Return
  list(data = data, # Returning the scaled version of the data
       Y = Y, # Return the response vector
       sigma.mat = sigma.mat, # Scaling factors
       J.draw = J.draw, A.draw = A.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
       VStar.draw = VStar.draw, # Components that predict Y,
       Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
       scores = scores, # Scores if provided by another method 
       ranks = c(r, r.vec), # Ranks
       tau2.draw = tau2.draw, beta.draw = beta.draw, # Regression parameters
       gamma.draw = gamma.draw, p.draw = p.draw) # Sparsity parameters
  
}

# Initialize with Y as a source and DO scale the data and Y
bpmf_test_scale <- function(data, Y, nninit = TRUE, model_params, ranks = NULL, scores = NULL, sparsity = FALSE, nsample, progress = TRUE, starting_values = NULL, rmt.data = TRUE, rmt.Y = TRUE) {
  # Gibbs sampling algorithm for sampling the underlying structure and the 
  # regression coefficient vector for a response vector. 
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # 
  # data = matrix with list entries corresponding to each data source
  # Y = column vector with outcome or NULL
  # nninit = should the model be initialized with a nuclear norm penalized objective? if FALSE, provide ranks
  # model_params = (error_vars, joint_vars, indiv_vars, beta_vars = NULL, response_vars)
  # ranks = joint and individual ranks (1st entry = joint rank, kth for k>1 is the individual rank for the k-1'st source)
  # scores = if using structure from another method to fit a linear model, provide joint and individual scores here.
  #   in this case, ranks should be provided and nninit = FALSE, Y != NULL
  # nsample = number of Gibbs sampling iterations
  # progress = should the progress of the sampler be displayed?
  # starting_values = list of starting values for V, U, W, Vs. If NULL and nninit = TRUE, init with BIDIFAC+,
  #    if NULL and nninit = FALSE, init from prior. If not NULL, will init with provided starting values unless 
  #    nninit. 
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  q <- nrow(data) # Number of sources
  p.vec <- apply(data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(data[[1,1]]) # Number of subjects
  
  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------
  
  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  # response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
  
  # ---------------------------------------------------------------------------
  # Check for missingness in data
  # ---------------------------------------------------------------------------
  
  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
  # If there is missingness, initialize the missing values with 0s
  for (s in 1:q) {
    data[[s,1]][missing_obs[[s]]] <- 0
  }
  
  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------
  
  response_given <- !is.null(Y[[1,1]]) 
  
  # If so, what kind of response is it?
  if (response_given) {
    Y <- matrix(unlist(Y))
    
    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"
    
    # If there is a response, is there missingness in the outcome?
    missingness_in_response <- any(is.na(Y))
    
    # Which entries are missing?
    missing_obs_Y <- which(is.na(Y))
    
    # If there are missing entries, initialize them with 0s
    Y[missing_obs_Y,] <- 0
  }
  
  # ---------------------------------------------------------------------------
  # Obtaining the ranks 
  # ---------------------------------------------------------------------------
  
  scaled_data <- data
  scaled_Y <- Y
  
  # Initializing with BIDIFAC+
  if (nninit) {
    
    # Create a matrix to store the estimated error sd
    if (!response_given) {
      sigma.mat <- matrix(1, nrow = q, ncol = 1)
    }
    
    if (response_given) {
      sigma.mat <- matrix(1, nrow = q+1, ncol = 1)
    }
    
    # Scale the data to have error variance 1
    if (rmt.data) {
      
      # Iterate through the sources and estimate the error variance
      for (s in 1:q) {
        
        # Save the estimated error variance
        sigma.mat[s,] <- sigma.rmt(data[[s,1]])
        
        # Scale the data
        scaled_data[[s,1]] <- data[[s,1]]/sigma.mat[s,]
      }
      
    }
    
    # If a response vector is given, scale Y
    if (response_given & rmt.Y) {
      
      # Save the estimated error variance using the natural lasso
      design_lasso <- do.call(rbind, scaled_data)
      nl_cv <- nlasso_cv(x = t(design_lasso)[-c(missing_obs_Y),], y = Y[-c(missing_obs_Y),]) # Remove missing observations
      sigma.mat[s+1,] <- nl_cv$sig_obj
      
      # Scale the response
      scaled_Y <- Y/sigma.mat[s+1,]
    }
    
    # Combine the data together 
    if (!response_given) {
      data_combined <- matrix(list(), nrow = q, ncol = 1)
    }
    
    if (response_given) {
      data_combined <- matrix(list(), nrow = (q+1), ncol = 1)
    }
    
    for (s in 1:q) {
      # Append the scaled data
      data_combined[[s,1]] <- scaled_data[[s,1]]
      
      # Add back in missing data if any
      data_combined[[s,1]][missing_obs[[1]]] <- NA
    }
    
    # Initialize the indices of observations in each source 
    p.ind <- lapply(1:q, function(s) {
      if (s == 1) {
        1:p.vec[s]
      } else {
        (p.vec[s-1] + 1):cumsum(p.vec)[s]
      }
    })
    
    # Save indices of samples per source
    n.ind <- list(1:n)
    
    # If no response is given, p.ind.list should only include the sources
    if (!response_given) {
      
      # Save the indices for features for each identified structure
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- p.ind[[s]]
      }
    }
    
    # Include the response as a source if a response is given and add it to p.ind, p.ind.list
    if (response_given) {
      
      # Append the response to the sources
      data_combined[[q+1,1]] <- t(scaled_Y)
      
      # Add back in any missingness if exists
      data_combined[[q+1,1]][missing_obs_Y] <- NA
      
      # Add the response indices to p.ind
      p.ind[[q+1]] <- cumsum(p.vec)[q] + 1 # For the response
      
      # Save the indices for features for each identified structure
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- c(p.ind[[s]], p.ind[[q+1]])
      }
      
    }
    
    # Collapse the data into one matrix
    data_combined <- do.call(rbind, data_combined)
    
    # Save the indices for samples in each identified structure
    if (!response_given) {
      n.ind.list <- lapply(1:q, function(s) c(1:n))
    }
    
    if (response_given) {
      n.ind.list <- lapply(1:(q+1), function(s) c(1:n))
    }
    
    # Run BIDIFAC+ --
    
    # If there is no missing data
    if (!missingness_in_data & !missingness_in_response) {
      rank_init <- bidifac.plus.given(data_combined, p.ind = p.ind, n.ind = n.ind,
                                      p.ind.list = p.ind.list, n.ind.list = n.ind.list)
    }
    
    # If there is missing data
    if (missingness_in_data | missingness_in_response) {
      # Save the indices of the missing values
      all.miss <- which(is.na(data_combined))
      
      # Initialize
      rank_init <- bidifac.plus.impute(data_combined, p.ind = p.ind, n.ind = n.ind,
                                       p.ind.list = p.ind.list, n.ind.list = n.ind.list,
                                       all.miss = all.miss)
    }
    
    # Print when finished
    print("Posterior mode obtained, ranks determined.")
    
    # Obtain the structures and their ranks
    S <- rank_init$S
    
    # Joint rank
    r <- Matrix::rankMatrix(S[[1]])[[1]]
    
    # Individual ranks if no response or pairwise-shared ranks of each source with Y if response given
    r.vec <- sapply(2:(q+1), function(s) Matrix::rankMatrix(S[[s]])[[1]])
    
  }
  
  if (!nninit) {
    r <- ranks[1]
    r.vec <- ranks[-1]
    sigma.mat <- matrix(nrow = q, 1)
    for (s in 1:q) {
      sigma.mat[s,1] <- 1
    }
  }
  
  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total
  
  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    beta_vars <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    diag(Sigma_beta) <- beta_vars
  }
  
  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  U.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  W.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = q))
  
  if (!response_given) {
    beta.draw <- Z.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (!missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  if (!sparsity) {
    gamma.draw <- p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (response_given) {
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    if (sparsity) {
      gamma.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
      p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
    }
  }
  
  if (missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  # ---------------------------------------------------------------------------
  # Initialize V, U, V, W
  # ---------------------------------------------------------------------------
  
  # If initializing with nuclear norm, initialize sampling at posterior mode
  if (nninit) {
    
    # Initialize joint scores, V
    
    V0 <- matrix(list(), nrow = 1, ncol = 1)
    beta_joint0 <- matrix(list(), nrow = 1, ncol = 1)
    
    # If there is joint structure
    if (r > 0) {
      svd.joint <- svd(S[[1]])
      V0[[1,1]] <- (svd.joint$v[,1:r, drop = FALSE]) %*% diag(svd.joint$d[1:r], nrow = r)
      
      # Save beta_joint if response is given
      if (response_given) {
        beta_joint0[[1,1]] <- t(svd(S[[1]])$u[p.ind[[q+1]],1:r, drop = FALSE])
      }
    } 
    
    # If there is no joint structure
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      
      # Save beta joint if response is given
      if (response_given) {
        beta_joint0[[1,1]] <- matrix(0, nrow = r, ncol = 1)
      }
    }
    
    U0 <- matrix(list(), nrow = q, ncol = 1)
    Vs0 <- matrix(list(), nrow = 1, ncol = q)
    W0 <- matrix(list(), nrow = q, ncol = q)
    beta_indiv0 <- matrix(list(), nrow = q, ncol = 1)
    
    for (s in 1:q) {
      
      # Initialize joint loadings, U 
      if (r > 0) {
        U0[[s,1]] <- svd(S[[1]])$u[p.ind[[s]],1:r, drop = FALSE]
      } 
      if (r == 0) {
        U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
      }
      
      # Initialize individual loadings, W, and individual scores, Vs, and beta_indiv if response is given
      if (r.vec[s] > 0) {
        
        # Compute SVD
        svd.indiv.s <- svd(S[[s+1]])
        
        # Save scores and loadings
        Vs0[[1,s]] <- (svd.indiv.s$v[,1:r.vec[s], drop = FALSE]) %*% diag(svd.indiv.s$d[1:r.vec[s]], nrow = r.vec[s])
        W0[[s,s]] <- svd.indiv.s$u[p.ind[[s]],1:r.vec[s], drop = FALSE]
        
        # Save beta_indiv
        if (response_given) {
          beta_indiv0[[s,1]] <- t(svd.indiv.s$u[p.ind[[q+1]],1:r.vec[s], drop = FALSE])
        }
        
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      } 
      
      # If there is no individual/pairwise-shared structure, set to 0
      if (r.vec[s] == 0) {
        
        # Saved scores and laodings
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
        W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        
        # Save beta_indiv if response is given
        if (response_given) {
          beta_indiv0[[s,1]] <- matrix(0, nrow = r.vec[s], ncol = 1)
        }
        
        # Fill in off-diagonal W with 0s
        for (ss in 1:q) {
          if (ss != s) {
            if (r.vec[ss] > 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
            }
            
            if (r.vec[ss] == 0) {
              W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
            }
          }
        }
      }
      
    }
    
    # Combine the coefficients
    beta0 <- rbind(rnorm(1, 0, Sigma_beta[1,1]), beta_joint0[[1,1]], do.call(rbind, beta_indiv0))
    
    # Combining the scores together 
    V0.star <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) V0.star[[1,1]] <- V0[[1,1]] else V0.star[[1,1]] <- matrix(nrow = n, ncol = r)
    
    Vs0.star <- Vs0
    for (s in 1:q) {
      if (r.vec[s] > 0) Vs0.star[[1,s]] <- Vs0[[1,s]] else Vs0.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }
    VStar0 <- cbind(1, do.call(cbind, V0.star), do.call(cbind, Vs0.star))
    
    # Initialize the latent variable for a binary outcome
    if (response_given) {
      Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
    } 
    
  }
  
  # If ranks provided, initialize with prior or use given starting values
  if (!nninit) {
    
    # If no starting values were provided, initialize from priors
    if (is.null(starting_values)) {
      V0 <- matrix(list(), nrow = 1, ncol = 1)
      if (r > 0) {
        V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
      } 
      if (r == 0) {
        V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      U0 <- matrix(list(), nrow = q, ncol = 1)
      Vs0 <- matrix(list(), nrow = 1, ncol = q)
      W0 <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        
        # Initialize U
        if (r > 0) {
          U0[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
        } 
        if (r == 0) {
          U0[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
        
        # Initialize W and V
        if (r.vec[s] > 0) {
          Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
          W0[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        } 
        if (r.vec[s] == 0) {
          Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
          W0[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
      }
    }
    
    # If starting values were provided, use them as initial values
    if (!is.null(starting_values)) {
      V0 <- starting_values$V
      U0 <- starting_values$U
      W0 <- starting_values$W
      Vs0 <- starting_values$Vs
    }
    
  }
  
  # If imputing missingness in Y
  if (response_given) {
    tau20 <- matrix(1) # matrix(1/rgamma(1, shape = shape, rate = rate))  
    
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20[[1,1]])))[missing_obs_Y,, drop = FALSE]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }
  
  # If there is missingness in the data, generate starting values for the missing entries
  if (missingness_in_data) {
    Xm0 <- matrix(list(), ncol = 1, nrow = q)
    for (s in 1:q) {
      Xm0[[s,1]] <- rep(0, length(missing_obs[s]))
    }
  }
  
  # If there is missingness in Y, generate starting values for the missing entries
  if (response_given) {
    if (missingness_in_response) {
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20)))[missing_obs_Y,, drop = FALSE]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,, drop = FALSE]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Storing the initial values 
  # ---------------------------------------------------------------------------
  
  V.draw[[1]] <- V0
  U.draw[[1]] <- U0
  Vs.draw[[1]] <- Vs0
  W.draw[[1]] <- W0
  
  if (response_given) {
    beta.draw[[1]][[1,1]] <- beta0
    Z.draw[[1]][[1,1]] <- Z0
    tau2.draw[[1]][[1,1]] <- tau20
    VStar.draw[[1]][[1,1]] <- VStar0
    
    if (missingness_in_response) {
      Ym.draw[[1]][[1,1]] <- Ym0
    }
    
    if (sparsity) {
      gamma.draw[[1]][[1,1]] <- gamma0
      p.draw[[1]][[1,1]] <- p0
    }
  }
  
  if (missingness_in_data) {
    Xm.draw[[1]] <- Xm0
  }
  
  # ---------------------------------------------------------------------------
  # Computing the inverses that don't change from iteration to iteration
  # ---------------------------------------------------------------------------
  
  if (!response_given) {
    # Error variance for X. 
    SigmaVInv <- diag(rep(1/error_vars, p.vec))
  }
  
  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined precisions between data and Z
      SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))
      
      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
      }
    } 
    
    # For beta - Combined precisions between intercept and all betas
    SigmaBetaInv <- solve(Sigma_beta)
  }
  
  # ---------------------------------------------------------------------------
  # If structure from another method is given, save the scores as VStar
  # ---------------------------------------------------------------------------
  
  if (!is.null(scores)) {
    VStar <- cbind(1, scores)
  }
  
  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------
  
  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))
    
    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------
    
    V.iter <- V.draw[[iter]]
    U.iter <- U.draw[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.draw[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.draw[[iter]][[1,1]] 
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), nrow = q, ncol = 1)
      
      # Breaking beta down into the intercept,
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE] 
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
      
      if (response_type == "binary") {
        Z.iter <- Z.draw[[iter]][[1,1]]
      }
      
      if (response_type == "continuous") {
        tau2.iter <- tau2.draw[[iter]][[1,1]]
      }
      
      if (missingness_in_response) {
        # Save the current imputations for the missing values
        Ym.iter <- Ym.draw[[iter]][[1,1]]
        
        # Creating the completed outcome vector
        Y_complete <- Y
        
        # Filling in the missing entries for R1 and R2. 
        Y_complete[missing_obs_Y,] <- Ym.iter
      }
      
      if (!missingness_in_response) {
        Y_complete <- Y
      }
      
      if (sparsity) {
        gamma.iter <- gamma.draw[[iter]][[1,1]]
        p.iter <- p.draw[[iter]][[1,1]]
      }
    }
    
    if (missingness_in_data) {
      # Creating the completed matrices. 
      X_complete <- data
      
      # Fill in the completed matrices with the imputed values
      for (s in 1:q) {
        X_complete[[s,1]][missing_obs[[s]]] <- Xm.draw[[iter]][[s,1]]
      }
    }
    
    if (!missingness_in_data) {
      X_complete <- data
    }
    
    # -------------------------------------------------------------------------
    # Computing the inverse that changes with tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        # For V - Combined error variances between X1, X2, and Y
        SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1/tau2.iter[[1,1]]))
        
        # For Vs
        SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
        
        for (s in 1:q) {
          SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1/tau2.iter[[1,1]]))
        }
      }
    }
    
    # If estimating the underlying structure
    if (is.null(scores)) {
      # -------------------------------------------------------------------------
      # Posterior sample for V
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        if (!response_given) {
          # Concatenating Ui's together
          U.iter.combined <- do.call(rbind, U.iter)
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv)
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          # The combined centered Xis with the latent response vector
          X.iter <- do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <-  tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
        if (response_given) {
          # Concatenating Ui's together
          U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))
          
          # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
          tU_Sigma <- crossprod(U.iter.combined, SigmaVInv) 
          tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
          
          Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
          
          if (response_type == "binary") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          if (response_type == "continuous") {
            # The combined centered Xis with the latent response vector
            X.iter <- rbind(do.call(rbind, X_complete) - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                            t(Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter)))
          }
          
          V.draw[[iter+1]][[1,1]] <- t(matrix(sapply(1:n, function(i) {
            bv <- tU_Sigma %*% X.iter[,i]
            
            Vi <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
            Vi
          }), nrow = r))
        }
        
      }
      
      if (r == 0) {
        V.draw[[iter+1]][[1,1]] <- matrix(0, nrow = n, ncol = 1)
      }
      
      # Updating the value of V
      V.iter <- V.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Us
      # -------------------------------------------------------------------------
      
      if (r > 0) {
        for (s in 1:q) {
          Xs.iter <- X_complete[[s,1]] - W.iter[[s,s]] %*% t(Vs.iter[[1,s]])
          Bu <- solve((1/error_vars[s]) * t(V.iter[[1,1]]) %*% V.iter[[1,1]] + (1/sigma2_joint) * diag(r))
          U.draw[[iter+1]][[s,1]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bu <- (1/error_vars[s]) * t(V.iter[[1,1]]) %*% Xs.iter[j, ]
            
            U1j <- mvrnorm(1, mu = Bu %*% bu, Sigma = Bu)
            U1j
          }), nrow = r))
        }
      }
      
      if (r == 0) {
        for (s in 1:q) {
          U.draw[[iter+1]][[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
      }
      
      U.iter <- U.draw[[iter+1]]
      
      # -------------------------------------------------------------------------
      # Posterior sample for Vs, s=1,...,q
      # -------------------------------------------------------------------------
      
      if (!response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
            Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      if (response_given) {
        for (s in 1:q) {
          if (r.vec[s] > 0) {
            # Combined Ws and beta
            W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))
            
            tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
            tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)
            
            Bvs <- solve(tW_Sigma_W + (1/indiv_vars[s]) * diag(r.vec[s]))
            
            if (response_type == "binary") {
              # Combined centered Xs and Z
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            if (response_type == "continuous") {
              # Combined centered Xs and Y
              Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                               t(Y_complete - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                                   do.call(cbind, Vs.iter[1, !(1:q %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(1:q %in% s), 1])))
            }
            
            Vs.draw[[iter+1]][[1,s]] <- t(matrix(sapply(1:n, function(i) {
              bvs <- tW_Sigma %*% Xs.iter[, i]
              
              Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
              Vsi
            }), nrow = r.vec[s]))
          }
          
          if (r.vec[s] == 0) {
            Vs.draw[[iter+1]][[1,s]] <- matrix(0, nrow = n, ncol = 1)
          }
        }
      }
      
      # Update the current value of V
      Vs.iter <- Vs.draw[[iter+1]]
      
      # Combine current values of V and V.
      V.iter.star.joint <- V.iter
      if (r == 0) {
        V.iter.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
      } 
      
      Vs.iter.star <- Vs.iter
      for (s in 1:q) {
        if (r.vec[s] == 0) {
          Vs.iter.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
        } 
      }
      
      VStar.iter <- cbind(1, do.call(cbind, V.iter.star.joint), do.call(cbind, Vs.iter.star))
      
      # Save the current VStar
      VStar.draw[[iter+1]][[1,1]] <- VStar.iter
      
      # -------------------------------------------------------------------------
      # Posterior sample for W
      # -------------------------------------------------------------------------
      
      for (s in 1:q) {
        if (r.vec[s] > 0) {
          Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
          Bws <- solve((1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Vs.iter[[1,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
          
          W.draw[[iter+1]][[s,s]] <- t(matrix(sapply(1:p.vec[s], function(j) {
            bws <- (1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Xs.iter[j,] 
            
            Wsj <- mvrnorm(1, mu = Bws %*% bws, Sigma = Bws)
            Wsj
          }), nrow = r.vec[s]))
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
        
        if (r.vec[s] == 0) {
          W.draw[[iter+1]][[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
          
          for (ss in 1:q) {
            if (ss != s) {
              if (r.vec[ss] > 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
              }
              
              if (r.vec[ss] == 0) {
                W.draw[[iter+1]][[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
              }
            }
          }
        }
      }
      
      # Update the current value of W
      W.iter <- W.draw[[iter+1]]
      
    }
    
    # If structure from another method is provided
    if (!is.null(scores)) {
      VStar.iter <- VStar
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------
    
    if (response_given) {
      
      if (sparsity) {
        # Change the precision for the intercept
        diag(SigmaBetaInv)[1] <- 1/beta_vars[1]
        
        # Change the precision of those betas under the slab, excluding the intercept
        diag(SigmaBetaInv)[-1][gamma.iter[-1] == 1] <- 1/(beta_vars[-1][gamma.iter[-1] == 1])
        
        # Change the precision of those betas under the spike
        diag(SigmaBetaInv)[-1][gamma.iter[-1] == 0] <- 1000 # Can change this precision
      }
      
      if (response_type == "binary") {
        Bbeta <- solve(t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- t(VStar.iter) %*% Z.iter
      }
      
      if (response_type == "continuous") {
        Bbeta <- solve((1/tau2.iter[[1,1]]) * t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- (1/tau2.iter[[1,1]]) * t(VStar.iter) %*% Y_complete
      }
      
      beta.draw[[iter+1]][[1,1]] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      
      # Update the current value of beta
      beta.iter <- beta.draw[[iter+1]][[1,1]]
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      
      # Joint effect
      if (r != 0) beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE] else beta_joint.iter <- matrix(0)
      
      # Individual effects
      if (sum(r.vec) > 0) beta_indiv.iter.temp <- beta.iter[(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is no individual effect
        if (r.vec[s] == 0) beta_indiv.iter[[s, 1]] <- matrix(0)
        
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE] 
          if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for spike-and-slab indicators and probabilities
    # -------------------------------------------------------------------------
    
    if (response_given & sparsity) {
      # Calculate the likelihood for beta under the spike and slab
      like_spike <- dnorm(beta.iter, mean = 0, sd = sqrt(1/1000), log = TRUE)
      like_slab <- dnorm(beta.iter, mean = 0, sd = sqrt(beta_vars), log = TRUE)
      
      # Calculating the probability that each gamma equals 1
      prob = sapply(1:n_beta, function(rs) {
        x = log(p.iter) + like_slab[rs,]
        y = log(1 - p.iter) + like_spike[rs,]
        exp(x - logSum(c(x,y)))
      })
      prob[1] <- 1 # Always include the intercept
      
      # Generating the gammas
      gamma.draw[[iter+1]][[1,1]] <- matrix(rbinom(n_beta, size = 1, prob = prob), ncol = 1)
      
      # Saving the current value of gamma
      gamma.iter <- gamma.draw[[iter+1]][[1,1]]
      
      # Generating the prior probabilities (excluding the intercept, hence n_beta - 1 and sum(gamma[-1]))
      p.draw[[iter+1]][[1,1]] <- rbeta(1, 1 + sum(gamma.iter[-1,]), 1 + (n_beta-1) - sum(gamma.iter[-1,]))
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        tau2.draw[[iter+1]][[1,1]] <- matrix(1) # matrix(1/rgamma(1, shape = shape + (n/2), rate = rate + 0.5 * sum((Y_complete - VStar.iter %*% beta.iter)^2)))
        
        # Update the current value of tau2
        tau2.iter <- tau2.draw[[iter+1]][[1,1]]
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Z.draw[[iter+1]][[1,1]] <- matrix(sapply(1:n, function(i) {
          if (Y_complete[i,] == 1) {
            rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          } else {
            rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          }
        }), ncol = 1)
      }
    }
    
    # -------------------------------------------------------------------------
    # Impute missing data
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (missingness_in_response) {
        if (response_type == "continuous") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter[[1,1]])), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
        
        if (response_type == "binary") {
          Ym.draw[[iter+1]][[1,1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)[missing_obs_Y,, drop = FALSE]
        }
      }
    }
    
    if (missingness_in_data) {
      for (s in 1:q) {
        Es <-  matrix(rnorm(p.vec[s]*n, 0, sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
        Xm.draw[[iter+1]][[s,1]] <- matrix((U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]])
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the joint and individual structure, scaled to the data
  # ---------------------------------------------------------------------------
  
  # Storing the structures at each Gibbs sampling iteration
  J.draw <- A.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  # Calculating the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  if (is.null(scores)) {
    for (iter in 1:nsample) {
      for (s in 1:q) {
        # Calculating the joint structure and scaling by sigma.mat
        J.draw[[iter]][[s,1]] <- (U.draw[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]]))
        
        # Calculating the individual structure and scaling by sigma.mat
        A.draw[[iter]][[s,1]] <- (W.draw[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) 
      }
      
      # Calculate the structure for Y
      if (response_given) {
        EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]]
      }
    }
  }
  
  # Return
  list(data = data, # Returning the scaled version of the data
       Y = Y, # Return the response vector
       sigma.mat = sigma.mat, # Scaling factors
       J.draw = J.draw, A.draw = A.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
       VStar.draw = VStar.draw, # Components that predict Y,
       Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
       scores = scores, # Scores if provided by another method 
       ranks = c(r, r.vec), # Ranks
       tau2.draw = tau2.draw, beta.draw = beta.draw, # Regression parameters
       gamma.draw = gamma.draw, p.draw = p.draw) # Sparsity parameters
  
}

# -----------------------------------------------------------------------------
# Validation Simulation
# -----------------------------------------------------------------------------

bpmf_sim <- function(nsample, n_clust, p.vec, n, true_params, model_params, nsim = 1000, s2nX = NULL, s2nY = NULL, center = FALSE, nninit, ranks, 
                     response = NULL, missingness = NULL, prop_missing = NULL, entrywise = NULL, sparsity = FALSE) {
  
  # ---------------------------------------------------------------------------
  # Check availability of parameters
  # ---------------------------------------------------------------------------
  results_available <- c(TRUE, TRUE, !is.null(response), check_availability(response, "continuous"), 
                         check_availability(missingness, "missingness_in_data") | check_availability(missingness, "both"), 
                         check_availability(missingness, "missingness_in_response") | check_availability(missingness, "both"))
  names(results_available) <- c("joint structure", "indiv structure", "EY", "tau2", "Xm", "Ym")
  
  # sim_results <- lapply(1:nsim, function(i) list())
  cl <- makeCluster(n_clust)
  registerDoParallel(cl)
  funcs <- c("bpmf_data", "center_data", "bpmf", "get_results", "BIDIFAC",
             "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
             "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum")
  packs <- c("Matrix", "MASS", "truncnorm")
  sim_results <- foreach (sim_iter = 1:nsim, .packages = packs, .export = funcs, .verbose = TRUE) %dopar% {
    # for (sim_iter in 1:nsim) {
    # svMisc::progress(sim_iter/(nsim/100))
    
    # -------------------------------------------------------------------------
    # Generating the data
    # -------------------------------------------------------------------------
    
    sim_data <- bpmf_data(p.vec, n, ranks, true_params, s2nX, s2nY, response, missingness, entrywise, prop_missing, sparsity)
    
    # Saving the data
    true_data <- sim_data$data
    q <- length(p.vec)
    joint.structure <- sim_data$joint.structure
    indiv.structure <- sim_data$indiv.structure
    
    missing_data <- sim_data$missing_data
    missing_obs <- sim_data$missing_obs
    
    # The response
    Y <- sim_data$Y
    Y_missing <- sim_data$Y_missing
    missing_obs_Y <- sim_data$missing_obs_Y
    
    # The standardizing coefficients
    s2nX_coef <- sim_data$s2nX_coef
    s2nY_coef <- sim_data$s2nY_coef
    
    # The response parameters
    beta <- sim_data$beta
    tau2 <- sim_data$tau2
    EY <- sim_data$EY
    
    # Setting the ranks
    r <- ranks[1]
    r.vec <- ranks[-1]
    
    # -------------------------------------------------------------------------
    # Setting the proper variable names for the sampler
    # -------------------------------------------------------------------------
    
    if (is.null(missingness)) {
      observed_data <- true_data
      Y_observed <- Y
    }
    
    if (!is.null(missingness)) {
      if (missingness == "missingness_in_data") {
        observed_data <- missing_data
        Y_observed <- Y
      }
      
      if (missingness == "missingness_in_response") {
        observed_data <- true_data
        Y_observed <- Y_missing
      }
      
      if (missingness == "both") {
        observed_data <- missing_data
        Y_observed <- Y_missing
      }
    }
    
    # -------------------------------------------------------------------------
    # Center the data
    # -------------------------------------------------------------------------
    
    if (center) {
      center_the_data <- center_data(data = observed_data, structure = structure)
      
      # Saving the centered data
      observed_data <- center_the_data$data_centered
      means_for_centering <- center_the_data$means_for_centering
      structure_centered <- center_the_data$structure_centered
    }
    
    # -------------------------------------------------------------------------
    # Running the Gibbs sampling algorithm
    # -------------------------------------------------------------------------
    
    # Gibbs sampling
    res <- bpmf(data = observed_data, Y = Y_observed, nninit = nninit, model_params = model_params, ranks = ranks, scores = NULL, sparsity = sparsity, nsample, progress = TRUE)
    
    # -------------------------------------------------------------------------
    # Extracting the results for each of decomposition matrices
    # -------------------------------------------------------------------------
    U.draw <- res$U.draw
    V.draw <- res$V.draw
    W.draw <- res$W.draw
    Vs.draw <- res$Vs.draw
    beta.draw <- res$beta.draw
    Z.draw <- res$Z.draw
    tau2.draw <- res$tau2.draw
    Xm.draw <- res$Xm.draw
    Ym.draw <- res$Ym.draw
    
    # -------------------------------------------------------------------------
    # Adding a burn-in
    # -------------------------------------------------------------------------
    burnin <- nsample/2
    
    U.burnin <- U.draw[burnin:nsample]
    V.burnin <- V.draw[burnin:nsample]
    W.burnin <- W.draw[burnin:nsample]
    Vs.burnin <- Vs.draw[burnin:nsample]
    beta.burnin <- beta.draw[burnin:nsample]
    Z.burnin <- Z.draw[burnin:nsample]
    tau2.burnin <- tau2.draw[burnin:nsample]
    Xm.burnin <- Xm.draw[burnin:nsample]
    Ym.burnin <- Ym.draw[burnin:nsample]
    
    # -------------------------------------------------------------------------
    # Computing the underlying structure for each Xs after burn-in from the sampler
    # -------------------------------------------------------------------------
    joint.structure.burnin <- lapply(1:(burnin+1), function(res) {
      mat <- matrix(list(), nrow = q, ncol = 1)
      for (s in 1:q) {
        mat[[s,1]] <- U.burnin[[res]][[s,1]] %*% t(V.burnin[[res]][[1,1]])
      }
      mat
    })
    
    indiv.structure.burnin <- lapply(1:(burnin+1), function(res) {
      mat <- matrix(list(), nrow = q, ncol = 1)
      for (s in 1:q) {
        mat[[s,1]] <- W.burnin[[res]][[s,s]] %*% t(Vs.burnin[[res]][[1,s]])
      }
      mat
    })
    
    # -------------------------------------------------------------------------
    # Calculating the estimate for E(Y) at each iteration
    # -------------------------------------------------------------------------
    if (is.null(response)) {
      EY.burnin <- lapply(1:(burnin+1), function(res) matrix(list(), nrow = 1, ncol = 1))
    }
    
    if (!is.null(response)) {
      if (response == "continuous") {
        EY.burnin <- lapply(1:(burnin+1), function(res) {
          mat <- matrix(list(), nrow = 1, ncol = 1)
          
          V.burnin.star <- V.burnin[[res]]
          if (ranks[1] == 0) V.burnin.star[[1,1]] <- matrix(nrow = n, ncol = r)
          
          Vs.burnin.star <- Vs.burnin[[res]]
          for (s in 1:q) {
            if (r.vec[s] == 0) Vs.burnin.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
          }
          
          VStar.iter <- cbind(1, do.call(cbind, V.burnin.star), do.call(cbind, Vs.burnin.star))
          mat[[1,1]] <- VStar.iter %*% beta.burnin[[res]][[1,1]]
          mat
        })
      }
      
      if (response == "binary") {
        EY.burnin <- lapply(1:(burnin+1), function(res) {
          mat <- matrix(list(), nrow = 1, ncol = 1)
          
          V.burnin.star <- V.burnin[[res]]
          if (ranks[1] == 0) V.burnin.star[[1,1]] <- matrix(nrow = n, ncol = r)
          
          Vs.burnin.star <- Vs.burnin[[res]]
          for (s in 1:q) {
            if (r.vec[s] == 0) Vs.burnin.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
          }
          
          VStar.iter <- cbind(1, do.call(cbind, V.burnin.star), do.call(cbind, Vs.burnin.star))
          
          mat[[1,1]] <- pnorm(VStar.iter %*% beta.burnin[[res]][[1,1]])
          mat
        })
      }
    }
    
    # -------------------------------------------------------------------------
    # Saving the results from the Gibbs sampler
    # -------------------------------------------------------------------------
    
    # Saving the draws from the sampler together
    draws <- list(joint.structure.burnin = joint.structure.burnin,
                  indiv.structure.burnin = indiv.structure.burnin,
                  EY = EY.burnin,
                  tau2 = tau2.burnin,
                  Xm = Xm.burnin,
                  Ym = Ym.burnin)
    
    # -------------------------------------------------------------------------
    # Saving the truth together for comparison
    # -------------------------------------------------------------------------
    
    # Returns the true values that were not observed in the data (X. and y)
    # If no missing data, Xm <- Ym <- NULL
    Xm <- return_missing(true_data, missing_obs)
    Ym <- return_missing(Y, missing_obs_Y)
    
    truth <- list(joint.structure = joint.structure,
                  indiv.structure = indiv.structure,
                  EY = EY,
                  tau2 = tau2,
                  Xm = Xm,
                  Ym = Ym)
    
    # -------------------------------------------------------------------------
    # Checking coverage, MSE, and CI width
    # -------------------------------------------------------------------------
    
    sim_iter_results <- get_results(truth, draws, burnin, results_available, missing_obs, missing_obs_Y)
    
    # ---------------------------------------------------------------------------
    # Returning the results at the end of the loop
    # ---------------------------------------------------------------------------
    
    # Add the indices for the missing values in this sim_iter
    sim_iter_results$any_missing <- list(missing_obs = missing_obs,
                                         missing_obs_Y = missing_obs_Y)
    
    # Share the ranks used
    sim_iter_results$ranks <- list(ranks = res$ranks)
    
    # Return 
    sim_iter_results
    # sim_results[[sim_iter]] <- sim_iter_results
  }
  stopCluster(cl)
  
  # ---------------------------------------------------------------------------
  # Averaging the results
  # ---------------------------------------------------------------------------
  
  # Calculate the denominator for the average
  denominator <- calculate_denominator(sim_results, q, p.vec, n, nsim, results_available)
  
  # Taking the averages of all the results
  sim_results_avg <- average_results(sim_results, denominator, p.vec, n, q, nsim, results_available, missingness)
  
  # Save all the ranks
  sim_ranks <- matrix(list(), nrow = 1, ncol = 1)
  sim_ranks[[1,1]] <- do.call(rbind, lapply(sim_results, function(sim_iter) sim_iter$ranks$ranks))
  sim_results_avg$ranks_results <- sim_ranks
  
  # ---------------------------------------------------------------------------
  # Returning the results
  # ---------------------------------------------------------------------------
  
  sim_results_avg
  
}

# -----------------------------------------------------------------------------
# BIDIFAC Functions
# -----------------------------------------------------------------------------

frob <- function(X){ sum(X^2,na.rm=T) }

diag2 <- function(x) {
  if (length(x)==1) return(as.matrix(x))
  else if (length(x)>1) return(diag(x))
}

sample2 <- function(x) {
  if (length(x)==1) return(x)
  else if (length(x)>1) return(sample(x))
}

sigma.rmt=function(X){ estim_sigma(X,method="MAD") }

softSVD=function(X, lambda){
  svdX=svd(X)
  nuc=pmax(svdX$d-lambda,0)
  out=tcrossprod(svdX$u, tcrossprod( svdX$v,diag(nuc) ))
  return(list(out=out, nuc=sum(nuc)))
}

data.rearrange=function(data,rmt=F,sigma=NULL){
  out=NULL
  p=nrow(data)
  q=ncol(data)
  
  m.vec=rep(NA,p)
  n.vec= ncol(data[[1,1]]) # do.call(c, lapply(data[1,], ncol))
  
  if (is.null(sigma)) sigma=matrix(1,p,q)
  
  for (i in 1:p){
    dimm=do.call(cbind, lapply(data[i,],dim))
    m1=unique(dimm[1,])
    if (length(m1)==1 ){m.vec[i]=m1 } 
    else{ stop("the number of rows do not match.") }
    if (!all(dimm[2,], n.vec)){ stop("the number of columns do not match")}
    
    for (j in 1:q){
      if (rmt) sigma[i,j]=sigma.rmt(data[[i,j]])
      data[[i,j]]=data[[i,j]]/sigma[i,j]
    }
    
    out=rbind(out,do.call(cbind,data[i,]))
  }
  
  return(list(out=out, nrows=m.vec, ncols=n.vec, sigma.mat=sigma))
}

estim_sigma <- function (X, k = NA, method = c("LN", "MAD"), center = "TRUE") {
  method <- match.arg(method, c("LN", "MAD", "ln", "mad", "Ln", 
                                "Mad"), several.ok = T)[1]
  method <- tolower(method)
  if (inherits(X, "data.frame")) {
    X <- as.matrix(X)
  }
  if (sum(sapply(X, is.numeric)) < ncol(X)) {
    stop("all the variables must be numeric")
  }
  if (center == "TRUE") {
    X <- scale(X, scale = F)
  }
  n = nrow(X)
  p = ncol(X)
  svdX = svd(X)
  if (method == "ln" & is.na(k)) {
    warning("Since you did not specify k, k was estimated using the FactoMineR estim_ncp function")
    k <- estim_ncp(X, scale = F)$ncp
    print(paste("k = ", k))
  }
  if (center == "TRUE") {
    N <- (n - 1)
  }
  else {
    N <- n
  }
  if ((k >= min(N, p)) & (method == "ln")) {
    stop("the number k specified has to be smaller than the minimum of the number of rows or columns")
  }
  if (method == "ln") {
    if (k == 0) {
      sigma = sqrt(sum(svdX$d^2)/(N * p))
    }
    else {
      sigma <- sqrt(sum(svdX$d[-c(1:k)]^2)/(N * p - N * 
                                              k - p * k + k^2))
    }
  }
  else {
    beta <- min(n, p)/max(n, p)
    lambdastar <- sqrt(2 * (beta + 1) + 8 * beta/((beta + 
                                                     1 + (sqrt(beta^2 + 14 * beta + 1)))))
    wbstar <- 0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta + 
      1.43
    sigma <- median(svdX$d)/(sqrt(max(n, p)) * (lambdastar/wbstar))
  }
  return(sigma)
}

BIDIFAC=function(data,rmt=T, sigma=NULL,
                 start=NULL, out=FALSE,
                 eps=1e-3, max.iter=1000, pbar=TRUE, seed=NULL, scale_back = TRUE, ...){
  if (!is.null(seed)){set.seed(seed)}
  if (!rmt & class(sigma)!="matrix") stop("sigma must be a matrix.")
  
  fit=data.rearrange(data, rmt, sigma)
  sigma.mat=fit$sigma.mat
  X00=fit$out
  
  mvec=fit$nrows; nvec=fit$ncols
  p=length(mvec); q=length(nvec)
  rm(fit)
  
  start.ind.m=c(1, cumsum(mvec)[1:(p-1)]+1)
  end.ind.m=cumsum(mvec)
  
  start.ind.n=c(1, cumsum(nvec)[1:(q-1)]+1)
  end.ind.n=cumsum(nvec)
  
  lambda.G=sqrt(sum(mvec))+sqrt(sum(nvec))
  lambda.R=sqrt(mvec)+sqrt(sum(nvec))
  lambda.C=sqrt(sum(mvec))+sqrt(nvec)
  lambda.I=tcrossprod(sqrt(mvec), rep(1, length(nvec)))+
    tcrossprod(rep(1, length(mvec)),sqrt(nvec))
  
  if (!is.null(start)){
    G00=start[[1]]; R00=start[[2]]
    C00=start[[3]]; I00=start[[4]]
  } else {
    G00= matrix(0, nrow = sum(mvec), ncol = sum(nvec)) # replicate(sum(nvec),rnorm(sum(mvec)))
    R00= matrix(0, nrow = sum(mvec), ncol = sum(nvec)) # replicate(sum(nvec),rnorm(sum(mvec)))
    C00=replicate(sum(nvec),rnorm(sum(mvec)))
    I00=replicate(sum(nvec),rnorm(sum(mvec)))
  }
  
  G00.nuc=NA; R00.nuc=rep(NA, p)
  C00.nuc=rep(NA, q); I00.nuc=matrix(NA,p,q)
  
  bool=TRUE
  count=1;crit0=0
  if (pbar) pb = txtProgressBar(min = 0, max=max.iter, initial=0, char="-", style = 3)
  while (bool){
    if (pbar){  setTxtProgressBar(pb, count)  }
    crit0.old = crit0
    
    #Update G to 0
    fit1=softSVD(X00-R00-C00-I00,lambda.G)
    G00 <- matrix(0, nrow = sum(mvec), ncol = sum(nvec)) # G00=fit1$out; 
    G00.nuc=fit1$nuc
    
    #update R to 0
    for (i in 1:p){
      ind=start.ind.m[i]:end.ind.m[i]
      fit1=softSVD(X00[ind,]-G00[ind,]-C00[ind,]-I00[ind,], lambda.R[i])
      # R00[ind,]=fit1$out; 
      R00 <- matrix(0, nrow = sum(mvec), ncol = sum(nvec))
      R00.nuc[i]=fit1$nuc
    }
    
    for (j in 1:q){
      ind=start.ind.n[j]:end.ind.n[j]
      fit1=softSVD(X00[,ind]-G00[,ind]-R00[,ind]-I00[,ind], lambda.C[j])
      C00[,ind]=fit1$out; C00.nuc[j]=fit1$nuc
    }
    
    for (i in 1:p){
      for (j in 1:q){
        ind1= start.ind.m[i]:end.ind.m[i]
        ind2=start.ind.n[j]:end.ind.n[j]
        fit1=softSVD(X00[ind1,ind2]-G00[ind1,ind2]-R00[ind1,ind2]-C00[ind1,ind2], lambda.I[i,j])
        I00[ind1,ind2]=fit1$out; I00.nuc[i,j]=fit1$nuc
      }
    }
    
    crit0 = frob(X00-G00-R00-C00-I00)+
      2*lambda.G*G00.nuc+2*sum(lambda.R*R00.nuc)+
      2*sum(lambda.C*C00.nuc)+2*sum(lambda.I*I00.nuc)
    
    if (abs(crit0.old-crit0)<eps){ bool=FALSE }
    else if (count==max.iter){ bool=FALSE}
    else{ count = count+1 }
  }
  
  if (scale_back) {
    S00.mat=G00.mat=R00.mat=C00.mat=I00.mat=data
    for (i in 1:p){
      ind1= start.ind.m[i]:end.ind.m[i]
      for (j in 1:q){
        ind2=start.ind.n[j]:end.ind.n[j]
        G00.mat[[i,j]]=G00[ind1,ind2] *sigma.mat[i,j] 
        R00.mat[[i,j]]=R00[ind1,ind2] *sigma.mat[i,j]
        C00.mat[[i,j]]=C00[ind1,ind2] *sigma.mat[i,j]
        I00.mat[[i,j]]=I00[ind1,ind2] *sigma.mat[i,j]
        S00.mat[[i,j]]=G00.mat[[i,j]]+R00.mat[[i,j]]+C00.mat[[i,j]]+I00.mat[[i,j]]
      }
    }
  }
  
  if (!scale_back) {
    S00.mat=G00.mat=R00.mat=C00.mat=I00.mat=data
    for (i in 1:p){
      ind1= start.ind.m[i]:end.ind.m[i]
      for (j in 1:q){
        ind2=start.ind.n[j]:end.ind.n[j]
        G00.mat[[i,j]]=G00[ind1,ind2] #*sigma.mat[i,j] Remove the scaling back by the error variance
        R00.mat[[i,j]]=R00[ind1,ind2] #*sigma.mat[i,j]
        C00.mat[[i,j]]=C00[ind1,ind2] #*sigma.mat[i,j]
        I00.mat[[i,j]]=I00[ind1,ind2] #*sigma.mat[i,j]
        S00.mat[[i,j]]=G00.mat[[i,j]]+R00.mat[[i,j]]+C00.mat[[i,j]]+I00.mat[[i,j]]
      }
    }
  }
  
  return(list(X=data, S=S00.mat,
              G=G00.mat, R=R00.mat, C=C00.mat, I=I00.mat,
              sigma.mat=sigma.mat, n.vec=nvec,m.vec=mvec))
}

bidifac.plus.given <- function(X0,p.ind,n.ind,p.ind.list,n.ind.list, S = list(), pen=c(), given.inits = FALSE,max.iter=500,conv.thresh=0.001){
  library(RSpectra)
  
  n.source <- length(p.ind)
  n.type <- length(n.ind)
  max.comb=length(p.ind.list)
  lambda=1
  if(given.inits==FALSE){
    for(i in 1:max.comb){
      S[[i]]=array(rep(0,prod(dim(X0))),dim=dim(X0))
      pen[i]=0}
  }
  obj.vec <- c(sum(X0^2))
  max.comb=length(p.ind.list)
  X0.resid = X0-Reduce('+',S)
  obj.prev <- sum(X0.resid^2)+2*sum(pen)
  
  conv <- FALSE
  while(!conv) {
    for(jj in 1:max.iter){
      print(jj)
      for(k in 1:max.comb){
        print(paste('*',k))
        X0.temp <- X0.resid+S[[k]]
        cur.p.ind <- p.ind.list[[k]] 
        cur.n.ind <- n.ind.list[[k]] 
        a <- svd(X0.temp[cur.p.ind,cur.n.ind],nu=0,nv=0)$d
        nc <- sum(a>(lambda*(sqrt(length(cur.p.ind))+sqrt(length(cur.n.ind)))))
        if(nc>0){
          SVD <- svd(X0.temp[cur.p.ind,cur.n.ind], nu=nc,nv=nc)
          s.vals <- pmax(SVD$d[1:nc]-lambda*(sqrt(length(cur.p.ind))+sqrt(length(cur.n.ind))),0)
          Diag <- diag(s.vals,nrow=nc,ncol=nc)
          Est <- SVD$u%*%Diag%*%t(SVD$v)
          S[[k]] <- array(rep(0,prod(dim(X0))),dim=dim(X0))
          S[[k]][cur.p.ind,cur.n.ind] <- Est
          X0.resid[cur.p.ind,cur.n.ind] <- X0.temp[cur.p.ind,cur.n.ind]-Est
          X0.resid = X0-Reduce('+',S)
          pen[k] <- lambda*(sqrt(length(cur.n.ind))+sqrt(length(cur.p.ind)))*(sum(s.vals))
          obj.vec <- c(obj.vec,sum(X0.resid^2)+2*sum(pen))
        }
      }
      obj.cur <- sum(X0.resid^2)+2*sum(pen)
      conv_check <- abs(obj.cur-obj.prev)<conv.thresh
      if(conv_check) break
      obj.prev <- sum(X0.resid^2)+2*sum(pen)
    }
    
    # If model converged
    if (conv_check) {
      conv <- TRUE
    }
    
    # If model did not converge
    if (!conv_check) {
      conv <- FALSE
      max.iter <- 2*max.iter
    }
    
    # Break after awhile
    if (max.iter > 10000) break
  }
  
  Sums <- array(dim=c(max.comb,n.source,n.type))
  for(kk in 1:max.comb){for(j in 1:n.source){ for(i in 1:n.type){
    Sums[kk,j,i] = sum(S[[kk]][p.ind[[j]],n.ind[[i]]]^2)
  }}}
  
  return(list(S=S,Sums=Sums,obj.vec=obj.vec, conv=conv))
}

bidifac.plus.impute <- function(X0,p.ind,n.ind,p.ind.list,n.ind.list, all.miss, S = list(), pen=c(), given.inits = FALSE,max.iter=500,conv.thresh=0.001,temp.iter=2){
  X0[all.miss]=0  
  S <- list()
  pen <- c()
  temp.fac <- svd(X0,nu=0,nv=0)$d[1]/(sum(sqrt(dim(X0))))-1
  n.source=length(p.ind)
  n.type=length(n.ind)
  max.comb=length(p.ind.list)
  for(i in 1:max.comb){
    S[[i]]=array(rep(0,prod(dim(X0))),dim=dim(X0))
    pen[i]=0}
  obj.vec <- c(sum(X0^2))
  X0.resid = X0-Reduce('+',S)
  obj.prev <- sum(X0.resid^2)+2*sum(pen)
  
  conv <- FALSE
  while (!conv) {
    for(jj in 1:max.iter){
      print(jj)
      if(jj < temp.iter){
        lambda <- 1+(temp.iter-1-jj)/temp.iter*temp.fac
      }
      #for(k in 1:max.comb){
      for(k in c(1:max.comb)){
        #     print(paste('*',k))
        X0.temp <- X0.resid+S[[k]]
        cur.p.ind <- p.ind.list[[k]] 
        cur.n.ind <- n.ind.list[[k]] 
        a <- svd(X0.temp[cur.p.ind,cur.n.ind],nu=0,nv=0)$d
        nc <- sum(a>(lambda*(sqrt(length(cur.p.ind))+sqrt(length(cur.n.ind)))))
        if(nc>0){
          SVD <- svd(X0.temp[cur.p.ind,cur.n.ind], nu=nc,nv=nc)
          s.vals <- pmax(SVD$d[1:nc]-lambda*(sqrt(length(cur.p.ind))+sqrt(length(cur.n.ind))),0)
          Diag <- diag(s.vals,nrow=nc,ncol=nc)
          Est <- SVD$u%*%Diag%*%t(SVD$v)
          S[[k]] <- array(rep(0,prod(dim(X0))),dim=dim(X0))
          S[[k]][cur.p.ind,cur.n.ind] <- Est
          X0.resid[cur.p.ind,cur.n.ind] <- X0.temp[cur.p.ind,cur.n.ind]-Est
          X0.resid = X0-Reduce('+',S)
          pen[k] <- lambda*(sqrt(length(cur.n.ind))+sqrt(length(cur.p.ind)))*(sum(s.vals))
          obj.vec <- c(obj.vec,sum(X0.resid^2)+2*sum(pen))
        }
      }
      Sig=Reduce('+',S)
      X0[all.miss]=Sig[all.miss]
      X0.resid = X0-Reduce('+',S)
      obj.cur <- sum(X0.resid^2)+2*sum(pen)
      conv_check <- abs(obj.cur-obj.prev)<conv.thresh
      if(conv_check) break
      obj.prev <- sum(X0.resid^2)+2*sum(pen)
    }
    
    if (conv_check) {
      conv <- TRUE
    }
    
    if (!conv_check) {
      conv <- FALSE
      max.iter <- 2*max.iter
    }
    
    if (max.iter > 10000) break
  }
  
  Sums <- array(dim=c(max.comb,n.source,n.type))
  for(kk in 1:max.comb){for(j in 1:n.source){ for(i in 1:n.type){
    Sums[kk,j,i] = sum(S[[kk]][p.ind[[j]],n.ind[[i]]]^2)
  }}}
  
  return(list(S=S,Sums=Sums,obj.vec=obj.vec,Sig=Sig, conv = conv))
}

# -----------------------------------------------------------------------------
# Validation simulation functions
# -----------------------------------------------------------------------------

logSum <- function(l){
  #given l=log(v), where v is a vector,
  #computes log(sum(v))
  max(l) + log(sum(exp(l-max(l))))
}

# Checks if a response is given (and what type) and if there is missingness (and what type)
check_availability <- function(param, compare) {
  if (is.null(param)) return(FALSE)
  if (!is.null(param)) return(param == compare)
}

# Generate fake data depending on conditions
bpmf_data <- function(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response, missingness, entrywise, prop_missing, sparsity, identically_zero = FALSE, num_in_spike = NULL) {
  # Generates fake data depending on the dims provided in p.vec, n, and ranks
  # and the true parameters provided in `true_params`
  
  # ---------------------------------------------------------------------------
  # Arguments:
  # 
  # p.vec = vector with the number of predictors in each source
  # n = the sample size for each source (number of columns)
  # true_params = list of parameters for the underlying data
  # s2nX = signal-to-noise ratio in the data
  # s2nY = signal-to-noise ratio in the response (if continuous)
  # response = NULL if no response is desired, "continuous", or "binary" 
  # missingness = NULL, "missingness_in_data", "missingness_in_response"
  # entrywise = NULL if missingness = NULL, TRUE if entrywise missingness, FALSE if columnwise missingness
  # prop_missing = NULL if missingness = NULL, otherwise the proportion of entries missingness
  # sparsity = TRUE if generate regression coefficients under spike-and-slab prior, FALSE otherwise
  # identically_zero = TRUE if generating response with sparsity and want spike coefficients to be exactly 0
  # num_in_spike (vec of ints) = number of coefficients to be generated from spike if generating response with sparsity.
  #     Should have an integer for joint and each individual structure. May be NULL if want to generate data under prior. 
  # ---------------------------------------------------------------------------
  
  # -------------------------------------------------------------------------
  # Setting the dimensions and latent components
  # -------------------------------------------------------------------------
  
  q <- length(p.vec)
  n_beta <- 1 + sum(ranks)
  r <- ranks[1]
  r.vec <- ranks[-1]
  
  # Setting the true parameters
  error_vars <- true_params$error_vars
  sigma2_joint <- joint_var <- true_params$joint_var
  sigma2_indiv <- indiv_vars <- true_params$indiv_vars
  beta_vars <- true_params$beta_vars
  
  # -------------------------------------------------------------------------
  # Generating the underlying structure
  # -------------------------------------------------------------------------
  
  joint.structure <- indiv.structure <- matrix(list(), ncol = 1, nrow = q)
  
  V <- matrix(list(), nrow = 1, ncol = 1)
  
  if (r > 0) V[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
  if (r == 0) V[[1,1]] <- matrix(0, nrow = n, ncol = 1)
  
  U <- matrix(list(), nrow = q, ncol = 1)
  Vs <- matrix(list(), nrow = 1, ncol = q)
  W <- matrix(list(), nrow = q, ncol = q)
  
  E <- matrix(list(), nrow = q, ncol = 1)
  
  for (s in 1:q) {
    if (r > 0) U[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
    if (r == 0) U[[s,1]] <- matrix(0, nrow = p.vec[s], ncol = 1)
    
    if (r.vec[s] > 0) {
      Vs[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
      W[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])
      
      for (ss in 1:q) {
        if (ss != s) {
          if (r.vec[ss] > 0) W[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
          
          if (r.vec[ss] == 0) W[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = 1)
        }
      }
    }
    if (r.vec[s] == 0) {
      Vs[[1,s]] <- matrix(0, nrow = n, ncol = 1)
      W[[s,s]] <- matrix(0, nrow = p.vec[s], ncol = 1)
      
      for (ss in 1:q) {
        if (ss != s) {
          if (r.vec[ss] > 0) W[[s,ss]] <- matrix(0, nrow = p.vec[s], ncol = r.vec[ss])
          
          if (r.vec[ss] == 0) W[[s,ss]] <- matrix(0, nrow = p.vec[s], ncol = 1)
        }
      }
    }
    
    E[[s,1]] <- matrix(rnorm(p.vec[s]*n, sd = sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
    
    joint.structure[[s,1]] <- U[[s,1]] %*% t(V[[1,1]])
    indiv.structure[[s,1]] <- W[[s,s]] %*% t(Vs[[1,s]])
  }
  
  # -------------------------------------------------------------------------
  # Standardizing the variance of the signal in the data
  # -------------------------------------------------------------------------
  
  if (is.null(s2nX)) {
    s2nX_coef <- NULL
  }
  
  if (!is.null(s2nX)) {
    # Calculating the scaling coefficient so that the variance of the underlying structure = s2nX * noise variance
    s2nX_coef <- rep(0, q)
    
    for (s in 1:q) {
      s2nX_coef[s] <- sqrt(s2nX) * sd(c(E[[s,1]]))/sd(c(joint.structure[[s,1]] + indiv.structure[[s,1]]))
      
      joint.structure[[s,1]] <- s2nX_coef[s] * joint.structure[[s,1]]
      indiv.structure[[s,1]] <- s2nX_coef[s] * indiv.structure[[s,1]]
    }
  }
  
  # -------------------------------------------------------------------------
  # Calculating the observed data
  # -------------------------------------------------------------------------
  
  data <- matrix(list(), nrow = q, ncol = 1)
  
  for (s in 1:q) {
    data[[s,1]] <- joint.structure[[s,1]] + indiv.structure[[s,1]] + E[[s,1]]
  }
  
  # -------------------------------------------------------------------------
  # Adding a response if desired
  # -------------------------------------------------------------------------
  
  if (is.null(response)) {
    Y <- EY <- Y_missing <- beta <- tau2 <- gamma <- p.prior <- matrix(list(), nrow = 1, ncol = 1)
    s2nY_coef <- NULL
  }
  
  if (!is.null(response)) {
    
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    diag(Sigma_beta) <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    
    if (!sparsity) {
      # Generate betas
      beta <- matrix(list(), nrow = 1, ncol = 1)
      beta[[1,1]] <- matrix(mvrnorm(1, mu = rep(0, n_beta), Sigma = Sigma_beta), ncol = 1)
      p.prior <- gamma <- matrix(list(), nrow = 1, ncol = 1)
    }
    
    if (sparsity) {
      
      p.prior <- matrix(rbeta(1, 1, 1), ncol = 1)
      
      # If no prior specification on number of coefficients in spike
      if (is.null(num_in_spike)) {
        gamma <- matrix(rbinom(n_beta, size = 1, prob = p.prior), ncol = 1)
        gamma[1,] <- 1 # Always include the intercept
      }
      
      # If user specifies number of coefficients in spike, randomly choose which
      if (!is.null(num_in_spike)) {
        gamma <- matrix(1, nrow = n_beta, ncol = 1) 
        
        # Randomly choose which components to be in the spike from each joint and individual structure
        spike_inds <- c()
        
        # Randomly choose joint components for spike
        joint_spike_inds <- sort(sample(x = c(2:r), size = num_in_spike[1], replace = FALSE))
        spike_inds <- c(spike_inds, joint_spike_inds) 
        
        # Randomly choose which individual components for spike
        for (s in 1:q) {
          if (s == 1) {
            indiv_spike_inds_s <- 1 + r + sort(sample(x = c(1:r.vec[s]), size = num_in_spike[2], replace = FALSE))
          }
          
          if (s != 1) {
            indiv_spike_inds_s <- 1 + r + sum(r.vec[1:(s-1)]) + sort(sample(x = c(1:r.vec[s]), size = num_in_spike[s], replace = FALSE))
          }
          spike_inds <- c(spike_inds, indiv_spike_inds_s)
        }
        
        gamma[spike_inds,] <- 0
      }
  
      diag(Sigma_beta)[gamma == 0] <- 1/1000
      beta <- matrix(list(), nrow = 1, ncol = 1)
      beta[[1,1]] <- matrix(mvrnorm(1, mu = rep(0, n_beta), Sigma = Sigma_beta), ncol = 1)
      
      # If desired, set spike coefficients to be identically 0
      if (identically_zero) {
        beta[[1,1]][gamma == 0] <- 0
      }
    }
    
    # Combine the Vs
    V.star.joint <- V
    if (r == 0) {
      V.star.joint[[1,1]] <- matrix(nrow = n, ncol = r)
    }
    
    Vs.star <- Vs
    for (s in 1:q) {
      if (r.vec[s] == 0) Vs.star[[1,s]] <- matrix(nrow = n, ncol = r.vec[s])
    }
    
    VStar <- cbind(1, do.call(cbind, V.star.joint), do.call(cbind, Vs.star))
    
    if (response == "binary") {
      Y <- EY <- matrix(list(), nrow = 1, ncol = 1)
      EY[[1,1]] <- pnorm(VStar %*% beta[[1,1]]) # True probability of being a case
      Y[[1,1]] <- matrix(rbinom(n, size = 1, prob = EY[[1,1]]), ncol = 1)
      tau2 <- matrix(list(), nrow = 1, ncol = 1)
      s2nY_coef <- NULL
    }
    
    if (response == "continuous") {
      Y <- EY <- matrix(list(), nrow = 1, ncol = 1)
      EY[[1,1]] <- VStar %*% beta[[1,1]]
      error_y <- matrix(rnorm(n, mean = 0, sd = sqrt(error_vars[q+1])), ncol = 1)
      
      # -------------------------------------------------------------------------
      # Standardizing the variance of the signal in the response
      # -------------------------------------------------------------------------
      
      if (is.null(s2nY)) {
        s2nY_coef <- NULL
      }
      
      if (!is.null(s2nY)) {
        # Calculating the scaling coefficient so that the variance of the response = s2nY * noise variance
        s2nY_coef <- sqrt(s2nY) * sd(error_y)/sd(EY[[1,1]])
        EY[[1,1]] <- s2nY_coef * EY[[1,1]]
      }
      
      Y[[1,1]] <- EY[[1,1]] + error_y
    }
  }
  
  # -------------------------------------------------------------------------
  # Adding missingness if desired
  # -------------------------------------------------------------------------
  
  if (is.null(missingness)) {
    missing_data <- missing_obs <-  matrix(list(), nrow = q, ncol = 1)
    missing_obs_Y <- Y_missing <- matrix(list(), nrow = 1, ncol = 1)
  }
  
  if (!is.null(missingness)) {
    if (missingness != "missingness_in_data" & missingness != "both") { # If missing in response
      missing_data <- missing_obs <- matrix(list(), nrow = q, ncol = 1)
    }
    
    if (missingness != "missingness_in_response" & missingness != "both") { # If missing in data
      Y_missing <- missing_obs_Y <- matrix(list(), nrow = 1, ncol = 1)
    }
    
    if (missingness == "missingness_in_response" | missingness == "both") { # If missing in response
      missing_obs_Y <- matrix(list(), nrow = 1, ncol = 1)
      missing_obs_Y[[1,1]] <- sort(sample(1:n, size = prop_missing * n, replace = FALSE))
      
      Y_missing <- Y
      Y_missing[[1,1]][missing_obs_Y[[1,1]],] <- NA
    }
    
    if (missingness == "missingness_in_data" | missingness == "both") { # If missing in data
      
      missing_obs <- missing_data <- missing_cols <- matrix(list(), nrow = q, ncol = 1)
      
      if (entrywise) { # if removing observations entrywise
        for (s in 1:q) {
          # these are counters going down the columns of R. So 9 would be the 9th entry counting down. 
          missing_obs[[s,1]] <- sort(sample(x = 1:length(data[[s,1]]), size = prop_missing*length(data[[s,1]]), replace = FALSE))
          
          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][missing_obs[[s,1]]] <- NA
        }
      } else { # if removing entire columns
        # Gives the column indices to remove
        
        for (s in 1:q) {
          # These are counters going down the columns of X. So 9 would be the 9th entry counting down. 
          missing_cols[[s,1]] <- sort(sample(x=1:n, size = n*prop_missing, replace = FALSE))
          
          if (s != 1) {
            avail_obs <- c(1:n)[!(c(1:n) %in% unlist(missing_cols[1:(s-1),]))]
            missing_cols[[s,1]] <- sample(x=avail_obs, size = n*prop_missing, replace = FALSE)
          }
          
          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][,missing_cols[[s,1]]] <- NA
          missing_obs[[s,1]] <- sort(which(is.na(missing_data[[s,1]])))
        }
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # Return
  # -------------------------------------------------------------------------
  
  list(data = data, # The "observed data"
       Y = Y, # The "observed outcome"
       missing_data = missing_data, # Missing data 
       missing_obs = missing_obs, # Missing data 
       Y_missing = Y_missing, missing_obs_Y = missing_obs_Y, # Missing data 
       s2nX = s2nX, s2nX_coef = s2nX_coef, # Scaling for s2n in data
       s2nY = s2nY, s2nY_coef = s2nY_coef, # Scaling for the response
       joint.structure = joint.structure, # Joint structure
       indiv.structure = indiv.structure, # Individual structure
       V = V, U = U, Vs = Vs, W = W, # Components of the structure
       beta = beta, EY = EY, gamma = gamma, p.prior = p.prior)
}

# Returns the entries in the true data that were not observed
return_missing <- function(true_param, missing_obs_inds) {
  # Returns the values in X or Y that were missing in the true data
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  #
  # true_param = the true value of X or Y
  # missing_obs_inds = the indices of the missing values in X or Y
  # ---------------------------------------------------------------------------
  
  any_missing <- !is.null(unlist(missing_obs_inds))
  
  if (!any_missing) {
    dm <- missing_obs_inds
  }
  
  if (any_missing) {
    dim_param <- nrow(missing_obs_inds)
    dm <- matrix(list(), nrow = dim_param, ncol = 1)
    
    for (s in 1:dim_param) {
      if (dim_param == 1) dm[[s,1]] <- true_param[[s,1]][missing_obs_inds[[s,1]],, drop = FALSE]
      if (dim_param != 1) dm[[s,1]] <- matrix(true_param[[s,1]][missing_obs_inds[[s,1]]])
    }
  }
  
  dm
}

# Center the data and the underlying structure 
center_data <- function(data, structure) {
  # Center each entry in data 
  q <- nrow(data)
  data_centered <- matrix(list(), nrow = q, ncol = 1)
  means_for_centering <- lapply(1:q, function(s) list())
  structure_centered <- list(joint = matrix(list(), nrow = q, ncol = 1),
                             indiv = matrix(list(), nrow = q, ncol = 1))
  
  for (s in 1:q) {
    # Center the data itself
    data_centered[[s,1]] <- scale(data[[s,1]], center = TRUE, scale = FALSE)
    
    # Use the means for centering to center the structure
    means_for_centering[[s]] <- attr(data_centered[[s,1]], "scaled:center")
    
    # Remove the attribute for the column means
    attr(data_centered[[s,1]], "scaled:center") <- NULL
    
    # Scale the structure
    structure_centered$joint[[s,1]] <- sweep(structure$joint[[s,1]], 2, means_for_centering[[s]])
    structure_centered$indiv[[s,1]] <- sweep(structure$indiv[[s,1]], 2, means_for_centering[[s]])
  }
  
  # Return
  list(data_centered = data_centered,
       means_for_centering = means_for_centering,
       structure_centered = structure_centered)
}

# Checking coverage for coverage simulation
check_coverage <- function(truth, draws, burnin) {
  # Checks whether the elements of truth are contained in the credible intervals
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # 
  # truth = matrix with the true values
  # draws = list of matrix results from a Gibbs sampler after burn-in that the 
  #         credible intervals will be based on.
  # ---------------------------------------------------------------------------
  
  # The dimensions of the true matrix
  N <- nrow(truth)
  M <- ncol(truth)
  
  Contained_Matrix <- sapply(1:M, function(col) { # for every column entry
    sapply(1:N, function(row) { # for every row entry
      
      # for every R matrix from the sampler
      row.col.entry <- sapply(1:(burnin+1), function(iter) { 
        draws[[iter]][row,col]
      })
      
      # Compute credible interval
      ci.row.col <- c(quantile(row.col.entry, 0.025), quantile(row.col.entry, 0.975))
      
      # Check if the true value is contained in the credible interval
      truth[row, col] >= ci.row.col[1] & truth[row,col] <= ci.row.col[2]
      
    })
  })
  
  return(Contained_Matrix)
}

# Computing the MSE
mse <- function(truth, draws) {
  # Computes the MSE between Gibbs sampling draws after burn-in and the 
  # true component being estimated. Returns a single value. 
  
  # ---------------------------------------------------------------------------
  # Arguments:
  # 
  # truth = true value of the parameter
  # draws = list of Gibbs samples for paramter
  # ---------------------------------------------------------------------------
  
  # Check if binary response
  is_binary_response <- all(unique(truth) %in% c(0,1))
  
  # Calculate the posterior mean for the draws after burn-in
  avg_res <- Reduce("+", draws)/length(draws)
  
  # Calculate the relative MSE
  if (!is_binary_response) {
    rel_mse <- frob(truth - avg_res)/frob(truth)
  }
  
  if (is_binary_response) {
    rel_mse <- frob(truth - avg_res)/frob(truth - mean(truth))
  }
  
  # Return
  rel_mse
  
}

# Computing credible interval width
ci_width <- function(draws, burnin) {
  cols <- ncol(draws[[1]])
  rows <- nrow(draws[[1]])
  
  sapply(1:cols, function(col) {
    sapply(1:rows, function(row) {
      
      # For every output matrix from the sampler
      row.col.entry <- sapply(1:(burnin+1), function(res) {
        draws[[res]][row,col]
      })
      
      # Compute credible interval
      ci.row.col <- c(quantile(row.col.entry, 0.025), quantile(row.col.entry, 0.975))
      
      # Calculate the length
      diff(ci.row.col)
    })
  })
}

# Combining checking coverage, MSE, and CI width calculations in one function
get_results <- function(truth, draws, burnin, results_available, missing_obs, missing_obs_Y) {
  # Combining checking coverage, MSE, and CI width calculations in one function
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # 
  # truth = list of true values for each parameter
  # draws = list of Gibbs samples for each parameter
  # burnin = number of samples used in burn-in
  # results_available = vector that checks which parameters are present under
  #     current simulation condition
  # missing_obs = list with the indices of the missing entries in X.
  # ---------------------------------------------------------------------------
  
  # Save the number of model parameters to check coverage for
  n_param <- length(draws)
  
  # Create a list to save the results 
  results <- lapply(1:n_param, function(i) list())
  names(results) <- names(truth)
  
  # Iterate through the parameters, checking the coverage, MSE, and CI width
  for (param in 1:n_param) {
    # Check the dimension of the current parameter
    dim_param <- nrow(truth[[param]])
    
    # Creating a matrix to store the results. 
    # (This is what will be returned if no results are available) 
    results[[param]] <- matrix(list(), nrow = dim_param, ncol = 1) 
    
    # Only calculate results if there are results
    if (results_available[param]) {
      
      # Iterate through the sources
      for (s in 1:dim_param) {
        # Save the samples for this source
        current_draws <- lapply(1:(burnin+1), function(iter) draws[[param]][[iter]][[s,1]])
        
        # Calculate coverage
        coverage_s <- check_coverage(truth[[param]][[s,1]], current_draws, burnin = burnin)
        
        # Calculate CI width
        ci_width_s <- ci_width(current_draws, burnin = burnin)
        
        # Calculate MSE
        # (Calculate the MSE separately for missing and not missing entries in the structure)
        if (param %in% c(1,2)) {
          if (!results_available[5]) { # If there is NO missing data
            mse_s_observed <- mse(truth[[param]][[s,1]], current_draws)
            mse_s_missing <- NULL
            mse_s <- list(observed = mse_s_observed, missing = mse_s_missing)
          }
          
          if (results_available[5]) { # If there IS missing data
            # For the entries in the structure that ARE observed -- 
            obs_inds <- 1:length(truth[[param]][[s,1]])
            
            # Save the true values that correspond to observed entries
            truth_observed <- truth[[param]][[s,1]][obs_inds[!(obs_inds %in% missing_obs[[s,1]])]]
            
            # Save the Gibbs samples corresponding to observed entries
            draws_observed <- lapply(current_draws, function(sim_iter) sim_iter[obs_inds[!(obs_inds %in% missing_obs[[s,1]])]])
            
            # MSE
            mse_s_observed <- mse(truth_observed, draws_observed)
            
            # For the entries in the structure that are NOT observed --
            
            # Save the true values that correspond to missing entries
            truth_missing <- truth[[param]][[s,1]][missing_obs[[s,1]]]
            
            # Save the Gibbs samples corresponding to missing entries
            draws_missing <- lapply(current_draws, function(sim_iter) sim_iter[missing_obs[[s,1]]])
            
            # MSE
            mse_s_missing <- mse(truth_missing, draws_missing)
            
            # Save
            mse_s <- list(observed = mse_s_observed, missing = mse_s_missing)
          }
        }
        
        if (param == 3) {
          if (!results_available[6]) { # If there is NO missing data
            mse_s_observed <- mse(truth[[param]][[s,1]], current_draws)
            mse_s_missing <- NULL
            mse_s <- list(observed = mse_s_observed, missing = mse_s_missing)
          }
          
          if (results_available[6]) { # If there IS missing data
            # For the entries in the structure that ARE observed -- 
            obs_inds <- 1:length(truth[[param]][[s,1]])
            
            # Save the true values that correspond to observed entries
            truth_observed <- truth[[param]][[s,1]][obs_inds[!(obs_inds %in% missing_obs_Y[[s,1]])]]
            
            # Save the Gibbs samples corresponding to observed entries
            draws_observed <- lapply(current_draws, function(sim_iter) sim_iter[obs_inds[!(obs_inds %in% missing_obs_Y[[s,1]])]])
            
            # MSE
            mse_s_observed <- mse(truth_observed, draws_observed)
            
            # For the entries in the structure that are NOT observed --
            
            # Save the true values that correspond to missing entries
            truth_missing <- truth[[param]][[s,1]][missing_obs_Y[[s,1]]]
            
            # Save the Gibbs samples corresponding to missing entries
            draws_missing <- lapply(current_draws, function(sim_iter) sim_iter[missing_obs_Y[[s,1]]])
            
            # MSE
            mse_s_missing <- mse(truth_missing, draws_missing)
            
            # Save
            mse_s <- list(observed = mse_s_observed, missing = mse_s_missing)
          }
        }
        
        # For all other parameters
        if (param != 1 & param != 2 & param != 3) mse_s <- mse(truth[[param]][[s,1]], current_draws)
        
        # Save the results
        results[[param]][[s,1]] <- list(coverage_s, mse_s, ci_width_s)
      }
    }
  }
  
  # Return the results
  results
}


# Count the number of times each entry in each parameter was observed (i.e. not missing)
calculate_denominator <- function(sim_results, q, p.vec, n, nsim, results_available) {
  # Count the number of times each entry in each parameter was observed for the 
  # underlying structure and how many times each entry in X and Y were missing. 
  # For parameters that are not missing, return the number of simulation replications
  
  # ---------------------------------------------------------------------------
  # Arguments:
  #
  # sim_results = coverage, MSE, and CI width for each parameter at each sim replication
  # q = number of sources
  # p.vec = number of features per source
  # n = number of samples
  # nsim = number of times simulation was run
  # results_available = which parameters were available under this sim condition?
  # ---------------------------------------------------------------------------
  
  # sim_results (list) = list of results from the simulation
  n_param <- length(results_available)
  
  # Create a list to store the results by model parameter
  counts <- lapply(1:n_param, function(i) list())
  names(counts) <- names(sim_results[[1]])[1:n_param]
  
  for (param in 1:n_param) {
    if (results_available[param]) {
      
      # Counting how many times entries in the underlying structure corresponded to observed X.
      if (param %in% c(1, 2)) {
        counts[[param]] <- matrix(list(), nrow = q, ncol = 1)
        
        for (s in 1:q) {
          # Vector of indices for each entry
          obs_inds <- 1:(p.vec[s] * n) 
          
          # Initialize counter for how many times each entry was NOT missing
          counts[[param]][[s,1]] <- rep(0, length(obs_inds))
          
          for (sim_iter in 1:nsim) {
            # Saving the missing entries for sim_iter simulation replication
            current_missing_obs <- sim_results[[sim_iter]]$any_missing$missing_obs[[s,1]] 
            counts[[param]][[s,1]] <- counts[[param]][[s,1]] + !(obs_inds %in% current_missing_obs)
          }
        }
      }
      
      # For E(Y), count how many times each entry corresponded to an observed response value
      if (param == 3) {
        counts[[param]] <- matrix(list(), nrow = 1, ncol = 1)
        
        # Vector of indices for each entry
        obs_inds <- 1:n
        
        # Initialize counter for how many times each entry was NOT missing
        counts[[param]][[1,1]] <- rep(0, length(obs_inds))
        
        for (sim_iter in 1:nsim) {
          # Saving the missing entries for sim_iter simulation replication
          current_missing_obs <- sim_results[[sim_iter]]$any_missing$missing_obs_Y[[1,1]] 
          counts[[param]][[1,1]] <- counts[[param]][[1,1]] + !(obs_inds %in% current_missing_obs)
        }
      }
      
      # For tau2, return the number of simulation replications 
      if (param == 4) {
        counts[[param]] <- matrix(list(), nrow = 1, ncol = 1)
        counts[[param]][[1,1]] <- rep(nsim, length(sim_results[[1]][[param]][[1,1]][[1]]))
      }
      
      # For Xm, count how many times in X. each entry WAS missing
      if (param == 5) {
        counts[[param]] <- matrix(list(), nrow = q, ncol = 1)
        
        for (s in 1:q) {
          obs_inds <- 1:(p.vec[s] * n) 
          counts[[param]][[s,1]] <- rep(0, length(obs_inds))
          
          for (sim_iter in 1:nsim) {
            current_missing_obs <- sim_results[[sim_iter]]$any_missing$missing_obs[[s]]
            counts[[param]][[s,1]] <- counts[[param]][[s,1]] + (obs_inds %in% current_missing_obs)
          }
        }
      }
      
      # For Ym, count how many times in y each entry WAS missing
      if (param == 6) {
        counts[[param]] <- matrix(list(), nrow = 1, ncol = 1)
        obs_inds <- 1:n
        counts[[param]][[1,1]] <- rep(0, length(obs_inds))
        for (sim_iter in 1:nsim) {
          current_missing_obs <- sim_results[[sim_iter]]$any_missing$missing_obs_Y[[1,1]]
          counts[[param]][[1,1]] <- counts[[param]][[1,1]] + (obs_inds %in% current_missing_obs)
        }
      }
    }
    
    if (!results_available[param]) {
      counts[[param]] <- matrix(list(), nrow = 1, ncol = 1)
    }
  }
  
  # Return
  counts
}

# Returns a matrix with the coverage, MSE, and CI width for each entry when the entry was observed (observed = TRUE)
# or unobserved (observed = FALSE)
compile_missing_results <- function(param_by_sim_iter, s, dims, nsim, missing_obs_inds, type, observed) {
  # Create a matrix for the results from each metric
  coverage_mat <- ci_width_mat <- matrix(0, nrow = dims[1], ncol = dims[2])
  
  # Calculate one MSE for the entire parameter
  mse_mat <- matrix(0, nrow = 1, ncol = 1)
  
  for (sim_iter in 1:nsim) {
    # Save the missing observations
    missing_obs_sim_iter <- missing_obs_inds[[sim_iter]]
    
    if (type == "structure") {
      if (observed) {
        obs_inds <- 1:(dims[1] * dims[2])
        observed_obs_sim_iter <- obs_inds[!(obs_inds %in% missing_obs_sim_iter)]
        coverage_mat[observed_obs_sim_iter] <- coverage_mat[observed_obs_sim_iter] + param_by_sim_iter[[sim_iter]][[s,1]][[1]][observed_obs_sim_iter]
        mse_mat <- mse_mat + param_by_sim_iter[[sim_iter]][[s,1]][[2]]$observed
        ci_width_mat[observed_obs_sim_iter] <- ci_width_mat[observed_obs_sim_iter] + param_by_sim_iter[[sim_iter]][[s,1]][[3]][observed_obs_sim_iter]
      }
      
      if (!observed) {
        coverage_mat[missing_obs_sim_iter] <- coverage_mat[missing_obs_sim_iter] + param_by_sim_iter[[sim_iter]][[s,1]][[1]][missing_obs_sim_iter]
        mse_mat <- mse_mat + param_by_sim_iter[[sim_iter]][[s,1]][[2]]$missing
        ci_width_mat[missing_obs_sim_iter] <- ci_width_mat[missing_obs_sim_iter] + param_by_sim_iter[[sim_iter]][[s,1]][[3]][missing_obs_sim_iter]
      }
    }
    
    if (type != "structure") {
      coverage_mat[missing_obs_sim_iter] <- coverage_mat[missing_obs_sim_iter] + param_by_sim_iter[[sim_iter]][[s,1]][[1]]
      mse_mat <- mse_mat + param_by_sim_iter[[sim_iter]][[s,1]][[2]]
      ci_width_mat[missing_obs_sim_iter] <- ci_width_mat[missing_obs_sim_iter] + param_by_sim_iter[[sim_iter]][[s,1]][[3]]
    }
  }
  
  # Save the results
  compiled_results <- list(coverage = coverage_mat, mse = mse_mat, ci_width = ci_width_mat)
  
  # Return
  compiled_results
}

# Average the results before returning from the function
average_results <- function(sim_results, denominator, p.vec, n, q, nsim, results_available, missingness) {
  # results_list (list) = list of results to return that should be averaged
  # observed_counts (list) = list of how many times each observation in each dataset was missing
  # nsim (int) = the number of simulation iterations that were run
  
  # Save the number of parameters
  n_param <- length(results_available)
  
  # Initialize a list to contain the averages
  results_avg <- lapply(1:n_param, function(i) list())
  names(results_avg) <- names(results_available)
  
  for (param in 1:n_param) {
    # Save the coverage, mse, ci widths across the simulation replications for this parameter
    param_by_sim_iter <- lapply(sim_results, function(sim_iter) sim_iter[[param]])
    
    # Check if there are results
    results_present <- results_available[param]
    
    if (results_present) {
      # Initialize a list for the results
      dim_param <- nrow(param_by_sim_iter[[1]])
      results_for_param <- matrix(list(), nrow = dim_param, ncol = 1)
      
      # For the underlying structure
      if (param %in% c(1,2)) {
        for (s in 1:dim_param) {
          # Save the results for source s
          results_by_source <- lapply(param_by_sim_iter, function(sim_iter) sim_iter[[s,1]])
          
          # If there is NO missingness
          if (is.null(missingness)) {
            # Calculate the average coverage
            avg_coverage_source <- Reduce("+", lapply(results_by_source, function(sim_iter) sim_iter[[1]]))/denominator[[param]][[s,1]]
            
            # Calculate the average MSE
            avg_mse_source <- mean(sapply(results_by_source, function(sim_iter) sim_iter[[2]]$observed)) 
            
            # Calculate the average CI width
            avg_ci_width_source <- Reduce("+", lapply(results_by_source, function(sim_iter) sim_iter[[3]]))/denominator[[param]][[s,1]]
            
            # Save the results
            results_for_param[[s,1]] <- list(avg_coverage = mean(avg_coverage_source),
                                             avg_mse = mean(avg_mse_source),
                                             avg_ci_width = mean(avg_ci_width_source))
          }
          
          # If there IS missingness 
          if (!is.null(missingness)) {
            # If missingness is in the response, calculate avg as if there is no missingness
            if (missingness == "missingness_in_response") {
              # Calculate the average coverage
              avg_coverage_source <- Reduce("+", lapply(results_by_source, function(sim_iter) sim_iter[[1]]))/denominator[[param]][[s,1]]
              
              # Calculate the average MSE
              avg_mse_source <- mean(sapply(results_by_source, function(sim_iter) sim_iter[[2]]$observed)) 
              
              # Calculate the average CI width
              avg_ci_width_source <- Reduce("+", lapply(results_by_source, function(sim_iter) sim_iter[[3]]))/denominator[[param]][[s,1]]
              
              # Save the results
              results_for_param[[s,1]] <- list(avg_coverage = mean(avg_coverage_source),
                                               avg_mse = mean(avg_mse_source),
                                               avg_ci_width = mean(avg_ci_width_source))
            }
            
            # If missingness is in the data, calculate avg separately for missing and non-missing data
            if (missingness == "missingness_in_data" | missingness == "both") {
              
              # ---------------------------------------------------------------
              # Check the results for structure corresponding to missing entries in each dataset 
              # ---------------------------------------------------------------
              
              # Save the indices of underlying structure corresponding to missing X.
              missing_obs_inds <- lapply(sim_results, function(sim_iter) sim_iter$any_missing$missing_obs[[s,1]])
              
              # Save total coverage, MSE, and CI width across sim_iters for entries in structure corresponding to missing X.
              results_compiled_missing <- compile_missing_results(param_by_sim_iter, s, dims = c(p.vec[s], n), nsim, missing_obs_inds, type = "structure", observed = FALSE)       
              
              # Calculate the average coverage for missing entries
              avg_coverage_missing_source <- results_compiled_missing[[1]]/(nsim - denominator[[param]][[s,1]])
              
              # Calculate the average MSE for missing entries
              avg_mse_missing_source <- results_compiled_missing[[2]]/nsim
              
              # Calculate the average CI width for missing entries
              avg_ci_width_missing_source <- results_compiled_missing[[3]]/(nsim - denominator[[param]][[s,1]])
              
              # ---------------------------------------------------------------
              # Check the results for structure corresponding to non-missing entries in each dataset 
              # ---------------------------------------------------------------
              
              # Save total coverage, MSE, and CI width across sim_iters for entries in structure corresponding to observed X.
              results_compiled_observed <- compile_missing_results(param_by_sim_iter, s, dims = c(p.vec[s], n), nsim, missing_obs_inds, type = "structure", observed = TRUE)
             
               # Calculate the average for observed entries
              avg_coverage_observed_source <- results_compiled_observed[[1]]/(denominator[[param]][[s,1]])
              
              # Calculate the average MSE for observed entries
              avg_mse_observed_source <- results_compiled_observed[[2]]/nsim
              
              # Calculate the average CI width for observed entries
              avg_ci_width_observed_source <- results_compiled_observed[[3]]/(denominator[[param]][[s,1]])
              
              # Save the results
              results_for_param[[s,1]] <- list(observed = list(avg_coverage = mean(avg_coverage_observed_source[!is.nan(avg_coverage_observed_source)]),
                                                          avg_mse = mean(avg_mse_observed_source),
                                                          avg_ci_width = mean(avg_ci_width_observed_source[!is.nan(avg_ci_width_observed_source)])),
                                               missing = list(avg_coverage = mean(avg_coverage_missing_source[!is.nan(avg_coverage_missing_source)]),
                                                              avg_mse = mean(avg_mse_missing_source),
                                                              avg_ci_width = mean(avg_ci_width_missing_source[!is.nan(avg_ci_width_missing_source)])))
                                               
            }
          }
        }
      }
      
      # For E(Y)
      if (param == 3) {
        # Save the results for source s
        results_by_source <- lapply(param_by_sim_iter, function(sim_iter) sim_iter[[1,1]])
        
        # If there is NO missingness
        if (is.null(missingness)) {
          # Calculate the average coverage
          avg_coverage_source <- Reduce("+", lapply(results_by_source, function(sim_iter) sim_iter[[1]]))/denominator[[param]][[1,1]]
          
          # Calculate the average MSE
          avg_mse_source <- mean(sapply(results_by_source, function(sim_iter) sim_iter[[2]]$observed)) 
          
          # Calculate the average CI width
          avg_ci_width_source <- Reduce("+", lapply(results_by_source, function(sim_iter) sim_iter[[3]]))/denominator[[param]][[1,1]]
          
          # Save the results
          results_for_param[[1,1]] <- list(avg_coverage = mean(avg_coverage_source),
                                           avg_mse = mean(avg_mse_source),
                                           avg_ci_width = mean(avg_ci_width_source))
        }
        
        # If there IS missingness 
        if (!is.null(missingness)) {
          # If missingness is in the data, calculate avg as if there is no missingness
          if (missingness == "missingness_in_data") {
            # Calculate the average coverage
            avg_coverage_source <- Reduce("+", lapply(results_by_source, function(sim_iter) sim_iter[[1]]))/denominator[[param]][[1,1]]
            
            # Calculate the average MSE
            avg_mse_source <- mean(sapply(results_by_source, function(sim_iter) sim_iter[[2]]$observed)) 
            
            # Calculate the average CI width
            avg_ci_width_source <- Reduce("+", lapply(results_by_source, function(sim_iter) sim_iter[[3]]))/denominator[[param]][[1,1]]
            
            # Save the results
            results_for_param[[1,1]] <- list(avg_coverage = mean(avg_coverage_source),
                                             avg_mse = mean(avg_mse_source),
                                             avg_ci_width = mean(avg_ci_width_source))
          }
          
          # If missingness is in the response, calculate avg separately for missing and non-missing response
          if (missingness == "missingness_in_response" | missingness == "both") {
            
            # ---------------------------------------------------------------
            # Check the results for structure corresponding to missing entries in each dataset 
            # ---------------------------------------------------------------
            
            # Save the indices of underlying structure corresponding to missing X.
            missing_obs_inds <- lapply(sim_results, function(sim_iter) sim_iter$any_missing$missing_obs_Y[[1,1]])
            
            # Save total coverage, MSE, and CI width across sim_iters for entries in structure corresponding to missing Y
            results_compiled_missing <- compile_missing_results(param_by_sim_iter, s = 1, dims = c(n, 1), nsim, missing_obs_inds, type = "structure", observed = FALSE)       
            
            # Calculate the average coverage for missing entries
            avg_coverage_missing_source <- results_compiled_missing[[1]]/(nsim - denominator[[param]][[1,1]])
            
            # Calculate the average MSE for missing entries
            avg_mse_missing_source <- results_compiled_missing[[2]]/nsim
            
            # Calculate the average CI width for missing entries
            avg_ci_width_missing_source <- results_compiled_missing[[3]]/(nsim - denominator[[param]][[1,1]])
            
            # ---------------------------------------------------------------
            # Check the results for structure corresponding to non-missing entries in each dataset 
            # ---------------------------------------------------------------
            
            # Save total coverage, MSE, and CI width across sim_iters for entries in structure corresponding to observed X.
            results_compiled_observed <- compile_missing_results(param_by_sim_iter, s = 1, dims = c(n, 1), nsim, missing_obs_inds, type = "structure", observed = TRUE)
            
            # Calculate the average for observed entries
            avg_coverage_observed_source <- results_compiled_observed[[1]]/(denominator[[param]][[1,1]])
            
            # Calculate the average MSE for observed entries
            avg_mse_observed_source <- results_compiled_observed[[2]]/nsim
            
            # Calculate the average CI width for observed entries
            avg_ci_width_observed_source <- results_compiled_observed[[3]]/(denominator[[param]][[1,1]])
            
            # Save the results
            results_for_param[[1,1]] <- list(observed = list(avg_coverage = mean(avg_coverage_observed_source),
                                                             avg_mse = mean(avg_mse_observed_source),
                                                             avg_ci_width = mean(avg_ci_width_observed_source)),
                                             missing = list(avg_coverage = mean(avg_coverage_missing_source[!is.nan(avg_coverage_missing_source)]),
                                                            avg_mse = mean(avg_mse_missing_source),
                                                            avg_ci_width = mean(avg_ci_width_missing_source[!is.nan(avg_ci_width_missing_source)])))
            
          }
        }
      }
      
      # For tau2
      if (param == 4) {
        # Calculate the average
        avg_coverage <- Reduce("+", lapply(param_by_sim_iter, function(res) res[[1,1]][[1]]))/(denominator[[param]][[1,1]])

        # Calculate the average MSE
        avg_mse <- Reduce("+", lapply(param_by_sim_iter, function(res) res[[1,1]][[2]]))/nsim

        # Calculate the average CI width
        avg_ci_width <- Reduce("+", lapply(param_by_sim_iter, function(res) res[[1,1]][[3]]))/(denominator[[param]][[1,1]])

        # Save the results
        results_for_param[[1,1]] <- list(avg_coverage = mean(avg_coverage),
                                         avg_mse = mean(avg_mse),
                                         avg_ci_width = mean(avg_ci_width))
      }
      
      # For Xm
      if (param == 5) {
        for (s in 1:q) {
          # Save the indices for missing observations in X.
          missing_obs_inds <- lapply(sim_results, function(sim_iter) sim_iter$any_missing$missing_obs[[s,1]])
          
          # Save the total coverage, MSE, and CI width for missing observations in X.
          results_compiled <- compile_missing_results(param_by_sim_iter, s, dims = c(p.vec[s], n), nsim, missing_obs_inds, type = "observed_data")

          # Calculate the average
          avg_coverage_source <- results_compiled[[1]]/(denominator[[param]][[s,1]])
          
          # Calculate the average MSE
          avg_mse_source <- results_compiled[[2]]/nsim
          
          # Calculate the average CI width
          avg_ci_width_source <- results_compiled[[3]]/(denominator[[param]][[s,1]])
          
          # Save the results
          results_for_param[[s,1]] <- list(avg_coverage = mean(avg_coverage_source[!is.nan(avg_coverage_source)]),
                                           avg_mse = mean(avg_mse_source),
                                           avg_ci_width = mean(avg_ci_width_source[!is.nan(avg_ci_width_source)]))
        }
      }
        
      # For Ym
      if (param == 6) {
        missing_obs_inds <- lapply(sim_results, function(sim_iter) sim_iter$any_missing$missing_obs_Y[[1,1]])
        results_compiled <- compile_missing_results(param_by_sim_iter, s = 1, dims = c(n, 1), nsim, missing_obs_inds, type = "observed_data")
        
        # Calculate the average
        avg_coverage_source <- results_compiled[[1]]/(denominator[[param]][[1,1]])
        
        # Calculate the average MSE
        avg_mse_source <- results_compiled[[2]]/nsim
        
        # Calculate the average CI width
        avg_ci_width_source <- results_compiled[[3]]/(denominator[[param]][[1,1]])
        
        # Save the results
        results_for_param[[1,1]] <- list(avg_coverage = mean(avg_coverage_source[!is.nan(avg_coverage_source)]),
                                         avg_mse = mean(avg_mse_source),
                                         avg_ci_width = mean(avg_ci_width_source[!is.nan(avg_ci_width_source)]))
      }
    }
    
    if (!results_present) {
      results_for_param <- matrix(list(), nrow = nrow(param_by_sim_iter[[1]]), ncol = 1)
    }
    
    # Save the overall results
    results_avg[[param]] <- results_for_param
  }

  # Return
  results_avg
}

# Checking convergence
log_joint_density <- function(data, U.iter, V.iter, W.iter, Vs.iter, model_params, ranks, Y = NULL, beta.iter = NULL, tau2.iter = NULL, Xm.iter = NULL, Ym.iter = NULL, gamma.iter = NULL, p.iter = NULL) {
  
  # ---------------------------------------------------------------------------
  # Arguments:
  # 
  # data = matrix of lists, the observed data 
  # U.iter, V.iter, W.iter, Vs.iter = the values of the parameters at the 
  #   current Gibbs sampling iteration
  # model_params = the parameters in the model
  # ranks = the ranks of the joint and individual structures
  # Y = the response vector, if available
  # beta.iter = the value of the coefficient vector at the current Gibbs 
  #   sampling iteration
  # tau2.iter = the value of the variance of the response if a continuous
  #   response is given
  # Xm.draw = the imputed values at the current Gibbs sampling iteration
  # Ym.draw = the imputed values at the current Gibbs sampling iteration
  # ---------------------------------------------------------------------------
  
  library(invgamma)
  
  # How many sources are there?
  q <- nrow(data)
  
  # How many observations?
  n <- ncol(data[[1,1]])
  
  # Saving the model parameters
  error_vars <- model_params$error_vars # Error variances
  sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
  beta_vars <- model_params$beta_vars # Variances on betas
  response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
  
  r <- ranks[1]
  r.vec <- ranks[-1]
  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total
  
  Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
  diag(Sigma_beta) <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
  
  # Check if there is a response
  response_given <- !is.null(Y)
  
  # If there is a response, what type of response is it?
  if (response_given) {
    response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"
  }
  
  # Check if there is missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
  # If there is missingness, fill in the missing values with the current imputed values
  if (missingness_in_data) {
    # Creating the completed matrices. 
    X_complete <- data
    
    # Fill in the completed matrices with the imputed values
    for (s in 1:q) {
      X_complete[[s,1]][missing_obs[[s]]] <- Xm.iter[[s,1]]
    }
  }
  
  # If there is no missingness, rename the data
  if (!missingness_in_data) {
    X_complete <- data
  }
  
  # Check if there is missingness in the response
  missingness_in_response <- any(is.na(Y[[1,1]]))
  
  # Which entries are missing?
  missing_obs_Y <- which(is.na(Y[[1,1]])) 
  
  # If there is missingness in the response, fill in the missing values with the current imputed values
  if (missingness_in_response) {
    # Create a completed response
    Y_complete <- Y
    
    # Fill in the missing values
    Y_complete[[1,1]][missing_obs_Y] <- Ym.iter[[1,1]]
  }
  
  # If there is no missingness, rename the response
  if (!missingness_in_response) {
    Y_complete <- Y
  }
  
  # Check if the model was fit with sparsity
  sparsity <- !is.null(gamma.iter)
  
  # Saving the contributions of each term to the joint density
  like <- 0
  
  # Contribution of V to the joint density
  like <- like + sum(sapply(1:r, function(rs) {
    dnorm(V.iter[[1,1]][,rs], mean = 0, sd = sqrt(sigma2_joint), log = TRUE)
  }))
  
  for (s in 1:q) {
    # Contribution of the observed data to the joint density
    data_s <- X_complete[[s,1]]
    like <- like + sum(sapply(1:n, function(i) {
      dnorm(data_s[,i], mean = (U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]))[,i], sd = sqrt(error_vars[s]), log = TRUE)
    }))
    
    # Contribution of Us to the joint density
    like <- like + sum(sapply(1:r, function(rs) {
      dnorm(U.iter[[s,1]][,rs], mean = 0, sd = sqrt(sigma2_joint), log = TRUE)
    }))
    
    # Contribution of Ws to the joint density
    like <- like + sum(sapply(1:r.vec[s], function(rs) {
      dnorm(W.iter[[s,s]][,rs], mean = 0, sd = sqrt(sigma2_indiv[s]), log = TRUE)
    }))
    
    # Contribution of Vs to the joint density
    like <- like + sum(sapply(1:r.vec[s], function(rs) {
      dnorm(Vs.iter[[1,s]][,rs], mean = 0, sd = sqrt(sigma2_indiv[s]), log = TRUE)
    }))
    
    # If there is a response
    if (response_given) {
      VStar.iter <- cbind(1, do.call(cbind, V.iter), do.call(cbind, Vs.iter))
      
      # The contribution of beta to the joint density
      if (!sparsity) {
        like <- like + sum(sapply(1:n_beta, function(rs) {
          dnorm(beta.iter[[1,1]][rs,], mean = 0, sd = sqrt(Sigma_beta[rs,rs]), log = TRUE)
        }))
      }
      
      if (sparsity) {
        like <- like + sum(sapply(1:n_beta, function(rs) {
          # Update the var-covar of beta given the current value of gamma
          if (gamma.iter[[1,1]][rs,] == 1) {
            var_beta <- Sigma_beta[rs,rs]
          } 
          if (gamma.iter[[1,1]][rs,] == 0) {
            var_beta <- 1/1000
          }
          
          dnorm(beta.iter[[1,1]][rs,], mean = 0, sd = sqrt(var_beta), log = TRUE)
        }))
      }
      
      if (response_type == "continuous") {
        # The contribution of the observed response to the joint density
        like <- like + sum(sapply(1:n, function(i) {
          dnorm(Y_complete[[1,1]][i,], mean = (VStar.iter %*% beta.iter[[1,1]])[i,], sd = sqrt(tau2.iter[[1,1]]), log = TRUE)
        }))
        
        # The contribution of tau2 to the joint density
        like <- like + log(dinvgamma(tau2.iter[[1,1]], shape = shape, scale = 1/rate))
      }
      
      if (response_type == "binary") {
        # The contribution of the observed response to the joint density
        like <- like + sum(log(sapply(1:n, function(i) {
          dbinom(Y_complete[[1,1]][i,], size = 1, prob = pnorm(VStar.iter %*% beta.iter[[1,1]]))
        })))
      }
      
      # The contribution of the gammas and the prior probability of inclusion, p, to the likelihood
      if (sparsity) {
        # The gammas
        like <- like + sum(sapply(2:n_beta, function(rs) { # Ignore the intercept, which is always included
          dbinom(gamma.iter[[1,1]][rs,], size = 1, prob = p.iter[[1,1]], log = TRUE)
        }))
        
        # Prior inclusion probability, p
        like <- like + dbeta(p.iter[[1,1]], shape1 = 1, shape2 = 1, log = TRUE)
      }
    }
  }

  # Return
  like
}

# -----------------------------------------------------------------------------
# Model comparison simulations
# -----------------------------------------------------------------------------

# Calculate the predicted Y
Y_predicted <- function(scores.draw, beta.draw, nsample) {
  # Calculates the predicted Y value at each iteration of the Gibbs sampler
  
  # ---------------------------------------------------------------------------
  # Arguments:
  #
  # scores.draw = VStar (with intercept) at each Gibbs sampling iteration
  # beta.draw = coefficients at each Gibbs sampling iteration
  # nsample = number of Gibbs sampling iterations
  # ---------------------------------------------------------------------------
  
  # Is Y continuous or binary?
  response_type <- if (all(unique(Y) %in% c(0, 1, NA))) "binary" else "continuous"
  
  # If Y is continuous:
  if (response_type == "continuous") {
    pred_y <- lapply(1:nsample, function(iter) {
      scores.draw[[iter]] %*% beta.draw[[iter]][[1,1]]
    })
  }
  
  # If Y is binary:
  if (response_type == "binary") {
    pred_y <- sapply(1:nsample, function(iter) {
      pnorm(scores.draw[[iter]] %*% beta.draw[[iter]][[1,1]])
    })
  }
  
  # Return
  pred_y
}

# Run each model being compared in model comparison simulation study
model_comparison <- function(mod, p.vec, n, ranks, response, true_params, model_params, s2nX, s2nY, sparsity, nsim, nsample = 2000, n_clust, estim_ranks = TRUE) {
  
  # ---------------------------------------------------------------------------
  # Arguments:
  #
  # mod = string in c("sJIVE", "BIDIFAC+", "JIVE", "MOFA", "BPMF")
  # p.vec = number of features per source
  # n = sample size
  # ranks = vector of joint and individual ranks = c(joint rank, indiv rank 1, indiv rank 2, ...)
  # response = string in c(NULL, "continuous", "binary")
  # true_params = the list of true parameters under which to generate data
  # s2nX = signal-to-noise ratio in the data, X.
  # s2nY = signal-to-noise ratio in the response, Y. NULL if response == "binary"
  # nsim = number of simulations to run
  # nsample = number of Gibbs sampling iterations to draw for the linear model
  # n_clust = how many clusters to run simulation in parallel?
  # estim_ranks = should the ranks be estimated? Default is TRUE
  # ---------------------------------------------------------------------------
  
  # Loading in the packages
  # library(r.jive)
  # library(MOFA2)
  library(doParallel)
  library(foreach)
  # library(sup.r.jive)
  # library(natural)
  # library(RSpectra)
  library(BIPnet)
  
  # The model options
  models <- c("sJIVE", "BIDIFAC", "JIVE", "MOFA", "BPMF_Data_Mode", "BIP", "BPMF_test", "BPMF_test_scale")
  
  cl <- makeCluster(n_clust)
  registerDoParallel(cl)
  funcs <- c("bpmf_data", "center_data", "bpmf_full_mode", "bpmf_data_mode", "bpmf_test", "bpmf_test_scale", "get_results", "BIDIFAC",
             "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
             "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum",
             "bidifac.plus.impute", "bidifac.plus.given")
  packs <- c("Matrix", "MASS", "truncnorm", "BIPnet")
             # , "r.jive", "sup.r.jive", "natural", "RSpectra", "MOFA2")
  sim_results <- foreach (sim_iter = 1:nsim, .packages = packs, .export = funcs, .verbose = TRUE, .combine = rbind) %dopar% {
    # Set seed
    if (mod == "test") {
      seed <- sim_iter + which(models %in% "BIDIFAC+") * nsim
    }
    
    if (mod != "test") {
      seed <- sim_iter + which(models %in% mod) * nsim
    }
    set.seed(seed)
    
    # -------------------------------------------------------------------------
    # Generate data 
    # -------------------------------------------------------------------------
    
    # Generate 2*n samples to split into equally-sized training and test datasets
    sim_data <- bpmf_data(p.vec, 2*n, ranks, true_params, s2nX, s2nY, response, missingness = NULL, entrywise = NULL, prop_missing = NULL, sparsity = sparsity)
    
    # Saving the data
    true_data <- sim_data$data
    true_data_list <- lapply(true_data, function(s) s)
    
    # Saving the underlying structure
    joint.structure <- sim_data$joint.structure
    indiv.structure <- sim_data$indiv.structure
    EY <- sim_data$EY
    
    # Save the number of sources (not including Y)
    q <- length(p.vec)
    
    # The response
    Y <- sim_data$Y
    
    # The standardizing coefficients
    s2nX_coef <- sim_data$s2nX_coef
    s2nY_coef <- sim_data$s2nY_coef
    
    # The response parameters
    beta <- sim_data$beta
    
    # Save a burn-in
    burnin <- nsample/2
    
    # Save the error variance estimates as NULL
    sigma.mat <- matrix(nrow = q+1, ncol = 1)
    
    # Set the indices of the sources
    p.ind <- lapply(1:q, function(s) {
      if (s == 1) {
        1:p.vec[s]
      } else {
        (p.vec[s-1] + 1):cumsum(p.vec)[s]
      }
    })
    p.ind[[q+1]] <- cumsum(p.vec)[q] + 1 # For the response
    
    # Set the indices of the samples
    n.ind <- list(1:n)
    
    # -------------------------------------------------------------------------
    # Split the data into a training and test set
    # -------------------------------------------------------------------------
    
    # Training data
    training_data <- matrix(list(), nrow = q, ncol = 1)
    training_data_list <- lapply(1:q, function(s) list())
    
    joint.structure_train <- matrix(list(), nrow = q, ncol = 1)
    indiv.structure_train <- matrix(list(), nrow = q, ncol = 1)
    
    # Test data
    test_data <- matrix(list(), nrow = q, ncol = 1)
    test_data_list <- lapply(1:q, function(s) list())
    
    joint.structure_test <- matrix(list(), nrow = q, ncol = 1)
    indiv.structure_test <- matrix(list(), nrow = q, ncol = 1)
    
    for (s in 1:q) {
      training_data[[s,1]] <- true_data[[s,1]][,1:n]
      training_data_list[[s]] <- training_data[[s,1]]
      
      joint.structure_train[[s,1]] <- joint.structure[[s,1]][,1:n]
      indiv.structure_train[[s,1]] <- indiv.structure[[s,1]][,1:n]
      
      test_data[[s,1]] <- true_data[[s,1]][,(n+1):(2*n)]
      test_data_list[[s]] <- test_data[[s,1]]
      
      joint.structure_test[[s,1]] <- joint.structure[[s,1]][,(n+1):(2*n)]
      indiv.structure_test[[s,1]] <- indiv.structure[[s,1]][,(n+1):(2*n)]
    }
    
    Y_train <- matrix(list(), nrow = 1, ncol = 1)
    Y_train[[1,1]] <- Y[[1,1]][1:n,, drop = FALSE]
    
    EY_train <- matrix(list(), nrow = 1, ncol = 1)
    EY_train[[1,1]] <- EY[[1,1]][1:n,, drop = FALSE]
    
    Y_test <- matrix(list(), nrow = 1, ncol = 1)
    Y_test[[1,1]] <- Y[[1,1]][(n+1):(2*n),, drop = FALSE]
    
    EY_test <- matrix(list(), nrow = 1, ncol = 1)
    EY_test[[1,1]] <- EY[[1,1]][(n+1):(2*n),, drop = FALSE]
    
    # -------------------------------------------------------------------------
    # Center X and y to have mean 0
    # -------------------------------------------------------------------------
    
    for (s in 1:q) {
      # Center the full true data (transpose to center the cols (features) than transpose back)
      true_data[[s,1]] <- t(scale(t(true_data[[s,1]]), center = TRUE, scale = FALSE))
      true_data_list[[s]] <- t(scale(t(true_data_list[[s]]), center = TRUE, scale = FALSE))
      
      # Save the means
      row_means_true <- attr(true_data[[s,1]], "scaled:center")
      
      # Subtract the row means from the structure
      joint.structure[[s,1]] <- sweep(joint.structure[[s,1]], 1, row_means_true)
      indiv.structure[[s,1]] <- sweep(indiv.structure[[s,1]], 1, row_means_true)
      
      # Center the observed training data 
      training_data[[s,1]] <- t(scale(t(training_data[[s,1]]), center = TRUE, scale = FALSE))
      training_data_list[[s]] <- t(scale(t(training_data_list[[s]]), center = TRUE, scale = FALSE))

      # Save the means 
      row_means_train <- attr(training_data[[s,1]], "scaled:center")

      # Subtract each row by the same mean
      joint.structure_train[[s,1]] <- sweep(joint.structure_train[[s,1]], 1, row_means_train)
      indiv.structure_train[[s,1]] <- sweep(indiv.structure_train[[s,1]], 1, row_means_train)

      # Center the observed test data
      test_data[[s,1]] <- t(scale(t(test_data[[s,1]]), center = TRUE, scale = FALSE))
      test_data_list[[s]] <- t(scale(t(test_data_list[[s]]), center = TRUE, scale = FALSE))

      # Save the means and variances
      row_means_test <- attr(test_data[[s,1]], "scaled:center")

      # Subtract each column by the same mean and divide by the same sd above
      joint.structure_test[[s,1]] <- sweep(joint.structure_test[[s,1]], 1, row_means_test)
      indiv.structure_test[[s,1]] <- sweep(indiv.structure_test[[s,1]], 1, row_means_test)
    }

    Y_train[[1,1]] <- scale(Y_train[[1,1]], center = TRUE, scale = FALSE)
    EY_train[[1,1]] <- (EY_train[[1,1]] - attr(Y_train[[1,1]], "scaled:center"))

    Y_test[[1,1]] <- scale(Y_test[[1,1]], center = TRUE, scale = FALSE)
    EY_test[[1,1]] <- (EY_test[[1,1]] - attr(Y_test[[1,1]], "scaled:center"))

    Y[[1,1]] <- scale(Y[[1,1]], center = TRUE, scale = FALSE)
    EY[[1,1]] <- (EY[[1,1]] - attr(Y[[1,1]], "scaled:center"))
    
    # -------------------------------------------------------------------------
    # Fit each model on generated data to obtain estimate of underlying structure
    # -------------------------------------------------------------------------
    
    if (mod == "sJIVE") {
      # Setting the number of tuning parameters to compare
      eta <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
      
      # Running sJIVE with rank estimation
      if (estim_ranks) {
        mod.out <- sJIVE(X = training_data_list, Y = c(Y_train[[1,1]]), eta = eta, rankA = NULL, rankJ = NULL, method = "permute", threshold = 0.001, center.scale = FALSE, reduce.dim = TRUE)
      }
      
      # Running sJIVE with fixed ranks
      if (!estim_ranks) {
        mod.out <- sJIVE(X = training_data_list, Y = c(Y_train[[1,1]]), eta = eta, rankA = ranks[-1], rankJ = ranks[1], threshold = 0.001, center.scale = FALSE, reduce.dim = TRUE)
      }

      # Saving the joint structure
      mod.joint <- lapply(1:q, function(s) {
        t(t(mod.out$U_I[[s]])) %*% mod.out$S_J
      })
      
      # Saving the individual structure
      mod.individual <- lapply(1:q, function(s) {
        mod.out$W_I[[s]] %*% mod.out$S_I[[s]]
      })
      
      # Saving the joint rank
      joint.rank <- mod.out$rankJ
      
      # Saving the individual ranks
      indiv.rank <- mod.out$rankA
      
      # Saving the joint scores
      if (joint.rank != 0) {
        joint.scores <- t(mod.out$S_J)
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Saving the individual scores
      indiv.scores <- lapply(1:q, function(s) {
        if (indiv.rank[s] != 0) {
          t(mod.out$S_I[[s]])
        }
      })
      
      # Saving the E(Y)
      EY.fit <- t(mod.out$fittedY)
      
      # Save the ranks
      mod.ranks <- c(joint.rank, indiv.rank)
      
      # Combining all scores together
      all.scores <- cbind(Y_train[[1,1]], joint.scores, do.call(cbind, indiv.scores))
      colnames(all.scores) <- c("y", rep("joint", joint.rank), rep("indiv", sum(indiv.rank)))
      
      # Do not compute results for coverage
      coverage_EY_train <- coverage_EY_test <- NA
      coverage_Y_train <- coverage_Y_test <- NA
      
      # Do not compute results for test data
      joint.recovery.structure.test <- indiv.recovery.structure.test <- mse_EY_test <- NA
    }
    
    if (mod == "BIDIFAC") {
      # Run model
      mod.out <- BIDIFAC(true_data, rmt = TRUE, pbar = FALSE, scale_back = TRUE)
      
      # Saving the column structure (the joint structure)
      mod.joint <- lapply(1:q, function(s) {
        mod.out$C[[s,1]]
      })
      
      # Saving the individual structure 
      mod.individual <- lapply(1:q, function(s) {
        mod.out$I[[s,1]]
      })
      
      # Saving the joint rank
      joint.rank <- rankMatrix(mod.out$C[[1,1]])[1]
      
      # Saving the individual ranks
      indiv.rank <- sapply(mod.out$I, function(s) {
        rankMatrix(s)[1]
      }) 
      
      # Obtaining the joint scores
      if (joint.rank != 0)  {
        svd.joint <- svd(mod.joint[[1]])
        joint.scores <- (svd.joint$v[,1:joint.rank,drop=FALSE]) %*% diag(svd.joint$d[1:joint.rank], nrow = joint.rank)
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Obtaining the individual scores
      indiv.scores <- lapply(1:q, function(s) {
        if (indiv.rank[s] != 0) {
          svd.source <- svd(mod.individual[[s]])
          (svd.source$v[,1:indiv.rank[s],drop=FALSE]) %*% diag(svd.source$d[1:indiv.rank[s]], nrow = indiv.rank[s])
        }
      })
    }
    
    if (mod == "BIDIFAC+") {
      # Scale the data to have error variance 1
      scaled_true_data <- true_data
      sigma.mat <- matrix(nrow = q+1, ncol = 1)
      
      # Iterate through the sources and estimate the error variance
      for (s in 1:q) {
        
        # Save the estimated error variance
        sigma.mat[s,] <- sigma.rmt(true_data[[s,1]])
        
        # Scale the data
        scaled_true_data[[s,1]] <- true_data[[s,1]]/sigma.mat[s,]
      }
      
      # Setting the test response to 0 while estimating the error variance
      Y_NA_for_test <- Y
      Y_NA_for_test[[1,1]][(n+1):(2*n),] <- 0
      
      # Scaling the response
      sigma.mat[q+1,] <- sigma.rmt(Y_NA_for_test[[1,1]]) 
      Y_NA_for_test[[1,1]] <- Y_NA_for_test[[1,1]]/sigma.mat[q+1,]
      
      # Now set the test Y to missing
      Y_NA_for_test[[1,1]][(n+1):(2*n),] <- NA
      
      # Combine the data into one matrix
      true_data_y_combined <- rbind(do.call(rbind, scaled_true_data), t(Y_NA_for_test[[1,1]]))
      
      # Set the indices of the structures to be estimated
      p.ind.list <- list(c(unlist(p.ind))) # Joint structure
      
      for (s in 1:q) {
        p.ind.list[[s+1]] <- c(p.ind[[s]], p.ind[[q+1]])
      }
      
      # Save the indices for samples in each identified structure
      n.ind.list <- lapply(1:(q+1), function(s) c(1:(2*n)))
      
      # Run model
      mod.out <- bidifac.plus.impute(X0 = true_data_y_combined, p.ind = p.ind, n.ind = n.ind,
                                     p.ind.list = p.ind.list, n.ind.list = n.ind.list,
                                     all.miss = which(is.na(true_data_y_combined)))
      
      # Saving the joint structure
      mod.joint <- lapply(1:q, function(s) {
        mod.out$S[[1]][p.ind[[s]],]
      })
      
      # Saving the shared structure with y
      mod.individual <- lapply(1:q, function(s) {
        mod.out$S[[s+1]][p.ind[[s]],]
      })
      
      # Saving the joint rank
      joint.rank <- rankMatrix(mod.out$S[[1]])[1]
      
      # Saving the shared rank structure with y
      indiv.rank <- sapply(1:q, function(s) {
        rankMatrix(mod.out$S[[q+1]])[1]
      }) 
      
      # Obtaining the joint scores and predicted E(Y) from joint
      if (joint.rank != 0)  {
        # Obtaining the joint scores
        svd.joint <- svd(mod.joint[[1]])
        joint.scores <- (svd.joint$v[,1:joint.rank,drop=FALSE]) %*% diag(svd.joint$d[1:joint.rank], nrow = joint.rank)
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Obtaining the individual scores
      indiv.scores <- lapply(1:q, function(s) {
        if (indiv.rank[s] != 0) {
          svd.source <- svd(mod.individual[[s]])
          (svd.source$v[,1:indiv.rank[s]]) %*% diag(svd.source$d[1:indiv.rank[s]], nrow = indiv.rank[s])
        }
      })
      
      # Saving the E(Y) from the joint structure
      joint.EY <- t(mod.out$S[[1]][p.ind[[q+1]],,drop=FALSE])
      indiv.EY <- lapply(1:q, function(s) {
        t(mod.out$S[[s+1]][p.ind[[q+1]],,drop=FALSE])
      })
    
      Y.fit <- (joint.EY + Reduce("+", indiv.EY)) * sigma.mat[q+1,]
      
      # Do not compute results for coverage
      coverage_EY_train <- coverage_EY_test <- NA
    }
    
    if (mod == "JIVE") {
      # Running JIVE with rank estimation
      if (estim_ranks) {
        mod.out <- jive(true_data_list, center = FALSE, scale = FALSE, method = "perm")
      }
      
      # Running JIVE with fixed ranks
      if (!estim_ranks) {
        mod.out <- jive(true_data_list, rankJ = ranks[1], rankA = ranks[-1], center = FALSE, scale = FALSE, method = "given")
      }
      
      # Saving the joint structure
      mod.joint <- mod.out$joint
      
      # Saving the individual structure
      mod.individual <- mod.out$individual
      
      # Saving the joint rank
      joint.rank <- mod.out$rankJ
      
      # Saving the individual ranks
      indiv.rank <- mod.out$rankA
      
      # Obtaining the joint scores
      if (joint.rank != 0) { 
        svd.joint <- svd(mod.joint[[1]])
        joint.scores <- (svd.joint$v[,1:joint.rank,drop=FALSE]) %*% diag(svd.joint$d[1:joint.rank], nrow = joint.rank)
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Obtaining the individual scores
      indiv.scores <- lapply(1:q, function(s) {
        if (indiv.rank[s] != 0) {
          svd.source <- svd(mod.individual[[s]])
          (svd.source$v[,1:indiv.rank[s],drop=FALSE]) %*% diag(svd.source$d[1:indiv.rank[s]], nrow = indiv.rank[s])
        }
      })
    }
    
    if (mod == "MOFA") {
      # Create an untrained MOFA object
      mofa_pre_train <- create_mofa(true_data_list)
      
      # Set the data options so that the data is not additionally centered and fix the ranks if desired
      data_opts <- get_default_data_options(mofa_pre_train)
      model_opts <- get_default_model_options(mofa_pre_train)
      
      # No centering
      data_opts$center_groups <- FALSE
      
      # If using fixed ranks
      if (!estim_ranks) {
        model_opts$num_factors <- sum(ranks)
      }
      
      # model_opts$spikeslab_factors <- TRUE
      # model_opts$ard_factors <- TRUE
      
      # Create the MOFA object
      MOFAobject <- prepare_mofa(
        object = mofa_pre_train,
        model_options = model_opts
      )
      
      # Train the MOFA model
      mod.out <- run_mofa(MOFAobject)
      
      # Getting the variance explained in each source by each factor
      mod.var.exp <- get_variance_explained(mod.out)$r2_per_factor$group1 # variance explained by factor per view
      
      # OLD --
      
      # Set threshold for determining factor inclusion for joint/individual structures
      # threshold <- 0.1
      
      # Determining for which views each factor is active
      # joint.factors <- which(apply(mod.var.exp, 1, function(factor) all(factor > threshold)))
      # indiv.factors <- lapply(1:q, function(s) {
      #   which(apply(mod.var.exp, 1, function(factor) factor[s] > threshold & factor[c(1:q)[!c(1:q) %in% s]] < threshold))
      # })
      
      # OLD -- 
      
      # Save which views each factor applies to
      joint_or_individual <- lapply(1:nrow(mod.var.exp), function(factor) {
        # Save variance explained for current factor
        row <- mod.var.exp[factor,]
        
        ind.max <- which.max(row) # Highest var explained
        ind.min <- which.min(row) # Minimum var explained
        
        # First, check that factor explains at least 1% of variation in each source
        greater_than_1p <- any(row > 1)
        
        # If explains more than 1% variation in at least one source
        if (greater_than_1p) {
          
          # If factor explains substantial (max var is less than 2 * min var) variation in both sources, it is joint
          if (row[ind.max] < 2*row[-ind.max]) { 
            c(1:q)
          } else { # If factor explains substantial variation in just one source, it is individual
            ind.max
          }
    
        }
        
      })
      
      # Save the joint factors
      joint.factors <- c(1:length(joint_or_individual))[sapply(joint_or_individual, function(factor) length(factor) == q)]
      
      # Save the individual factors
      indiv.factors <- lapply(1:q, function(s) {
        c(1:length(joint_or_individual))[sapply(joint_or_individual, function(factor) all(factor %in% s) & length(factor) > 0)]
      })
      
      # Saving the joint rank
      joint.rank <- length(joint.factors)
      
      # Saving the individual rank
      indiv.rank <- sapply(indiv.factors, length)
      
      # Save the underlying structure
      mod.joint <- lapply(1:q, function(s) {
        t(mod.out@expectations$Z$group1[, joint.factors, drop = FALSE] %*% t(mod.out@expectations$W[[s]][, joint.factors, drop = FALSE]))
      })
      mod.individual <- lapply(1:q, function(s) {
        t(mod.out@expectations$Z$group1[, indiv.factors[[s]], drop = FALSE] %*% t(mod.out@expectations$W[[s]][, indiv.factors[[s]], drop = FALSE]))
      })
      
      # Save the MOFA scores (all.equal(get_factors(mod.out)$group1, mod.out@expectations$Z$group1[, joint.factors, drop = FALSE]) # TRUE!)
      mofa.scores <- get_factors(mod.out)$group1
      
      # Save the joint scores
      if (joint.rank != 0) {
        joint.scores <- mofa.scores[,joint.factors, drop = FALSE]
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Save the individual scores
      indiv.scores <- lapply(1:q, function(s) {
        if (indiv.rank[s] != 0) {
          indiv.scores <- mofa.scores[,unlist(indiv.factors[[s]]), drop = FALSE]
        }
      })
    }
    
    if (mod == "BPMF_Data_Mode") {
      # Setting the test response to NA
      Y_NA_for_test <- Y
      Y_NA_for_test[[1,1]][(n+1):(2*n),] <- NA
      
      # Running BPMF with rank estimation
      if (estim_ranks) {
        mod.out <- bpmf_data_mode(data = true_data, Y = Y_NA_for_test, nninit = TRUE, model_params = model_params, nsample = nsample)
      }
      
      # Running BPMF with fixed ranks
      if (!estim_ranks) {
        mod.out <- bpmf_data_mode(data = true_data, Y = Y_NA_for_test, ranks = ranks, nninit = FALSE, model_params = model_params, nsample = nsample)
      }
      
      
      # Saving the joint structure
      mod.joint.iter <- lapply(1:nsample, function(iter) {
        lapply(1:q, function(s) {
          mod.out$J.draw[[iter]][[s,1]]
        })
      })
      
      # Taking the posterior mean
      mod.joint <- lapply(1:q, function(s) {
        # Save the joint structure at each iteration for source
        joint.source <- lapply((burnin+1):nsample, function(iter) {
          mod.joint.iter[[iter]][[s]]
        })
        # Take the mean
        Reduce("+", joint.source)/length(joint.source)
      })
      
      # Saving the individual structure
      mod.individual.iter <- lapply(1:nsample, function(iter) {
        lapply(1:q, function(s) {
          mod.out$A.draw[[iter]][[s,1]]
        })
      })
      
      # Taking the posterior mean
      mod.individual <- lapply(1:q, function(s) {
        # Save the joint structure at each iteration for source
        indiv.source <- lapply((burnin+1):nsample, function(iter) {
          mod.individual.iter[[iter]][[s]]
        })
        # Take the mean
        Reduce("+", indiv.source)/length(indiv.source)
      })
      
      # Saving the joint rank
      joint.rank <- mod.out$ranks[1]
      
      # Saving the individual ranks
      indiv.rank <- mod.out$ranks[-1]
      
      # Saving the joint scores
      if (joint.rank != 0) {
        joint.scores <- Reduce("+", lapply((burnin+1):nsample, function(iter) {
          mod.out$V.draw[[iter]][[1,1]]
        }))/(nsample-burnin)
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Saving the individual scores
      indiv.scores <- lapply(1:q, function(source) {
        if (indiv.rank[source] != 0) {
          indiv.scores <- Reduce("+", lapply((burnin+1):nsample, function(iter) {
            mod.out$Vs.draw[[iter]][[1,source]]
          }))/(nsample-burnin)
        }
      })
      
      # Save the estimated error standard deviation for X and y
      sigma.mat <- matrix(nrow = 3, ncol = 1)
      sigma.mat[1:2,] <- mod.out$sigma.mat
      sigma.mat[3,] <- mean(sqrt(unlist(mod.out$tau2.draw))) 
    }
    
    if (mod == "BIP") {
      # Transpose the data sets so they are in nxp orientation
      training_data_list_bip <- lapply(training_data_list, function(s) t(s))
      
      # Add the response to the training data list
      training_data_list_bip[[3]] <- Y_train[[1,1]]
      
      # If estimating the ranks, set to a very high value 
      if (estim_ranks) {
        bip_num_ranks <- 10
      }
      
      # If not estimating the ranks, fix at the true ranks
      if (!estim_ranks) {
        bip_num_ranks <- sum(ranks)
      }
      
      # Running the model
      mod.out <- BIP(dataList = training_data_list_bip, IndicVar = c(0,0,1), Method = "BIP", sample = 5000, burnin = 2500, nbrcomp = bip_num_ranks)
      
      # Saving the results
      bip.scores <- mod.out$EstU
      bip.loadings <- mod.out$EstLoad
      
      # Saving the feature sds to scale the estimated structure
      bip.est.sd <- mod.out$SDData
      
      # Save the component selection
      factor_mpp_by_source <- mod.out$CompoSelMean[1:q,]
      factor_mpp_by_source_list <- split(factor_mpp_by_source, rep(1:ncol(factor_mpp_by_source), each = nrow(factor_mpp_by_source)))
      bip.joint.individual <- lapply(factor_mpp_by_source_list, function(comp) {
        # Save lowest MPP
        ind.min <- which.min(comp)
        
        # Save the highest MPP
        ind.max <- which.max(comp)
        
        # If the smallest MPP is still greater than 0.5 then it is joint
        if (comp[ind.min] >= 0.5) {
          1:q
        } 
        
        # If the smallest MPP is less than 0.5, check the max
        else if (comp[ind.min] < 0.5) {
          if (comp[ind.max] > 0.5) {
            ind.max
          } 
        }
      })
      
      # Save the joint factor indices
      joint.factors <- sapply(bip.joint.individual, function(comp) {
        length(comp) == q
      })
      
      # Save the individual factor indices
      indiv.factors <- lapply(1:q, function(s) {
        sapply(bip.joint.individual, function(comp) {
          all(comp %in% s) & length(comp) > 0
        })
      })
      
      # Save the joint rank
      joint.rank <- sum(joint.factors)
      
      # Save the individual rank
      indiv.rank <- sapply(1:q, function(s) {
        sum(indiv.factors[[s]])
      })
      
      # Save the ranks
      mod.ranks <- c(joint.rank, indiv.rank)
      
      # Save the posterior mean of the joint structure
      mod.joint <- lapply(1:q, function(s) {
        t(bip.scores[,joint.factors,drop=FALSE] %*% bip.loadings[[s]][joint.factors,,drop=FALSE] %*% diag(bip.est.sd[[s]]))
      })
      
      # Save the posterior mean of the individual structure
      mod.individual <- lapply(1:q, function(s) {
        t(bip.scores[,indiv.factors[[s]],drop=FALSE] %*% bip.loadings[[s]][indiv.factors[[s]],,drop=FALSE] %*% diag(bip.est.sd[[s]]))
      })
      
      # Save the predicted Y
      comp.inc.y <- sapply(mod.out$CompoSelMean[q+1,], function(comp) comp > 0.5)
      EY.fit <- mod.out$EstIntcp + bip.scores[,comp.inc.y,drop=FALSE] %*% bip.loadings[[q+1]][comp.inc.y,,drop=FALSE]
      
      # Save the estimated error sd for Y
      sigma.mat <- matrix(nrow = q+1, ncol = 1)
      sigma.mat[q+1,] <- sqrt(mod.out$EstSig2[[q+1]])
      
      # Do not compute results for coverage
      coverage_EY_train <- coverage_EY_test <- NA
      coverage_Y_train <- coverage_Y_test <- NA
      
      # Do not compute results for test data
      joint.recovery.structure.test <- indiv.recovery.structure.test <- mse_EY_test <- NA
    }
    
    if (mod == "BPMF_test") {
      # Setting the test response to NA
      Y_NA_for_test <- Y
      Y_NA_for_test[[1,1]][(n+1):(2*n),] <- NA
      
      # Running BPMF
      mod.out <- bpmf_test(data = true_data, Y = Y_NA_for_test, nninit = TRUE, model_params = model_params, nsample = nsample)
      
      # Saving the joint structure
      mod.joint.iter <- lapply(1:nsample, function(iter) {
        lapply(1:q, function(s) {
          mod.out$J.draw[[iter]][[s,1]]
        })
      })
      
      # Taking the posterior mean
      mod.joint <- lapply(1:q, function(s) {
        # Save the joint structure at each iteration for source
        joint.source <- lapply((burnin+1):nsample, function(iter) {
          mod.joint.iter[[iter]][[s]]
        })
        # Take the mean
        Reduce("+", joint.source)/length(joint.source)
      })
      
      # Saving the individual structure
      mod.individual.iter <- lapply(1:nsample, function(iter) {
        lapply(1:q, function(s) {
          mod.out$A.draw[[iter]][[s,1]]
        })
      })
      
      # Taking the posterior mean
      mod.individual <- lapply(1:q, function(s) {
        # Save the joint structure at each iteration for source
        indiv.source <- lapply((burnin+1):nsample, function(iter) {
          mod.individual.iter[[iter]][[s]]
        })
        # Take the mean
        Reduce("+", indiv.source)/length(indiv.source)
      })
      
      # Saving the joint rank
      joint.rank <- mod.out$ranks[1]
      
      # Saving the individual ranks
      indiv.rank <- mod.out$ranks[-1]
      
      # Saving the joint scores
      if (joint.rank != 0) {
        joint.scores <- Reduce("+", lapply((burnin+1):nsample, function(iter) {
          mod.out$V.draw[[iter]][[1,1]]
        }))/(nsample-burnin)
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Saving the individual scores
      indiv.scores <- lapply(1:q, function(source) {
        if (indiv.rank[source] != 0) {
          indiv.scores <- Reduce("+", lapply((burnin+1):nsample, function(iter) {
            mod.out$Vs.draw[[iter]][[1,source]]
          }))/(nsample-burnin)
        }
      })
      
      # Save the error standard deviation estimates
      sigma.mat <- mod.out$sigma.mat
    }
    
    if (mod == "BPMF_test_scale") {
      # Setting the test response to NA
      Y_NA_for_test <- Y
      Y_NA_for_test[[1,1]][(n+1):(2*n),] <- NA
      
      # Running BPMF
      mod.out <- bpmf_test_scale(data = true_data, Y = Y_NA_for_test, nninit = TRUE, model_params = model_params, nsample = nsample)
      
      # Saving the joint structure
      mod.joint.iter <- lapply(1:nsample, function(iter) {
        lapply(1:q, function(s) {
          mod.out$J.draw[[iter]][[s,1]]
        })
      })
      
      # Taking the posterior mean
      mod.joint <- lapply(1:q, function(s) {
        # Save the joint structure at each iteration for source
        joint.source <- lapply((burnin+1):nsample, function(iter) {
          mod.joint.iter[[iter]][[s]]
        })
        # Take the mean
        Reduce("+", joint.source)/length(joint.source)
      })
      
      # Saving the individual structure
      mod.individual.iter <- lapply(1:nsample, function(iter) {
        lapply(1:q, function(s) {
          mod.out$A.draw[[iter]][[s,1]]
        })
      })
      
      # Taking the posterior mean
      mod.individual <- lapply(1:q, function(s) {
        # Save the joint structure at each iteration for source
        indiv.source <- lapply((burnin+1):nsample, function(iter) {
          mod.individual.iter[[iter]][[s]]
        })
        # Take the mean
        Reduce("+", indiv.source)/length(indiv.source)
      })
      
      # Saving the joint rank
      joint.rank <- mod.out$ranks[1]
      
      # Saving the individual ranks
      indiv.rank <- mod.out$ranks[-1]
      
      # Saving the joint scores
      if (joint.rank != 0) {
        joint.scores <- Reduce("+", lapply((burnin+1):nsample, function(iter) {
          mod.out$V.draw[[iter]][[1,1]]
        }))/(nsample-burnin)
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Saving the individual scores
      indiv.scores <- lapply(1:q, function(source) {
        if (indiv.rank[source] != 0) {
          indiv.scores <- Reduce("+", lapply((burnin+1):nsample, function(iter) {
            mod.out$Vs.draw[[iter]][[1,source]]
          }))/(nsample-burnin)
        }
      })
      
      # Save the error standard deviation estimates
      sigma.mat <- mod.out$sigma.mat
    }
    
    if (mod != "test" & mod != "sJIVE" & mod != "BIP") {
      # Combining the ranks
      mod.ranks <- c(joint.rank, indiv.rank)
      
      # Combining all scores together
      all.scores <- cbind(Y[[1,1]], joint.scores, do.call(cbind, indiv.scores))
      colnames(all.scores) <- c("y", rep("joint", joint.rank), rep("indiv", sum(indiv.rank)))
    }
    
    # -------------------------------------------------------------------------
    # As applicable, use structure in Bayesian linear model
    # -------------------------------------------------------------------------
    
    if (mod %in% c("BIDIFAC", "JIVE", "MOFA")) {
      # Subset the scores to just the training data
      all.scores.train <- all.scores[1:n,,drop=FALSE]
      
      # Fitting the Bayesian linear model
      mod.bayes <- bpmf_data_mode(data = training_data, Y = Y_train, nninit = FALSE, model_params = model_params, 
                                  ranks = mod.ranks, scores = all.scores.train[,-1], nsample = nsample)
      
      # Calculate the predicted E(Y) at each Gibbs sampling iteration using training and testing scores
      EY.fit.iter <- lapply((burnin+1):nsample, function(iter) {
        VStar.iter <- cbind(1, all.scores[,-1])
        beta.iter <- mod.bayes$beta.draw[[iter]][[1,1]] # Beta depends only on training data
        VStar.iter %*% beta.iter
      })
      
      EY.fit <- Reduce("+", EY.fit.iter)/length(EY.fit.iter)
      
      # Calculate the coverage of training E(Y) and test E(Y)
      EY.fit.iter <- do.call(cbind, EY.fit.iter)
      ci_by_EY <- apply(EY.fit.iter, 1, function(subj) c(quantile(subj, 0.025), quantile(subj, 0.975)))
      
      coverage_EY_train <- mean(sapply(1:n, function(i) {
        (EY[[1,1]][i,] >= ci_by_EY[1,i]) & (EY[[1,1]][i,] <= ci_by_EY[2,i])
      }))
      
      coverage_EY_test <- mean(sapply((n+1):(2*n), function(i) {
        (EY[[1,1]][i,] >= ci_by_EY[1,i]) & (EY[[1,1]][i,] <= ci_by_EY[2,i])
      }))
      
      # Calculate the predicted Y at each Gibbs sampling iteration using training and testing scores
      Y.fit.iter <- lapply((burnin+1):nsample, function(iter) {
        VStar.iter <- cbind(1, all.scores[,-1])
        beta.iter <- mod.bayes$beta.draw[[iter]][[1,1]] # Beta depends only on training data
        VStar.iter %*% beta.iter + rnorm(2*n, mean = 0, sd = sqrt(unlist(mod.bayes$tau2.draw[[iter]])))
      })
      Y.fit.iter <- do.call(cbind, Y.fit.iter)
      ci_by_Y <- apply(Y.fit.iter, 1, function(subj) c(quantile(subj, 0.025), quantile(subj, 0.975)))
      
      coverage_Y_train <- mean(sapply(1:n, function(i) {
        (Y[[1,1]][i,] >= ci_by_Y[1,i]) & (Y[[1,1]][i,] <= ci_by_Y[2,i])
      }))
      
      coverage_Y_test <- mean(sapply((n+1):(2*n), function(i) {
        (Y[[1,1]][i,] >= ci_by_Y[1,i]) & (Y[[1,1]][i,] <= ci_by_Y[2,i])
      }))
      
      # Save the posterior mean of the estimated error standard deviation for Y
      sigma.mat <- matrix(nrow = 3, ncol = 1)
      sigma.mat[3,] <- mean(sqrt(unlist(mod.bayes$tau2.draw)))
    }
    
    if (mod == "BPMF_Data_Mode" | mod == "BPMF_test" | mod == "BPMF_test_scale") {
      # Calculate the predicted E(Y) at each Gibbs sampling iteration
      EY.fit.iter <- lapply((burnin+1):nsample, function(iter) {
        mod.out$EY.draw[[iter]][[1,1]]
      })
      
      # Calculate the mean estimated E(Y)
      EY.fit <- Reduce("+", EY.fit.iter)/length(EY.fit.iter)
      
      # Calculate the coverage of training E(Y) and test E(Y)
      EY.fit.iter <- do.call(cbind, EY.fit.iter)
      ci_by_EY <- apply(EY.fit.iter, 1, function(subj) c(quantile(subj, 0.025), quantile(subj, 0.975)))
      
      # Calculate the coverage for true E(Y)
      coverage_EY_train <- mean(sapply(1:n, function(i) {
        (EY[[1,1]][i,] >= ci_by_EY[1,i]) & (EY[[1,1]][i,] <= ci_by_EY[2,i])
      }))
      
      coverage_EY_test <- mean(sapply((n+1):(2*n), function(i) {
        (EY[[1,1]][i,] >= ci_by_EY[1,i]) & (EY[[1,1]][i,] <= ci_by_EY[2,i])
      }))
      
      # Calculate the predicted Y at each Gibbs sampling iteration using training and testing scores
      Y.fit.iter <- lapply((burnin+1):nsample, function(iter) {
        mod.out$EY.draw[[iter]][[1,1]] + rnorm(2*n, mean = 0, sd = sqrt(unlist(mod.out$tau2.draw[[iter]])))
      })
      
      # Calculate the coverage of Y
      Y.fit.iter <- do.call(cbind, Y.fit.iter)
      ci_by_Y <- apply(Y.fit.iter, 1, function(subj) c(quantile(subj, 0.025), quantile(subj, 0.975)))
      
      coverage_Y_train <- mean(sapply(1:n, function(i) {
        (Y[[1,1]][i,] >= ci_by_Y[1,i]) & (Y[[1,1]][i,] <= ci_by_Y[2,i])
      }))
      
      coverage_Y_test <- mean(sapply((n+1):(2*n), function(i) {
        (Y[[1,1]][i,] >= ci_by_Y[1,i]) & (Y[[1,1]][i,] <= ci_by_Y[2,i])
      }))
    }
    
    # -------------------------------------------------------------------------
    # Testing the model fit on the true scores (to be removed later)
    # -------------------------------------------------------------------------
    
    if (mod == "test") {
      # Set the true ranks to be mod.ranks
      mod.ranks <- ranks
      joint.rank <- mod.ranks[1]
      indiv.rank <- mod.ranks[-1]
      
      # Set the true scores to be the all.scores.train
      true.scores <- cbind(Y[[1,1]], sim_data$V[[1,1]], do.call(cbind, sim_data$Vs))
      colnames(true.scores) <- c("y", rep("joint", joint.rank), rep("indiv", sum(indiv.rank)))
      
      # Subset the scores to just the training data
      true.scores.train <- true.scores[1:n,,drop=FALSE]
      
      # Fitting the Bayesian linear model
      mod.bayes <- bpmf(data = training_data, Y = Y_train, nninit = FALSE, model_params = model_params, 
                        ranks = mod.ranks, scores = true.scores.train[,-1], nsample = nsample)
      
      # Calculate the predicted E(Y) at each Gibbs sampling iteration using training and testing scores
      Y.fit.iter <- lapply((burnin+1):nsample, function(iter) {
        VStar.iter <- cbind(1, true.scores[,-1])
        beta.iter <- mod.bayes$beta.draw[[iter]][[1,1]] # Beta depends only on training data
        VStar.iter %*% beta.iter
      })
      
      Y.fit <- Reduce("+", Y.fit.iter)/length(Y.fit.iter)
      
      # Calculate the coverage of training E(Y) and test E(Y)
      Y.fit.iter <- do.call(cbind, Y.fit.iter)
      ci_by_Y <- apply(Y.fit.iter, 1, function(subj) c(quantile(subj, 0.025), quantile(subj, 0.975)))
      
      coverage_EY_train <- mean(sapply(1:n, function(i) {
        (EY[[1,1]][i,] >= ci_by_Y[1,i]) & (EY[[1,1]][i,] <= ci_by_Y[2,i])
      }))
      
      coverage_EY_test <- mean(sapply((n+1):(2*n), function(i) {
        (EY[[1,1]][i,] >= ci_by_Y[1,i]) & (EY[[1,1]][i,] <= ci_by_Y[2,i])
      }))
    }
    
    # -------------------------------------------------------------------------
    # Assess recovery of underlying joint and individual structure
    # -------------------------------------------------------------------------
    
    if (mod == "test") {
      joint.recovery.structure.train <- joint.recovery.structure.test <- indiv.recovery.structure.train <- indiv.recovery.structure.test <- NA
    }
    
    if (mod != "test") {
      # Joint structure
      joint.recovery.structure.train <- mean(sapply(1:q, function(s) {
        frob(mod.joint[[s]][,1:n] - joint.structure_train[[s,1]])/frob(joint.structure_train[[s,1]])
      }))
      
      # Individual structure
      indiv.recovery.structure.train <- mean(sapply(1:q, function(s) {
        frob(mod.individual[[s]][,1:n] - indiv.structure_train[[s,1]])/frob(indiv.structure_train[[s]])
      }))
    }
     
    if (mod != "test" & mod != "sJIVE" & mod != "BIP") {
      # Joint structure
      joint.recovery.structure.test <- mean(sapply(1:q, function(s) {
        frob(mod.joint[[s]][,(n+1):(2*n)] - joint.structure_test[[s,1]])/frob(joint.structure_test[[s,1]])
      }))
      
      # Individual structure
      indiv.recovery.structure.test <- mean(sapply(1:q, function(s) {
        frob(mod.individual[[s]][,(n+1):(2*n)] - indiv.structure_test[[s,1]])/frob(indiv.structure_test[[s,1]])
      }))
    }
    
    # -------------------------------------------------------------------------
    # Calculate prediction error for test data
    # -------------------------------------------------------------------------
    
    # Comparing the predicted Y to the training and test Y
    mse_EY_train <- frob(EY.fit[1:n,] - EY_train[[1,1]])/frob(EY_train[[1,1]])
    
    if (mod != "sJIVE" & mod != "BIP") {
      mse_EY_test <- frob(EY.fit[(n+1):(2*n),] - EY_test[[1,1]])/frob(EY_test[[1,1]])
    }
    
    # Save 
    if (estim_ranks) {
      file_name <- paste0("~/BayesianPMF/03Simulations/", mod, "/", mod, "_sim_", sim_iter, "_s2nX_", s2nX, "_s2nY_", s2nY, ".rda")
    }
    
    if (!estim_ranks) {
      file_name <- paste0("~/BayesianPMF/03Simulations/", mod, "_Fixed_Ranks", "/", mod,"_sim_", sim_iter, "_s2nX_", s2nX, "_s2nY_", s2nY, "_fixed_ranks.rda")
    }
    
    save(joint.recovery.structure.train, joint.recovery.structure.test,
         indiv.recovery.structure.train, indiv.recovery.structure.test,
         mse_EY_train, mse_EY_test, coverage_EY_train, coverage_EY_test,
         coverage_Y_train, coverage_Y_test, mod.ranks, sigma.mat,
         file = file_name)
    
    res <- c(joint.recovery.structure.train, joint.recovery.structure.test, indiv.recovery.structure.train, indiv.recovery.structure.test, mse_EY_train, mse_EY_test, coverage_EY_train, coverage_EY_test, coverage_Y_train, coverage_Y_test, mod.ranks)
    names(res) <- c("joint mse (train)", "joint mse (test)", "indiv mse (train)", "indiv mse (test)", "E(Y) mse (train)", "E(Y) mse (test)", "coverage E(y) (train)", "coverage E(y) (test)", "coverage y (train)", "coverage y (test)", "joint rank", paste("indiv rank", 1:q))
    
    res
  }
  stopCluster(cl)
  
  # Return
  sim_results
  
}

# Simulation study for assessing adjustment of label switching (permutation invariance)
identifiability_sim <- function(p.vec, n, ranks, response, true_params, model_params, sparsity = TRUE, identically_zero = TRUE, s2nX, s2nY, init_at_truth, num_in_spike = NULL, nsim, nsample, n_clust = 10) {
  
  # ---------------------------------------------------------------------------
  # Arguments:
  #
  # p.vec = number of features per source
  # n = sample size
  # ranks = vector of joint and individual ranks = c(joint rank, indiv rank 1, indiv rank 2, ...)
  # response = string in c(NULL, "continuous", "binary")
  # true_params = the list of true parameters under which to generate data
  # sparsity = should the data be generated with sparsity in the response? (Boolean)
  # identically_zero = should the coefficients in the spike under sparsity be identically 0?
  # s2nX, s2nY = signal-to-noise ratio in X and Y
  # nsim = number of simulations to run
  # nsample = number of Gibbs sampling iterations to draw for the linear model
  # n_clust = how many clusters to run simulation in parallel?
  # ---------------------------------------------------------------------------
  
  # Loading in the packages
  library(doParallel)
  library(foreach)
  
  cl <- makeCluster(n_clust)
  registerDoParallel(cl)
  funcs <- c("bpmf_data", "center_data", "bpmf", "get_results", "BIDIFAC",
             "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
             "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum", "factor_switching", 
             "factor_switching_permutation_plus_rotation", "erosheva_curtis_fs")
  packs <- c("Matrix", "MASS", "truncnorm")
  sim_results <- foreach (sim_iter = 1:nsim, .packages = packs, .export = funcs, .verbose = TRUE, .combine = rbind) %dopar% {
    
    # Set seed
    seed <- sim_iter + 200
    set.seed(seed)
    
    # -------------------------------------------------------------------------
    # Generate data 
    # -------------------------------------------------------------------------
    
    # Generate n samples
    sim_data <- bpmf_data(p.vec, n, ranks, true_params, s2nX = s2nX, s2nY = s2nY, response, missingness = NULL, entrywise = NULL, prop_missing = NULL, sparsity = sparsity, identically_zero = identically_zero, num_in_spike = num_in_spike)
    
    # Saving the data
    true_data <- sim_data$data
    true_data_list <- lapply(true_data, function(s) s)
    q <- length(p.vec)
    joint.structure <- sim_data$joint.structure
    indiv.structure <- sim_data$indiv.structure
    
    # The response
    Y <- sim_data$Y
    
    # Saving the ranks
    r <- ranks[1]
    r.vec <- ranks[-1]
    n_beta <- 1 + r + sum(r.vec)
    
    # The inclusion indicators
    true.gammas <- sim_data$gamma
    
    # Separating the inclusion indicators by joint and individual structure
    true.joint.gammas <- true.gammas[2:(r+1),,drop=FALSE]
    true.indiv.gammas.temp <- true.gammas[(2+r):n_beta,,drop=FALSE]
    
    # Initializing a list to store the individual gammas for each source
    true.indiv.gammas <- lapply(1:q, function(iter) list())
    
    for (s in 1:q) {
      if (s == 1) {
        true.indiv.gammas[[s]] <- true.indiv.gammas.temp[1:r.vec[s],,drop=FALSE]
      }
      
      if (s != 1) {
        true.indiv.gammas[[s]] <- true.indiv.gammas.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
      }
    }
    
    # Save the true components of the factorization
    true.V <- sim_data$V
    true.Vs <- sim_data$Vs
    true.U <- sim_data$U
    true.W <- sim_data$W

    # Save a burn-in
    burnin <- nsample/2
    thinned_iters <- seq(1, nsample, by = 10)
    thinned_iters_burnin <- seq(burnin, nsample, by = 10)
    
    # Save the true values to use as the initial values
    starting_values <- list(V = true.V, U = true.U, Vs = true.Vs, W = true.W)
    
    # -------------------------------------------------------------------------
    # Run the model and extract results
    # -------------------------------------------------------------------------
    
    # If initializing with the priors
    if (!init_at_truth) {
      res <- bpmf(data = true_data, Y = Y, nninit = FALSE, model_params = model_params, ranks = ranks, scores = NULL, sparsity = sparsity, nsample, progress = TRUE)
      
    }
    
    # If initializing at the true values
    if (init_at_truth) {
      res <- bpmf(data = true_data, Y = Y, nninit = FALSE, model_params = model_params, ranks = ranks, scores = NULL, sparsity = sparsity, nsample, progress = TRUE, starting_values = starting_values)
    }
  
    # Extract the posterior samples
    U.draw <- res$U.draw
    V.draw <- res$V.draw
    W.draw <- res$W.draw
    Vs.draw <- res$Vs.draw
    betas.draw <- res$beta.draw
    gammas.draw <- res$gamma.draw
    
    # Calculate the non-label swapped posterior inclusion indicators
    gammas.draw.non.ls <- do.call(cbind, lapply(gammas.draw, function(iter) iter[[1,1]]))
    
    # Calculate the posterior inclusion indicators
    post.non.ls <- rowMeans(gammas.draw.non.ls[,thinned_iters_burnin])
    post.non.ls[post.non.ls >= 0.5] <- 1
    post.non.ls[post.non.ls < 0.5] <- 0
    
    # -------------------------------------------------------------------------
    # Apply the label switching algorithm using correlation
    # -------------------------------------------------------------------------
    
    # Apply algorithm
    res.ls <- factor_switching(U.draw, V.draw, W.draw, Vs.draw, betas = betas.draw, 
                              gammas = gammas.draw, r = ranks[1], r.vec = ranks[-1],
                              nsample = nsample, thinned_iters_burnin = thinned_iters_burnin,
                              nninit = FALSE, pivot = list(true.V, true.Vs))
    
    # Save gammas
    gammas.draw.ls <- do.call(cbind, lapply(res.ls$swapped_gammas, function(iter) iter[[1,1]]))
    
    # Calculate the posterior inclusion indicators
    post.ls <- rowMeans(gammas.draw.ls[,thinned_iters_burnin])
    post.ls[post.ls >= 0.5] <- 1
    post.ls[post.ls < 0.5] <- 0
    
    # -------------------------------------------------------------------------
    # Apply the label switching algorithm using rotation correlation
    # -------------------------------------------------------------------------
    
    # Apply algorithm
    res.ls.rotate <- factor_switching_permutation_plus_rotation(U.draw, V.draw, W.draw, Vs.draw, betas.draw = betas.draw, 
                                                                gammas.draw = gammas.draw, r = ranks[1], r.vec = ranks[-1],
                                                                nsample = nsample, thinned_iters_burnin = thinned_iters_burnin,
                                                                nninit = FALSE, pivot = list(true.V, true.Vs))
    
    # Save gammas
    gammas.draw.ls.rotate <- do.call(cbind, lapply(res.ls.rotate$swapped_gammas, function(iter) iter[[1,1]]))
    
    # Calculate the posterior inclusion indicators
    post.ls.rotate <- rowMeans(gammas.draw.ls.rotate[,thinned_iters_burnin])
    post.ls.rotate[post.ls.rotate >= 0.5] <- 1
    post.ls.rotate[post.ls.rotate < 0.5] <- 0
    
    # -------------------------------------------------------------------------
    # Apply the label switching algorithm using the Erosheva and Curtis approach
    # -------------------------------------------------------------------------
    
    # Apply algorithm
    res.es <- erosheva_curtis_fs(U.draw, V.draw, W.draw, Vs.draw, betas = betas.draw, 
                                   gammas = gammas.draw, r = ranks[1], r.vec = ranks[-1],
                                   n, p.vec, nsample = nsample, thinned_iters_burnin = thinned_iters_burnin,
                                   BIDIFAC_solution = list(V = true.V, U = true.U, Vs = true.Vs, W = true.W))
    
    # Save gammas
    gammas.draw.es <- do.call(cbind, lapply(res.es$swapped_gammas, function(iter) iter[[1,1]]))
    
    # Calculate the posterior inclusion indicators
    post.es <- rowMeans(gammas.draw.es[,thinned_iters_burnin])
    post.es[post.es >= 0.5] <- 1
    post.es[post.es < 0.5] <- 0
    
    # -------------------------------------------------------------------------
    # Compare the true gammas indicator vector to the label swapped and non-label
    # swapped versions
    # -------------------------------------------------------------------------
    
    # Label swapped with correlation
    ssd.ls <- frob(true.gammas - post.ls)/length(true.gammas)
    
    # Label swapped with correlation and rotate
    ssd.ls.rotate <- frob(true.gammas - post.ls.rotate)/length(true.gammas)
    
    # Label swapped with log density
    ssd.es <- frob(true.gammas - post.es)/length(true.gammas)
    
    # Non-label swapped
    ssd.non.ls <- frob(true.gammas - post.non.ls)/length(true.gammas)
    
    # -------------------------------------------------------------------------
    # Save the joint and individual structures that explain variability in Y in 
    # the true data
    # -------------------------------------------------------------------------
    
    # Initialize variables to store the structure that explains variability in Y and structure that does not 
    true.joint.structure.slab <- true.joint.structure.spike <- true.indiv.structure.slab <- true.indiv.structure.spike <- matrix(list(), nrow = q, ncol = 1)
    
    # Save the columns that correspond to the slab
    true_joint_slab_inds <- true.joint.gammas == 1
    
    for (s in 1:q) {
      # Save the joint structure that explains variability
      true.joint.structure.slab[[s,1]] <- (sim_data$U[[s,1]][,true_joint_slab_inds, drop = FALSE]) %*% t(sim_data$V[[1,1]][,true_joint_slab_inds, drop = FALSE])
      
      # Save the joint structure that does not explain variability
      true.joint.structure.spike[[s,1]] <- (sim_data$U[[s,1]][,!true_joint_slab_inds, drop = FALSE]) %*% t(sim_data$V[[1,1]][,!true_joint_slab_inds, drop = FALSE])
      
      # Save the columns for each individual structure that corresponds to the slab
      true_indiv_slab_inds <- true.indiv.gammas[[s]] == 1
      
      # Save the individual structure that explains variability
      true.indiv.structure.slab[[s,1]] <- (sim_data$W[[s,s]][,true_indiv_slab_inds, drop = FALSE]) %*% t(sim_data$Vs[[1,s]][,true_indiv_slab_inds, drop = FALSE])
      
      # Save the individual structure that does not explain variability
      true.indiv.structure.spike[[s,1]] <- (sim_data$W[[s,s]][,!true_indiv_slab_inds, drop = FALSE]) %*% t(sim_data$Vs[[1,s]][,!true_indiv_slab_inds, drop = FALSE])
    }
    
    # -------------------------------------------------------------------------
    # Save the joint and individual structures that explain variability in Y 
    # at each Gibbs sampling iteration
    # ------------------------------------------------------------------------
    
    # Initialize the variables to store the structure that explains variability in Y and structure that does not at each iteration
    joint.structure.slab <- joint.structure.spike <- indiv.structure.slab <- indiv.structure.spike <- lapply(1:length(thinned_iters_burnin), function(iter) matrix(list(), nrow = q, ncol = 1))
    
    # Initialize index to store results
    ind <- 1
    
    # Go through the iterations after burn-in and calculate these structures
    for (iter in thinned_iters_burnin) {
      # Save the current Gibbs sampling iteration
      U.iter <- U.draw[[iter]]
      V.iter <- V.draw[[iter]]
      Vs.iter <- Vs.draw[[iter]]
      W.iter <- W.draw[[iter]]
      gamma.iter <- gammas.draw[[iter]][[1,1]]
      
      # Decompose the gammas into joint and individual
      indiv.gammas.temp.iter <- gamma.iter[(2+r):n_beta,,drop=FALSE]
      
      # Separating the inclusion indicators by joint and individual structure
      joint.gammas.iter <- gamma.iter[2:(r+1),,drop=FALSE]
      indiv.gammas.temp <- lapply(1:q, function(iter) list())
      
      for (s in 1:q) {
        if (s == 1) {
          indiv.gammas.temp[[s]] <- indiv.gammas.temp.iter[1:r.vec[s],,drop=FALSE]
        }
        
        if (s != 1) {
          indiv.gammas.temp[[s]] <- indiv.gammas.temp.iter[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
      }
      
      # Save the columns in the joint structure that correspond to the slab
      joint_slab_inds <- joint.gammas.iter == 1
      
      # Save the structures that do and don't explain variability in joint and individual structures
      for (s in 1:q) {
        # Save the joint structure that explains variability
        joint.structure.slab[[ind]][[s,1]] <- (U.iter[[s,1]][,joint_slab_inds, drop = FALSE]) %*% t(V.iter[[1,1]][,joint_slab_inds, drop = FALSE])
        
        # Save the joint structure that does not explain variability
        joint.structure.spike[[ind]][[s,1]] <- (U.iter[[s,1]][,!joint_slab_inds, drop = FALSE]) %*% t(V.iter[[1,1]][,!joint_slab_inds, drop = FALSE])
        
        # Save the columns for each indivdidual structure that corresponds to the slab
        indiv_slab_inds <- indiv.gammas.temp[[s]] == 1
        
        # Save the individual structure that explains variability
        indiv.structure.slab[[ind]][[s,1]] <- (W.iter[[s,s]][,indiv_slab_inds, drop = FALSE]) %*% t(Vs.iter[[1,s]][,indiv_slab_inds, drop = FALSE])
        
        # Save the individual structure that does not explain variability
        indiv.structure.spike[[ind]][[s,1]] <- (W.iter[[s,s]][,!indiv_slab_inds, drop = FALSE]) %*% t(Vs.iter[[1,s]][,!indiv_slab_inds, drop = FALSE])
      }
      
      # Update index
      ind <- ind + 1
    } 
    
    # Calculate the posterior mean
    joint.structure.slab.mean <- joint.structure.spike.mean <- indiv.structure.slab.mean <- indiv.structure.spike.mean <-  matrix(list(), nrow = q, ncol = 1)
    
    for (s in 1:q) {
      joint.structure.slab.mean[[s,1]] <- Reduce("+", lapply(joint.structure.slab, function(iter) iter[[s,1]]))/length(thinned_iters_burnin)
      joint.structure.spike.mean[[s,1]] <- Reduce("+", lapply(joint.structure.spike, function(iter) iter[[s,1]]))/length(thinned_iters_burnin)
      indiv.structure.slab.mean[[s,1]] <- Reduce("+", lapply(indiv.structure.slab, function(iter) iter[[s,1]]))/length(thinned_iters_burnin)
      indiv.structure.spike.mean[[s,1]] <- Reduce("+", lapply(indiv.structure.spike, function(iter) iter[[s,1]]))/length(thinned_iters_burnin)
    }
    
    # -------------------------------------------------------------------------
    # Compare the structures that do and don't explain variability 
    # -------------------------------------------------------------------------
    
    joint.slab.mse <- mean(sapply(1:q, function(s) {
      (frob(true.joint.structure.slab[[s,1]] - joint.structure.slab.mean[[s,1]]))/frob(true.joint.structure.slab[[s,1]])
    }))
    
    joint.spike.mse <- mean(sapply(1:q, function(s) {
      (frob(true.joint.structure.spike[[s,1]] - joint.structure.spike.mean[[s,1]]))/frob(true.joint.structure.spike[[s,1]])
    }))
    
    indiv.slab.mse <- mean(sapply(1:q, function(s) {
      (frob(true.indiv.structure.slab[[s,1]] - indiv.structure.slab.mean[[s,1]]))/frob(true.indiv.structure.slab[[s,1]])
    }))
    
    indiv.spike.mse <- mean(sapply(1:q, function(s) {
      (frob(true.indiv.structure.spike[[s,1]] - indiv.structure.spike.mean[[s,1]]))/frob(true.indiv.structure.spike[[s,1]])
    }))

    # -------------------------------------------------------------------------
    # Return
    # -------------------------------------------------------------------------
    
    ssd_and_mse <- c(ssd.ls, ssd.ls.rotate, ssd.es, ssd.non.ls, joint.slab.mse, joint.spike.mse, indiv.slab.mse, indiv.spike.mse)
    names(ssd_and_mse) <- c("Corrected (Corr)", "Corrected (Corr + Rotate)", "Corrected (E&S)", "Non-Corrected", "Joint Slab MSE", "Joint Spike MSE", "Indiv Slab MSE", "Indiv Spike MSE")
    ssd_and_mse
  }
  stopCluster(cl)
  
  # Return the results
  sim_results
}

# -----------------------------------------------------------------------------
# Create tables for results
# -----------------------------------------------------------------------------

create_validation_table <- function(results_list, condition) {
  # Create a set of rows in the validation table for the validation simulation
  
  # -----------------------------------------------------------------------------
  # Arguments:
  # 
  # results_list = contains list of results, some of which may be NULL
  # condition (string) = name of the validation condition
  # -----------------------------------------------------------------------------
  
  # Create a 3-row dataframe to return
  dt <- data.frame(Condition = rep(condition, 3),
                   Metric = c("Coverage", "MSE", "CI Width"),
                   Joint_Obs = numeric(3),
                   Indiv_Obs = numeric(3),
                   Joint_Mis = numeric(3),
                   Indiv_Mis = numeric(3),
                   EY_Obs = numeric(3),
                   EY_Mis = numeric(3),
                   tau2 = numeric(3))
  
  # Fill in the table --
  
  # Joint structure
  if (length(results_list$`joint structure`[[1,1]]) != 2) { # If there is NOT missingness
    dt$Joint_Obs <- c(mean(sapply(results_list$`joint structure`, function(source) source$avg_coverage)),
                      mean(sapply(results_list$`joint structure`, function(source) source$avg_mse)),
                      mean(sapply(results_list$`joint structure`, function(source) source$avg_ci_width)))
    
    # Individual structure
    dt$Indiv_Obs <- c(mean(sapply(results_list$`indiv structure`, function(source) source$avg_coverage)),
                      mean(sapply(results_list$`indiv structure`, function(source) source$avg_mse)),
                      mean(sapply(results_list$`indiv structure`, function(source) source$avg_ci_width)))
    
    dt$Joint_Mis <- NA
    
    dt$Indiv_Mis <- NA
  }
  
  if (length(results_list$`joint structure`[[1,1]]) == 2) { # If there IS missingness
    # Observed joint structure
    dt$Joint_Obs <- c(mean(sapply(results_list$`joint structure`, function(source) source$observed$avg_coverage)),
                      mean(sapply(results_list$`joint structure`, function(source) source$observed$avg_mse)),
                      mean(sapply(results_list$`joint structure`, function(source) source$observed$avg_ci_width)))
    
    # Observed individual structure
    dt$Indiv_Obs <- c(mean(sapply(results_list$`indiv structure`, function(source) source$observed$avg_coverage)),
                      mean(sapply(results_list$`indiv structure`, function(source) source$observed$avg_mse)),
                      mean(sapply(results_list$`indiv structure`, function(source) source$observed$avg_ci_width)))
    
    # Missing joint structure
    dt$Joint_Mis <- c(mean(sapply(results_list$`joint structure`, function(source) source$missing$avg_coverage)),
                      mean(sapply(results_list$`joint structure`, function(source) source$missing$avg_mse)),
                      mean(sapply(results_list$`joint structure`, function(source) source$missing$avg_ci_width)))
    
    # Missing individual structure
    dt$Indiv_Mis <- c(mean(sapply(results_list$`indiv structure`, function(source) source$missing$avg_coverage)),
                      mean(sapply(results_list$`indiv structure`, function(source) source$missing$avg_mse)),
                      mean(sapply(results_list$`indiv structure`, function(source) source$missing$avg_ci_width)))
  }
  
  
  # E(Y)
  if (is.null(results_list$EY[[1,1]])) {
    dt$EY_Obs <- NA
    dt$EY_Mis <- NA
  }
  
  if (!is.null(results_list$EY[[1,1]])) {
    if (length(results_list$EY[[1,1]]) != 2) { # If there is NO missingness in response
      dt$EY_Obs <- c(mean(sapply(results_list$EY, function(source) source$avg_coverage)),
                     mean(sapply(results_list$EY, function(source) source$avg_mse)),
                     mean(sapply(results_list$EY, function(source) source$avg_ci_width)))
      
      dt$EY_Mis <- NA
    }
    
    if (length(results_list$EY[[1,1]]) == 2) { # If there IS missingness in response
      # Observed joint structure
      dt$EY_Obs <- c(mean(sapply(results_list$EY, function(source) source$observed$avg_coverage)),
                     mean(sapply(results_list$EY, function(source) source$observed$avg_mse)),
                     mean(sapply(results_list$EY, function(source) source$observed$avg_ci_width)))
      
      # Missing joint structure
      dt$EY_Mis <- c(mean(sapply(results_list$EY, function(source) source$missing$avg_coverage)),
                        mean(sapply(results_list$EY, function(source) source$missing$avg_mse)),
                        mean(sapply(results_list$EY, function(source) source$missing$avg_ci_width)))
    }
  }
  
  # tau2
  if (is.null(results_list$tau2[[1,1]])) {
    dt$tau2 <- NA
  }
  
  if (!is.null(results_list$tau2[[1,1]])) {
    dt$tau2 <- c(results_list$tau2[[1,1]]$avg_coverage, results_list$tau2[[1,1]]$avg_mse, results_list$tau2[[1,1]]$avg_ci_width)
  }
  
  # Return
  dt
}

create_simulation_table <- function(simulation_results, mod.list, path.list, s2nX, s2nY) {
  
  # ---------------------------------------------------------------------------
  # Arguments:
  #
  # mod.list = vector with model names desired
  # path.list = list of strings for file directory containing results from 
  #   each model. List entries should be named by model
  # s2nX = current s2n for data to load in
  # s2nY = current s2n for response to load in
  # ---------------------------------------------------------------------------
  
  # Set up the results table
  simulation_results$Metric <- c("MSE Joint (train)", "MSE Joint (test)", "MSE Indiv (train)",
                                 "MSE Indiv (test)", "MSE EY (train)", "MSE EY (test)", 
                                 "Coverage EY (train)", "Coverage EY (test)", "Joint Rank", 
                                 "Indiv Rank 1", "Indiv Rank 2", "Estim Y Err Var")
  simulation_results$s2nX <- s2nX
  simulation_results$s2nY <- s2nY
  
  # Iterate through the models
  for (mod in mod.list) {
    # Load in the files for this model
    path <- path.list[[mod]]
    all_files <- list.files(path)
    all_files_split <- strsplit(all_files, split = "_")
    
    # Save the indices of the current s2n
    s2nX_ind <- which(all_files_split[[1]] == "s2nX") + 1
    s2nY_ind <- which(all_files_split[[1]] == "s2nY") + 1
    
    # Save the names of the current results
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[s2nX_ind] == s2nX) & (file[s2nY_ind] == paste0(s2nY, ".rda")))]
    
    # Create a temporary table to load the results into
    simulation_results_temp <- data.frame(joint.structure.train = numeric(length(files_for_s2nX_s2nY)),
                                          joint.structure.test = numeric(length(files_for_s2nX_s2nY)),
                                          indiv.structure.train = numeric(length(files_for_s2nX_s2nY)),
                                          indiv.structure.test = numeric(length(files_for_s2nX_s2nY)),
                                          mse_ey_train = numeric(length(files_for_s2nX_s2nY)),
                                          mse_ey_test = numeric(length(files_for_s2nX_s2nY)),
                                          coverage_EY_train = numeric(length(files_for_s2nX_s2nY)),
                                          coverage_EY_test = numeric(length(files_for_s2nX_s2nY)),
                                          joint.rank = numeric(length(files_for_s2nX_s2nY)),
                                          indiv.rank1 = numeric(length(files_for_s2nX_s2nY)),
                                          indiv.rank2 = numeric(length(files_for_s2nX_s2nY)),
                                          est.y.err.var = numeric(length(files_for_s2nX_s2nY)))
    
    # Iteratively load in the results
    for (ind in 1:length(files_for_s2nX_s2nY)) {
      # Load in the results
      file <- files_for_s2nX_s2nY[ind]
      res <- load(paste0(path, file))
      
      # If no error variance for Y was estimated
      if (is.null(sigma.mat)) {
        err.y <- NA 
      }
      if (!is.null(sigma.mat)) {
        err.y <- sigma.mat[3]
      }
      
      simulation_results_temp[ind,] <- c(unlist(sapply(res, function(met) get(met)))[1:11], err.y)
    }
    
    # Take the mean for the results
    simulation_results_temp_mean <- colMeans(simulation_results_temp)
    
    # Save the results
    simulation_results[,mod] <- simulation_results_temp_mean
  }
  
  # Return 
  simulation_results
}

# -----------------------------------------------------------------------------
# Functions for data application
# -----------------------------------------------------------------------------

run_model_with_cv <- function(mod, hiv_copd_data, outcome, outcome_name, ind_of_pairs, model_params, nsample, results_wd = "~/BayesianPMF/04DataApplication/") {
  
  # ---------------------------------------------------------------------------
  # For a given model, run cross validation in parallel
  #
  # Arguments:
  # mod (character): in c(JIVE, BIDIFAC, MOFA)
  # hiv_copd_data (list): data for the model
  # outcome (double): response vector for Bayesian modeling
  # outcome_name (character): name of the outcome for saving files
  # ind_of_pairs (list): indices to cross validate over
  # model_params (list): model parameters for the Bayesian model
  # nsample (int): how many Gibbs samples to generate
  # ---------------------------------------------------------------------------
  
  # Set the functions and packages for the parallel computation
  funcs <- c("bpmf_data", "center_data", "bpmf_full_mode", "bpmf_data_mode", "bpmf_test", "bpmf_test_scale", "get_results", "BIDIFAC",
             "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
             "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum",
             "bidifac.plus.impute", "bidifac.plus.given")
  packs <- c("Matrix", "MASS", "truncnorm", "r.jive", "sup.r.jive", "natural", "RSpectra", "MOFA2")
  
  # Load in the full training data fit
  results_path <- paste0("~/BayesianPMF/04DataApplication/", mod, "/", mod, "_training_data_fit.rda") 
  out <- load(results_path, verbose = TRUE)
  
  # For JIVE
  if (mod == "JIVE") {
    # Saving the joint structure
    mod.joint <- get(out)$joint
    
    # Saving the individual structure
    mod.individual <- get(out)$individual
    
    # Saving the joint rank
    joint.rank <- get(out)$rankJ
    
    # Saving the individual ranks
    indiv.rank <- get(out)$rankA
    
    # Obtaining the joint scores
    if (joint.rank != 0)  {
      svd.joint <- svd(mod.joint[[1]])
      joint.scores <- (svd.joint$v[,1:joint.rank,drop=FALSE]) %*% diag(svd.joint$d[1:joint.rank], nrow = joint.rank)
    }
    if (joint.rank == 0) {
      joint.scores <- NULL
    }
    
    # Obtaining the individual scores
    indiv.scores <- lapply(1:q, function(s) {
      if (indiv.rank[s] != 0) {
        svd.source <- svd(mod.individual[[s]])
        (svd.source$v[,1:indiv.rank[s],drop=FALSE]) %*% diag(svd.source$d[1:indiv.rank[s]], nrow = indiv.rank[s])
      }
    })
  }
  
  # For BIDIFAC
  if (mod == "BIDIFAC") {
    # Saving the column structure (the joint structure)
    mod.joint <- lapply(1:q, function(s) {
      get(out)$C[[s,1]]
    })
    
    # Saving the individual structure 
    mod.individual <- lapply(1:q, function(s) {
      get(out)$I[[s,1]]
    })
    
    # Saving the joint rank
    joint.rank <- rankMatrix(get(out)$C[[1,1]])[1]
    
    # Saving the individual ranks
    indiv.rank <- sapply(get(out)$I, function(s) {
      rankMatrix(s)[1]
    }) 
    
    # Obtaining the joint scores
    if (joint.rank != 0)  {
      svd.joint <- svd(mod.joint[[1]])
      joint.scores <- (svd.joint$v[,1:joint.rank,drop=FALSE]) %*% diag(svd.joint$d[1:joint.rank], nrow = joint.rank)
    }
    if (joint.rank == 0) {
      joint.scores <- NULL
    }
    
    # Obtaining the individual scores
    indiv.scores <- lapply(1:q, function(s) {
      if (indiv.rank[s] != 0) {
        svd.source <- svd(mod.individual[[s]])
        (svd.source$v[,1:indiv.rank[s],drop=FALSE]) %*% diag(svd.source$d[1:indiv.rank[s]], nrow = indiv.rank[s])
      }
    })
  }
  
  # For MOFA
  if (mod == "MOFA") {
    # Getting the variance explained in each source by each factor
    mod.var.exp <- get_variance_explained(get(out))$r2_per_factor$group1 # variance explained by factor per view
    
    # Save which views each factor applies to
    joint_or_individual <- lapply(1:nrow(mod.var.exp), function(factor) {
      # Save variance explained for current factor
      row <- mod.var.exp[factor,]
      
      ind.max <- which.max(row) # Highest var explained
      ind.min <- which.min(row) # Minimum var explained
      
      # First, check that factor explains at least 1% of variation in each source
      greater_than_1p <- any(row > 1)
      
      # If explains more than 1% variation in at least one source
      if (greater_than_1p) {
        
        # If factor explains substantial (max var is less than 2 * min var) variation in both sources, it is joint
        if (row[ind.max] < 2*row[-ind.max]) { 
          c(1:q)
        } else { # If factor explains substantial variation in just one source, it is individual
          ind.max
        }
        
      }
      
    })
    
    # Save the joint factors
    joint.factors <- c(1:length(joint_or_individual))[sapply(joint_or_individual, function(factor) length(factor) == q)]
    
    # Save the individual factors
    indiv.factors <- lapply(1:q, function(s) {
      c(1:length(joint_or_individual))[sapply(joint_or_individual, function(factor) all(factor %in% s) & length(factor) > 0)]
    })
    
    # Saving the joint rank
    joint.rank <- length(joint.factors)
    
    # Saving the individual rank
    indiv.rank <- sapply(indiv.factors, length)
    
    # Save the underlying structure
    mod.joint <- lapply(1:q, function(s) {
      t(get(out)@expectations$Z$group1[, joint.factors, drop = FALSE] %*% t(get(out)@expectations$W[[s]][, joint.factors, drop = FALSE]))
    })
    mod.individual <- lapply(1:q, function(s) {
      t(get(out)@expectations$Z$group1[, indiv.factors[[s]], drop = FALSE] %*% t(get(out)@expectations$W[[s]][, indiv.factors[[s]], drop = FALSE]))
    })
    
    # Save the MOFA scores (all.equal(get_factors(mod.out)$group1, mod.out@expectations$Z$group1[, joint.factors, drop = FALSE]) # TRUE!)
    mofa.scores <- get_factors(get(out))$group1
    
    # Save the joint scores
    if (joint.rank != 0) {
      joint.scores <- mofa.scores[,joint.factors, drop = FALSE]
    }
    if (joint.rank == 0) {
      joint.scores <- NULL
    }
    
    # Save the individual scores
    indiv.scores <- lapply(1:q, function(s) {
      if (indiv.rank[s] != 0) {
        indiv.scores <- mofa.scores[,unlist(indiv.factors[[s]]), drop = FALSE]
      }
    })
  }
  
  # Combine the scores together
  all.scores <- cbind(joint.scores, do.call(cbind, indiv.scores))
  colnames(all.scores) <- c(rep("joint", joint.rank), rep("indiv", sum(indiv.rank)))
  
  # Running in parallel
  cl <- makeCluster(10)
  registerDoParallel(cl)
  fev1pp_cv <- foreach(pair = ind_of_pairs, .packages = packs, .export = funcs, .verbose = TRUE) %dopar% {
    # Create a new vector of the outcome with the current pair set to NA
    outcome_cv <- outcome
    outcome_cv[[1,1]][pair:(pair+1),] <- NA
    
    # Fitting the Bayesian linear model
    mod.bayes <- bpmf_data_mode(data = hiv_copd_data, Y = outcome_cv, nninit = FALSE, model_params = model_params, 
                                ranks = c(joint.rank, indiv.rank), scores = all.scores, nsample = nsample)
    
    # Save the imputed outcomes
    Ym.draw_pair <- mod.bayes$Ym.draw
    
    # Save just the relevant output
    ranks <- c(joint.rank, indiv.rank)
    save(Ym.draw_pair, ranks, file = paste0(results_wd, mod,"/", outcome_name, "_CV_", mod, "_Pair_", pair, ".rda"))
    
    # Remove large objects
    rm(mod.bayes)
    
    # Garbage collection
    gc()
  }
  stopCluster(cl)
  
}

# -----------------------------------------------------------------------------
# Addressing the factor switching problem
# -----------------------------------------------------------------------------

# The factor switching algorithm: matches results to the posterior mode
factor_switching <- function(U.draw, V.draw, W.draw, Vs.draw, betas = NULL, gammas = NULL, r, r.vec, nsample, thinned_iters_burnin, nninit = TRUE, pivot = NULL) {
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # draws (list): the draws to swap according to the pivot
  # rank (int): number of columns to iterate through
  # loadings (list): the corresponding loadings for each draw of the matrix to swap
  # betas (double vec): list of vectors of the coefficients drawn at each iteration
  # gammas (int vec): list of vectors of the inclusion indicators drawn at each iteration
  # r (int) = joint rank
  # r.vec (c(ints)) = vector of individual ranks
  # nsample (int) = number of posterior samples
  # thinned_iters_burnin (c(ints)) = vector of indices for samples after burn-in + thinning if desired
  # nninit (boolean) = TRUE if initialized at theoretical posterior mode or FALSE if use posterior mean as pivot
  # pivot (list) = list of V and Vs to use as pivots. If NULL, must provide nninit. 
  # ---------------------------------------------------------------------------
  
  # Saving the number of sources
  q <- length(r.vec)
  
  # Returning the swapped matrices
  swapped_U.draw <- lapply(1:length(U.draw), function(iter) list())
  swapped_V.draw <- lapply(1:length(V.draw), function(iter) list())
  swapped_W.draw <- lapply(1:length(W.draw), function(iter) list())
  swapped_Vs.draw <- lapply(1:length(Vs.draw), function(iter) list())
  swapped_betas <- lapply(1:length(betas), function(iter) list())
  swapped_gammas <- lapply(1:length(gammas), function(iter) list())
  
  # Returning the swaps made so they can be undone
  swaps_made <- lapply(1:length(swapped_U.draw), function(iter) list())
  
  # Storing the sign changes so they can be undone
  signs_changed <- lapply(1:length(swapped_U.draw), function(iter) list())
  
  # Setting the pivots to the posterior mode (result from BIDIFAC) or posterior mean
  if (is.null(pivot)) {
    if (nninit) {
      pivot_V <- V.draw[[1]]
      pivot_Vs <- Vs.draw[[1]]
    }
    
    if (!nninit) {
      V.draw.thinned.burnin <- lapply(V.draw[thinned_iters_burnin], function(iter) iter[[1,1]])
      pivot_V <- matrix(list(), nrow = 1, ncol = 1)
      pivot_V[[1,1]] <- Reduce("+", V.draw.thinned.burnin)/length(V.draw.thinned.burnin)
      
      Vs.draw.thinned.burnin <- Vs.draw[thinned_iters_burnin]
      pivot_Vs <- matrix(list(), nrow = 1, ncol = q)
      
      for (s in 1:q) {
        pivot_Vs[[1,s]] <- Reduce("+", lapply(Vs.draw.thinned.burnin, function(iter) iter[[1,s]]))/length(Vs.draw.thinned.burnin)
      }
    }
  }
  
  # Setting the pivot to be user-chosen
  if (!is.null(pivot)) {
    pivot_V <- pivot[[1]]
    pivot_Vs <- pivot[[2]]
  }
  
  for (iter in 1:length(U.draw)) {
    # -------------------------------------------------------------------------
    # Initializing for iter
    # 
    # For each iteration, arrange the columns of each draw to match that 
    # of the pivot. 
    # -------------------------------------------------------------------------
    
    # Initializing the rearranged version
    tilde_V <- matrix(list(), nrow = 1, ncol = 1)
    tilde_U <- matrix(list(), nrow = q, ncol = 1)
    tilde_W <- matrix(list(), nrow = q, ncol = q)
    tilde_Vs <- matrix(list(), nrow = 1, ncol = q)
    
    tilde_beta <- matrix(list(), nrow = 1, ncol = 1)
    tilde_beta_joint <- matrix(list(), nrow = 1, ncol = 1)
    tilde_beta_indiv <- matrix(list(), nrow = q, ncol = 1)
    
    tilde_gamma <- matrix(list(), nrow = 1, ncol = 1)
    tilde_gamma_joint <- matrix(list(), nrow = 1, ncol = 1)
    tilde_gamma_indiv <- matrix(list(), nrow = q, ncol = 1)
    
    tilde_V[[1,1]] <- matrix(nrow = nrow(V.draw[[1]][[1,1]]), ncol = ncol(V.draw[[1]][[1,1]]))
    
    for (s in 1:q) {
      tilde_U[[s,1]] <- matrix(nrow = nrow(U.draw[[1]][[s,1]]), ncol = ncol(U.draw[[1]][[s,1]]))
      tilde_W[[s,s]] <- matrix(nrow = nrow(W.draw[[1]][[s,s]]), ncol = ncol(W.draw[[1]][[s,s]]))
      tilde_Vs[[1,s]] <- matrix(nrow = nrow(Vs.draw[[1]][[1,s]]), ncol = ncol(Vs.draw[[1]][[1,s]]))
      
      for (ss in 1:q) {
        if (ss != s) {
          tilde_W[[s,ss]] <- matrix(0, nrow = nrow(W.draw[[1]][[s,s]]), ncol = ncol(W.draw[[1]][[ss,ss]]))
        }
      }
    }
    
    n_beta <- 1 + r + sum(r.vec)
    tilde_beta[[1,1]] <- matrix(nrow = n_beta, ncol = 1)
    tilde_beta_joint[[1,1]] <- matrix(nrow = r, ncol = 1)
    
    tilde_gamma[[1,1]] <- matrix(nrow = n_beta, ncol = 1)
    tilde_gamma_joint[[1,1]] <- matrix(nrow = r, ncol = 1)
    
    for (s in 1:q) {
      tilde_beta_indiv[[s,1]] <- matrix(nrow = r.vec[s], ncol = 1)
      tilde_gamma_indiv[[s,1]] <- matrix(nrow = r.vec[s], ncol = 1)
    }
    
    # -------------------------------------------------------------------------
    # Storing the current V, U, beta, and gamma
    # -------------------------------------------------------------------------
    
    current_V <- V.draw[[iter]]
    current_U <- U.draw[[iter]]
    current_Vs <- Vs.draw[[iter]]
    current_W <- W.draw[[iter]]
    
    if (!is.null(betas)) {
      current_beta <- betas[[iter]]
      current_beta_indiv <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        current_beta_joint <- current_beta[[1,1]][2:(r+1),, drop = FALSE] 
      } else {
        current_beta_joint <- NULL
      }
      
      # Individual effects
      current_beta_indiv.temp <- current_beta[[1,1]][(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) {
            current_beta_indiv[[s, 1]] <- current_beta_indiv.temp[1:r.vec[s],, drop = FALSE] 
          }
          if (s != 1) {
            current_beta_indiv[[s, 1]] <- current_beta_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
      }
    }
    
    if (!is.null(gammas)) {
      current_gamma <- gammas[[iter]] 
      current_gamma_indiv <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        current_gamma_joint <- current_gamma[[1,1]][2:(r+1),, drop = FALSE]
      } else {
        current_gamma_joint <- NULL
      }
      
      # Individual effects
      current_gamma_indiv.temp <- current_gamma[[1,1]][(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) {
            current_gamma_indiv[[s,1]] <- current_gamma_indiv.temp[1:r.vec[s],, drop = FALSE] 
          }
          if (s != 1) {
            current_gamma_indiv[[s,1]] <- current_gamma_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Finding the optimal order
    # -------------------------------------------------------------------------
    
    # Storing the swaps and signs
    current_swaps_joint <- current_signs_joint <- c()
    current_swaps_indiv <- current_signs_indiv <- lapply(1:q, function(s) c())
    
    # Start with the joint structure
    
    for (k in 1:r) {
      # for each column in the pivot, rearrange the columns in draws to match the order
      corrs <- apply(current_V[[1,1]], 2, function(current_col) cor(current_col, pivot_V[[1,1]][,k])) # compute correlation between each 
      jk <- which.max(abs(corrs)) # storing the index of the highest correlation
      ak <- sign(corrs[jk]) # storing the sign of that correlation
      
      # Setting the kth column of tilde_V to be the ak*current_V[,jk] column. 
      tilde_V[[1,1]][,k] <- ak*current_V[[1,1]][,jk]
      
      # Rearranging the corresponding column in U:
      for (s in 1:q) {
        tilde_U[[s,1]][,k] <- ak*current_U[[s,1]][,jk]
      }
      
      # Rearranging the corresponding rows in betas:
      if (!is.null(betas)) {
        tilde_beta_joint[[1,1]][k,] <- ak*current_beta_joint[jk,]
      }
      
      # Rearranging the corresponding rows in the gammas:
      if (!is.null(gammas)) {
        tilde_gamma_joint[[1,1]][k,] <- current_gamma_joint[jk,]
      }
      
      # Recording what swap was made. The first entry would be the column index
      # of current_V that was moved to the first spot. The second index would be 
      # the column index of current_V that was moved to the second column, etc. 
      current_swaps_joint[k] <- jk
      
      # Recording the sign change (if any) that was made
      current_signs_joint[k] <- ak
      
      # Removing the already-swapped columns in current_V from the running:
      current_V[[1,1]][,jk] <- NA
    }
    
    # Then the individual structures
    
    for (s in 1:q) {
      for (k in 1:r.vec[s]) {
        # for each column in the pivot, rearrange the columns in draws to match the order
        corrs <- apply(current_Vs[[1,s]], 2, function(current_col) cor(current_col, pivot_Vs[[1,s]][,k])) # compute correlation between each 
        jk <- which.max(abs(corrs)) # storing the index of the highest correlation
        ak <- sign(corrs[jk]) # storing the sign of that correlation
        
        # Setting the kth column of \tilde_Vs to be the ak*current_Vs[,jk] column. 
        tilde_Vs[[1,s]][,k] <- ak*current_Vs[[1,s]][,jk]
        
        # Rearranging the corresponding column in W:
        tilde_W[[s,s]][,k] <- ak*current_W[[s,s]][,jk]
        
        # Rearranging the corresponding rows in betas:
        if (!is.null(betas)) {
          tilde_beta_indiv[[s,1]][k,] <- ak*current_beta_indiv[[s,1]][jk,]
        }
        
        # Rearranging the corresponding rows in the gammas:
        if (!is.null(gammas)) {
          tilde_gamma_indiv[[s,1]][k,] <- current_gamma_indiv[[s,1]][jk,]
        }
        
        # Recording what swap was made. The first entry would be the column index
        # of current_V that was moved to the first spot. The second index would be 
        # the column index of current_V that was moved to the second column, etc. 
        current_swaps_indiv[[s]][k] <- jk
        
        # Recording the sign change (if any) that was made
        current_signs_indiv[[s]][k] <- ak
        
        # Removing the already-swapped columns in current_V from the running:
        current_Vs[[1,s]][,jk] <- NA
      }
    }
    
    # Combine the betas 
    
    if (!is.null(betas)) {
      tilde_beta[[1,1]] <- rbind(current_beta[[1,1]][1], tilde_beta_joint[[1,1]], do.call(rbind, tilde_beta_indiv))
    }
    
    # Combine the gammas
    if (!is.null(gammas)) {
      tilde_gamma[[1,1]] <- rbind(current_gamma[[1,1]][1], tilde_gamma_joint[[1,1]], do.call(rbind, tilde_gamma_indiv))
    }
    
    # Storing the swapped samples
    swapped_V.draw[[iter]] <- tilde_V
    swapped_U.draw[[iter]] <- tilde_U
    swapped_W.draw[[iter]] <- tilde_W
    swapped_Vs.draw[[iter]] <- tilde_Vs
    
    if (!is.null(betas)) swapped_betas[[iter]] <- tilde_beta
    
    if (!is.null(gammas)) swapped_gammas[[iter]] <- tilde_gamma
    
    # Storing the swaps made
    swaps_made[[iter]] <- list(current_swaps_joint = current_swaps_joint, current_swaps_indiv = current_swaps_indiv)
    
    # Storing the signs changed
    signs_changed[[iter]] <- list(current_signs_joint = current_signs_joint, current_signs_indiv = current_signs_indiv)
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the posterior mean for each parameter
  # ---------------------------------------------------------------------------
  
  V.mean <- matrix(list(), nrow = 1, ncol = 1)
  V.mean[[1,1]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
    swapped_V.draw[[iter]][[1,1]]
  }))/length(thinned_iters_burnin)
  
  U.mean <- matrix(list(), nrow = q, ncol = 1)
  W.mean <- matrix(list(), nrow = q, ncol = q)
  Vs.mean <- matrix(list(), nrow = 1, ncol = q)
  
  for (s in 1:q) {
    U.mean[[s,1]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_U.draw[[iter]][[s,1]]
    }))/length(thinned_iters_burnin)
    
    Vs.mean[[1,s]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_Vs.draw[[iter]][[1,s]]
    }))/length(thinned_iters_burnin)
    
    W.mean[[s,s]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_W.draw[[iter]][[s,s]]
    }))/length(thinned_iters_burnin)
    
    for (ss in 1:q) {
      if (ss != s) {
        W.mean[[s,ss]] <- swapped_W.draw[[1]][[s,ss]]
      }
    }
  }
  
  if (!is.null(betas)) {
    betas.mean <- matrix(list(), nrow = 1, ncol = 1)
    
    betas.mean[[1,1]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_betas[[iter]][[1,1]]
    }))/length(thinned_iters_burnin)
  }
  
  if (!is.null(gammas)) {
    gammas.mean <- matrix(list(), nrow = 1, ncol = 1)
    
    gammas.mean[[1,1]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_gammas[[iter]][[1,1]]
    }))/length(thinned_iters_burnin)
  }

  # ---------------------------------------------------------------------------
  # Ordering the results in the posterior mean to be ordered most-to-least
  # variance explained
  # ---------------------------------------------------------------------------
  
  # Joint structure
  
  # Calculate the variance explained by each of the r factors (using Frobenius norm)
  var_exp_joint <- c()
  for (i in 1:r) {
    var_exp_joint[i] <- frob(U.mean[[1,1]][,i, drop = FALSE] %*% t(V.mean[[1,1]])[i,, drop = FALSE])
  }
  
  # Order the factors 
  order_joint <- order(var_exp_joint, decreasing = TRUE)
  
  # Reorder the joint and individual factors
  V.mean.reorder <- V.mean; U.mean.reorder <- U.mean
  V.mean.reorder[[1,1]] <- V.mean[[1,1]][,order_joint]
  
  for (s in 1:q) {
    U.mean.reorder[[s,1]] <- U.mean[[s,1]][,order_joint]
  }
  
  # Individual structure
  var_exp_indiv <- order_indiv <- matrix(list(), nrow = q, ncol = 1)
  W.mean.reorder <- W.mean; Vs.mean.reorder <- Vs.mean
  
  for (s in 1:q) {
    for (i in 1:r.vec[s]) {
      var_exp_indiv[[s,1]][i] <- frob(W.mean[[s,s]][,i, drop = FALSE] %*% t(Vs.mean[[1,s]])[i,, drop = FALSE])
    }
    order_indiv[[s,1]] <- order(var_exp_indiv[[s,1]], decreasing = TRUE)
    
    # Order the factors
    W.mean.reorder[[s,s]] <- W.mean[[s,s]][,order_indiv[[s,1]]]
    Vs.mean.reorder[[1,s]] <- Vs.mean[[1,s]][,order_indiv[[s,1]]]
  }
  
  # Betas
  if (!is.null(betas)) {
    betas.mean.reorder <- betas.mean
    betas.mean.indiv <- betas.mean.indiv.reorder <- matrix(list(), nrow = q, ncol = 1)
    
    # Joint effect
    if (r != 0) {
      betas.mean.joint <- betas.mean.reorder[[1,1]][2:(r+1),, drop = FALSE] 
      betas.mean.joint.reorder <- betas.mean.joint[order_joint,, drop = FALSE]
    } else {
      betas.mean.joint <- betas.mean.joint.reorder <-  NULL
    }
    
    # Individual effects
    betas.mean.indiv.temp <- betas.mean.reorder[[1,1]][(r+2):n_beta,, drop = FALSE]
    
    for (s in 1:q) {
      # If there is an individual effect
      if (r.vec[s] != 0) {
        if (s == 1) {
          betas.mean.indiv[[s, 1]] <- betas.mean.indiv.temp[1:r.vec[s],, drop = FALSE] 
        }
        if (s != 1) {
          betas.mean.indiv[[s, 1]] <- betas.mean.indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
        betas.mean.indiv.reorder[[s,1]] <- betas.mean.indiv[[s, 1]][order_indiv[[s,1]],, drop = FALSE]
      }
    }
    
    betas.mean.reorder <- rbind(betas.mean[[1,1]][1,], betas.mean.joint.reorder, do.call(rbind, betas.mean.indiv.reorder))
  }
  
  if (is.null(betas)) {
    betas.mean.reorder <- NULL
  }
  
  # Gammas
  
  if (!is.null(gammas)) {
    gammas.mean.reorder <- gammas.mean
    gammas.mean.indiv <- gammas.mean.indiv.reorder <- matrix(list(), nrow = q, ncol = 1)
    
    # Joint effect
    if (r != 0) {
      gammas.mean.joint <- gammas.mean[[1,1]][2:(r+1),, drop = FALSE]
      gammas.mean.joint.reorder <- gammas.mean.joint[order_joint,, drop = FALSE]
    } else {
      gammas.mean.joint <- gammas.mean.joint.reorder <- NULL
    }
    
    # Individual effects
    gammas.mean.indiv.temp <- gammas.mean[[1,1]][(r+2):n_beta,, drop = FALSE]
    
    for (s in 1:q) {
      # If there is an individual effect
      if (r.vec[s] != 0) {
        if (s == 1) {
          gammas.mean.indiv[[s,1]] <- gammas.mean.indiv.temp[1:r.vec[s],, drop = FALSE] 
        }
        if (s != 1) {
          gammas.mean.indiv[[s,1]] <- gammas.mean.indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
        gammas.mean.indiv.reorder[[s,1]] <- gammas.mean.indiv[[s,1]][order_indiv[[s,1]],, drop = FALSE]
      }
    }
    
    gammas.mean.reorder <- rbind(gammas.mean[[1,1]][1,], gammas.mean.joint.reorder, do.call(rbind, gammas.mean.indiv.reorder))
  }
  
  if (is.null(gammas)) {
    gammas.mean.reorder <- NULL
  }
  
  # Return
  list(swapped_V.draw = swapped_V.draw, 
       swapped_U.draw = swapped_U.draw, 
       swapped_Vs.draw = swapped_Vs.draw, 
       swapped_W.draw = swapped_W.draw,
       swapped_betas = swapped_betas,
       swapped_gammas = swapped_gammas,
       swaps_made = swaps_made, 
       signs_changed = signs_changed,
       V.mean.reorder = V.mean.reorder,
       U.mean.reorder = U.mean.reorder,
       Vs.mean.reorder = Vs.mean.reorder,
       W.mean.reorder = W.mean.reorder,
       betas.mean.reorder = betas.mean.reorder,
       gammas.mean.reorder = gammas.mean.reorder,
       var_exp_joint = var_exp_joint, 
       order_joint = order_joint, 
       var_exp_indiv = var_exp_indiv,
       order_indiv = order_indiv)
}

# Creating an alternative factor switching algorithm
erosheva_curtis_fs <- function(U.draw, V.draw, W.draw, Vs.draw, betas = NULL, gammas = NULL, r, r.vec, n, p.vec, nsample, thinned_iters_burnin, BIDIFAC_solution) {
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # draws (list): the draws to swap according to the pivot
  # rank (int): number of columns to iterate through
  # loadings (list): the corresponding loadings for each draw of the matrix to swap
  # betas (double vec): list of vectors of the coefficients drawn at each iteration
  # gammas (int vec): list of vectors of the inclusion indicators drawn at each iteration
  # r (int) = joint rank
  # r.vec (c(ints)) = vector of individual ranks
  # nsample (int) = number of posterior samples
  # thinned_iters_burnin (c(ints)) = vector of indices for samples after burn-in + thinning if desired
  # nninit (boolean) = TRUE if initialized at theoretical posterior mode or FALSE if use posterior mean as pivot
  # pivot (list) = list of V and Vs to use as pivots. If NULL, must provide nninit. 
  # ---------------------------------------------------------------------------
  
  # Saving the number of sources
  q <- length(r.vec)
  
  # Returning the swapped matrices
  swapped_U.draw <- lapply(1:length(U.draw), function(iter) list())
  swapped_V.draw <- lapply(1:length(V.draw), function(iter) list())
  swapped_W.draw <- lapply(1:length(W.draw), function(iter) list())
  swapped_Vs.draw <- lapply(1:length(Vs.draw), function(iter) list())
  swapped_betas <- lapply(1:length(betas), function(iter) list())
  swapped_gammas <- lapply(1:length(gammas), function(iter) list())
  
  # Returning the swaps made so they can be undone
  swaps_made <- lapply(1:length(swapped_U.draw), function(iter) list())
  
  # Storing the sign changes so they can be undone
  signs_changed <- lapply(1:length(swapped_U.draw), function(iter) list())
  
  # Storing the BIDIFAC+ solutions
  V.bid <- BIDIFAC_solution$V
  U.bid <- BIDIFAC_solution$U
  Vs.bid <- BIDIFAC_solution$Vs
  W.bid <- BIDIFAC_solution$W
  
  # Save the posterior variances
  sd_joint <- sqrt(1/(sqrt(n) + sqrt(sum(p.vec))))
  sd_indiv <- sqrt(sapply(p.vec, function(ps) 1/(sqrt(n) + sqrt(ps))))
  
  for (iter in 1:length(U.draw)) {
    # -------------------------------------------------------------------------
    # Initializing for iter
    # 
    # For each iteration, arrange the columns of each draw to match that 
    # of the pivot. 
    # -------------------------------------------------------------------------
    
    # Initializing the rearranged version
    tilde_V <- matrix(list(), nrow = 1, ncol = 1)
    tilde_U <- matrix(list(), nrow = q, ncol = 1)
    tilde_W <- matrix(list(), nrow = q, ncol = q)
    tilde_Vs <- matrix(list(), nrow = 1, ncol = q)
    
    tilde_beta <- matrix(list(), nrow = 1, ncol = 1)
    tilde_beta_joint <- matrix(list(), nrow = 1, ncol = 1)
    tilde_beta_indiv <- matrix(list(), nrow = q, ncol = 1)
    
    tilde_gamma <- matrix(list(), nrow = 1, ncol = 1)
    tilde_gamma_joint <- matrix(list(), nrow = 1, ncol = 1)
    tilde_gamma_indiv <- matrix(list(), nrow = q, ncol = 1)
    
    tilde_V[[1,1]] <- matrix(nrow = nrow(V.draw[[1]][[1,1]]), ncol = ncol(V.draw[[1]][[1,1]]))
    
    for (s in 1:q) {
      tilde_U[[s,1]] <- matrix(nrow = nrow(U.draw[[1]][[s,1]]), ncol = ncol(U.draw[[1]][[s,1]]))
      tilde_W[[s,s]] <- matrix(nrow = nrow(W.draw[[1]][[s,s]]), ncol = ncol(W.draw[[1]][[s,s]]))
      tilde_Vs[[1,s]] <- matrix(nrow = nrow(Vs.draw[[1]][[1,s]]), ncol = ncol(Vs.draw[[1]][[1,s]]))
      
      for (ss in 1:q) {
        if (ss != s) {
          tilde_W[[s,ss]] <- matrix(0, nrow = nrow(W.draw[[1]][[s,s]]), ncol = ncol(W.draw[[1]][[ss,ss]]))
        }
      }
    }
    
    n_beta <- 1 + r + sum(r.vec)
    tilde_beta[[1,1]] <- matrix(nrow = n_beta, ncol = 1)
    tilde_beta_joint[[1,1]] <- matrix(nrow = r, ncol = 1)
    
    tilde_gamma[[1,1]] <- matrix(nrow = n_beta, ncol = 1)
    tilde_gamma_joint[[1,1]] <- matrix(nrow = r, ncol = 1)
    
    for (s in 1:q) {
      tilde_beta_indiv[[s,1]] <- matrix(nrow = r.vec[s], ncol = 1)
      tilde_gamma_indiv[[s,1]] <- matrix(nrow = r.vec[s], ncol = 1)
    }
    
    # -------------------------------------------------------------------------
    # Storing the current V, U, beta, and gamma
    # -------------------------------------------------------------------------
    
    current_V <- V.draw[[iter]]
    current_U <- U.draw[[iter]]
    current_Vs <- Vs.draw[[iter]]
    current_W <- W.draw[[iter]]
    
    if (!is.null(betas)) {
      current_beta <- betas[[iter]]
      current_beta_indiv <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        current_beta_joint <- current_beta[[1,1]][2:(r+1),, drop = FALSE] 
      } else {
        current_beta_joint <- NULL
      }
      
      # Individual effects
      current_beta_indiv.temp <- current_beta[[1,1]][(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) {
            current_beta_indiv[[s, 1]] <- current_beta_indiv.temp[1:r.vec[s],, drop = FALSE] 
          }
          if (s != 1) {
            current_beta_indiv[[s, 1]] <- current_beta_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
      }
    }
    
    if (!is.null(gammas)) {
      current_gamma <- gammas[[iter]] 
      current_gamma_indiv <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        current_gamma_joint <- current_gamma[[1,1]][2:(r+1),, drop = FALSE]
      } else {
        current_gamma_joint <- NULL
      }
      
      # Individual effects
      current_gamma_indiv.temp <- current_gamma[[1,1]][(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) {
            current_gamma_indiv[[s,1]] <- current_gamma_indiv.temp[1:r.vec[s],, drop = FALSE] 
          }
          if (s != 1) {
            current_gamma_indiv[[s,1]] <- current_gamma_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Finding the optimal order
    # -------------------------------------------------------------------------
    
    # Storing the swaps and signs
    current_swaps_joint <- current_signs_joint <- c()
    current_swaps_indiv <- current_signs_indiv <- lapply(1:q, function(s) c())
    
    # Start with the joint structure
    
    for (k in 1:r) {
      # For each factor, calculate the log-density given the BIDIFAC+ solution for both signs
      dens_pos <- apply(current_V[[1,1]], 2, function(current_col) -sum(dnorm(current_col, mean = V.bid[[1,1]][,k], sd = rep(sd_joint, n), log = TRUE))) # corrs <- apply(current_V[[1,1]], 2, function(current_col) cor(current_col, pivot_V[[1,1]][,k]))
      dens_neg <- apply(current_V[[1,1]], 2, function(current_col) -sum(dnorm(-current_col, mean = V.bid[[1,1]][,k], sd = rep(sd_joint, n), log = TRUE))) # corrs <- apply(current_V[[1,1]], 2, function(current_col) cor(current_col, pivot_V[[1,1]][,k]))
      
      # Which sign had the higher negative log density?
      ak <- if (max(dens_pos, dens_neg, na.rm = TRUE) %in% dens_pos) 1 else -1
      
      # Save the index if the sign change was pos or neg
      if (ak == 1) {
        jk <- which.max(abs(dens_pos)) 
      }
      
      if (ak == -1) {
        jk <- which.max(abs(dens_neg)) 
      }
      
      # Setting the kth column of tilde_V to be the ak*current_V[,jk] column. 
      tilde_V[[1,1]][,k] <- ak*current_V[[1,1]][,jk]
      
      # Rearranging the corresponding column in U:
      for (s in 1:q) {
        tilde_U[[s,1]][,k] <- ak*current_U[[s,1]][,jk]
      }
      
      # Rearranging the corresponding rows in betas:
      if (!is.null(betas)) {
        tilde_beta_joint[[1,1]][k,] <- ak*current_beta_joint[jk,]
      }
      
      # Rearranging the corresponding rows in the gammas:
      if (!is.null(gammas)) {
        tilde_gamma_joint[[1,1]][k,] <- current_gamma_joint[jk,]
      }
      
      # Recording what swap was made. The first entry would be the column index
      # of current_V that was moved to the first spot. The second index would be 
      # the column index of current_V that was moved to the second column, etc. 
      current_swaps_joint[k] <- jk
      
      # Recording the sign change (if any) that was made
      current_signs_joint[k] <- ak
      
      # Removing the already-swapped columns in current_V from the running:
      current_V[[1,1]][,jk] <- NA
    }
    
    # Then the individual structures
    
    for (s in 1:q) {
      for (k in 1:r.vec[s]) {
        # For each factor, calculate the log-density given the BIDIFAC+ solution for both signs
        dens_pos <- apply(current_Vs[[1,s]], 2, function(current_col) -sum(dnorm(current_col, mean = Vs.bid[[1,s]][,k], sd = rep(sd_indiv[s], n), log = TRUE))) # corrs <- apply(current_Vs[[1,s]], 2, function(current_col) cor(current_col, pivot_Vs[[1,s]][,k])) 
        dens_neg <- apply(current_Vs[[1,s]], 2, function(current_col) -sum(dnorm(-current_col, mean = Vs.bid[[1,s]][,k], sd = rep(sd_indiv[s], n), log = TRUE))) # corrs <- apply(current_Vs[[1,s]], 2, function(current_col) cor(current_col, pivot_Vs[[1,s]][,k])) 
        
        # Which sign had the higher negative log density?
        ak <- if (max(dens_pos, dens_neg, na.rm = TRUE) %in% dens_pos) 1 else -1
        
        # Save the index if the sign change was pos or neg
        if (ak == 1) {
          jk <- which.max(abs(dens_pos)) 
        }
        
        if (ak == -1) {
          jk <- which.max(abs(dens_neg)) 
        }
        
        # Setting the kth column of \tilde_Vs to be the ak*current_Vs[,jk] column. 
        tilde_Vs[[1,s]][,k] <- ak*current_Vs[[1,s]][,jk]
        
        # Rearranging the corresponding column in W:
        tilde_W[[s,s]][,k] <- ak*current_W[[s,s]][,jk]
        
        # Rearranging the corresponding rows in betas:
        if (!is.null(betas)) {
          tilde_beta_indiv[[s,1]][k,] <- ak*current_beta_indiv[[s,1]][jk,]
        }
        
        # Rearranging the corresponding rows in the gammas:
        if (!is.null(gammas)) {
          tilde_gamma_indiv[[s,1]][k,] <- current_gamma_indiv[[s,1]][jk,]
        }
        
        # Recording what swap was made. The first entry would be the column index
        # of current_V that was moved to the first spot. The second index would be 
        # the column index of current_V that was moved to the second column, etc. 
        current_swaps_indiv[[s]][k] <- jk
        
        # Recording the sign change (if any) that was made
        current_signs_indiv[[s]][k] <- ak
        
        # Removing the already-swapped columns in current_V from the running:
        current_Vs[[1,s]][,jk] <- NA
      }
    }
    
    # Combine the betas 
    
    if (!is.null(betas)) {
      tilde_beta[[1,1]] <- rbind(current_beta[[1,1]][1], tilde_beta_joint[[1,1]], do.call(rbind, tilde_beta_indiv))
    }
    
    # Combine the gammas
    if (!is.null(gammas)) {
      tilde_gamma[[1,1]] <- rbind(current_gamma[[1,1]][1], tilde_gamma_joint[[1,1]], do.call(rbind, tilde_gamma_indiv))
    }
    
    # Storing the swapped samples
    swapped_V.draw[[iter]] <- tilde_V
    swapped_U.draw[[iter]] <- tilde_U
    swapped_W.draw[[iter]] <- tilde_W
    swapped_Vs.draw[[iter]] <- tilde_Vs
    
    if (!is.null(betas)) swapped_betas[[iter]] <- tilde_beta
    
    if (!is.null(gammas)) swapped_gammas[[iter]] <- tilde_gamma
    
    # Storing the swaps made
    swaps_made[[iter]] <- list(current_swaps_joint = current_swaps_joint, current_swaps_indiv = current_swaps_indiv)
    
    # Storing the signs changed
    signs_changed[[iter]] <- list(current_signs_joint = current_signs_joint, current_signs_indiv = current_signs_indiv)
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the posterior mean for each parameter
  # ---------------------------------------------------------------------------
  
  V.mean <- matrix(list(), nrow = 1, ncol = 1)
  V.mean[[1,1]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
    swapped_V.draw[[iter]][[1,1]]
  }))/length(thinned_iters_burnin)
  
  U.mean <- matrix(list(), nrow = q, ncol = 1)
  W.mean <- matrix(list(), nrow = q, ncol = q)
  Vs.mean <- matrix(list(), nrow = 1, ncol = q)
  
  for (s in 1:q) {
    U.mean[[s,1]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_U.draw[[iter]][[s,1]]
    }))/length(thinned_iters_burnin)
    
    Vs.mean[[1,s]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_Vs.draw[[iter]][[1,s]]
    }))/length(thinned_iters_burnin)
    
    W.mean[[s,s]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_W.draw[[iter]][[s,s]]
    }))/length(thinned_iters_burnin)
    
    for (ss in 1:q) {
      if (ss != s) {
        W.mean[[s,ss]] <- swapped_W.draw[[1]][[s,ss]]
      }
    }
  }
  
  if (!is.null(betas)) {
    betas.mean <- matrix(list(), nrow = 1, ncol = 1)
    
    betas.mean[[1,1]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_betas[[iter]][[1,1]]
    }))/length(thinned_iters_burnin)
  }
  
  if (!is.null(gammas)) {
    gammas.mean <- matrix(list(), nrow = 1, ncol = 1)
    
    gammas.mean[[1,1]] <- Reduce("+", lapply(thinned_iters_burnin, function(iter) {
      swapped_gammas[[iter]][[1,1]]
    }))/length(thinned_iters_burnin)
  }
  
  # ---------------------------------------------------------------------------
  # Ordering the results in the posterior mean to be ordered most-to-least
  # variance explained
  # ---------------------------------------------------------------------------
  
  # Joint structure
  
  # Calculate the variance explained by each of the r factors (using Frobenius norm)
  var_exp_joint <- c()
  for (i in 1:r) {
    var_exp_joint[i] <- frob(U.mean[[1,1]][,i, drop = FALSE] %*% t(V.mean[[1,1]])[i,, drop = FALSE])
  }
  
  # Order the factors 
  order_joint <- order(var_exp_joint, decreasing = TRUE)
  
  # Reorder the joint and individual factors
  V.mean.reorder <- V.mean; U.mean.reorder <- U.mean
  V.mean.reorder[[1,1]] <- V.mean[[1,1]][,order_joint]
  
  for (s in 1:q) {
    U.mean.reorder[[s,1]] <- U.mean[[s,1]][,order_joint]
  }
  
  # Individual structure
  var_exp_indiv <- order_indiv <- matrix(list(), nrow = q, ncol = 1)
  W.mean.reorder <- W.mean; Vs.mean.reorder <- Vs.mean
  
  for (s in 1:q) {
    for (i in 1:r.vec[s]) {
      var_exp_indiv[[s,1]][i] <- frob(W.mean[[s,s]][,i, drop = FALSE] %*% t(Vs.mean[[1,s]])[i,, drop = FALSE])
    }
    order_indiv[[s,1]] <- order(var_exp_indiv[[s,1]], decreasing = TRUE)
    
    # Order the factors
    W.mean.reorder[[s,s]] <- W.mean[[s,s]][,order_indiv[[s,1]]]
    Vs.mean.reorder[[1,s]] <- Vs.mean[[1,s]][,order_indiv[[s,1]]]
  }
  
  # Betas
  if (!is.null(betas)) {
    betas.mean.reorder <- betas.mean
    betas.mean.indiv <- betas.mean.indiv.reorder <- matrix(list(), nrow = q, ncol = 1)
    
    # Joint effect
    if (r != 0) {
      betas.mean.joint <- betas.mean.reorder[[1,1]][2:(r+1),, drop = FALSE] 
      betas.mean.joint.reorder <- betas.mean.joint[order_joint,, drop = FALSE]
    } else {
      betas.mean.joint <- betas.mean.joint.reorder <-  NULL
    }
    
    # Individual effects
    betas.mean.indiv.temp <- betas.mean.reorder[[1,1]][(r+2):n_beta,, drop = FALSE]
    
    for (s in 1:q) {
      # If there is an individual effect
      if (r.vec[s] != 0) {
        if (s == 1) {
          betas.mean.indiv[[s, 1]] <- betas.mean.indiv.temp[1:r.vec[s],, drop = FALSE] 
        }
        if (s != 1) {
          betas.mean.indiv[[s, 1]] <- betas.mean.indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
        betas.mean.indiv.reorder[[s,1]] <- betas.mean.indiv[[s, 1]][order_indiv[[s,1]],, drop = FALSE]
      }
    }
    
    betas.mean.reorder <- rbind(betas.mean[[1,1]][1,], betas.mean.joint.reorder, do.call(rbind, betas.mean.indiv.reorder))
  }
  
  if (is.null(betas)) {
    betas.mean.reorder <- NULL
  }
  
  # Gammas
  
  if (!is.null(gammas)) {
    gammas.mean.reorder <- gammas.mean
    gammas.mean.indiv <- gammas.mean.indiv.reorder <- matrix(list(), nrow = q, ncol = 1)
    
    # Joint effect
    if (r != 0) {
      gammas.mean.joint <- gammas.mean[[1,1]][2:(r+1),, drop = FALSE]
      gammas.mean.joint.reorder <- gammas.mean.joint[order_joint,, drop = FALSE]
    } else {
      gammas.mean.joint <- gammas.mean.joint.reorder <- NULL
    }
    
    # Individual effects
    gammas.mean.indiv.temp <- gammas.mean[[1,1]][(r+2):n_beta,, drop = FALSE]
    
    for (s in 1:q) {
      # If there is an individual effect
      if (r.vec[s] != 0) {
        if (s == 1) {
          gammas.mean.indiv[[s,1]] <- gammas.mean.indiv.temp[1:r.vec[s],, drop = FALSE] 
        }
        if (s != 1) {
          gammas.mean.indiv[[s,1]] <- gammas.mean.indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
        }
        gammas.mean.indiv.reorder[[s,1]] <- gammas.mean.indiv[[s,1]][order_indiv[[s,1]],, drop = FALSE]
      }
    }
    
    gammas.mean.reorder <- rbind(gammas.mean[[1,1]][1,], gammas.mean.joint.reorder, do.call(rbind, gammas.mean.indiv.reorder))
  }
  
  if (is.null(gammas)) {
    gammas.mean.reorder <- NULL
  }
  
  # Return
  list(swapped_V.draw = swapped_V.draw, 
       swapped_U.draw = swapped_U.draw, 
       swapped_Vs.draw = swapped_Vs.draw, 
       swapped_W.draw = swapped_W.draw,
       swapped_betas = swapped_betas,
       swapped_gammas = swapped_gammas,
       swaps_made = swaps_made, 
       signs_changed = signs_changed,
       V.mean.reorder = V.mean.reorder,
       U.mean.reorder = U.mean.reorder,
       Vs.mean.reorder = Vs.mean.reorder,
       W.mean.reorder = W.mean.reorder,
       betas.mean.reorder = betas.mean.reorder,
       gammas.mean.reorder = gammas.mean.reorder,
       var_exp_joint = var_exp_joint, 
       order_joint = order_joint, 
       var_exp_indiv = var_exp_indiv,
       order_indiv = order_indiv)
}

# Papastamoulis and Ntzoufras algorithm for factor switching (from the factor.switching package)
rsp_exact <- function (lambda_mcmc, maxIter = 100, threshold = 1e-06, verbose = TRUE, 
                       rotate = TRUE, printIter = 1000) {
  lambda_mcmc <- as.matrix(lambda_mcmc)
  cnames <- colnames(lambda_mcmc)
  d <- dim(lambda_mcmc)[2]
  if (sum(grepl("LambdaV", cnames)) < d) {
    stop("Column names of lambda_mcmc should be formatted as: LambdaV1_1 ... LambdaV1_q   ... LambdaVp_1 ... LambdaVp_q, where p and q denote the number of variables and factors, respectively.")
  }
  q <- as.numeric(strsplit(cnames[d], split = "_")[[1]][2])
  p <- d/q
  mcmcIterations <- dim(lambda_mcmc)[1]
  threshold <- threshold * mcmcIterations * p * q
  checkNames <- paste0(rep(paste0("LambdaV", 1:p, "_"), each = q), 
                       1:q)
  if (sum(checkNames == cnames) != d) {
    stop("Column names of lambda_mcmc are wrong, please check input.")
  }
  lambda_mcmc_varimax <- lambda_mcmc
  if (rotate) {
    if (q > 1) {
      for (iter in 1:mcmcIterations) {
        lambda <- matrix(lambda_mcmc[iter, ], nrow = p, 
                         byrow = T)
        v <- varimax(lambda, normalize = F)
        rotated_lambda <- v$loadings
        class(rotated_lambda) <- "matrix"
        lambda_mcmc_varimax[iter, ] <- c(t(rotated_lambda))
      }
    }
  }
  all_c <- array(data = 1, dim = c(2^q, q))
  l <- 1
  if (q > 1) {
    for (i in 1:(q - 1)) {
      pos <- combn(1:q, i)
      j <- dim(pos)[2]
      for (k in 1:j) {
        l <- l + 1
        all_c[l, pos[, k]] <- rep(-1, i)
      }
    }
  }
  all_c[l + 1, ] <- rep(-1, q)
  lambda_hat <- matrix(colMeans(lambda_mcmc_varimax), ncol = q, 
                       nrow = p, byrow = TRUE)
  lambda_hat_zero <- lambda_hat
  c_vectors <- matrix(rep(1, q), ncol = q, nrow = mcmcIterations, 
                      byrow = TRUE)
  v_vectors <- matrix(1:q, ncol = q, nrow = mcmcIterations, 
                      byrow = TRUE)
  st <- 1:q
  dim_all_c <- 2^q
  dim_all_v <- factorial(q)
  perm <- matrix(1:q, ncol = q, nrow = dim_all_c, byrow = TRUE)
  costs <- numeric(dim_all_c)
  f <- numeric(maxIter)
  lambda_hat_values <- array(data = NA, dim = c(p, q, maxIter))
  totalIterations <- 0
  criterion = TRUE
  cost.matrix <- matrix(numeric(q * q), nrow = q, ncol = q)
  t1 <- numeric(maxIter)
  start_time <- Sys.time()
  while ((criterion == TRUE) & (totalIterations < maxIter)) {
    totalIterations <- totalIterations + 1
    cat(paste0("* iteration: ", totalIterations), "\n")
    lambda_hat_new <- 0 * lambda_hat
    objective <- 0
    for (iter in 1:mcmcIterations) {
      lambda <- matrix(lambda_mcmc_varimax[iter, ], ncol = q, 
                       nrow = p, byrow = TRUE)
      for (i in 1:dim_all_c) {
        c_vec <- all_c[i, ]
        lambda_switch <- matrix(c_vec, ncol = q, nrow = p, 
                                byrow = T) * lambda_hat
        for (j in 1:q) {
          temp <- (lambda - lambda_switch[, j])^2
          cost.matrix[j, ] <- colSums(temp)
        }
        matr <- lp.assign(cost.matrix)$solution
        for (j in 1:q) {
          perm[i, j] <- st[matr[, j] > 0]
        }
        perm[i, ] <- order(perm[i, ])
        costs[i] <- sum(cost.matrix * matr)
      }
      minIndex <- order(costs)[1]
      v_vectors[iter, ] <- perm[minIndex, ]
      c_vectors[iter, ] <- all_c[minIndex, ]
      switchedMatrix <- matrix(c_vectors[iter, ], ncol = q, 
                               nrow = p, byrow = T) * lambda[, v_vectors[iter, 
                               ]]
      objective <- objective + costs[minIndex]
      lambda_hat_new <- lambda_hat_new + switchedMatrix
      if ((iter%%printIter == 0) && (verbose == TRUE)) {
        cat(paste0("          mcmc draw = ", iter, ":  sum f = ", 
                   round(objective, 3)), "\r")
      }
    }
    f[totalIterations] <- objective
    if (totalIterations > 1) {
      if (f[totalIterations - 1] - f[totalIterations] < 
          threshold) {
        criterion = FALSE
      }
    }
    if (verbose == TRUE) {
      cat(paste0("   -----  objective function = ", round(objective, 
                                                          3), " -----"), "\n")
      cat("\n")
    }
    lambda_hat_new <- lambda_hat_new/mcmcIterations
    lambda_hat <- lambda_hat_new
    lambda_hat_values[, , totalIterations] <- lambda_hat
    end_time <- Sys.time()
    t1[totalIterations] <- as.numeric(difftime(end_time, 
                                               start_time, units = "min"))
  }
  c_vec <- rep(1, q)
  v_vec <- 1:q
  f_zero <- 0
  for (i in 1:mcmcIterations) {
    lambda <- matrix(lambda_mcmc_varimax[iter, ], ncol = q, 
                     nrow = p, byrow = TRUE)
    switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow = T) * 
      lambda[, v_vec]
    f_zero <- f_zero + sum((switchedMatrix - lambda_hat)^2)
  }
  t_exact <- c(0, t1[1:totalIterations])
  f_exact <- c(f_zero, f[1:totalIterations])
  objective_function <- data.frame(time = t_exact, value = f_exact)
  lambda_reordered_mcmc <- lambda_mcmc
  for (i in 1:mcmcIterations) {
    lambda <- matrix(lambda_mcmc_varimax[i, ], ncol = q, 
                     nrow = p, byrow = TRUE)
    switchedMatrix <- matrix(c_vectors[i, ], ncol = q, nrow = p, 
                             byrow = T) * lambda[, v_vectors[i, ]]
    lambda_reordered_mcmc[i, ] <- c(t(switchedMatrix))
  }
  lambda_reordered_mcmc <- as.mcmc(lambda_reordered_mcmc)
  result <- vector("list", length = 5)
  result[[1]] <- lambda_reordered_mcmc
  result[[2]] <- c_vectors
  result[[3]] <- v_vectors
  result[[4]] <- lambda_hat
  result[[5]] <- objective_function
  names(result) <- c("lambda_reordered_mcmc", "sign_vectors", 
                     "permute_vectors", "lambda_hat", "objective_function")
  class(result) <- c("list", "rsp")
  return(result)
}

# Modying the Papastamoulis and Ntzoufras algorithm
rsp_test <- function(U.draw, V.draw, W.draw, Vs.draw, betas.draw = NULL, gammas.draw = NULL, r, r.vec, BIDIFAC_solution, maxIter = 10) {
  
  # ---------------------------------------------------------------------------
  # Arguments:
  # ---------------------------------------------------------------------------
  
  # Load in solver for identifying best sign swapping solution
  library(lpSolve)
  
  # Save basic information
  nsample <- length(U.draw)
  n <- nrow(V.draw[[1]][[1,1]])
  p.vec <- sapply(U.draw[[1]], nrow)
  q <- length(p.vec)
  n_beta <- 1 + r + sum(r.vec)
  
  # Save the BIDIFAC solution
  V.bid <- BIDIFAC_solution$V
  U.bid <- BIDIFAC_solution$U
  W.bid <- BIDIFAC_solution$W
  Vs.bid <- BIDIFAC_solution$Vs
  
  # Create lists for swapped factors and loadings
  U.draw.swapped <- U.draw
  V.draw.swapped <- V.draw
  W.draw.swapped <- W.draw
  Vs.draw.swapped <- Vs.draw
  
  # Create lists for the swapped betas and gammas if provided
  betas.swapped <- betas.draw
  gammas.swapped <- gammas.draw
  
  betas.joint.swapped <- gammas.joint.swapped <- lapply(1:nsample, function(iter) matrix(list(), nrow = 1, ncol = 1))
  betas.indiv.swapped <- gammas.indiv.swapped <- lapply(1:nsample, function(iter) matrix(list(), nrow = q, ncol = 1))
  
  # ---------------------------------------------------------------------------
  # Permute the joint factors
  # ---------------------------------------------------------------------------
  
  # Enumerate all possible sign changes
  all_signs <- array(data = 1, dim = c(2^r, r))
  l <- 1
  if (r > 1) {
    for (i in 1:(r - 1)) {
      pos <- combn(1:r, i)
      j <- dim(pos)[2]
      for (k in 1:j) {
        l <- l + 1
        all_signs[l, pos[, k]] <- rep(-1, i)
      }
    }
  }
  all_signs[l + 1, ] <- rep(-1, r)
  
  # Create a cost matrix for V and Vs
  cost.matrix.V <- matrix(numeric(r*r), nrow = r, ncol = r)
  cost.matrix.Vs <- lapply(1:q, function(s) matrix(numeric(r.vec[s]*r.vec[s]), nrow = r.vec[s], ncol = r.vec[s]))
  
  # Number of possible sign combinations
  dim_all_r_for_V <- 2^r
  dim_all_r_for_Vs <- 2^r.vec
  
  # Number of column permutations
  dim_all_perm_for_V <- factorial(r)
  dim_all_perm_for_Vs <- factorial(r.vec)
  
  # Index for each factor
  ind_V <- 1:r
  ind_Vs <- lapply(1:q, function(s) 1:r.vec[s])
  
  # Matrix of possible sign changes
  perm <- matrix(1:r, ncol = r, nrow = dim_all_r_for_V, byrow = TRUE) 
  
  # Cost for a particular sign change
  costs <- numeric(dim_all_r_for_V) 
  
  # Sign changes for each iteration and for each factor
  sign_vectors <- matrix(rep(1, r), ncol = r, nrow = nsample, byrow = TRUE) # c_vectors
  
  # Permutations for each iteration and for each factor
  perm_vectors <- matrix(1:r, ncol = r, nrow = nsample, byrow = TRUE) # v_vectors
  
  # Vector to monitor convergence
  objective_vec <- numeric(maxIter)
  
  # Set up for while loop
  criterion <- TRUE
  total.iter <- 1
  
  # Iterate until convergence
  while (criterion & (total.iter <= maxIter)) {
    
    # Initialize the value of our objective function to monitor convergence
    objective <- 0
    
    # Iterate through MCMC iterations 
    for (iter in 1:nsample) {
      
      # Save current scores and loadings
      current.V <- V.draw[[iter]]
      current.U <- U.draw[[iter]]
      
      # Initialize beta.joint and gamma.joint
      if (!is.null(betas.draw)) {
        beta.joint <- betas.draw[[iter]][[1,1]][2:(r+1),,drop = FALSE]
      }
      
      if (!is.null(gammas.draw)) {
        gamma.joint <- gammas.draw[[iter]][[1,1]][2:(r+1),,drop = FALSE]
      }
      
      # Iterate through all possible sign changes
      for (i in 1:dim_all_r_for_V) {
        # Save the current sign 
        sign_vec <- all_signs[i, ]
        
        # Try a new sign
        V.bid_switch <- matrix(sign_vec, ncol = r, nrow = n, byrow = T) * V.bid[[1,1]]
        
        # Calculate the cost of this sign change
        for (j in 1:r) {
          temp <- (current.V[[1,1]] - V.bid_switch[, j])^2
          cost.matrix.V[j, ] <- colSums(temp)
        }
        
        # Obtain the optimal sign
        matr <- lp.assign(cost.matrix.V)$solution
        for (j in 1:q) {
          perm[i, j] <- ind_V[matr[, j] > 0]
        }
        # Now we have the best permutation for the ith sign change
        
        # Now we find the index that corresponds to this permutation
        perm[i, ] <- order(perm[i, ])
        
        # Then we compute the total cost of this permutation for this sign change?
        costs[i] <- sum(cost.matrix.V * matr)
      }
      
      # Find the index of the minimum cost for the sign changes
      minIndex <- order(costs)[1]
      
      # Find the permutation that corresponds to this minimum cost sign change for this iteration
      perm_vectors[iter, ] <- perm[minIndex, ]
      
      # Find the sign change that corresponds to this minimum cost sign change for this iteration
      sign_vectors[iter, ] <- all_signs[minIndex, ]
      
      # Swap the joint factors
      V.draw.swapped[[iter]][[1,1]] <- matrix(sign_vectors[iter, ], ncol = r, nrow = n, byrow = T) * current.V[[1,1]][, perm_vectors[iter,]]
      
      for (s in 1:q) {
        U.draw.swapped[[iter]][[s,1]] <- matrix(sign_vectors[iter, ], ncol = r, nrow = p.vec[s], byrow = T) * current.U[[s,1]][, perm_vectors[iter,]]
      }
      
      # Swap the betas for joint factors
      if (!is.null(betas.draw)) {
        betas.joint.swapped[[iter]][[1,1]] <- sign_vectors[iter, ] * beta.joint[perm_vectors[iter,],,drop=FALSE]
      }
      
      # Swap the gammas for joint factors
      if (!is.null(gammas.draw)) {
        gammas.joint.swapped[[iter]][[1,1]] <- gamma.joint[perm_vectors[iter,],,drop=FALSE]
      }
      
      # Update the objective
      objective <- objective + costs[minIndex]
    }
    
    # Check for convergence
    objective_vec[total.iter] <- objective
    if (total.iter > 1) {
      if (objective_vec[total.iter - 1] - objective_vec[total.iter] < threshold) criterion = FALSE
    }

    # Tabulating the total number of iterations
    total.iter <- total.iter + 1
  }
  
  # ---------------------------------------------------------------------------
  # Permute the individual factors
  # ---------------------------------------------------------------------------
  
  # Create a cost matrix for V and Vs
  cost.matrix.Vs <- lapply(1:q, function(s) matrix(numeric(r.vec[s]*r.vec[s]), nrow = r.vec[s], ncol = r.vec[s]))
  
  # Number of possible sign combinations
  dim_all_r_for_Vs <- 2^r.vec
  
  # Number of column permutations
  dim_all_perm_for_Vs <- factorial(r.vec)
  
  # Index for each factor
  ind_Vs <- lapply(1:q, function(s) 1:r.vec[s])
  
  # Matrix of possible sign changes
  perm_Vs <- lapply(1:q, function(s) matrix(1:r.vec[s], ncol = r.vec[s], nrow = dim_all_r_for_Vs[s], byrow = TRUE) )
  
  # Cost for a particular sign change
  costs_Vs <- lapply(1:q, function(s) numeric(dim_all_r_for_Vs[s]) )
  
  # Sign changes for each iteration and for each factor
  sign_vectors_Vs <- lapply(1:q, function(s) matrix(rep(1, r.vec[s]), ncol = r.vec[s], nrow = nsample, byrow = TRUE)) # c_vectors
  
  # Permutations for each iteration and for each factor
  perm_vectors_Vs <- lapply(1:q, function(s) matrix(1:r.vec[s], ncol = r.vec[s], nrow = nsample, byrow = TRUE)) # v_vectors
  
  # For each source, iterate until convergence
  for (s in 1:q) {
    # Enumerate all possible sign changes
    all_signs_s <- array(data = 1, dim = c(2^r.vec[s], r.vec[s]))
    l <- 1
    if (r.vec[s] > 1) {
      for (i in 1:(r.vec[s] - 1)) {
        pos <- combn(1:r.vec[s], i)
        j <- dim(pos)[2]
        for (k in 1:j) {
          l <- l + 1
          all_signs_s[l, pos[, k]] <- rep(-1, i)
        }
      }
    }
    all_signs_s[l + 1, ] <- rep(-1, r.vec[s])
    
    # Vector to monitor convergence
    objective_vec <- numeric(maxIter)
    
    # Set up for while loop
    criterion <- TRUE
    total.iter <- 1
    
    while (criterion & (total.iter <= maxIter)) {
      
      # Initialize the value of our objective function to monitor convergence
      objective <- 0
      
      # Iterate through MCMC iterations 
      for (iter in 1:nsample) {
        
        # Save current scores and loadings
        current.Vs <- Vs.draw[[iter]][[1,s]]
        current.W <- W.draw[[iter]][[s,s]]
        
        # Initialize beta.indiv.s
        if (!is.null(betas.draw)) {
          # Select all the betas for individual factors
          current_betas_indiv.temp <- betas.draw[[iter]][[1,1]][(r+2):n_beta,, drop = FALSE]
          
          # If we are at the first source
          if (s == 1) {
            beta.indiv.s <- current_betas_indiv.temp[1:r.vec[s],, drop = FALSE]
          }
          
          # If we are not at the first source
          if (s != 1) {
            beta.indiv.s <- current_betas_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
        
        # Initialize gamma.indiv.s
        if (!is.null(gammas.draw)) {
          # Select all the gammas for individual factors
          current_gammas_indiv.temp <- gammas.draw[[iter]][[1,1]][(r+2):n_beta,, drop = FALSE]
          
          # If we are at the first source
          if (s == 1) {
            gamma.indiv.s <- current_gammas_indiv.temp[1:r.vec[s],, drop = FALSE]
          }
          
          # If we are not at the first source
          if (s != 1) {
            gamma.indiv.s <- current_gammas_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
        
        # Iterate through all possible sign changes
        for (i in 1:dim_all_r_for_Vs[s]) {
          # Save the current sign 
          sign_vec <- all_signs_s[i, ]
          
          # Try a new sign
          Vs.bid_switch <- matrix(sign_vec, ncol = r.vec[s], nrow = n, byrow = T) * Vs.bid[[1,s]]
          
          # Calculate the cost of this sign change
          for (j in 1:r.vec[s]) {
            temp <- (current.Vs - Vs.bid_switch[, j])^2
            cost.matrix.Vs[[s]][j, ] <- colSums(temp)
          }
          
          # Obtain the optimal sign
          matr <- lp.assign(cost.matrix.Vs[[s]])$solution
          for (j in 1:q) {
            perm_Vs[[s]][i, j] <- ind_Vs[[s]][matr[, j] > 0]
          }
          # Now we have the best permutation for the ith sign change
          
          # Now we find the index that corresponds to this permutation
          perm_Vs[[s]][i, ] <- order(perm_Vs[[s]][i, ])
          
          # Then we compute the total cost of this permutation for this sign change?
          costs_Vs[[s]][i] <- sum(cost.matrix.Vs[[s]] * matr)
        }
        
        # Find the index of the minimum cost for the sign changes
        minIndex <- order(costs_Vs[[s]])[1]
        
        # Find the permutation that corresponds to this minimum cost sign change for this iteration
        perm_vectors_Vs[[s]][iter, ] <- perm_Vs[[s]][minIndex, ]
        
        # Find the sign change that corresponds to this minimum cost sign change for this iteration
        sign_vectors_Vs[[s]][iter, ] <- all_signs_s[minIndex, ]
        
        # Swap the individual factors
        Vs.draw.swapped[[iter]][[1,s]] <- matrix(sign_vectors_Vs[[s]][iter, ], ncol = r.vec[s], nrow = n, byrow = T) * current.Vs[, perm_vectors_Vs[[s]][iter,]]
        W.draw.swapped[[iter]][[s,s]] <- matrix(sign_vectors_Vs[[s]][iter, ], ncol = r.vec[s], nrow = p.vec[s], byrow = T) * current.W[, perm_vectors_Vs[[s]][iter,]]
        
        # Swap the betas for the individual factors
        if (!is.null(betas.draw)) {
          betas.indiv.swapped[[iter]][[s,1]] <- sign_vectors_Vs[[s]][iter, ] * beta.indiv.s[perm_vectors[iter,],,drop=FALSE]
        }
        
        # Swap the gammas for the individual factors
        if (!is.null(gammas.draw)) {
          gammas.indiv.swapped[[iter]][[s,1]] <- gamma.indiv.s[perm_vectors[iter,],,drop=FALSE]
        }
        
        # Update the objective
        objective <- objective + costs[minIndex]
      }
      
      # Check for convergence
      objective_vec[total.iter] <- objective
      if (total.iter > 1) {
        if (objective_vec[total.iter - 1] - objective_vec[total.iter] < threshold) criterion = FALSE
      }
      
      # Tabulating the total number of iterations
      total.iter <- total.iter + 1
    }
  }
  
  # Combine the betas and gammas together
  for (iter in 1:nsample) {
    betas.swapped[[iter]][[1,1]] <- rbind(betas.swapped[[iter]][[1,1]][1,], betas.joint.swapped[[iter]][[1,1]], do.call(rbind, betas.indiv.swapped[[iter]]))
    gammas.swapped[[iter]][[1,1]] <- rbind(gammas.swapped[[iter]][[1,1]][1,], gammas.joint.swapped[[iter]][[1,1]], do.call(rbind, gammas.indiv.swapped[[iter]]))
  }
  
  # Return the swapped matrices
  list(V.draw.swapped = V.draw.swapped,
       U.draw.swapped = U.draw.swapped,
       W.draw.swapped = W.draw.swapped,
       Vs.draw.swapped = Vs.draw.swapped,
       betas.swapped = betas.swapped,
       gammas.swapped = gammas.swapped)
}

# Creating a new rotation-sign-permutation algorithm based on the Kabsch Algorithm
factor_switching_plus_rotation <- function(U.draw, V.draw, W.draw, Vs.draw, betas = NULL, gammas = NULL, r, r.vec, nsample, thinned_iters_burnin, nninit = TRUE, pivot = NULL) {
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # draws (list): the draws to swap according to the pivot
  # rank (int): number of columns to iterate through
  # loadings (list): the corresponding loadings for each draw of the matrix to swap
  # betas (double vec): list of vectors of the coefficients drawn at each iteration
  # gammas (int vec): list of vectors of the inclusion indicators drawn at each iteration
  # r (int) = joint rank
  # r.vec (c(ints)) = vector of individual ranks
  # nsample (int) = number of posterior samples
  # thinned_iters_burnin (c(ints)) = vector of indices for samples after burn-in + thinning if desired
  # nninit (boolean) = TRUE if initialized at theoretical posterior mode or FALSE if use posterior mean as pivot
  # pivot (list) = list of V and Vs to use as pivots. If NULL, must provide nninit. 
  # ---------------------------------------------------------------------------
  
  # Load in library for Kabsch and Procrustes algorithms
  library(pracma)
  
  # Saving the number of sources
  q <- length(r.vec)
  
  # Returning the swapped matrices
  swapped_U.draw <- lapply(1:length(U.draw), function(iter) list())
  swapped_V.draw <- lapply(1:length(V.draw), function(iter) list())
  swapped_W.draw <- lapply(1:length(W.draw), function(iter) list())
  swapped_Vs.draw <- lapply(1:length(Vs.draw), function(iter) list())
  swapped_betas <- lapply(1:length(betas), function(iter) list())
  swapped_gammas <- lapply(1:length(gammas), function(iter) list())
  
  # Returning the swaps made so they can be undone
  swaps_made <- lapply(1:length(swapped_U.draw), function(iter) list())
  
  # Storing the sign changes so they can be undone
  signs_changed <- lapply(1:length(swapped_U.draw), function(iter) list())
  
  # Setting the pivots to the posterior mode (result from BIDIFAC) or posterior mean
  if (is.null(pivot)) {
    if (nninit) {
      pivot_V <- V.draw[[1]]
      pivot_Vs <- Vs.draw[[1]]
    }
    
    if (!nninit) {
      V.draw.thinned.burnin <- lapply(V.draw[thinned_iters_burnin], function(iter) iter[[1,1]])
      pivot_V <- matrix(list(), nrow = 1, ncol = 1)
      pivot_V[[1,1]] <- Reduce("+", V.draw.thinned.burnin)/length(V.draw.thinned.burnin)
      
      Vs.draw.thinned.burnin <- Vs.draw[thinned_iters_burnin]
      pivot_Vs <- matrix(list(), nrow = 1, ncol = q)
      
      for (s in 1:q) {
        pivot_Vs[[1,s]] <- Reduce("+", lapply(Vs.draw.thinned.burnin, function(iter) iter[[1,s]]))/length(Vs.draw.thinned.burnin)
      }
    }
  }
  
  # Setting the pivot to be user-chosen
  if (!is.null(pivot)) {
    pivot_V <- pivot[[1]]
    pivot_Vs <- pivot[[2]]
  }
  
  for (iter in 1:length(U.draw)) {
    # -------------------------------------------------------------------------
    # Initializing for iter
    # 
    # For each iteration, arrange the columns of each draw to match that 
    # of the pivot. 
    # -------------------------------------------------------------------------
    
    # Initializing the rearranged version
    tilde_V <- matrix(list(), nrow = 1, ncol = 1)
    tilde_U <- matrix(list(), nrow = q, ncol = 1)
    tilde_W <- matrix(list(), nrow = q, ncol = q)
    tilde_Vs <- matrix(list(), nrow = 1, ncol = q)
    
    tilde_beta <- matrix(list(), nrow = 1, ncol = 1)
    tilde_beta_joint <- matrix(list(), nrow = 1, ncol = 1)
    tilde_beta_indiv <- matrix(list(), nrow = q, ncol = 1)
    
    tilde_gamma <- matrix(list(), nrow = 1, ncol = 1)
    tilde_gamma_joint <- matrix(list(), nrow = 1, ncol = 1)
    tilde_gamma_indiv <- matrix(list(), nrow = q, ncol = 1)
    
    tilde_V[[1,1]] <- matrix(nrow = nrow(V.draw[[1]][[1,1]]), ncol = ncol(V.draw[[1]][[1,1]]))
    
    for (s in 1:q) {
      tilde_U[[s,1]] <- matrix(nrow = nrow(U.draw[[1]][[s,1]]), ncol = ncol(U.draw[[1]][[s,1]]))
      tilde_W[[s,s]] <- matrix(nrow = nrow(W.draw[[1]][[s,s]]), ncol = ncol(W.draw[[1]][[s,s]]))
      tilde_Vs[[1,s]] <- matrix(nrow = nrow(Vs.draw[[1]][[1,s]]), ncol = ncol(Vs.draw[[1]][[1,s]]))
      
      for (ss in 1:q) {
        if (ss != s) {
          tilde_W[[s,ss]] <- matrix(0, nrow = nrow(W.draw[[1]][[s,s]]), ncol = ncol(W.draw[[1]][[ss,ss]]))
        }
      }
    }
    
    n_beta <- 1 + r + sum(r.vec)
    tilde_beta[[1,1]] <- matrix(nrow = n_beta, ncol = 1)
    tilde_beta_joint[[1,1]] <- matrix(nrow = r, ncol = 1)
    
    tilde_gamma[[1,1]] <- matrix(nrow = n_beta, ncol = 1)
    tilde_gamma_joint[[1,1]] <- matrix(nrow = r, ncol = 1)
    
    for (s in 1:q) {
      tilde_beta_indiv[[s,1]] <- matrix(nrow = r.vec[s], ncol = 1)
      tilde_gamma_indiv[[s,1]] <- matrix(nrow = r.vec[s], ncol = 1)
    }
    
    # -------------------------------------------------------------------------
    # Storing the current V, U, beta, and gamma
    # -------------------------------------------------------------------------
    
    current_V <- V.draw[[iter]]
    current_U <- U.draw[[iter]]
    current_Vs <- Vs.draw[[iter]]
    current_W <- W.draw[[iter]]
    
    if (!is.null(betas)) {
      current_beta <- betas[[iter]]
      current_beta_indiv <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        current_beta_joint <- current_beta[[1,1]][2:(r+1),, drop = FALSE] 
      } else {
        current_beta_joint <- NULL
      }
      
      # Individual effects
      current_beta_indiv.temp <- current_beta[[1,1]][(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) {
            current_beta_indiv[[s, 1]] <- current_beta_indiv.temp[1:r.vec[s],, drop = FALSE] 
          }
          if (s != 1) {
            current_beta_indiv[[s, 1]] <- current_beta_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
      }
    }
    
    if (!is.null(gammas)) {
      current_gamma <- gammas[[iter]] 
      current_gamma_indiv <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        current_gamma_joint <- current_gamma[[1,1]][2:(r+1),, drop = FALSE]
      } else {
        current_gamma_joint <- NULL
      }
      
      # Individual effects
      current_gamma_indiv.temp <- current_gamma[[1,1]][(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) {
            current_gamma_indiv[[s,1]] <- current_gamma_indiv.temp[1:r.vec[s],, drop = FALSE] 
          }
          if (s != 1) {
            current_gamma_indiv[[s,1]] <- current_gamma_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Finding the optimal rotation
    # -------------------------------------------------------------------------
    
    # Solve the Procrustes problem for the joint structure
    joint_rotate <- procrustes(pivot_V[[1,1]], current_V[[1,1]])
    
    # Save the rotation matrix
    joint_rotate_mat <- joint_rotate$Q
    
    # Save the rotated V
    current_V_rotate <- current_V
    current_V_rotate[[1,1]] <- current_V[[1,1]] %*% t(joint_rotate_mat)
    
    # Rotate the loadings accordingly
    current_U_rotate <- current_U
    
    for (s in 1:q) {
      current_U_rotate[[s,1]] <- current_U[[s,1]] %*% t(joint_rotate_mat)
    }
    
    # Follow the same procedure for the individual structure
    current_Vs_rotate <- current_Vs
    current_W_rotate <- current_W
    
    for (s in 1:q) {
      # Solve the Procrustes problem for the individual structure
      indiv_rotate <- procrustes(pivot_Vs[[1,s]], current_Vs[[1,s]])
      
      # Save the rotation matrix
      indiv_rotate_mat <- indiv_rotate$Q
      
      # Save the rotated Vs
      current_Vs_rotate[[1,s]] <- current_Vs[[1,s]] %*% t(indiv_rotate_mat)
      
      # Rotate the loadings accordingly
      current_W_rotate[[s,s]] <- current_W[[s,s]] %*% t(indiv_rotate_mat)
    }
    
    # -------------------------------------------------------------------------
    # Finding the optimal order
    # -------------------------------------------------------------------------
    
    # Storing the swaps and signs
    current_swaps_joint <- current_signs_joint <- c()
    current_swaps_indiv <- current_signs_indiv <- lapply(1:q, function(s) c())
    
    # Start with the joint structure
    
    for (k in 1:r) {
      # for each column in the pivot, rearrange the columns in draws to match the order
      corrs <- apply(current_V_rotate[[1,1]], 2, function(current_col) cor(current_col, pivot_V[[1,1]][,k])) # compute correlation between each 
      jk <- which.max(abs(corrs)) # storing the index of the highest correlation
      ak <- sign(corrs[jk]) # storing the sign of that correlation
      
      # Setting the kth column of tilde_V to be the ak*current_V[,jk] column. 
      tilde_V[[1,1]][,k] <- ak*current_V_rotate[[1,1]][,jk]
      
      # Rearranging the corresponding column in U:
      for (s in 1:q) {
        tilde_U[[s,1]][,k] <- ak*current_U_rotate[[s,1]][,jk]
      }
      
      # Rearranging the corresponding rows in betas:
      if (!is.null(betas)) {
        tilde_beta_joint[[1,1]][k,] <- ak*current_beta_joint[jk,]
      }
      
      # Rearranging the corresponding rows in the gammas:
      if (!is.null(gammas)) {
        tilde_gamma_joint[[1,1]][k,] <- current_gamma_joint[jk,]
      }
      
      # Recording what swap was made. The first entry would be the column index
      # of current_V that was moved to the first spot. The second index would be 
      # the column index of current_V that was moved to the second column, etc. 
      current_swaps_joint[k] <- jk
      
      # Recording the sign change (if any) that was made
      current_signs_joint[k] <- ak
      
      # Removing the already-swapped columns in current_V from the running:
      current_V_rotate[[1,1]][,jk] <- NA
    }
    
    # Then the individual structures
    
    for (s in 1:q) {
      for (k in 1:r.vec[s]) {
        # for each column in the pivot, rearrange the columns in draws to match the order
        corrs <- apply(current_Vs_rotate[[1,s]], 2, function(current_col) cor(current_col, pivot_Vs[[1,s]][,k])) # compute correlation between each 
        jk <- which.max(abs(corrs)) # storing the index of the highest correlation
        ak <- sign(corrs[jk]) # storing the sign of that correlation
        
        # Setting the kth column of \tilde_Vs to be the ak*current_Vs[,jk] column. 
        tilde_Vs[[1,s]][,k] <- ak*current_Vs_rotate[[1,s]][,jk]
        
        # Rearranging the corresponding column in W:
        tilde_W[[s,s]][,k] <- ak*current_W_rotate[[s,s]][,jk]
        
        # Rearranging the corresponding rows in betas:
        if (!is.null(betas)) {
          tilde_beta_indiv[[s,1]][k,] <- ak*current_beta_indiv[[s,1]][jk,]
        }
        
        # Rearranging the corresponding rows in the gammas:
        if (!is.null(gammas)) {
          tilde_gamma_indiv[[s,1]][k,] <- current_gamma_indiv[[s,1]][jk,]
        }
        
        # Recording what swap was made. The first entry would be the column index
        # of current_V that was moved to the first spot. The second index would be 
        # the column index of current_V that was moved to the second column, etc. 
        current_swaps_indiv[[s]][k] <- jk
        
        # Recording the sign change (if any) that was made
        current_signs_indiv[[s]][k] <- ak
        
        # Removing the already-swapped columns in current_V from the running:
        current_Vs_rotate[[1,s]][,jk] <- NA
      }
    }
    
    # Combine the betas 
    
    if (!is.null(betas)) {
      tilde_beta[[1,1]] <- rbind(current_beta[[1,1]][1], tilde_beta_joint[[1,1]], do.call(rbind, tilde_beta_indiv))
    }
    
    # Combine the gammas
    if (!is.null(gammas)) {
      tilde_gamma[[1,1]] <- rbind(current_gamma[[1,1]][1], tilde_gamma_joint[[1,1]], do.call(rbind, tilde_gamma_indiv))
    }
    
    # Storing the swapped samples
    swapped_V.draw[[iter]] <- tilde_V
    swapped_U.draw[[iter]] <- tilde_U
    swapped_W.draw[[iter]] <- tilde_W
    swapped_Vs.draw[[iter]] <- tilde_Vs
    
    if (!is.null(betas)) swapped_betas[[iter]] <- tilde_beta
    
    if (!is.null(gammas)) swapped_gammas[[iter]] <- tilde_gamma
    
    # Storing the swaps made
    swaps_made[[iter]] <- list(current_swaps_joint = current_swaps_joint, current_swaps_indiv = current_swaps_indiv)
    
    # Storing the signs changed
    signs_changed[[iter]] <- list(current_signs_joint = current_signs_joint, current_signs_indiv = current_signs_indiv)
  }
  
  # Return
  list(swapped_V.draw = swapped_V.draw, 
       swapped_U.draw = swapped_U.draw, 
       swapped_Vs.draw = swapped_Vs.draw, 
       swapped_W.draw = swapped_W.draw,
       swapped_betas = swapped_betas,
       swapped_gammas = swapped_gammas,
       swaps_made = swaps_made, 
       signs_changed = signs_changed)
}

# Creating yet another rotation-sign-permutation algorithm
factor_switching_permutation_plus_rotation <- function(U.draw, V.draw, W.draw, Vs.draw, betas.draw = NULL, gammas.draw = NULL, r, r.vec, nsample, thinned_iters_burnin, nninit = TRUE, pivot = NULL) {
  
  # ---------------------------------------------------------------------------
  # Arguments: 
  # draws (list): the draws to swap according to the pivot
  # rank (int): number of columns to iterate through
  # loadings (list): the corresponding loadings for each draw of the matrix to swap
  # betas (double vec): list of vectors of the coefficients drawn at each iteration
  # gammas (int vec): list of vectors of the inclusion indicators drawn at each iteration
  # r (int) = joint rank
  # r.vec (c(ints)) = vector of individual ranks
  # nsample (int) = number of posterior samples
  # thinned_iters_burnin (c(ints)) = vector of indices for samples after burn-in + thinning if desired
  # nninit (boolean) = TRUE if initialized at theoretical posterior mode or FALSE if use posterior mean as pivot
  # pivot (list) = list of V and Vs to use as pivots. If NULL, must provide nninit. 
  # ---------------------------------------------------------------------------
  
  # Load in library for Procrustes algorithms and permutations
  library(pracma)
  library(combinat)
  library(svMisc)
  
  # Saving the number of sources
  q <- length(r.vec)
  
  # Returning the swapped matrices
  swapped_U.draw <- lapply(1:length(U.draw), function(iter) list())
  swapped_V.draw <- lapply(1:length(V.draw), function(iter) list())
  swapped_W.draw <- lapply(1:length(W.draw), function(iter) list())
  swapped_Vs.draw <- lapply(1:length(Vs.draw), function(iter) list())
  swapped_betas <- lapply(1:length(betas), function(iter) list())
  swapped_gammas <- lapply(1:length(gammas), function(iter) list())
  
  # Returning the swaps made so they can be undone
  swaps_made <- lapply(1:length(swapped_U.draw), function(iter) list())
  
  # Storing the sign changes so they can be undone
  signs_changed <- lapply(1:length(swapped_U.draw), function(iter) list())
  
  # Setting the pivots to the posterior mode (result from BIDIFAC) or posterior mean
  if (is.null(pivot)) {
    if (nninit) {
      pivot_V <- V.draw[[1]]
      pivot_Vs <- Vs.draw[[1]]
    }
    
    if (!nninit) {
      V.draw.thinned.burnin <- lapply(V.draw[thinned_iters_burnin], function(iter) iter[[1,1]])
      pivot_V <- matrix(list(), nrow = 1, ncol = 1)
      pivot_V[[1,1]] <- Reduce("+", V.draw.thinned.burnin)/length(V.draw.thinned.burnin)
      
      Vs.draw.thinned.burnin <- Vs.draw[thinned_iters_burnin]
      pivot_Vs <- matrix(list(), nrow = 1, ncol = q)
      
      for (s in 1:q) {
        pivot_Vs[[1,s]] <- Reduce("+", lapply(Vs.draw.thinned.burnin, function(iter) iter[[1,s]]))/length(Vs.draw.thinned.burnin)
      }
    }
  }
  
  # Setting the pivot to be user-chosen
  if (!is.null(pivot)) {
    pivot_V <- pivot[[1]]
    pivot_Vs <- pivot[[2]]
  }
  
  # Creating a matrix with all possible column permutations
  perm_V <- do.call(rbind, permn(r))
  perm_Vs <- lapply(1:q, function(s) do.call(rbind, permn(r.vec[s])))
  
  # How many permutations are there?
  n_perm_V <- factorial(r)
  n_perm_Vs <- factorial(r.vec)
  
  for (iter in 1:length(U.draw)) {
    progress(iter/(length(U.draw)/100))
    
    # -------------------------------------------------------------------------
    # Storing the current V, U, beta, and gamma
    # -------------------------------------------------------------------------
    
    current_V <- V.draw[[iter]]
    current_U <- U.draw[[iter]]
    current_Vs <- Vs.draw[[iter]]
    current_W <- W.draw[[iter]]
    
    if (!is.null(betas.draw)) {
      current_beta <- betas.draw[[iter]]
      current_beta_indiv <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        current_beta_joint <- current_beta[[1,1]][2:(r+1),, drop = FALSE] 
      } else {
        current_beta_joint <- NULL
      }
      
      # Individual effects
      current_beta_indiv.temp <- current_beta[[1,1]][(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) {
            current_beta_indiv[[s, 1]] <- current_beta_indiv.temp[1:r.vec[s],, drop = FALSE] 
          }
          if (s != 1) {
            current_beta_indiv[[s, 1]] <- current_beta_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
      }
    }
    
    if (!is.null(gammas.draw)) {
      current_gamma <- gammas.draw[[iter]] 
      current_gamma_indiv <- matrix(list(), nrow = q, ncol = 1)
      
      # Joint effect
      if (r != 0) {
        current_gamma_joint <- current_gamma[[1,1]][2:(r+1),, drop = FALSE]
      } else {
        current_gamma_joint <- NULL
      }
      
      # Individual effects
      current_gamma_indiv.temp <- current_gamma[[1,1]][(r+2):n_beta,, drop = FALSE]
      
      for (s in 1:q) {
        # If there is an individual effect
        if (r.vec[s] != 0) {
          if (s == 1) {
            current_gamma_indiv[[s,1]] <- current_gamma_indiv.temp[1:r.vec[s],, drop = FALSE] 
          }
          if (s != 1) {
            current_gamma_indiv[[s,1]] <- current_gamma_indiv.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
          }
        }
      }
    }
    
    # -------------------------------------------------------------------------
    # Finding the optimal rotation for each permutation
    # -------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------
    # Joint structure
    # -------------------------------------------------------------------------
    
    costs_V <- c()
    
    for (ind in 1:n_perm_V) {
      
      # Attempt current permutation
      current_perm <- perm_V[ind,]
      current_V_perm <- current_V 
      current_V_perm[[1,1]] <- current_V[[1,1]][,current_perm]
      
      # For this permutation, what is the optimal rotation?
      
      # Solve the Procrustes problem for the joint structure
      joint_rotate <- procrustes(pivot_V[[1,1]], current_V_perm[[1,1]])
      
      # Save the rotation matrix
      joint_rotate_mat <- joint_rotate$Q
      
      # Save the rotated V
      current_V_perm_rotate <- current_V_perm
      current_V_perm_rotate[[1,1]] <- current_V_perm[[1,1]] %*% t(joint_rotate_mat)
      
      # Calculate the Frob difference between pivot and current perm/rotate
      costs_V[ind] <- frob(pivot_V[[1,1]] - current_V_perm_rotate[[1,1]])
    }
    
    # For the smallest cost, permute and rotate V and U accordingly
    min_ind <- which.min(costs_V)
    best_perm <- perm_V[min_ind,]
    current_V_best_perm_rotate <- current_V
    
    # Solve the Procrustes problem for the joint structure
    joint_rotate_best <- procrustes(pivot_V[[1,1]], current_V_best_perm_rotate[[1,1]])
    
    # Save the rotation matrix
    joint_rotate_mat_best <- joint_rotate_best$Q
    
    # Save the rotated V
    current_V_best_perm_rotate[[1,1]] <- current_V[[1,1]][,best_perm]
    current_V_best_perm_rotate[[1,1]] <- current_V_best_perm_rotate[[1,1]] %*% t(joint_rotate_mat)
    
    # Rotate the loadings accordingly
    current_U_best_perm_rotate <- current_U
    
    for (s in 1:q) {
      current_U_best_perm_rotate[[s,1]] <- current_U_best_perm_rotate[[s,1]][,best_perm]
      current_U_best_perm_rotate[[s,1]] <- current_U_best_perm_rotate[[s,1]] %*% t(joint_rotate_mat)
    }
    
    # Rearranging the corresponding rows in betas:
    if (!is.null(betas.draw)) {
      best_beta_joint <- current_beta_joint
      best_beta_joint <- current_beta_joint[best_perm,,drop=FALSE]
    }
    
    # Rearranging the corresponding rows in the gammas:
    if (!is.null(gammas.draw)) {
      best_gamma_joint <- current_gamma_joint
      best_gamma_joint <- current_gamma_joint[best_perm,,drop=FALSE]
    }
    
    # -------------------------------------------------------------------------
    # Individual structure
    # ------------------------------------------------------------------------
    
    # Follow the same procedure for the individual structure
    costs_Vs <- lapply(1:q, function(s) c())
    
    for (s in 1:q) {
      for (ind in 1:n_perm_Vs[s]) {
      
        # Attempt current permutation
        current_perm <- perm_Vs[[s]][ind,]
        current_Vs_perm <- current_Vs
        current_Vs_perm[[1,s]] <- current_Vs[[1,s]][,current_perm]
        
        # For this permutation, what is the optimal rotation?
        
        # Solve the Procrustes problem for the joint structure
        indiv_rotate <- procrustes(pivot_Vs[[1,s]], current_Vs_perm[[1,s]])
        
        # Save the rotation matrix
        indiv_rotate_mat <- indiv_rotate$Q
        
        # Save the rotated V
        current_Vs_perm_rotate <- current_Vs_perm
        current_Vs_perm_rotate[[1,s]] <- current_Vs_perm[[1,1]] %*% t(indiv_rotate_mat)
        
        # Calculate the Frob difference between pivot and current perm/rotate
        costs_Vs[[s]][ind] <- frob(pivot_Vs[[1,s]] - current_Vs_perm_rotate[[1,s]])
      }
    }
    
    # For the smallest cost, permute and rotate Vs and W accordingly
    min_ind_Vs <- sapply(costs_Vs, function(costs) which.min(costs))
    best_perm_Vs <- lapply(1:q, function(s) perm_Vs[[s]][min_ind_Vs[s],])
    
    # Initialize the permutated and rotated results
    current_Vs_best_perm_rotate <- current_Vs
    current_W_best_perm_rotate <- current_W
    best_beta_indiv <- current_beta_indiv
    best_gamma_indiv <- current_gamma_indiv
    
    for (s in 1:q) {
      # Apply the best permutation
      current_Vs_best_perm_rotate[[1,s]] <- current_Vs[[1,s]][,best_perm_Vs[[s]]]
      
      # Obtain the corresponding rotation
      indiv_rotate <- procrustes(pivot_Vs[[1,s]], current_Vs_best_perm_rotate[[1,s]])
      
      # Save the rotation matrix
      indiv_rotate_mat <- indiv_rotate$Q

      # Rotate the scores
      current_Vs_best_perm_rotate[[1,s]] <- current_Vs_best_perm_rotate[[1,s]] %*% t(indiv_rotate_mat)
      
      # Permute the loadings
      current_W_best_perm_rotate[[s,s]] <- current_W[[s,s]][,best_perm_Vs[[s]]]
      
      # Rotate the loadings
      current_W_best_perm_rotate[[s,s]] <- current_W_best_perm_rotate[[s,s]] %*% t(indiv_rotate_mat)
      
      # Rearranging the corresponding rows in betas:
      if (!is.null(betas.draw)) {
        best_beta_indiv[[s,1]] <- current_beta_indiv[[s,1]][best_perm_Vs[[s]],,drop=FALSE]
      }
    
      # Rearranging the corresponding rows in the gammas:
      if (!is.null(gammas.draw)) {
        best_gamma_indiv[[s,1]] <- current_gamma_indiv[[s,1]][best_perm_Vs[[s]],,drop=FALSE]
      }
    }
    
    # -------------------------------------------------------------------------
    # Combine and save the results
    # -------------------------------------------------------------------------
    
    # Combine the betas 
    if (!is.null(betas.draw)) {
      current_best_beta <- current_beta
      current_best_beta[[1,1]] <- rbind(current_beta[[1,1]][1], best_beta_joint, do.call(rbind, best_beta_indiv))
    }
    
    # Combine the gammas
    if (!is.null(gammas.draw)) {
      current_best_gamma <- current_gamma
      current_best_gamma[[1,1]] <- rbind(current_gamma[[1,1]][1], best_gamma_joint, do.call(rbind, best_gamma_indiv))
    }
    
    # Save in the lists
    swapped_V.draw[[iter]] <- current_V_best_perm_rotate
    swapped_U.draw[[iter]] <- current_U_best_perm_rotate
    swapped_Vs.draw[[iter]] <- current_Vs_best_perm_rotate
    swapped_W.draw[[iter]] <- current_W_best_perm_rotate
    
    if (!is.null(betas.draw)) {
      swapped_betas[[iter]] <- current_best_beta
    }
    
    if (!is.null(gammas.draw)) {
      swapped_gammas[[iter]] <- current_best_gamma
    }
  }
  
  # Return
  list(swapped_V.draw = swapped_V.draw, 
       swapped_U.draw = swapped_U.draw, 
       swapped_Vs.draw = swapped_Vs.draw, 
       swapped_W.draw = swapped_W.draw,
       swapped_betas = swapped_betas,
       swapped_gammas = swapped_gammas)
}

# Edit the jointRot function from the infinite factor package to use the Varimax-rotated BIDIFAC solutions 
jointRot_multi <- function(lambda, eta, var_betas = NULL) {
  
  # ---------------------------------------------------------------------------
  # Arguments:
  #
  # lambda (list): list of loadings matrices. For BPMF, this will be combined
  #                loadings and betas. 
  # eta (list): list of scores matrices
  # var_betas (matrix): prior variance matrix for regression coefficients, if
  #                     considering an outcome
  # 
  # Acknowlegements: code adapted from infinitefactor package by Evan Poworoznek (2021)
  # ---------------------------------------------------------------------------
  
  # Apply the Varimax rotation to the loadings
  vari = lapply(lambda, varimax)
  
  # Save the resulting Varimax-rotated loadings
  loads = lapply(vari, `[[`, 1)
  
  # Save the resulting rotation matrix
  rots = lapply(vari, `[[`, 2)
  
  # Apply the rotation matrix to the corresponding scores
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)
  
  # Calculate the pivot to match to for the data
  pivot.data <- loads[[1]][1:(nrow(loads[[1]])-1),] # The rotated BIDIFAC solution
  
  # Create a pivot for the betas based on the Varimax-BIDIFAC solution and posterior mean for betas
  if (!is.null(var_betas)) {
    X <- rotfact[[1]] # Saving the Varimax-rotated scores at the BIDIFAC solutio
    pivot.betas <- solve(t(X) %*% X + solve(var_betas)) %*% (t(X) %&% y[[1,1]])
    
    # Combine the pivots together to create a fully joint pivot
    piv <- rbind(pivot.data, t(as.matrix(pivot.betas)))
  }
  
  if (is.null(var_betas)) {
    # Set the pivot to the data pivot
    piv <- pivot.data
  }

  # Match to the defined pivot
  matches = lapply(loads, msfOUT, piv)
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)
  
  # Return
  return(list(lambda = lamout, eta = etaout))
}

# Modifying the factor switching method from Poworoznek et al. (2021)
match_align_bpmf <- function(BPMF.fit, y, model_params, p.vec) {
  
  # ---------------------------------------------------------------------------
  # Arguments:
  # 
  # BPMF.fit (list): a list of results from the BPMF model fit 
  # y (matrix of lists): outcome vector
  # model_params (list): model parameters used in model fitting
  # p.vec (vector): vector with number of variables in each source
  # ---------------------------------------------------------------------------
  
  # Loading in the library for factor switching
  library(infinitefactor)
  
  # Save the ranks from the model fit
  ranks <- BPMF.fit$ranks
  
  # Save the number of datasets
  q <- nrow(BPMF.fit$data)
  
  # Save the number of posterior samples
  nsample <- length(BPMF.fit$V.draw)
  
  # Create a vector for the indices of each rank
  beta.ind <- lapply(1:length(ranks), function(i) {
    if (i == 1) { 
      (1:ranks[i]) + 1
    } else {
      ((sum(ranks[1:(i-1)])+1):sum(ranks[1:i])) + 1
    }
  })
  
  # Joint structure (plus joint regression coefficients) --
  
  # Save the overall number of features
  p <- sum(p.vec)
  
  # Combine joint loadings and betas
  joint.loadings <- lapply(1:nsample, function(iter) {
    rbind(do.call(rbind, BPMF.fit$U.draw[[iter]]), # Joint loadings  
          t(BPMF.fit$beta.draw[[iter]][[1,1]][beta.ind[[1]],])) # Joint regression coefficients
  })
  
  # Save joint scores
  joint.scores <- lapply(1:nsample, function(iter) BPMF.fit$V.draw[[iter]][[1,1]])
  
  # Save the prior variance on the betas for joint factors
  joint_var_betas <- diag(rep(model_params$beta_vars[2], ranks[1]))
  
  # Apply the factor switching method to the joint structure
  joint.results.rotate <- jointRot_multi(joint.loadings, joint.scores, var_betas = joint_var_betas)
  
  # Separate the joint loadings from the joint scores
  joint.loadings.final <- lapply(joint.results.rotate$lambda, function(iter) iter[1:p,])
  joint.betas.final <- lapply(joint.results.rotate$lambda, function(iter) t(iter[p+1,,drop=FALSE]))
  
  # Applying the permutation and sign switching to the scores
  joint.scores.final <- joint.results.rotate$eta
  
  # Individual structure (plus individual regression coefficients) --
  
  # Combine the individual loadings and betas
  individual.loadings <- lapply(1:q, function(s) lapply(1:nsample, function(iter) {
    rbind(BPMF.fit$W.draw[[iter]][[s,s]], 
          BPMF.fit$beta.draw[[iter]][[1,1]][beta.ind[[s+1]],])
  }))
  individual.scores <- lapply(1:q, function(s) lapply(1:nsample, function(iter) BPMF.fit$Vs.draw[[iter]][[1,s]]))
  
  # Saving the Varimax-rotated BIDIFAC solution for the data which serves as the pivot
  individual.pivot.data <- lapply(1:q, function(s) individual.loadings[[s]]$lambda[[1]][1:p.vec[s],])
  
  # Save the prior variance on the betas for the individual factors
  indiv_var_betas <- lapply(1:q, function(s) diag(rep(model_params$beta_vars[s+2], ranks[s+1])))
  
  # Trying it a different way
  individual.results.rotate <- lapply(1:q, function(s) {
    jointRot_multi(individual.loadings[[s]], individual.scores[[s]], var_betas = indiv_var_betas[[s]])
  })
  
  # Save the final individual loadings and betas
  individual.loadings.final <- lapply(1:q, function(s) lapply(individual.results.rotate[[s]]$lambda, function(iter) iter[1:p.vec[s],,drop=FALSE]))
  individual.betas.final <- lapply(1:q, function(s) lapply(individual.results.rotate[[s]]$lambda, function(iter) t(iter[p.vec[s]+1,,drop=FALSE])))
  
  # Applying the permutation and sign switching to the scores
  individual.scores.final <- lapply(1:q, function(s) individual.results.rotate[[s]]$eta)
  
  # Return the final loadings, scores, and betas
  list(joint.scores.final = joint.scores.final, 
       joint.loadings.final = joint.loadings.final, 
       joint.betas.final = joint.betas.final, 
       individual.scores.final = individual.scores.final, 
       individual.loadings.final = individual.loadings.final, 
       individual.betas.final = individual.betas.final)
  
  # Check the individual structure after the algorithm matches the original structure
  
  # Check the joint structure after the algorithm matches the original structure
  # joint.structure.original <- lapply(1:nsample, function(iter) {
  #   joint.loadings[[iter]] %*% t(joint.scores[[iter]])
  # })
  # 
  # joint.structure.new <- lapply(1:nsample, function(iter) {
  #   joint.results.rotate$lambda[[iter]] %*% t(joint.results.rotate$eta[[iter]])
  # })
  # 
  # # Check
  # all.equal(joint.structure.original, joint.structure.new)
  
  # Calculate the original structure
  # individual.structure.original <- lapply(1:q, function(s) {
  #   lapply(1:nsample, function(iter) {
  #     individual.loadings[[s]][[iter]] %*% t(individual.scores[[s]][[iter]])
  #   })
  # })
  # 
  # # Calculate the new structure
  # individual.structure.new <- lapply(1:q, function(s) {
  #   lapply(1:nsample, function(iter) {
  #     individual.results.rotate[[s]]$lambda[[iter]] %*% t(individual.results.rotate[[s]]$eta[[iter]])
  #   })
  # })
  # 
  # # Check
  # lapply(1:q, function(s) all.equal(individual.structure.original[[s]], individual.structure.new[[s]]))
  
}