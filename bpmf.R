# -----------------------------------------------------------------------------
# Helper functions for Bayesian PMF
# -----------------------------------------------------------------------------

# Packages
library(doParallel)
library(foreach)
library(Matrix)
library(MASS)
library(truncnorm)
# library(MCMCpack)

# -----------------------------------------------------------------------------
# Bayesian PMF functions
# -----------------------------------------------------------------------------

bpmf <- function(data, Y, nninit = TRUE, model_params, ranks = NULL, scores = NULL, sparsity = FALSE, nsample, progress = TRUE) {
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
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (!missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  if (!sparsity) {
    gamma.draw <- p.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (response_given) {
    beta.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    Z.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    tau2.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    Ym.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
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
        Vs0[[1,s]] <- (svd.indiv.s$v[,1:r.vec[s], drop = FALSE]) %*% diag(svd.joint$d[1:r.vec[s]], nrow = r.vec[s])
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
  
  # If ranks provided, initialize with prior
  if (!nninit) {
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
  
  # Storing the joint structure at each Gibbs sampling iteration
  Joint.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  # Storing the individual structure at each Gibbs sampling iteration
  Indiv.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  if (is.null(scores)) {
    for (iter in 1:nsample) {
      for (s in 1:q) {
        # Calculating the joint structure and scaling by sigma.mat
        Joint.draw[[iter]][[s,1]] <- (U.draw[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]])) * sigma.mat[s,1]
        
        # Calculating the individual structure and scaling by sigma.mat
        Indiv.draw[[iter]][[s,1]] <- (W.draw[[iter]][[s,1]] %*% t(Vs.draw[[iter]][[1,1]])) * sigma.mat[s,1]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the posterior mean
  # ---------------------------------------------------------------------------
  
  # Storing the posterior mean for the joint and individual structure
  Joint.mean <- matrix(list(), nrow = q, ncol = 1)
  Indiv.mean <- matrix(list(), nrow = q, ncol = 1)
  
  if (is.null(scores)) {
    for (s in 1:q) {
      Joint.mean[[s,1]] <- Reduce("+", lapply(1:nsample, function(iter) Joint.draw[[iter]][[s,1]]))/length(Joint.draw)
      Indiv.mean[[s,1]] <- Reduce("+", lapply(1:nsample, function(iter) Indiv.draw[[iter]][[s,1]]))/length(Indiv.draw)
    }
  }
  
  # Return
  list(data = data, # Returning the scaled version of the data
        Y = Y, # Return the response vector
        sigma.mat = sigma.mat, # Scaling factors
        Joint.draw = Joint.draw, Indiv.draw = Indiv.draw, # Underlying structure
        Joint.mean = Joint.mean, Indiv.mean = Indiv.mean, # Posterior mean of structures
        V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
        Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
        scores = scores, # Scores if provided by another method 
        ranks = c(r, r.vec), # Ranks
        tau2.draw = tau2.draw, beta.draw = beta.draw, # Regression parameters
        gamma.draw = gamma.draw, p.draw = p.draw) # Sparsity parameters

}

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
# Helper functions for initializing with BIDIFAC
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

# -----------------------------------------------------------------------------
# Helper functions for validation simulations
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
bpmf_data <- function(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response, missingness, entrywise, prop_missing, sparsity) {
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
  response_vars <- true_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2]
  
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
    data[s,1][[1]] <- joint.structure[[s,1]] + indiv.structure[[s,1]] + E[[s,1]]
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
      gamma <- matrix(rbinom(n_beta, size = 1, prob = p.prior), ncol = 1)
      gamma[1,] <- 1 # Always include the intercept
      diag(Sigma_beta)[gamma == 0] <- 1/1000
      beta <- matrix(list(), nrow = 1, ncol = 1)
      beta[[1,1]] <- matrix(mvrnorm(1, mu = rep(0, n_beta), Sigma = Sigma_beta), ncol = 1)
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
      tau2 <- matrix(list(), nrow = 1, ncol = 1)
      tau2[[1,1]] <- matrix(1/rgamma(1, shape = shape, rate = rate)) 
      error_y <- matrix(rnorm(n, mean = 0, sd = sqrt(tau2[[1,1]])), ncol = 1)
      
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
       beta = beta, tau2 = tau2, EY = EY, gamma = gamma, p.prior = p.prior)
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
              results_for_param[[s,1]] <- list(observed = list(avg_coverage = mean(avg_coverage_observed_source),
                                                          avg_mse = mean(avg_mse_observed_source),
                                                          avg_ci_width = mean(avg_ci_width_observed_source)),
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
# Helper functions for model simulations
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
run_each_mod <- function(mod, p.vec, n, ranks, response, true_params, model_params, s2nX, s2nY, sparsity, nsim, nsample = 2000, n_clust) {
  
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
  # ---------------------------------------------------------------------------
  
  # Loading in the packages
  library(r.jive)
  library(MOFA2)
  library(doParallel)
  library(foreach)
  
  # The model options
  models <- c("sJIVE", "BIDIFAC+", "JIVE", "MOFA", "BPMF")
  
  start <- Sys.time()
  cl <- makeCluster(n_clust)
  registerDoParallel(cl)
  funcs <- c("bpmf_data", "center_data", "bpmf", "get_results", "BIDIFAC",
             "check_coverage", "mse", "ci_width", "data.rearrange", "return_missing",
             "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum",
             "sJIVE", "sJIVE.converge", "sJIVE.predict", "sJIVE.ranks")
  packs <- c("Matrix", "MASS", "truncnorm", "r.jive")
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
    q <- length(p.vec)
    joint.structure <- sim_data$joint.structure
    indiv.structure <- sim_data$indiv.structure
    
    # The response
    Y <- sim_data$Y
    
    # The standardizing coefficients
    s2nX_coef <- sim_data$s2nX_coef
    s2nY_coef <- sim_data$s2nY_coef
    
    # The response parameters
    beta <- sim_data$beta
    tau2 <- sim_data$tau2
    EY <- sim_data$EY
    
    # Save a burn-in
    burnin <- nsample/2
    
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
    # Center and scale X and y to have variance 1
    # -------------------------------------------------------------------------
    
    # for (s in 1:q) {
    #   # Center and scale the observed training data (transpose to scale the cols (features) than transpose back)
    #   training_data[[s,1]] <- t(scale(t(training_data[[s,1]]), center = TRUE, scale = TRUE))
    #   training_data_list[[s]] <- t(scale(t(training_data_list[[s]]), center = TRUE, scale = TRUE))
    #   
    #   # Save the means and variances
    #   row_means_train <- attr(training_data[[s,1]], "scaled:center")
    #   row_sd_train <- attr(training_data[[s,1]], "scaled:scale")
    #   
    #   # Subtract each column by the same mean and divide by the same sd above
    #   joint.structure_train[[s,1]] <- sweep(joint.structure_train[[s,1]], 1, row_means_train)
    #   joint.structure_train[[s,1]] <- sweep(joint.structure_train[[s,1]], 1, row_sd_train, FUN = "/")
    #   
    #   indiv.structure_train[[s,1]] <- sweep(indiv.structure_train[[s,1]], 1, row_means_train)
    #   indiv.structure_train[[s,1]] <- sweep(indiv.structure_train[[s,1]], 1, row_sd_train, FUN = "/")
    #   
    #   # Center and scale the observed test data
    #   test_data[[s,1]] <- t(scale(t(test_data[[s,1]]), center = TRUE, scale = TRUE))
    #   test_data_list[[s]] <- t(scale(t(test_data_list[[s]]), center = TRUE, scale = TRUE))
    #   
    #   # Save the means and variances
    #   row_means_test <- attr(test_data[[s,1]], "scaled:center")
    #   row_sd_test <- attr(test_data[[s,1]], "scaled:scale")
    #   
    #   # Subtract each column by the same mean and divide by the same sd above
    #   joint.structure_test[[s,1]] <- sweep(joint.structure_test[[s,1]], 1, row_means_test)
    #   joint.structure_test[[s,1]] <- sweep(joint.structure_test[[s,1]], 1, row_sd_test, FUN = "/")
    #   
    #   indiv.structure_test[[s,1]] <- sweep(indiv.structure_test[[s,1]], 1, row_means_test)
    #   indiv.structure_test[[s,1]] <- sweep(indiv.structure_test[[s,1]], 1, row_sd_test, FUN = "/")
    # }
    # 
    # Y_train[[1,1]] <- scale(Y_train[[1,1]], center = TRUE, scale = TRUE)
    # EY_train[[1,1]] <- (EY_train[[1,1]] - attr(Y_train[[1,1]], "scaled:center"))/attr(Y_train[[1,1]], "scaled:scale")
    # 
    # Y_test[[1,1]] <- scale(Y_test[[1,1]], center = TRUE, scale = TRUE)
    # EY_test[[1,1]] <- (EY_test[[1,1]] - attr(Y_test[[1,1]], "scaled:center"))/attr(Y_test[[1,1]], "scaled:scale")
    
    # -------------------------------------------------------------------------
    # Fit each model on generated data to obtain estimate of underlying structure
    # -------------------------------------------------------------------------
    
    if (mod == "sJIVE") {
      # Running sJIVE
      mod.out <- sJIVE(training_data_list, Y_train[[1,1]], method = "CV")
      
      # Saving the joint structure
      mod.joint <- lapply(1:q, function(source) {
        t(t(mod.out$U_I[[source]])) %*% mod.out$S_J
      })
      
      # Saving the individual structure
      mod.individual <- lapply(1:q, function(source) {
        mod.out$W_I[[source]] %*% mod.out$S_I[[source]]
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
      indiv.scores <- lapply(1:q, function(source) {
        if (indiv.rank[source] != 0) {
          t(mod.out$S_I[[source]])
        }
      })
      
      # Saving the E(Y)
      Y.fit <- t(mod.out$fittedY)
      
      # Do not compute results for coverage
      coverage_EY_train <- coverage_EY_test <- NA
    }
    
    if (mod == "BIDIFAC+") {
      # Run model
      mod.out <- BIDIFAC(true_data, rmt = TRUE, pbar = FALSE, scale_back = TRUE)
      
      # Saving the column structure (the joint structure)
      mod.joint <- lapply(1:q, function(source) {
        mod.out$C[[source,1]]
      })
      
      # Saving the individual structure 
      mod.individual <- lapply(1:q, function(source) {
        mod.out$I[[source,1]]
      })
      
      # Saving the joint rank
      joint.rank <- rankMatrix(mod.out$C[[1,1]])[1]
      
      # Saving the individual ranks
      indiv.rank <- sapply(mod.out$I, function(source) {
        rankMatrix(source)[1]
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
      indiv.scores <- lapply(1:q, function(source) {
        if (indiv.rank[source] != 0) {
          svd.source <- svd(mod.individual[[source]])
          (svd.source$v[,1:indiv.rank[source]]) %*% diag(svd.source$d[1:indiv.rank[source]], nrow = indiv.rank[source])
        }
      })
    }
    
    if (mod == "JIVE") {
      # Running JIVE
      mod.out <- jive(true_data_list, scale = FALSE)
      
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
        joint.scores <- svd(mod.joint[[1]])$v[,1:joint.rank]
      }
      if (joint.rank == 0) {
        joint.scores <- NULL
      }
      
      # Obtaining the individual scores
      indiv.scores <- lapply(1:q, function(source) {
        if (indiv.rank[source] != 0) {
          svd.source <- svd(mod.individual[[source]])
          svd.source$v[,1:indiv.rank[source]]
        }
      })
    }
    
    if (mod == "MOFA") {
      # Create the MOFA object
      MOFAobject <- prepare_mofa(
        object = create_mofa(true_data_list)
      )
      
      # Train the MOFA model
      mod.out <- run_mofa(MOFAobject)
      
      # Determining for which views each factor is active
      mod.var.exp <- get_variance_explained(mod.out)$r2_per_factor$group1 # variance explained by factor per view
      joint.factors <- which(apply(mod.var.exp, 1, function(factor) all(factor > 0.002)))
      indiv.factors <- lapply(1:q, function(source) {
        which(apply(mod.var.exp, 1, function(factor) factor[source] > 0.002 & factor[c(1:q)[!c(1:q) %in% source]] < 0.002))
      })
      
      # Saving the joint rank
      joint.rank <- length(joint.factors)
      
      # Saving the individual rank
      indiv.rank <- sapply(indiv.factors, length)
      
      # Save the underlying structure
      mod.joint <- lapply(1:q, function(source) {
        t(mod.out@expectations$Z$group1[, joint.factors, drop = FALSE] %*% t(mod.out@expectations$W[[source]][, joint.factors, drop = FALSE]))
      })
      mod.individual <- lapply(1:q, function(source) {
        t(mod.out@expectations$Z$group1[, indiv.factors[[source]], drop = FALSE] %*% t(mod.out@expectations$W[[source]][, indiv.factors[[source]], drop = FALSE]))
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
      indiv.scores <- lapply(1:q, function(source) {
        if (indiv.rank[source] != 0) {
          indiv.scores <- mofa.scores[,unlist(indiv.factors[[source]]), drop = FALSE]
        }
      })
    }
    
    if (mod == "BPMF") {
      # Setting the test response to NA
      Y_NA_for_test <- Y
      Y_NA_for_test[[1,1]][(n+1):(2*n),] <- NA
      
      # Running BPMF
      mod.out <- bpmf(true_data, Y = Y_NA_for_test, nninit = TRUE, model_params = model_params, nsample = nsample)
      
      # Saving the joint structure
      mod.joint.iter <- lapply(1:nsample, function(iter) {
        lapply(1:q, function(source) {
          (mod.out$U.draw[[iter]][[source,1]] %*% t(mod.out$V.draw[[iter]][[1]])) * mod.out$sigma.mat[source,]
        })
      })
      
      # Taking the posterior mean
      mod.joint <- lapply(1:q, function(source) {
        # Save the joint structure at each iteration for source
        joint.source <- lapply((burnin+1):nsample, function(iter) {
          mod.joint.iter[[iter]][[source]]
        })
        # Take the mean
        Reduce("+", joint.source)/length(joint.source)
      })
      
      # Saving the individual structure
      mod.individual.iter <- lapply(1:nsample, function(iter) {
        lapply(1:q, function(source) {
          (mod.out$W.draw[[iter]][[source,source]] %*% t(mod.out$Vs.draw[[iter]][[1,source]])) * mod.out$sigma.mat[source,]
        })
      })
      
      # Taking the posterior mean
      mod.individual <- lapply(1:q, function(source) {
        # Save the joint structure at each iteration for source
        indiv.source <- lapply((burnin+1):nsample, function(iter) {
          mod.individual.iter[[iter]][[source]]
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
    }
    
    if (mod != "test") {
      # Combining the ranks
      mod.ranks <- c(joint.rank, indiv.rank)
      
      # Combining all scores together
      all.scores <- cbind(Y[[1,1]], joint.scores, do.call(cbind, indiv.scores))
      colnames(all.scores) <- c("y", rep("joint", joint.rank), rep("indiv", sum(indiv.rank)))
    }
    
    # -------------------------------------------------------------------------
    # As applicable, use structure in Bayesian linear model
    # -------------------------------------------------------------------------
    
    if (mod %in% c("BIDIFAC+", "JIVE", "MOFA")) {
      # Subset the scores to just the training data
      all.scores.train <- all.scores[1:n,,drop=FALSE]
      
      # Fitting the Bayesian linear model
      mod.bayes <- bpmf(data = training_data, Y = Y_train, nninit = FALSE, model_params = model_params, 
                        ranks = mod.ranks, scores = all.scores.train[,-1], nsample = nsample)
      
      # Calculate the predicted E(Y) at each Gibbs sampling iteration using training and testing scores
      Y.fit.iter <- lapply((burnin+1):nsample, function(iter) {
        VStar.iter <- cbind(1, all.scores[,-1])
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
    
    if (mod == "BPMF") {
      # Calculate the predicted E(Y) at each Gibbs sampling iteration
      Y.fit.iter <- lapply((burnin+1):nsample, function(iter) {
        VStar.iter <- cbind(1, all.scores[,-1])
        beta.iter <- mod.out$beta.draw[[iter]][[1,1]]
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
      joint.recovery.structure.train <- mean(sapply(1:q, function(source) {
        frob(mod.joint[[source]][,1:n] - joint.structure_train[[source,1]])/frob(joint.structure_train[[source,1]])
      }))
      
      joint.recovery.structure.test <- mean(sapply(1:q, function(source) {
        frob(mod.joint[[source]][,(n+1):(2*n)] - joint.structure_test[[source,1]])/frob(joint.structure_test[[source,1]])
      }))
      
      # Individual structure
      indiv.recovery.structure.train <- mean(sapply(1:q, function(source) {
        frob(mod.individual[[source]][,1:n] - indiv.structure_train[[source,1]])/frob(indiv.structure_train[[source,1]])
      }))
      
      indiv.recovery.structure.test <- mean(sapply(1:q, function(source) {
        frob(mod.individual[[source]][,(n+1):(2*n)] - indiv.structure_test[[source,1]])/frob(indiv.structure_test[[source,1]])
      }))
    }
    
    # -------------------------------------------------------------------------
    # Calculate prediction error for test data
    # -------------------------------------------------------------------------
    
    # Comparing the predicted Y to the training and test Y
    mse_y_train <- frob(Y.fit[1:n,] - Y_train[[1,1]])/frob(Y_train[[1,1]])
    mse_y_test <- frob(Y.fit[(n+1):(2*n),] - Y_test[[1,1]])/frob(Y_test[[1,1]])
    
    # Save 
    save(joint.recovery.structure.train, joint.recovery.structure.test,
         indiv.recovery.structure.train, indiv.recovery.structure.test,
         mse_y_train, mse_y_test, coverage_EY_train, coverage_EY_test, mod.ranks, 
         file = paste0("~/BayesianPMF/03Simulations/", mod, "/", mod, "_sim_", sim_iter, "_s2nX_", s2nX, "_s2nY_", s2nY, ".rda"))
    
    res <- c(joint.recovery.structure.train, joint.recovery.structure.test, indiv.recovery.structure.train, indiv.recovery.structure.test, mse_y_train, mse_y_test, coverage_EY_train, coverage_EY_test, mod.ranks)
    names(res) <- c("joint mse (train)", "joint mse (test)", "indiv mse (train)", "indiv mse (test)", "y mse (train)", "y mse (test)", "coverage y (train)", "coverage y (test)", "joint rank", paste("indiv rank", 1:q))
    
    res
  }
  stopCluster(cl)
  end <- Sys.time()
  end - start
  
  # Return
  sim_results
  
}

# Simulation study for assessing adjustment of label switching (permutation invariance)
identifiability_sim <- function(p.vec, n, ranks, response, true_params, model_params, sparsity = TRUE, nsim, nsample, n_clust = 10) {
  
  # ---------------------------------------------------------------------------
  # Arguments:
  #
  # p.vec = number of features per source
  # n = sample size
  # ranks = vector of joint and individual ranks = c(joint rank, indiv rank 1, indiv rank 2, ...)
  # response = string in c(NULL, "continuous", "binary")
  # true_params = the list of true parameters under which to generate data
  # sparsity = should the data be generated with sparsity in the response? (Boolean)
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
             "sigma.rmt", "estim_sigma", "softSVD", "frob", "sample2", "logSum", "label_switching")
  packs <- c("Matrix", "MASS", "truncnorm", "r.jive")
  sim_results <- foreach (sim_iter = 1:nsim, .packages = packs, .export = funcs, .verbose = TRUE, .combine = rbind) %dopar% {
    
    # Set seed
    seed <- sim_iter
    set.seed(seed)
    
    # -------------------------------------------------------------------------
    # Generate data 
    # -------------------------------------------------------------------------
    
    # Generate n samples
    sim_data <- bpmf_data(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response, missingness = NULL, entrywise = NULL, prop_missing = NULL, sparsity = sparsity)
    
    # Saving the data
    true_data <- sim_data$data
    true_data_list <- lapply(true_data, function(s) s)
    q <- length(p.vec)
    joint.structure <- sim_data$joint.structure
    indiv.structure <- sim_data$indiv.structure
    
    # The response
    Y <- sim_data$Y
    
    # The inclusion indicators
    true.gammas <- sim_data$gamma
    
    # Save the true V and Vs
    true.V <- sim_data$V
    true.Vs <- sim_data$Vs

    # Save a burn-in
    burnin <- nsample/2
    thinned_iters <- seq(1, nsample, by = 10)
    thinned_iters_burnin <- seq(burnin, nsample, by = 10)
    
    # -------------------------------------------------------------------------
    # Run the model and extract results
    # -------------------------------------------------------------------------
    
    # Gibbs sampling
    res <- bpmf(data = true_data, Y = Y, nninit = FALSE, model_params = model_params, ranks = ranks, scores = NULL, sparsity = sparsity, nsample, progress = TRUE)
    
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
    # Apply the label switching algorithm 
    # -------------------------------------------------------------------------
    
    # Apply algorithm
    res.ls <- label_switching(U.draw, V.draw, W.draw, Vs.draw, betas = betas.draw, 
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
    # Compare the true gammas indicator vector to the label swapped and non-label
    # swapped versions
    # -------------------------------------------------------------------------
    
    # Label swapped
    ssd.ls <- frob(true.gammas - post.ls)/length(true.gammas)
    
    # Non-label swapped
    ssd.non.ls <- frob(true.gammas - post.non.ls)/length(true.gammas)
    
    # -------------------------------------------------------------------------
    # Return
    # -------------------------------------------------------------------------
    
    ssd <- c(ssd.ls, ssd.non.ls)
    names(ssd) <- c("Corrected", "Non-Correceted")
    ssd
  }
  
  # Return the results
  sim_results
}

# -----------------------------------------------------------------------------
# sJIVE functions
# -----------------------------------------------------------------------------

# sJIVE
sJIVE.converge <- function(X, Y, eta=NULL, max.iter=1000, threshold = 0.001,
                           show.error =F, rankJ=NULL, rankA=NULL, 
                           show.message=T, reduce.dim=T, center.scale=T){
  #X = list(X_1  X_2  ... X_k) with each row centered and scaled
  #Y is continuous vector centered and scaled
  #eta is between 0 and 1, when eta=NULL, no weighting is done
  #rank is prespecified rank of svd approximation
  
  
  optim.error <- function(X.tilde, U, theta1, Sj, W, Si, theta2, k, obs){
    WS <- NULL; thetaS <- 0
    for(i in 1:k){
      temp <- W[[i]] %*% Si[[i]]
      WS <- rbind(WS, temp)
      thetaS <- thetaS + theta2[[i]] %*% Si[[i]]
    }
    WS.new <- rbind(WS, thetaS)
    error  <- norm(X.tilde - rbind(U, theta1) %*% Sj - WS.new, type = "F")^2
    return(error)
  }
  
  
  #Set Ranks 
  if(is.null(rankJ) | is.null(rankA)){
    temp <- sJIVE.ranks(X,Y, eta=eta, max.iter = max.iter)
    rankJ <- temp$rankJ
    rankA <- temp$rankA
    cat(paste0("Using rankJ= ", rankJ, " and rankA= ", paste(rankA, collapse = " "), "\n"))
  }
  k <- length(X)
  
  
  if(center.scale){
    Y <- as.numeric(scale(as.numeric(Y)))
    for(i in 1:k){
      for(j in 1:nrow(X[[i]])){
        X[[i]][j,] <- as.numeric(scale(as.numeric(X[[i]][j,])))
      }
    }
  }
  
  k <- length(X)
  n <- ncol(X[[1]])
  svd.bigX <- list()
  for(i in 1:k){
    if(ncol(X[[i]]) != n){
      stop("Number of columns differ between datasets")
    }
    if(nrow(X[[i]])>n & reduce.dim){
      svd.bigX[[i]]<- svd(X[[i]], nu=n)
      if(svd.bigX[[i]]$u[1,1]<0){
        svd.bigX[[i]]$u <- svd.bigX[[i]]$u *-1
        svd.bigX[[i]]$v <- svd.bigX[[i]]$v *-1
      }
      X[[i]] <- diag(svd.bigX[[i]]$d) %*% t(svd.bigX[[i]]$v)
    }else{
      svd.bigX[[i]]<-NULL
    }
  }
  svd.bigX[[k+1]] <- 1
  
  #Step 1: Initialize values
  k <- length(X)
  obs <- list(); temp <- 0
  for(i in 1:k){  
    max.obs <- max(temp)
    temp <- (max.obs+1):(max.obs+nrow(X[[i]]))
    obs[[i]] <- temp
  }
  X.tilde <- NULL
  if(is.null(eta)==T){
    for(i in 1:k){X.tilde <- rbind(X.tilde, X[[i]]) }
    X.tilde <- rbind(X.tilde, Y)
  }else{
    for(i in 1:k){X.tilde <- rbind(X.tilde, sqrt(eta) * X[[i]])}
    X.tilde <- rbind(X.tilde, sqrt(1-eta)* Y)
  }
  y <- nrow(X.tilde)
  n <- ncol(X.tilde)
  
  #Initialize U, theta1, and Sj 
  if(rankJ == 0){
    X.svd <- svd(X.tilde, nu=1, nv=1)
    U.new <- U.old <- as.matrix(X.svd$u[-y,]) * 0
    theta1.new <- theta1.old <- t(as.matrix(X.svd$u[y,])) * 0
    Sj.new <- Sj.old <- as.matrix(X.svd$d[1] * t(X.svd$v)) * 0
  }else{
    X.svd <- svd(X.tilde, nu=rankJ, nv=rankJ)
    U.new <- U.old <- as.matrix(X.svd$u[-y,])
    theta1.new <- theta1.old <- t(as.matrix(X.svd$u[y,]))
    if(rankJ==1){Sj.new <- Sj.old <- as.matrix(X.svd$d[1] * t(X.svd$v))
    }else{Sj.new <- Sj.old <- diag(X.svd$d[1:rankJ]) %*% t(X.svd$v) }
  }
  
  #Initialize W and S_i = 0
  W.old <- S.old <- theta2.old <- list() 
  WS <- NULL
  for(i in 1:k){  
    #Get Wi, Si, and theta2i
    X.tilde_i <- X.tilde[c(obs[[i]],y),]
    yi <- nrow(X.tilde_i)
    xi <- (X.tilde_i - rbind(as.matrix(U.new[obs[[i]],]), theta1.new) %*% Sj.new)
    if(rankJ==0){
      vi <- diag(rep(1,n)) 
    }else{
      vi <- diag(rep(1,n)) -  X.svd$v %*% t(X.svd$v)
    }
    if(rankA[i] == 0){
      X2.svd <- svd(xi, nu=1, nv=1)
      W.old[[i]] <- as.matrix(X2.svd$u[-yi,]) *0
      theta2.old[[i]] <- t(as.matrix(X2.svd$u[yi,])) *0
      S.old[[i]] <- as.matrix(X2.svd$d[1] * t(X2.svd$v)) *0
    }else{
      X2.svd <- svd(xi %*% vi, nu=rankA[i], nv=rankA[i])
      W.old[[i]] <- as.matrix(X2.svd$u[-yi,])
      theta2.old[[i]] <- t(as.matrix(X2.svd$u[yi,]))
      if(rankA[i]==1){S.old[[i]] <- as.matrix(X2.svd$d[1] * t(X2.svd$v))
      }else{S.old[[i]] <- diag(X2.svd$d[1:rankA[i]]) %*% t(X2.svd$v) }
    }
    WS <- rbind(WS, W.old[[i]] %*% S.old[[i]])
  }
  error.old <- optim.error(X.tilde, U.old, theta1.old, Sj.old, W.old, S.old, theta2.old, k, obs)
  if(show.error){cat(paste0("Iter: ", 0, "  Error: ", error.old, "\n"))}
  e.vec <- NULL
  S.new <- S.old
  theta2.new <- theta2.old
  W.new <- W.old
  thetaS.sum <- matrix(rep(0,n), ncol=n)
  
  #Step 2: Interatively optimize problem
  for(iter in 1:max.iter){
    #Get U, theta1, and Sj 
    if(rankJ >0){
      X.svd <- svd(X.tilde-rbind(WS, thetaS.sum), nu=rankJ, nv=rankJ)
      neg <- sign(X.svd$u[1,1])
      U.new <- as.matrix(X.svd$u[-y,]) * neg
      theta1.new <- t(as.matrix(X.svd$u[y,])) *neg
      if(rankJ==1){Sj.new <- X.svd$d[1] * t(X.svd$v) * neg
      }else{Sj.new <- diag(X.svd$d[1:rankJ]) %*% t(X.svd$v) * neg }
    }
    
    if(show.error){
      e <- optim.error(X.tilde, U.new, theta1.new, Sj.new, W.old, S.old, theta2.old, k, obs)
      cat(paste0("Iter: ", iter, "  Error: ", e, "  Updated: U, theta1, Sj \n"))
    }
    
    #Get Wi, Si, and theta2i
    thetaS.sum <- 0; WS <- NULL
    for(i in 1:k){
      thetaS <-0
      for(j in 1:k){if(j != i){ thetaS <- thetaS + theta2.new[[j]] %*% S.new[[j]]}}
      X.tilde_i <- X.tilde[c(obs[[i]],y),]
      yi <- nrow(X.tilde_i)
      X.tilde_i[yi,] <- X.tilde_i[yi,] - thetaS
      
      if(rankA[i]>0){
        #Get Wi, Si, and theta2i
        xi <- (X.tilde_i - rbind(as.matrix(U.new[obs[[i]],]), theta1.new) %*% Sj.new)
        if(rankJ==0){
          vi <- diag(rep(1,n))
        }else{
          vi <- diag(rep(1,n)) -  X.svd$v %*% t(X.svd$v) 
        }
        X2.svd <- svd(xi %*% vi, nu=rankA[i], nv=rankA[i])
        
        #X2.svd <- svd(X.tilde_i - rbind(as.matrix(U.new[obs[[i]],]), theta1.new) %*% Sj.new, nu=rankA[i], nv=rankA[i])
        neg2 <- sign(X2.svd$u[1,1])
        W.new[[i]] <- as.matrix(X2.svd$u[-yi,]) * neg2
        theta2.new[[i]] <- t(as.matrix(X2.svd$u[yi,])) * neg2
        if(rankA[i]==1){S.new[[i]] <- as.matrix(X2.svd$d[1] * t(X2.svd$v)) * neg2
        }else{S.new[[i]] <- diag(X2.svd$d[1:rankA[i]]) %*% t(X2.svd$v) * neg2}
      }
      
      #prep for next iteration
      thetaS.sum <- thetaS.sum + theta2.new[[i]] %*% S.new[[i]]
      WS <- rbind(WS, W.new[[i]] %*% S.new[[i]])
      
      if(show.error){
        e <- optim.error(X.tilde, U.new, theta1.new, Sj.new, W.new, S.new, theta2.new, k, obs)
        cat(paste0("Iter: ", iter, "  Error: ", e, "  Updated: Wi, theta2i, Si  i=", i, "\n"))
      }
    }
    
    #Figure out the error
    error <- optim.error(X.tilde, U.new, theta1.new, Sj.new, W.new, S.new, theta2.new, k, obs)
    if(abs(error.old-error) < threshold){
      #If converged, then stop loop
      if(show.message){cat(paste0("Converged in ", iter, " iterations \n"))}
      break
    }else if(iter == max.iter){
      if(show.message){cat(paste0("Warning: ", iter, " iterations reached \n"))}
      break
    }else{
      #If didn't converge, prep for another loop
      e.vec <- c(e.vec, error)
      U.old <- U.new
      W.old <- W.new
      theta2.old <- theta2.new 
      theta1.old <- theta1.new
      Sj.old <- Sj.new
      S.old <- S.new
      error.old <- error
    }
  }
  
  #Scale so first value in U and W are always positive
  if(U.new[1,1]<0){
    U.new <- -1 * U.new
    theta1.new <- -1 * theta1.new
    Sj.new <- -1 * Sj.new
  }
  for(i in 1:k){
    if(W.new[[i]][1,1]<0){
      W.new[[i]] <- W.new[[i]] * -1
      theta_2[[i]] <- theta_2[[i]] * -1
      S.new[[i]] <- S.new(i) * -1
    }
  }
  
  
  #Step 3: Export the results
  U_i <- W <- theta_2 <- list()
  if(is.null(eta)){
    for(i in 1:k){
      U_i[[i]] <- U.new[obs[[i]],]
    }
    theta_1 <- theta1.new
    W <- W.new
    theta_2 <- theta2.new
  }else{
    theta_1 <- (1/sqrt(1-eta)) * theta1.new
    U.new<- (1/sqrt(eta)) * U.new
    for(i in 1:k){ 
      W[[i]] <- (1/sqrt(eta)) * W.new[[i]]
      theta_2[[i]] <- (1/sqrt(1-eta)) * theta2.new[[i]]
    }
    #ReScale to be norm 1
    for(j in 1:ncol(U.new)){
      U.norm <-norm(rbind(as.matrix(U.new[,j]), theta_1[j]),type="F")^2 
      if(U.norm==0){
        U.new[,j] <- as.matrix(U.new[,j])
        theta_1[j] <- theta_1[j]
      }else{
        U.new[,j] <- as.matrix(U.new[,j])/sqrt(U.norm)
        theta_1[j] <- theta_1[j]/sqrt(U.norm)
      }
      Sj.new[j,] <- Sj.new[j,] * sqrt(U.norm)
    }
    for(i in 1:k){
      for(j in 1:ncol(W[[i]])){
        W.norm <-norm(rbind(as.matrix(W[[i]][,j]),theta_2[[i]][j]),type="F")^2 
        if(W.norm==0){
          W[[i]][,j] <- as.matrix(W[[i]][,j])
          theta_2[[i]][j] <- theta_2[[i]][j]
        }else{
          W[[i]][,j] <- as.matrix(W[[i]][,j])/sqrt(W.norm)
          theta_2[[i]][j] <- theta_2[[i]][j]/sqrt(W.norm)
        }
        S.new[[i]][j,] <- S.new[[i]][j,] * sqrt(W.norm)
        U_i[[i]] <- U.new[obs[[i]],]
      }
    }
    
  }
  
  Ypred <- theta_1 %*% Sj.new
  for(i in 1:k){
    Ypred <- Ypred + theta_2[[i]] %*% S.new[[i]]
  }
  
  
  #Map X back to original space
  for(i in 1:k){
    if(is.null(svd.bigX[[i]]) == FALSE){
      U_i[[i]] <- svd.bigX[[i]]$u %*% U_i[[i]]
      W[[i]] <- svd.bigX[[i]]$u %*% W[[i]]
    }
  }
  
  return(list(S_J=Sj.new, S_I=S.new, U_I=U_i, W_I=W,
              theta1=theta_1, theta2=theta_2, fittedY=Ypred,
              error=error, all.error=e.vec,
              iterations = iter, rankJ=rankJ, rankA=rankA, eta=eta))
}

sJIVE.predict <- function(sJIVE.fit, newdata, threshold = 0.001, max.iter=2000){
  ##############################################
  # sJIVE.fit is the output from sJIVE
  # newdata is list with the same predictors and
  #     number of datasets as used in sJIVE.fit
  ##############################################
  sJIVE.pred.err <- function(X.tilde, U, Sj, W, Si, k){
    J <- A <- NULL
    for(i in 1:k){
      J <- rbind(J, as.matrix(U[[i]]) %*% as.matrix(Sj))
      A <- rbind(A, as.matrix(W[[i]]) %*% as.matrix(Si[[i]]))
    }
    temp <- X.tilde - J - A
    error <- norm(temp, type = "F")^2
    return(error)
  }
  
  
  
  if(sJIVE.fit$rankJ==0 & sum(sJIVE.fit$rankA)==0){
    return(list(Ypred = 0,
                Sj = 0,
                Si = 0,
                iteration = 0,
                error = NA))
  }
  
  #Initalize values
  k <- length(newdata)
  n <- ncol(newdata[[1]])
  W <- sJIVE.fit$W_I
  U <- sJIVE.fit$U_I
  rankJ <- ncol(as.matrix(U[[1]]))
  Sj <- matrix(rep(0,rankJ*n), ncol = n)
  
  obs <- rankA <- Si <- list(); temp <- 0; X.tilde <- NULL
  for(i in 1:k){  
    max.obs <- max(temp)
    temp <- (max.obs+1):(max.obs+nrow(newdata[[i]]))
    obs[[i]] <- temp
    
    X.tilde <- rbind(X.tilde, newdata[[i]])
    
    rankA[[i]] <- ncol(as.matrix(W[[i]]))
    Si[[i]] <- matrix(rep(0, rankA[[i]]*n), ncol=n)
  }
  
  #Get Error
  error.old <- sJIVE.pred.err(X.tilde, U, Sj, W, Si, k)
  
  for(iter in 1:max.iter){
    
    #Update Sj
    U.mat <- A <- NULL
    for(i in 1:k){
      U.mat <- rbind(U.mat, as.matrix(U[[i]]))
      A <- rbind(A, as.matrix(W[[i]]) %*% as.matrix(Si[[i]]))
    }
    Sj <- t(U.mat) %*% (X.tilde - A)
    
    #Update Si
    for(i in 1:k){
      Si[[i]] <- t(W[[i]]) %*% (newdata[[i]] - U[[i]] %*% Sj)
    }
    
    #Get Error
    error.new <- sJIVE.pred.err(X.tilde, U, Sj, W, Si, k)
    
    #Check for Convergence
    if(abs(error.old - error.new) < threshold){
      break
    }else{
      error.old <- error.new
    }
  }
  
  Ypred <- sJIVE.fit$theta1 %*% Sj
  for(i in 1:k){
    Ypred <- Ypred + sJIVE.fit$theta2[[i]] %*% Si[[i]]
  }
  
  return(list(Ypred = Ypred,
              Sj = Sj,
              Si = Si,
              iteration = iter,
              error = error.new))
}

sJIVE.ranks <- function(X, Y, eta=NULL, max.iter=1000, threshold = 0.01,
                        max.rank=100, center.scale=T,
                        reduce.dim=T){
  cat("Estimating joint and individual ranks via cross-validation... \n")
  k <- length(X)
  n <- ncol(X[[1]])
  
  #get cv folds
  fold <- list()
  cutoff <- round(quantile(1:n, c(.2,.4,.6,.8)))
  fold[[1]] <- 1:cutoff[1]
  fold[[2]] <- (cutoff[1]+1):cutoff[2]
  fold[[3]] <- (cutoff[2]+1):cutoff[3]
  fold[[4]] <- (cutoff[3]+1):cutoff[4]
  fold[[5]] <- (cutoff[4]+1):n
  
  rankJ <- 0
  rankA <- rep(0,k)
  
  #initialize error
  error.old  <- NULL
  for(i in 1:5){
    #get X train, Y.train
    train.X <- test.X <- list()
    for(j in 1:k){
      train.X[[j]] <- X[[j]][,-fold[[i]]]
      test.X[[j]] <- X[[j]][,fold[[i]]] 
    }
    train.Y <- Y[-fold[[i]]]
    test.Y <- Y[fold[[i]]]
    fit.old <- sJIVE.converge(train.X, train.Y, eta=eta, max.iter = max.iter, 
                              rankJ = rankJ, rankA = rankA, show.message = F, 
                              center.scale=center.scale,
                              reduce.dim=reduce.dim)
    new.data <- sJIVE.predict(fit.old, test.X) 
    error.old <- c(error.old, sum((new.data$Ypred-test.Y)^2) )
  }
  err.old <- mean(error.old)
  
  #iteravely add ranks
  for(iter in 1:1000){
    error.j <- NULL; error.a <- list()
    
    for(i in 1:5){
      #get X train, Y.train
      train.X <- test.X <- list()
      for(j in 1:k){
        train.X[[j]] <- X[[j]][,-fold[[i]]]
        test.X[[j]] <- X[[j]][,fold[[i]]] 
      }
      train.Y <- Y[-fold[[i]]]
      test.Y <- Y[fold[[i]]]
      
      #Add rank to joint
      if(rankJ < max.rank){
        fit.j <- sJIVE.converge(train.X, train.Y, eta=eta, max.iter = max.iter, 
                                rankJ = rankJ+1, rankA = rankA, show.message=F,
                                center.scale=center.scale,
                                reduce.dim=reduce.dim)
        new.data <- sJIVE.predict(fit.j, test.X)
        error.j <- c(error.j, sum((new.data$Ypred-test.Y)^2) )
      }else{
        error.j <- c(error.j, 99999999)
      }
      
      #Add rank to individual
      for(j in 1:k){
        if(rankA[j] < max.rank){
          rankA.new <- rankA
          rankA.new[j] <- rankA.new[j]+1
          
          #Add rank to individual
          fit.a <- sJIVE.converge(train.X, train.Y, eta=eta, max.iter = max.iter, 
                                  rankJ = rankJ, rankA = rankA.new, show.message=F,
                                  center.scale=center.scale,
                                  reduce.dim=reduce.dim)
          new.data <- sJIVE.predict(fit.a, test.X)
          if(length(error.a)<j){error.a[[j]] <- NA}
          error.a[[j]] <- c(error.a[[j]], sum((new.data$Ypred-test.Y)^2) )
        }else{
          if(length(error.a)<j){error.a[[j]] <- NA}
          error.a[[j]] <- c(error.a[[j]], 99999999 )
        }
      }
      
    }
    
    #average over folds
    err.j <- mean(error.j)
    err.a <- lapply(error.a, function(x) mean(x, na.rm=T))
    #cat(error.j)
    #cat(error.a[[1]])
    #cat(error.a[[2]])
    
    #Determine which rank to increase
    asd <- c(err.old - err.j, err.old - as.vector(unlist(err.a)))
    if(max(asd) < threshold){
      break
    }else{
      temp <- which(asd == max(asd))
      if(temp==1){
        #cat(asd)
        rankJ <- rankJ+1
        err.old <- err.j
      }else{
        rankA[temp-1] <- rankA[temp-1]+1
        err.old <- err.a[[temp-1]]
      }
    }
    #cat(paste0("Joint Rank: ", rankJ, "\n"))
    #cat(paste0("Individual Ranks: ", rankA, "\n"))
  }
  
  return(list(rankJ = rankJ, 
              rankA = rankA,
              error = err.old,
              error.joint = err.j,
              error.individual = err.a))
}

sJIVE <- function(X, Y, rankJ = NULL, rankA=NULL,eta=NULL, max.iter=1000,
                  threshold = 0.001,  method="permute",
                  center.scale=TRUE, reduce.dim = TRUE){
  ############################################################################
  #X is a list of 2 or more datasets, each with dimensions p_i by n
  #Y is continuous vector length n
  #eta is a tuning parameter between 0 and 1. When eta=NULL, a gridsearch
  #   is conducted to tune eta. You can specify a value of eta to use, 
  #   or supply a vector of eta values for sJIVE to consider.
  #rankJ is a value for the low-rank of the joint component
  #rankA is a vector of the ranks for each X dataset. When rankJ or rankA
  #   are NULL, a rank selection method (see method) will choose ranks
  #max.iter specifies the maximum number of iterations that will run
  #threshold specifies the criteria to determine when algorithm has 
  #   converged
  #Method = c("permute", "CV"). When ranks are not specified, ranks will
  #   be determined by JIVE's permutation method, or sJIVE's 
  #   cross-validation method
  #Center.scale is a true/false indicator for whether or not to center and 
  #   scale the data prior to fitting.
  ############################################################################
  
  k <- length(X)
  n <- ncol(X[[1]])
  if(length(Y) != n){stop("Number of columns differ between datasets")}
  
  if(center.scale){
    Y <- as.numeric(scale(as.numeric(Y)))
    for(i in 1:k){
      for(j in 1:nrow(X[[i]])){
        X[[i]][j,] <- as.numeric(scale(as.numeric(X[[i]][j,])))
      }
    }
  }
  
  
  if(is.null(eta)){
    e.vec=c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  }else{
    e.vec=eta
  }
  
  if(is.null(rankJ) | is.null(rankA)){
    if(method=="permute"){
      temp <- jive(X, Y, center=F, scale=F,orthIndiv = F)
      rankJ <- temp$rankJ
      rankA <- temp$rankA
    }else if(method=="CV"){
      temp <- sJIVE.ranks(X,Y, eta=eta, max.iter = max.iter, center.scale=center.scale,
                          reduce.dim=reduce.dim)
      rankJ <- temp$rankJ
      rankA <- temp$rankA
    }else{
      errorCondition("Invalid method chosen")
    }
    cat(paste0("Using rankJ= ", rankJ, " and rankA= ", paste(rankA, collapse = " "), "\n"))
  }
  
  if(length(e.vec)>1){
    #get cv folds
    n <- length(Y)
    fold <- list()
    cutoff <- round(quantile(1:n, c(.2,.4,.6,.8)))
    fold[[1]] <- 1:cutoff[1]
    fold[[2]] <- (cutoff[1]+1):cutoff[2]
    fold[[3]] <- (cutoff[2]+1):cutoff[3]
    fold[[4]] <- (cutoff[3]+1):cutoff[4]
    fold[[5]] <- (cutoff[4]+1):n
    
    cat("Choosing Tuning Parameter: eta \n")
    err.test <- NA
    for(e in e.vec){
      err.fold <- NA
      for(i in 1:5){
        #Get train/test sets
        sub.train.x <- sub.test.x <- list()
        sub.train.y <- Y[-fold[[i]]]
        sub.test.y <- Y[fold[[i]]]
        for(j in 1:length(X)){
          sub.train.x[[j]] <- X[[j]][,-fold[[i]]]
          sub.test.x[[j]] <- X[[j]][,fold[[i]]]
        }
        fit1 <- sJIVE.converge(sub.train.x, sub.train.y, max.iter = max.iter, 
                               rankJ = rankJ, rankA = rankA, eta = e, 
                               show.message=F, center.scale=center.scale,
                               reduce.dim=reduce.dim)
        #Record Error for fold
        fit_test1 <- sJIVE.predict(fit1, sub.test.x)
        fit.mse <- sum((fit_test1$Ypred - sub.test.y)^2)/length(sub.test.y)
        err.fold <- c(err.fold, fit.mse)
      }
      
      #Record Test Error (using validation set)
      fit.mse <- mean(err.fold, na.rm = T)
      err.test <- c(err.test, fit.mse)
      if(min(err.test, na.rm = T) == fit.mse){
        best.eta <- e
      }
    }
    cat(paste0("Using eta= ", best.eta, "\n"))
    test.best <- sJIVE.converge(X, Y, max.iter = max.iter, 
                                rankJ = rankJ, rankA = rankA, eta = best.eta,
                                threshold = threshold, center.scale=center.scale,
                                reduce.dim=reduce.dim)
  }else{
    test.best <- sJIVE.converge(X, Y, max.iter = max.iter, 
                                rankJ = rankJ, rankA = rankA, eta = e.vec,
                                threshold = threshold, center.scale=center.scale,
                                reduce.dim=reduce.dim)
  }
  
  
  return(test.best)
  
}

# -----------------------------------------------------------------------------
# BIP Functions
# -----------------------------------------------------------------------------

BIP <- function(dataList=dataList,IndicVar=IndicVar, groupList=NULL,Method=Method,nbrcomp=4, sample=5000, burnin=1000,nbrmaxmodels=50,
                priorcompselv=c(1,1),priorcompselo=c(1,1),priorb0=c(2,2),priorb=c(1,1),priorgrpsel=c(1,1),probvarsel=0.05) {
  
  if (sample<=burnin){
    stop("Argument burnin must be smaller than sample: the number of MCMC iterations.")
  }
  if (sample<=20){
    stop("Please specify a larger number of MCMC iterations")
  }
  
  if (is.null(nbrcomp)){
    D=which(IndicVar==0);
    mysvd=lapply(D, function(i)  svd(dataList[[i]]))
    mysumsvd=lapply(1:length(D), function(i) cumsum(mysvd[[i]]$d)/max(cumsum(mysvd[[i]]$d))*100)
    KMax=max(unlist(lapply(1:length(D), function(i) min(which(mysumsvd[[i]]>=80, arr.ind = TRUE))))) #chooses maximum from D cumulative proportions
    nbrcomp=min(KMax+1,10);
  }
  
  
  #Method="GroupInfo"
  if (Method=="BIP"){
    meth=0;
  } else if (Method=="BIPnet") {
    meth=1;
  } else {
    stop("You must provide either BIP or BIPnet")
  }
  
  Np=length(dataList);P=NULL
  n=nrow(dataList[[1]])
  P=NULL;
  MeanData=list();SD=list();
  for (i in 1:Np){
    dataList[[i]]=as.matrix(dataList[[i]])
    P[i]=ncol(dataList[[i]]);
    if ((IndicVar[i]==0)||(IndicVar[i]==2)){
      MeanData[[i]]=apply(dataList[[i]],2,mean);
      SD[[i]]=apply(dataList[[i]],2,sd);
      dataList[[i]]=scale(as.matrix(dataList[[i]]),T,T)
    }
    dataList[[i]]=t(dataList[[i]]);
  }
  datasets=unlist(dataList)
  
  if (is.null(groupList) && (meth==1)){
    stop("You must provide a list of grouping information")
  } else if (is.null(groupList)){
    for (ll in 1:length(IndicVar)){
      groupList[[ll]]=matrix(1,P[ll],1)
    }
  }
  
  
  ll=1;
  K=NULL;
  for (i in 1:Np){
    if (IndicVar[i]!=0){
      K[i]=1;
    } else {
      K[i]=ncol(groupList[[ll]]);
      ll=ll+1;
    }
  }
  groupList=lapply(groupList,t);Paths=unlist(groupList);
  
  result <- .C("mainfunction",
               Method1=as.integer(meth),n1=as.integer(n),P=as.integer(P),r1=as.integer(nbrcomp),Np1=as.integer(Np), datasets=as.double(datasets),IndVar=as.integer(IndicVar), K=as.integer(K), Paths=as.integer(Paths),maxmodel1=as.integer(nbrmaxmodels),                nbrsample1=as.integer(sample), burninsample1=as.integer(burnin),
               CompoSelMean=as.double(rep(0,Np*nbrcomp)),VarSelMean=as.double(rep(0,nbrcomp*sum(P))),
               VarSelMeanGlobal=as.double(rep(0,sum(P))),GrpSelMean=as.double(rep(0,nbrcomp*sum(K))),
               GrpEffectMean=as.double(rep(0,nbrcomp*sum(K))),IntGrpMean=as.double(rep(0,nbrcomp*Np)),
               EstU=as.double(rep(0,n*nbrcomp)),EstSig2=as.double(rep(0,sum(P))),InterceptMean=as.double(rep(0,1)),
               EstLoadMod=as.double(rep(0,nbrmaxmodels*nbrcomp*sum(P))),EstLoad=as.double(rep(0,nbrcomp*sum(P))),
               nbrmodel1=as.integer(0),postgam=rep(0,nbrmaxmodels),priorcompsel=priorcompselv,
               priorcompselo=priorcompselo,priorb0=priorb0,priorb=as.double(priorb),priorgrpsel=priorgrpsel,probvarsel=as.double(probvarsel));
  
  
  
  reseffect=result$EstLoadMod;
  nbrmodel=result$nbrmodel1;
  EstLoadModel=rep( list(list()), nbrmodel ) 
  for (mo in 1:nbrmodel){
    for (m in 1:Np){
      x=sum(P[1:(m-1)]);y=sum(P[1:m]);
      if (m==1) {x=0}
      init=1+x*nbrcomp+(mo-1)*nbrcomp*sum(P);
      final=y*nbrcomp+(mo-1)*nbrcomp*sum(P);
      EstLoadModel[[mo]][[m]]=matrix(reseffect[init:final],nbrcomp,byrow=T);
    }
  }
  
  resvarsel1=result$VarSelMeanGlobal;
  resvarsel2=result$VarSelMean;
  resvarsel3=result$GrpSelMean;
  resvarsel4=result$GrpEffectMean
  resvarsel5=result$EstLoad;
  EstimateU=matrix(result$EstU,n,byrow=T);
  CompoSelMean=matrix(result$CompoSelMean,Np,byrow=T);
  IntGrpMean=matrix(result$IntGrpMean,Np,byrow=T);
  VarSelMeanGlobal=list();
  VarSelMean=list();
  GrpSelMean=list();
  GrpEffectMean=list();
  EstLoad=list();
  EstSig2=list();
  m1=m2=m3=1;
  for (m in 1:Np){
    VarSelMeanGlobal[[m]]=resvarsel1[m1:(m1-1+P[m])]
    VarSelMean[[m]]=matrix(resvarsel2[m2:(m2-1+P[m]*nbrcomp)],P[m],byrow=T)
    GrpSelMean[[m]]=matrix(resvarsel3[m3:(m3-1+K[m]*nbrcomp)],nbrcomp,byrow=T)
    GrpEffectMean[[m]]=matrix(resvarsel4[m3:(m3-1+K[m]*nbrcomp)],nbrcomp,byrow=T);
    EstLoad[[m]]=matrix(resvarsel5[m2:(m2-1+P[m]*nbrcomp)],nbrcomp,byrow=T);
    EstSig2[[m]]=result$EstSig2[m1:(m1-1+P[m])]
    m1=m1+P[m];
    m2=m2+P[m]*nbrcomp;
    m3=m3+K[m]*nbrcomp;
    
  }
  
  return (list(EstU=EstimateU,VarSelMean=VarSelMean,VarSelMeanGlobal=VarSelMeanGlobal,CompoSelMean=CompoSelMean,GrpSelMean=GrpSelMean,GrpEffectMean=GrpEffectMean,IntGrpMean=IntGrpMean,EstLoad=EstLoad,EstLoadModel=EstLoadModel,nbrmodel=result$nbrmodel1,EstSig2=EstSig2,EstIntcp=result$InterceptMean,PostGam=result$postgam,IndicVar=IndicVar,nbrcomp=nbrcomp,MeanData=MeanData,SDData=SD))
}

BIPpredict <- function(dataListNew=dataListNew, Result=Result, meth="BMA"){
  IndicVar=Result$IndicVar;nbrcomp=Result$nbrcomp;
  X_new=list();
  np=which(IndicVar==1)
  Np=length(IndicVar);
  j=1;
  for (i in  1:Np){
    X_new[[i]]=NULL;
    if (IndicVar[i]!=1){ ## We do not center the response variable
      X_new[[i]]=dataListNew[[j]];
      X_new[[i]]=t((t(X_new[[i]])-Result$MeanData[[i]])/Result$SDData[[i]])
      j=j+1;
    }
  }
  
  nbrmodel=Result$nbrmodel;
  EstLoad=Result$EstLoadModel
  Upredict=matrix(0,nrow(X_new[[1]]),nbrcomp)
  #Upredict=matrix(0,nrow(X_new[[1]]),ncol(Latent_new))
  ypredict=rep(0,nrow(X_new[[1]]));
  nbrmodel1=nbrmodel;
  if (meth!="BMA"){
    nbrmodel1=1;
  }
  
  for (nb in 1:nbrmodel1){
    SelLoadPacked=NULL;SelLoadPackedFull=NULL
    SelVarXnew=NULL;SelVarXnewFull=NULL;
    Sig2=NULL;Sig2Full=NULL;
    if (meth=="BMA"){
      EstL=EstLoad[[nb]];
      pb=Result$PostGam[nb];
    } else {
      EstL=Result$EstLoad; pb=1;
    }
    for (m in 1:Np){
      
      #if (m<=Np-1){
      if (IndicVar[m]!=1){ # We exclude the response variable
        nc=nrow(EstL[[m]]);
        selvar=apply(EstL[[m]],2, function(x) length((which(x==0)))!=nc)
        SelLoadPacked=cbind(SelLoadPacked,EstL[[m]][,selvar]);
        Sig2=c(Sig2,Result$EstSig2[[m]][selvar])
        SelVarXnew=cbind(SelVarXnew,X_new[[m]][,selvar])
        
        
        #SelLoadPackedFull=cbind(SelLoadPackedFull,EstL[[m]][,selvar]);
        #Sig2Full=c(Sig2Full,Result$EstSig2[[m]][selvar])
        #SelVarXnewFull=cbind(SelVarXnewFull,X_new[[m]][,selvar])
      }
    }
    #SelLoadPacked=as.matrix(SelLoadPacked);
    #Sigma22=t(SelLoadPacked)%*%SelLoadPacked+diag(Sig2);
    #AoA=t(Result$EstLoad[[Np]])%*%SelLoadPacked
    #ypredict=(AoA%*%solve(Sigma22))%*%t(SelVarXnew)
    #mse=mean((as.vector(ypredict)-as.vector(y_new))^2);
    
    SelLoadPackedFull=SelLoadPacked;## I added this
    Sig2Full=Sig2;
    SelVarXnewFull=SelVarXnew;
    SelLoadPackedFull=as.matrix(SelLoadPackedFull);
    Sigma2Full=solve(SelLoadPackedFull%*%diag(1/Sig2Full)%*%t(SelLoadPackedFull)+diag(nrow(SelLoadPackedFull)));
    Upredict1=t(Sigma2Full%*%SelLoadPackedFull%*%diag(1/Sig2Full)%*%t(SelVarXnewFull))
    Upredict=Upredict+Upredict1*pb;
    ypredict=ypredict+as.matrix(Upredict1%*%EstL[[np]])*pb;
  }
  
  
  #mse=mean((as.vector(ypredict)+Result$EstIntcp-as.vector(y_new))^2);
  return (list(ypredict=as.vector(ypredict)+Result$EstIntcp,Upredtest=Upredict))
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

create_simulation_table <- function(mod.list, path.list, s2nX, s2nY) {
  
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
  simulation_results <- data.frame(s2nX = rep(s2nX, 11),
                                   s2nY = rep(s2nY, 11),
                                   Metric = numeric(11),
                                   BPMF = numeric(11),
                                   JIVE = numeric(11),
                                   `BIDIFAC+` = numeric(11),
                                   check.names = FALSE)
  simulation_results$Metric <- c("MSE Joint (train)", "MSE Joint (test)", "MSE Indiv (train)",
                                 "MSE Indiv (test)", "MSE Y (train)", "MSE Y (test)", 
                                 "Coverage EY (train)", "Coverage EY (test)", "Joint Rank", 
                                 "Indiv Rank 1", "Indiv Rank 2")
  
  # Iterate through the models
  for (mod in mod.list) {
    # Load in the files for this model
    path <- path.list[[mod]]
    all_files <- list.files(path)
    all_files_split <- strsplit(all_files, split = "_")
    
    # Save the names of the current results
    files_for_s2nX_s2nY <- all_files[sapply(all_files_split, function(file) (file[5] == s2nX) & (file[7] == paste0(s2nY, ".rda")))]
    
    # Create a temporary table to load the results into
    simulation_results_temp <- data.frame(joint.structure.train = numeric(length(files_for_s2nX_s2nY)),
                                          joint.structure.test = numeric(length(files_for_s2nX_s2nY)),
                                          indiv.structure.train = numeric(length(files_for_s2nX_s2nY)),
                                          indiv.structure.test = numeric(length(files_for_s2nX_s2nY)),
                                          mse_y_train = numeric(length(files_for_s2nX_s2nY)),
                                          mse_y_test = numeric(length(files_for_s2nX_s2nY)),
                                          coverage_EY_train = numeric(length(files_for_s2nX_s2nY)),
                                          coverage_EY_test = numeric(length(files_for_s2nX_s2nY)),
                                          joint.rank = numeric(length(files_for_s2nX_s2nY)),
                                          indiv.rank1 = numeric(length(files_for_s2nX_s2nY)),
                                          indiv.rank2 = numeric(length(files_for_s2nX_s2nY)))
    
    # Iteratively load in the results
    for (ind in 1:length(files_for_s2nX_s2nY)) {
      # Load in the results
      file <- files_for_s2nX_s2nY[ind]
      res <- load(paste0(path, file))
      simulation_results_temp[ind,] <- unlist(sapply(res, function(met) get(met)))
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
# Addressing the factor switching problem
# -----------------------------------------------------------------------------

# The label switching algorithm: matches results to the posterior mode
label_switching <- function(U.draw, V.draw, W.draw, Vs.draw, betas = NULL, gammas = NULL, r, r.vec, nsample, thinned_iters_burnin, nninit = TRUE, pivot = NULL) {
  
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
    
    # Storing the swaps and signs
    current_swaps_joint <- current_signs_joint <- c()
    current_swaps_indiv <- current_signs_indiv <- lapply(1:q, function(s) c())
    
    for (k in 1:r) {
      # for each column in the pivot, rearrange the columns in draws to match the order
      corrs <- apply(current_V[[1,1]], 2, function(current_col) cor(current_col, pivot_V[[1,1]][,k])) # compute correlation between each 
      jk <- which.max(abs(corrs)) # storing the index of the highest correlation
      ak <- sign(corrs[jk]) # storing the sign of that correlation
      
      # Setting the kth column of \tilde_V to be the ak*current_V[,jk] column. 
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

# Papastamoulis and Ntzoufras simulated annealing approach to address factor switching (from the factor.switching package)
rsp_full_sa <- function (lambda_mcmc, maxIter = 1000, threshold = 1e-06, verbose = TRUE, 
          sa_loops, rotate = TRUE, increaseIter = FALSE, temperatureSchedule = NULL, 
          printIter = 1000) {
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
    stop("Something is wrong, please check input.")
  }
  if (is.null(temperatureSchedule)) {
    temperatureSchedule <- function(sim_ann_iter) {
      1/log(sim_ann_iter + 1)
    }
  }
  objective_function <- function(lambda, lambda_hat, c_vec, 
                                 v_vec, q, p) {
    switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow = T) * 
      lambda[, v_vec]
    f <- sum((switchedMatrix - lambda_hat)^2)
    result <- vector("list", length = 2)
    result[[1]] <- f
    result[[2]] <- switchedMatrix
    names(result) <- c("value", "matrix")
    return(result)
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
  lambda_hat <- matrix(colMeans(lambda_mcmc_varimax), ncol = q, 
                       nrow = p, byrow = TRUE)
  lambda_hat_zero <- lambda_hat
  st <- 1:q
  c_vectors <- matrix(rep(1, q), ncol = q, nrow = mcmcIterations, 
                      byrow = TRUE)
  v_vectors <- matrix(1:q, ncol = q, nrow = mcmcIterations, 
                      byrow = TRUE)
  dim_all_c <- 2^q
  dim_all_v <- factorial(q)
  perm <- 1:q
  sim_an_iters <- sa_loops
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
      f_values <- numeric(sim_an_iters)
      c_vec <- c_vectors[iter, ]
      v_vec <- v_vectors[iter, ]
      obj <- objective_function(lambda = lambda, lambda_hat = lambda_hat, 
                                c_vec = c_vec, v_vec = v_vec, q = q, p = p)
      f_old <- obj$value
      f_values[1] <- f_old
      alpha <- 0.2
      temperature <- 1000
      accept_rate <- 1
      for (sim_ann_iter in 2:sim_an_iters) {
        temperature <- temperatureSchedule(sim_ann_iter)
        uR <- runif(1)
        if (uR < 0.8) {
          c_proposed <- c_vec
          rIndex <- sample(q, 1)
          c_proposed[rIndex] <- -c_vec[rIndex]
          random_pair <- sample(q, 2, replace = F)
          v_proposed <- v_vec
          v_proposed[random_pair] <- v_vec[random_pair[c(2, 
                                                         1)]]
        }
        else {
          if (uR < 0.9) {
            c_proposed <- c_vec
            rIndex <- sample(q, 1)
            c_proposed[rIndex] <- -c_vec[rIndex]
            v_proposed <- v_vec
          }
          else {
            c_proposed <- c_vec
            random_pair <- sample(q, 2, replace = F)
            v_proposed <- v_vec
            v_proposed[random_pair] <- v_vec[random_pair[c(2, 
                                                           1)]]
          }
        }
        f_new <- objective_function(lambda = lambda, 
                                    lambda_hat = lambda_hat, c_vec = c_proposed, 
                                    v_vec = v_proposed, q = q, p = p)$value
        if (f_new < f_old) {
          c_vec <- c_proposed
          v_vec <- v_proposed
          f_old <- f_new
          accept_rate <- accept_rate + 1
        }
        else {
          u <- runif(1)
          ar <- -(f_new - f_old)/temperature
          if (log(u) < ar) {
            c_vec <- c_proposed
            v_vec <- v_proposed
            f_old <- f_new
            accept_rate <- accept_rate + 1
          }
        }
        f_values[sim_ann_iter] <- f_old
      }
      v_vectors[iter, ] <- v_vec
      c_vectors[iter, ] <- c_vec
      switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, 
                               byrow = T) * lambda[, v_vec]
      objective <- objective + f_old
      lambda_hat_new <- lambda_hat_new + switchedMatrix
      if ((iter%%printIter == 0) && (verbose == TRUE)) {
        cat(paste0("          mcmc draw = ", iter, ":  sum f = ", 
                   round(objective, 3)), "\r")
      }
    }
    f[totalIterations] <- objective
    if (totalIterations > 10) {
      f_diff <- f[totalIterations - 1] - f[totalIterations]
      if ((f_diff < 0) && (increaseIter)) {
        sim_an_iters <- sim_an_iters + sa_loops
        cat(paste0("  WARNING: objective function increased."), 
            "\n")
        cat(paste0("           Simulated Annealing loops increased to ", 
                   sim_an_iters, "."), "\n")
        f_diff <- threshold + 1
      }
      if (abs(f_diff) < threshold) {
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

# Papastamoulis and Ntzoufras partial simulated annealing approach
rsp_partial_sa <- function (lambda_mcmc, maxIter = 1000, threshold = 1e-06, verbose = TRUE, 
                            sa_loops, rotate = TRUE, increaseIter = FALSE, temperatureSchedule = NULL, 
                            printIter = 1000) {
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
    stop("Something is wrong, please check input.")
  }
  if (is.null(temperatureSchedule)) {
    temperatureSchedule <- function(sim_ann_iter) {
      1/log(sim_ann_iter + 1)
    }
  }
  objective_function <- function(lambda, lambda_hat, c_vec, 
                                 v_vec, q, p) {
    switchedMatrix <- matrix(c_vec, ncol = q, nrow = p, byrow = T) * 
      lambda[, v_vec]
    f <- sum((switchedMatrix - lambda_hat)^2)
    result <- vector("list", length = 2)
    result[[1]] <- f
    result[[2]] <- switchedMatrix
    names(result) <- c("value", "matrix")
    return(result)
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
  lambda_hat <- matrix(colMeans(lambda_mcmc_varimax), ncol = q, 
                       nrow = p, byrow = TRUE)
  lambda_hat_zero <- lambda_hat
  st <- 1:q
  c_vectors <- matrix(rep(1, q), ncol = q, nrow = mcmcIterations, 
                      byrow = TRUE)
  v_vectors <- matrix(1:q, ncol = q, nrow = mcmcIterations, 
                      byrow = TRUE)
  dim_all_c <- 2^q
  dim_all_v <- factorial(q)
  perm <- 1:q
  sim_an_iters <- sa_loops
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
      f_values <- numeric(sim_an_iters)
      c_vec <- c_vectors[iter, ]
      v_vec <- v_vectors[iter, ]
      obj <- objective_function(lambda = lambda, lambda_hat = lambda_hat, 
                                c_vec = c_vec, v_vec = v_vec, q = q, p = p)
      f_old <- obj$value
      f_values[1] <- f_old
      alpha <- 0.2
      temperature <- 1000
      accept_rate <- 1
      for (sim_ann_iter in 2:sim_an_iters) {
        temperature <- temperatureSchedule(sim_ann_iter)
        if (runif(1) < 0.9) {
          c_proposed <- c_vec
          rIndex <- sample(q, 1)
          c_proposed[rIndex] <- -c_vec[rIndex]
        }
        else {
          c_proposed <- sample(c(-1, 1), q, replace = TRUE)
        }
        lambda_switch <- matrix(c_proposed, ncol = q, 
                                nrow = p, byrow = T) * lambda_hat
        for (j in 1:q) {
          temp <- (lambda - lambda_switch[, j])^2
          cost.matrix[j, ] <- colSums(temp)
        }
        matr <- lp.assign(cost.matrix)$solution
        for (j in 1:q) {
          perm[j] <- st[matr[, j] > 0]
        }
        v_proposed <- order(perm)
        cost_proposed <- sum(cost.matrix * matr)
        obj <- cost_proposed
        f_new <- obj
        if (f_new < f_old) {
          c_vec <- c_proposed
          v_vec <- v_proposed
          f_old <- f_new
          accept_rate <- accept_rate + 1
        }
        else {
          u <- runif(1)
          ar <- -(f_new - f_old)/temperature
          if (log(u) < ar) {
            c_vec <- c_proposed
            v_vec <- v_proposed
            f_old <- f_new
            accept_rate <- accept_rate + 1
          }
        }
        f_values[sim_ann_iter] <- f_old
      }
      v_vectors[iter, ] <- v_vec
      c_vectors[iter, ] <- c_vec
      switchedMatrix <- switchedMatrix <- matrix(c_vec, 
                                                 ncol = q, nrow = p, byrow = T) * lambda[, v_vec]
      objective <- objective + f_old
      lambda_hat_new <- lambda_hat_new + switchedMatrix
      if ((iter%%printIter == 0) && (verbose == TRUE)) {
        cat(paste0("          mcmc draw = ", iter, ":  sum f = ", 
                   round(objective, 3)), "\r")
      }
    }
    f[totalIterations] <- objective
    if (totalIterations > 10) {
      f_diff <- f[totalIterations - 1] - f[totalIterations]
      if ((f_diff < 0) && (increaseIter)) {
        sim_an_iters <- sim_an_iters + sa_loops
        cat(paste0("  WARNING: objective function increased."), 
            "\n")
        cat(paste0("           Simulated Annealing loops increased to ", 
                   sim_an_iters, "."), "\n")
        f_diff <- threshold + 1
      }
      if (abs(f_diff) < threshold) {
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
