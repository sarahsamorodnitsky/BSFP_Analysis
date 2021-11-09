# -----------------------------------------------------------------------------
# Helper functions for Bayesian PMF
# -----------------------------------------------------------------------------

# Packages
library(doParallel)
library(foreach)
library(Matrix)
library(MASS)
library(truncnorm)

# -----------------------------------------------------------------------------
# Bayesian PMF functions
# -----------------------------------------------------------------------------

bpmf <- function(data, Y, nninit = TRUE, model_params, ranks = NULL, nsample, progress = TRUE) {
  # Gibbs sampling algorithm for sampling the underlying structure and the 
  # regression coefficient vector for a response vector. 
  # 
  # Arguments: 
  # data = matrix with list entries corresponding to each data source
  # Y = column vector with outcome or NULL
  # nninit = should the model be initialized with a nuclear norm penalized objective? if FALSE, provide ranks
  # model_params = (error_vars, joint_vars, indiv_vars, beta_vars = NULL, response_vars)
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  q <- nrow(data)
  p.vec <- apply(data, 1, function(source) nrow(source[[1]]))
  p <- sum(p.vec)
  n <- ncol(data[[1,1]])
  
  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------
  
  error_vars <- model_params$error_vars
  sigma2_joint <- joint_var <- model_params$joint_var
  sigma2_indiv <- indiv_vars <- model_params$indiv_vars
  beta_vars <- model_params$beta_vars
  response_vars <- model_params$response_vars; a <- response_vars[1]; b <- response_vars[2]
  
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
  
  response_given <- !is.null(Y) 
  
  # If so, what kind of response is it?
  if (response_given) {
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
    rank_init <- BIDIFAC(data, rmt = TRUE, pbar = FALSE)
    
    # Print when finished
    print("Posterior mode obtained: joint and individual ranks determined.")
    
    # Saving the results
    sigma.mat <- rank_init$sigma.mat
    C <- rank_init$C
    r <- rankMatrix(C[[1,1]]) # Joint rank
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
  }
  
  r_total <- r + sum(r.vec)
  n_beta <- 1 + r_total
  
  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    diag(Sigma_beta) <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
  }
  
  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  U.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  W.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = q))
  
  if (!response_given) {
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- matrix(list(), nrow = 1, ncol = 1)
  }
  
  if (!missingness_in_data) {
    Xm.draw <- matrix(list(), nrow = q, ncol = 1)
  }
  
  if (response_given) {
    beta.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    Z.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
    if (response_type == "continuous") {
      tau2.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    }
    
    if (missingness_in_response) {
      Ym.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    }
  }
  
  if (missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  # ---------------------------------------------------------------------------
  # Initialize V, U, V, W
  # ---------------------------------------------------------------------------
  
  V0 <- matrix(list(), nrow = 1, ncol = 1)
  V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
  
  U0 <- matrix(list(), nrow = q, ncol = 1)
  Vs0 <- matrix(list(), nrow = 1, ncol = q)
  W0 <- matrix(list(), nrow = q, ncol = q)
  
  for (s in 1:q) {
    U0[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p.vec[s], ncol = r)
    
    Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
    
    W0[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = p.vec[s], ncol = r.vec[s])
    
    for (ss in 1:q) {
      if (ss != s) {
        W0[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
      }
    }
  }
  
  if (response_given) {
    # Combining the scores together
    VStar0 <- cbind(1, do.call(cbind, V0), do.call(cbind, Vs0))
    
    beta0 <- matrix(mvrnorm(1, mu = c(rep(0, n_beta)), Sigma = Sigma_beta))
    Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
    tau20 <- 1/rgamma(1, shape = a, rate = b)
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
        Ym0 <- matrix(list(), nrow = 1, ncol = 1)
        Ym0[[1,1]] <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20)))[missing_obs_Y,]
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(list(), nrow = 1, ncol = 1)
        Ym0[[1,1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))[missing_obs_Y,]
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
  }
  
  if (missingness_in_data) {
    Xm.draw[[1]] <- Xm0
  }
  
  # ---------------------------------------------------------------------------
  # Computing the inverses that don't change from iteration to iteration
  # ---------------------------------------------------------------------------
  
  if (!response_given) {
    # Error variance for X. 
    SigmaInv <- matrix(list(), nrow = q, ncol = q)
    
    for (s in 1:q) {
      SigmaInv[[s,s]] <- diag(rep(1/error_vars[s], p.vec[s]))
      
      # Fill in the off-diagonal entries
      for (ss in 1:q) {
        if (ss != s) {
          SigmaInv[[s, ss]] <- matrix(0, nrow = p.vec[s], ncol = p.vec[ss])
        }
      }
    }
    
  }
  
  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined error variances between data and Z
      SigmaVInv <- diag(c(rep(1/error_vars, p.vec), 1))
      
      # For Vs
      SigmaVsInv <- matrix(list(), nrow = q, ncol = q)
      
      for (s in 1:q) {
        SigmaVsInv[[s,s]] <- diag(c(rep(1/error_vars[s], p.vec[s]), 1))
      }
    } 
    
    # For beta - Combined error variances between intercept and all betas
    SigmaBetaInv <- solve(Sigma_beta)
  }
  
  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------
  
  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/(nsample/100))
    
    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------
    
    V.iter <- V.draw[[iter]]
    U.iter <- U.draw[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.draw[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- t(beta.draw[iter,, drop = FALSE])
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept, joint, and individual effects
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE]
      beta_indiv.iter.temp <- beta.iter[(r+2):nrow(beta.iter),, drop = FALSE]
      
      for (s in 1:q) {
        if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
        if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
      }
      
      if (response_type == "binary") {
        Z.iter <- Z.draw[iter,]
      }
      
      if (response_type == "continuous") {
        tau2.iter <- tau2.draw[iter]
      }
      
      if (missingness_in_response) {
        # Save the current imputations for the missing values
        Ym.iter <- Ym.draw[[iter]]
        
        # Creating the completed outcome vector
        Y_complete <- Y
        
        # Filling in the missing entries for R1 and R2. 
        Y_complete[missing_obs_Y] <- Ym.iter
      }
      
      if (!missingness_in_response) {
        Y_complete <- Y
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
    # Posterior sample for V
    # -------------------------------------------------------------------------
    
    if (!response_given) {
      # Concatenating Ui's together
      U.iter.combined <- data.rearrange(U.iter)$out
      SigmaInvV <- data.rearrange(SigmaInv)$out
      tU_Sigma <- crossprod(U.iter.combined, SigmaInvV)
      
      # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
      tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
      
      # The combined centered Xis with the latent response vector
      X.iter <- data.rearrange(X_complete)$out - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t))
      Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
      
      V.draw[[iter+1]][[1]] <- t(sapply(1:n, function(i) {
        bv <-  tU_Sigma %*% X.iter[,i]
        
        Vj <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
        Vj
      }))
    }
    
    if (response_given) {
      # Concatenating Ui's together
      U.iter.combined <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))
      
      # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
      tU_Sigma <- crossprod(U.iter.combined, SigmaVInv) 
      tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter.combined)
      
      if (response_type == "binary") {
        Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
        V.draw[[iter+1]][[1]] <- t(sapply(1:n, function(i) {
          # The combined centered Xis with the latent response vector
          X.iter <- rbind(data.rearrange(X_complete)$out - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                          (Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter))[i,])
          
          bv <- tU_Sigma %*% X.iter[,i]
          
          Vj <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
          Vj
        }))
      }
      
      if (response_type == "continuous") {
        Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
        V.draw[[iter+1]][[1]] <- t(sapply(1:n, function(i) {
          # The combined centered Xis with the latent response vector
          X.iter <- rbind(data.rearrange(X_complete)$out - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                          (Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% do.call(rbind, beta_indiv.iter))[i,])
          
          bv <- tU_Sigma %*% X.iter[,i]
          
          Vj <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
          Vj
        }))
      }
    }
    
    # Updating the value of V
    V.iter <- V.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for Us
    # -------------------------------------------------------------------------
    
    for (s in 1:q) {
      Xs.iter <- X_complete[[s,1]] - W.iter[[s,s]] %*% t(Vs.iter[[1,s]])
      Bu <- solve((1/error_vars[s]) * t(V.iter[[1,1]]) %*% V.iter[[1,1]] + (1/sigma2_joint) * diag(r))
      U.draw[[iter+1]][[s,1]] <- t(sapply(1:p.vec[s], function(j) {
        bu <- (1/error_vars[s]) * t(V.iter[[1,1]]) %*% Xs.iter[j, ]
        
        U1j <- mvrnorm(1, mu = Bu %*% bu, Sigma = Bu)
        U1j
      }))
    }

    U.iter <- U.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for Vs, s=1,...,q
    # -------------------------------------------------------------------------
    
    if (!response_given) {
      for (s in 1:q) {
        Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
        Bvs <- solve((1/error_vars[s]) * t(W.iter[[s,s]]) %*% W.iter[[s,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
        
        Vs.draw[[iter+1]][[1,s]] <- t(sapply(1:n, function(i) {
          bvs <- (1/error_vars[s]) * t(W.iter[[s,s]]) %*% Xs.iter[, i]
          
          V1j <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
          V1j
        }))
      }
    }
    
    if (response_given) {
      for (s in 1:q) {
        # Combined Ws and beta
        W.iter.combined <- rbind(W.iter[[s,s]], t(beta_indiv.iter[[s,1]]))
        
        tW_Sigma <- crossprod(W.iter.combined, SigmaVsInv[[s,s]])
        tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter.combined)
        
        if (response_type == "binary") {
          # Combined centered Xs and Z
          Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                          t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                              do.call(cbind, Vs.iter[1, !(1:ncol(Vs.iter) %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(c(1:q) %in% s), 1])))
          
          Bvs <- solve(tW_Sigma_W + (1/error_vars[s]) * diag(r.vec[s]))
          Vs.draw[[iter+1]][[1,s]] <- t(sapply(1:n, function(i) {
            bvs <- tW_Sigma %*% Xs.iter[, i]
            
            Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
            Vsi
          }))
        }
        
        if (response_type == "continuous") {
          # Combined centered Xs and Y
          Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                           t(Y_complete - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                               do.call(cbind, Vs.iter[1, !(1:ncol(Vs.iter) %in% s)]) %*% do.call(rbind, beta_indiv.iter[!(c(1:q) %in% s), 1])))
          
          Bvs <- solve(tW_Sigma_W + (1/error_vars[s]) * diag(r.vec[s]))
          Vs.draw[[iter+1]][[1,s]] <- t(sapply(1:n, function(i) {
            bvs <- tW_Sigma %*% Xs.iter[, i]
            
            Vsi <- mvrnorm(1, mu = Bvs %*% bvs, Sigma = Bvs)
            Vsi
          }))
        }
      }
    }
    
    # Update the current value of V
    Vs.iter <- Vs.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for W
    # -------------------------------------------------------------------------
    
    for (s in 1:q) {
      Xs.iter <- X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]])
      Bws <- solve((1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Vs.iter[[1,s]] + (1/indiv_vars[s]) * diag(r.vec[s]))
      
      W.draw[[iter+1]][[s,s]] <- t(sapply(1:p.vec[s], function(j) {
        bws <- (1/error_vars[s]) * t(Vs.iter[[1,s]]) %*% Xs.iter[j,] 
        
        Wsj <- mvrnorm(1, mu = Bws %*% bws, Sigma = Bws)
        Wsj
      }))
      
      for (ss in 1:q) {
        if (ss != s) W.draw[[iter+1]][[s, ss]] <- matrix(0, nrow = p.vec[s], ncol = r.vec[ss])
      }
    }
    
    # Update the current value of W
    W.iter <- W.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      # Combine current values of V and V_\cdot
      VStar.iter <- cbind(1, do.call(cbind, V.iter), do.call(cbind, Vs.iter))
      
      if (response_type == "continuous") {
        tau2.draw[iter+1] <- 1/rgamma(1, shape = a + (n/2), rate = b + 0.5 * sum((Y_complete - VStar.iter %*% beta.iter)^2))
        
        # Update the current value of tau2
        tau2.iter <- tau2.draw[iter+1]
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Bbeta <- solve(t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- t(VStar.iter) %*% Z.iter
        beta.draw[iter+1,] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      }
      
      if (response_type == "continuous") {
        Bbeta <- solve((1/tau2.iter) * t(VStar.iter) %*% VStar.iter + SigmaBetaInv)
        bbeta <- (1/tau2.iter) * t(VStar.iter) %*% Y_complete
        beta.draw[iter+1,] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      }
      
      # Update the current value of beta
      beta.iter <- t(beta.draw[iter+1,, drop = FALSE])
      
      # Creating a matrix of the joint and individual effects
      beta_indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept, joint, and individual effects
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE]
      beta_indiv.iter.temp <- beta.iter[(r+2):nrow(beta.iter),, drop = FALSE]
      
      for (s in 1:q) {
        if (s == 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
        if (s != 1) beta_indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Z.draw[iter+1,] <- matrix(sapply(1:n, function(i) {
          if (Y_complete[i,] == 1) z <- rtruncnorm(1, a = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          if (Y_complete[i,] == 0) z <- rtruncnorm(1, b = 0, mean = (VStar.iter %*% beta.iter)[i,], sd = 1)
          z
        }), ncol = 1)
      }
    }
    
    # -------------------------------------------------------------------------
    # Impute missing data
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (missingness_in_response) {
        if (response_type == "continuous") {
          Ym.draw[[iter+1]] <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter)), ncol = 1)[missing_obs_Y,]
        }
        
        if (response_type == "binary") {
          Ym.draw[[iter+1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)[missing_obs_Y,]
        }
      }
    }
    
    if (missingness_in_data) {
      for (s in 1:q) {
        Es <-  matrix(rnorm(p.vec[s]*n, 0, sqrt(error_vars[s])), nrow = p.vec[s], ncol = n)
        Xm.draw[[iter+1]][[s,1]] <- (U.iter[[s,1]] %*% t(V.iter[[1,1]]) + W.iter[[s,s]] %*% t(Vs.iter[[1,s]]) + Es)[missing_obs[[s]]]
      }
    }
  }
  
  # Return
  list(data = data, # Returning the scaled version of the data
        Y = Y, # Return the response vector
        sigma.mat = sigma.mat, # Returning the scaling factors
        V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw,
        Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw,
        tau2.draw = tau2.draw, beta.draw = beta.draw)

}

bpmf_sim <- function(nsample, p.vec, n, true_params, model_params, nsim = 1000, s2n, ranks, 
                     response = NULL, missingness = NULL, prop_missing = NULL, entrywise = NULL, nninit = TRUE) {
  
  # ---------------------------------------------------------------------------
  # Change the below to be fed in by the parameters and hyperparameters arguments above 
  # ---------------------------------------------------------------------------
  
  cl <- makeCluster(n_clust)
  registerDoParallel(cl)
  sim_results <- foreach (sim_iter = 1:nsim, .packages = c("Matrix", "MASS")) %dopar% {
    
    # -------------------------------------------------------------------------
    # Generating the data
    # -------------------------------------------------------------------------
  
    sim_data <- bpmf_data(p.vec, n, ranks, true_params, s2n, response, missingness, entrywise, prop_missing)
    
    # Saving the data
    data <- sim_data$data
    structure <- list(joint = sim_data$joint.structure.scale,
                      indiv = sim_data$indiv.structure.scale)
    
    missing_data <- sim_data$missing_data
    missing_obs <- sim_data$missing_obs
    
    # The response
    Y <- sim_data$Y
    Y_missing <- sim_data$Y_missing
    missing_obs_Y <- sim_data$missing_obs_Y
    
    # The standardizing coefficients
    s2n_coef <- sim_data$s2n_coef
    
    # The response parameters
    beta <- sim_data$beta
    tau2 <- sim_data$tau2
    
    # -------------------------------------------------------------------------
    # Center the data
    # -------------------------------------------------------------------------
    
    center_the_data <- center_data(data = data, structure = structure)
    
    # Saving the centered data
    data_centered <- center_the_data$data_centered
    means_for_centering <- center_the_data$means_for_centering
    structure_centered <- center_the_data$structure_centered
    
    # -------------------------------------------------------------------------
    # Running the Gibbs sampling algorithm
    # -------------------------------------------------------------------------
    
    # Gibbs sampling
    res <- bpmf(data_centered, Y, nninit = nninit, model_params, ranks = NULL, nsample, progress = TRUE)
    
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
    
    # Computing the underlying structure for each Xs after burn-in from the sampler
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
    # Checking the coverage
    # -------------------------------------------------------------------------
    
    # Saving the draws from the sampler together
    draws <- list(joint.structure.burnin = joint.structure.burnin,
                  indiv.structure.burnin = indiv.structure.burnin,
                  beta.burnin = beta.burnin,
                  tau2.burnin = tau2.burnin,
                  Xm.burnin = Xm.burnin,
                  Ym.burnin = Ym.burnin)
    
    # Saving the truth together for comparison
    truth <- list(joint.structure.scale = joint.structure.scale,
                  indiv.structure.scale = indiv.structure.scale,
                  beta = beta,
                  tau2 = tau2,
                  Xm = missing_data,
                  Ym = Y_missing)
    
    # Adding up the coverage, MSE, and CI width from this iteration
    sim_iter_results <- get_results(truth, draws, burnin)
    
    # ---------------------------------------------------------------------------
    # Returning the results at the end of the loop
    # ---------------------------------------------------------------------------
    
    # Add the indices for the missing values in this sim_iter
    sim_iter_results$any_missing <- list(missing_obs = missing_obs,
                                         missing_obs_Y = missing_obs_Y)
    
    # Return 
    sim_iter_results
  }
  stopCluster(cl)
  
  # ---------------------------------------------------------------------------
  # Averaging the results
  # ---------------------------------------------------------------------------
  
  # Calculating how many times each entry in the data was observed
  observed_counts <- calculate_observed(sim_results, p.vec, n, response)
  
  # Taking the averages of all the results
  sim_results_avg <- average_results(sim_results, observed_counts, nsim)
  
  # ---------------------------------------------------------------------------
  # Returning the results
  # ---------------------------------------------------------------------------
  
  sim_results_avg
  
}

# -----------------------------------------------------------------------------
# Helper functions for initializing with BIDIFAC
# -----------------------------------------------------------------------------

# functions
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

BIDsim=function(m.vec, n.vec, 
                rkG=1, rkC=1, rkR=1, rkI=1, SNR=1){
  p=length(m.vec); q=length(n.vec)
  
  m.vec.end=cumsum(m.vec)
  m.vec.start=c(1,m.vec.end[-p]+1)
  n.vec.end=cumsum(n.vec)
  n.vec.start=c(1,n.vec.end[-q]+1)
  
  X=S=G=R=C=I=matrix(list(), p,q)
  
  if (length(rkR)==1){rkR=rep(rkR, p)}
  if (length(rkC)==1){rkC=rep(rkC, q)}
  if (length(rkI)==1){rkI=matrix(rkI, p,q)}
  vec.rk=c(rkG, rkR, rkC, rkI)
  sum.rk=sum(vec.rk)
  
  rk.ind.end=cumsum(vec.rk)
  rk.ind.start=c(1,rk.ind.end[-((p+1)*(q+1))]+1)
  
  sigma.mat=1/sqrt(tcrossprod(m.vec,n.vec))/SNR
  
  m=sum(m.vec); n=sum(n.vec)
  svdX=svd(matrix(rnorm(m*n),m,n))
  min.mn=min(m,n)
  
  if (sum.rk>min.mn) stop("the rank of the signal is too high")
  
  ind=sample(min.mn, sum.rk)
  svdX.u=svdX$u[,ind]
  svdX.d=svdX$d[1:sum.rk]/sqrt(frob(svdX$d[1:sum.rk]))
  svdX.v=svdX$v[,ind]
  
  #generate G
  ind.G=rk.ind.start[1]:rk.ind.end[1]
  for (i in 1:p){
    U.g=svdX.u[m.vec.start[i]:m.vec.end[i],ind.G]
    for (j in 1:q){
      D.g=diag2(sample2(svdX$d[ind.G]))
      V.g=svdX.v[n.vec.start[j]:n.vec.end[j], ind.G]
      G[[i,j]]=tcrossprod(U.g, tcrossprod(V.g,D.g))
    }
  }
  rk.ind.start=rk.ind.start[-1]
  rk.ind.end=rk.ind.end[-1]
  
  #generate R
  for (i in 1:p){
    ind.R=rk.ind.start[i]:rk.ind.end[i]
    U.r=svdX.u[m.vec.start[i]:m.vec.end[i],ind.R]
    for (j in 1:q){
      D.r=diag2(sample2(svdX$d[ind.R]))
      V.r=svdX.v[n.vec.start[j]:n.vec.end[j], ind.R]
      R[[i,j]]=tcrossprod(U.r, tcrossprod(V.r,D.r))
    }
  }
  rk.ind.start=rk.ind.start[-(1:p)]
  rk.ind.end=rk.ind.end[-(1:p)]
  
  #generate C
  for (j in 1:q){
    ind.C=rk.ind.start[j]:rk.ind.end[j]
    V.c=svdX.v[n.vec.start[j]:n.vec.end[j], ind.C]
    for (i in 1:p){
      D.c=diag2(sample2(svdX$d[ind.C]))
      U.c=svdX.u[m.vec.start[i]:m.vec.end[i],ind.C]
      C[[i,j]]=tcrossprod(U.c, tcrossprod(V.c,D.c))
    }
  }
  rk.ind.start=matrix(rk.ind.start[-(1:q)],p)
  rk.ind.end=matrix(rk.ind.end[-(1:q)],p)
  
  #generate I
  for (i in 1:p){
    for (j in 1:q){
      ind.I=rk.ind.start[i,j]:rk.ind.end[i,j]
      U.i=svdX.u[m.vec.start[i]:m.vec.end[i],ind.I]
      D.i=diag2(sample2(svdX$d[ind.I]))
      V.i=svdX.v[n.vec.start[j]:n.vec.end[j], ind.I]
      I[[i,j]]=tcrossprod(U.i, tcrossprod(V.i,D.i))
    }
  }
  
  for (i in 1:p){
    for (j in 1:q){
      d=svd(G[[i,j]]+R[[i,j]]+C[[i,j]]+I[[i,j]])
      G[[i,j]]=G[[i,j]]/sqrt(frob(d$d))
      R[[i,j]]=R[[i,j]]/sqrt(frob(d$d))
      C[[i,j]]=C[[i,j]]/sqrt(frob(d$d))
      I[[i,j]]=I[[i,j]]/sqrt(frob(d$d))
      
      d$d=d$d/sqrt(frob(d$d))
      S[[i,j]]=d$u%*%diag(d$d)%*%t(d$v)
      X[[i,j]]=S[[i,j]]+replicate(n.vec[j],rnorm(m.vec[i],0,sigma.mat[i,j] ))
    }
  }
  
  return(list(X=X,S=S,G=G,R=R,C=C,I=I,sigma=sigma.mat))
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
                 eps=1e-3, max.iter=1000, pbar=TRUE, seed=NULL, ...){
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
    G00=replicate(sum(nvec),rnorm(sum(mvec)))
    R00=replicate(sum(nvec),rnorm(sum(mvec)))
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
    
    #Update G
    fit1=softSVD(X00-R00-C00-I00,lambda.G)
    G00=fit1$out; G00.nuc=fit1$nuc
    
    #update R
    for (i in 1:p){
      ind=start.ind.m[i]:end.ind.m[i]
      fit1=softSVD(X00[ind,]-G00[ind,]-C00[ind,]-I00[ind,], lambda.R[i])
      R00[ind,]=fit1$out; R00.nuc[i]=fit1$nuc
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
  
  S00.mat=G00.mat=R00.mat=C00.mat=I00.mat=data
  for (i in 1:p){
    ind1= start.ind.m[i]:end.ind.m[i]
    for (j in 1:q){
      ind2=start.ind.n[j]:end.ind.n[j]
      G00.mat[[i,j]]=G00[ind1,ind2]*sigma.mat[i,j]
      R00.mat[[i,j]]=R00[ind1,ind2]*sigma.mat[i,j]
      C00.mat[[i,j]]=C00[ind1,ind2]*sigma.mat[i,j]
      I00.mat[[i,j]]=I00[ind1,ind2]*sigma.mat[i,j]
      S00.mat[[i,j]]=G00.mat[[i,j]]+R00.mat[[i,j]]+C00.mat[[i,j]]+I00.mat[[i,j]]
    }
  }
  
  return(list(X=data, S=S00.mat,
              G=G00.mat, R=R00.mat, C=C00.mat, I=I00.mat,
              sigma.mat=sigma.mat, n.vec=nvec,m.vec=mvec))
}

# -----------------------------------------------------------------------------
# Helper functions for simulations
# -----------------------------------------------------------------------------

# Generate fake data depending on conditions
bpmf_data <- function(p.vec, n, ranks, true_params, s2n, response, missingness, entrywise, prop_missing) {
  # Generates fake data depending on the dims provided in p.vec, n, and ranks
  # and the true parameters provided in `true_params`
  # s2n = desired signal-to-noise ratio 
  
  # -------------------------------------------------------------------------
  # Setting the dimensions and latent components
  # -------------------------------------------------------------------------
  
  q <- length(p.vec)
  n_beta <- 1 + sum(ranks)
  
  # Setting the hyperparameters
  error_vars <- true_params$error_vars
  sigma2_joint <- joint_var <- true_params$joint_var
  sigma2_indiv <- indiv_vars <- true_params$indiv_vars
  beta_vars <- true_params$beta_vars
  response_vars <- true_params$response_vars; a <- response_vars[1]; b <- response_vars[2]
  
  # -------------------------------------------------------------------------
  # Generating the underlying structure
  # -------------------------------------------------------------------------
  
  joint.structure <- indiv.structure <- matrix(list(), ncol = 1, nrow = q)
  
  V <- matrix(list(), nrow = 1, ncol = 1)
  V[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
  
  U <- matrix(list(), nrow = q, ncol = 1)
  Vs <- matrix(list(), nrow = 1, ncol = q)
  W <- matrix(list(), nrow = q, ncol = q)
  
  E <- matrix(list(), nrow = q, ncol = 1)
  
  for (s in 1:q) {
    U[[s,1]] <- matrix(rnorm(p.vec[s]*r, mean = 0, sd = sqrt(10)), nrow = p.vec[s], ncol = r)
    
    Vs[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(5)), nrow = n, ncol = r.vec[s])
    
    W[[s,s]] <- matrix(rnorm(p.vec[s]*r.vec[s], mean = 0, sd = sqrt(5)), nrow = p.vec[s], ncol = r.vec[s])
    
    E[[s,1]] <- matrix(rnorm(p.vec[s]*n), nrow = p.vec[s], ncol = n)
    
    for (ss in 1:q) {
      if (ss != s) {
        W[[s,ss]] <- matrix(0, nrow = p.vec[[s]], ncol = r.vec[ss])
      }
    }
    
    joint.structure[[s,1]] <- U[[s,1]] %*% t(V[[1,1]])
    indiv.structure[[s,1]] <- W[[s,s]] %*% t(Vs[[1,s]])
  }
  
  # -------------------------------------------------------------------------
  # Standardizing the variance of the signal
  # -------------------------------------------------------------------------
  
  # Calculating the scaling coefficient so that the variance of the underlying structure = s2n * noise variance
  s2n_coef <- rep(0, q)
  joint.structure.scale <- joint.structure
  indiv.structure.scale <- indiv.structure
  
  for (s in 1:q) {
    s2n_coef[s] <- s2n * sd(c(E[[s,1]]))/sd(c(joint.structure[[s,1]] + indiv.structure[[s,1]]))
    
    joint.structure.scale[[s,1]] <- s2n_coef[s] * joint.structure[[s,1]]
    indiv.structure.scale[[s,1]] <- s2n_coef[s] * indiv.structure[[s,1]]
  }
  
  # -------------------------------------------------------------------------
  # Calculating the observed data
  # -------------------------------------------------------------------------
  
  data <- matrix(list(), nrow = q, ncol = 1)
  
  for (s in 1:q) {
    data[s,1][[1]] <- joint.structure.scale[[s,1]] + indiv.structure[[s,1]] + E[[s,1]]
  }
  
  # -------------------------------------------------------------------------
  # Adding a response if desired
  # -------------------------------------------------------------------------
  
  if (is.null(response)) {
    Y <- NULL
    Y_missing <- NULL
    beta <- NULL
    tau2 <- NULL
  }
  
  if (!is.null(response)) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    diag(Sigma_beta) <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
    
    # Generate betas
    beta <- matrix(mvrnorm(1, mu = rep(0, n_beta), Sigma = Sigma_beta), ncol = 1)
    
    # Combine the Vs
    VStar <- cbind(1, do.call(cbind, V), do.call(cbind, Vs))
    
    # True probability of being a case
    Prob.VStar.beta <- pnorm(VStar %*% beta)
    
    if (response == "binary") {
      Y <- matrix(rbinom(n, size = 1, prob = pnorm(VStar %*% beta)), ncol = 1)
      tau2 <- NULL
    }
    
    if (response == "continuous") {
      tau2 <- 1/rgamma(1, shape = a, rate = b)
      Y <- VStar %*% beta + rnorm(n, mean = 0, sd = sqrt(tau2))
    }
  }
  
  # -------------------------------------------------------------------------
  # Adding missingness if desired
  # -------------------------------------------------------------------------
  
  if (is.null(missingness)) {
    missing_obs_X1 <- NULL
    missing_obs_X2 <- NULL
    
    X1_missing <- NULL
    X2_missing <- NULL
    
    Y_missing <- NULL
    missing_obs_Y <- NULL
    
  }
  
  if (!is.null(missingness)) {
    if (missingness != "missingness_in_data" | missingness != "both") {
      missing_obs <- lapply(1:q, function(s) NULL)
      missing_data <- matrix(list(), nrow = q, ncol = 1);
    }
    
    if (missingness != "missingness_in_response" | missingness != "both") {
      Y_missing <- NULL
      missing_obs_Y <- NULL
    }
    
    if (missingness == "missingness_in_response" | missingness == "both") {
      missing_obs_Y <- sample(1:n, size = prop_missing * n, replace = FALSE)
      Y_missing <- Y
      Y_missing[missing_obs_Y] <- NA
    }
    
    if (missingness == "missingness_in_data" | missingness == "both") {
      
      missing_obs <- lapply(1:q, function(s) list())
      missing_data <- matrix(list(), nrow = q, ncol = 1)
      
      if (entrywise) { # if removing observations entrywise
        for (s in 1:q) {
          # these are counters going down the columns of R. So 9 would be the 9th entry counting down. 
          missing_obs[[s]] <- sample(x = 1:length(data[[s,1]]), size = prop_missing*length(data[[s,1]]), replace = FALSE) 

          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][missing_obs[[s]]] <- NA
        }
      } else { # if removing entire columns
        # Gives the column indices to remove
        
        for (s in 1:q) {
          # These are counters going down the columns of X. So 9 would be the 9th entry counting down. 
          missing_obs[[s]] <- sample(x=1:n, size = n*prop_missing, replace = FALSE)
          
          if (s != 1) {
            avail_obs <- c(1:n)[!(c(1:n) %in% unlist(missing_obs[1:(s-1)]))]
            missing_obs[[s]] <- sample(x=avail_obs, size = n*prop_missing, replace = FALSE)
          }
      
          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][,missing_obs[[s]]] <- NA
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
       s2n = s2n, s2n_coef = s2n_coef, # Scaling for s2n
       joint.structure.scale = joint.structure.scale, # Joint structure
       indiv.structure.scale = indiv.structure.scale, # Individual structure
       beta = beta, tau2 = tau2)
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
  # based on res.
  # truth = matrix with the true values
  # draws = list of matrix results from a Gibbs sampler after burn-in that the 
  #         credible intervals will be based on.
  
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
  # true component being estimated. 
  Reduce('+', lapply(draws, function(iter.val) {
    (iter.val - truth)^2
  }))/length(draws)
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
get_results <- function(truth, draws, burnin) {
  # Save the number of model parameters to check coverage for
  n_param <- length(draws)
  
  # Save the results 
  sim_results <- lapply(1:n_param, function(i) list())
  names(sim_results) <- names(draws)
  
  # Iterate through the parameters, checking the coverage, MSE, and CI width
  for (param in 1:n_param) {
    # Check the dimension of the current parameter
    dim_param <- dim(truth[[param]])
    q <- dim_param[1]
    
    # If there are results available
    if (!is.null(truth[[param]][[1,1]])) {
      
      # If the results are a matrix then they correspond to multiple datasets
      if (!is.null(dim_param)) {
        sim_results[[param]] <- matrix(list(), nrow = q, ncol = 1)
        
        for (s in 1:q) {
          current_draws <- lapply(1:(burnin+1), function(iter) draws[[param]][[iter]][[s,1]])
          sim_results[[param]][[s,1]] <- list(check_coverage(truth[[param]][[s,1]], current_draws, burnin = burnin),
                                              mse(truth[[param]][[s,1]], current_draws),
                                              ci_width(current_draws, burnin = burnin))
        }
      }
      
      # if the results do not correspond to a matrix
      if (is.null(dim_param)) {
        sim_results[[param]] <- list(check_coverage(truth[[param]], draws[[param]], burnin),
                                     mse(truth[[param]], draws[[param]]),
                                     ci_width(draws[[param]], burnin = burnin))
      }
    }
  }
  
  # Return the results
  sim_results
}

# Count the number of times each observation in each dataset was observed
calculate_observed <- function(sim_results, p.vec, n, response) {
  # sim_results (list) = list of results from the simulation

  # Count how many datasources there were (add 1 for the response)
  q <- nrow(sim_results[[1]]$joint.structure.burnin)
  
  # Create a list of indices for observation
  obs_inds <- lapply(1:q, function(i) 1:(p.vec[i] * n))
  
  # Create a list to store the results
  observed_counts <- lapply(1:q, function(i) rep(0, length(obs_inds[[i]])))
  
  # If there is a response
  if (response) {
    obs_inds[[q+1]] <- 1:n
    observed_counts[[q+1]] <- rep(0, length(obs_inds[[q+1]]))
  }
  
  for (sim_iter in 1:length(sim_results)) {
    for (s in 1:q) {
      current_missing_obs <- sim_results[[sim_iter]]$any_missing$missing_obs[[s]]
      observed_counts[[s]] <- observed_counts[[s]] + !(obs_inds[[s]] %in% current_missing_obs)
    }
    
    if (response) {
      current_missing_obs <- sim_results[[sim_iter]]$any_missing$missing_obs_Y
      observed_counts[[q+1]] <- observed_counts[[q+1]] + !(obs_inds[[q+1]] %in% current_missing_obs)
    }
  }

  # Return
  observed_counts
}

# Average the results before returning from the function
average_results <- function(sim_results, observed_counts, nsim) {
  # results_list (list) = list of results to return that should be averaged
  # observed_counts (list) = list of how many times each observation in each dataset was missing
  # nsim (int) = the number of simulation iterations that were run
  
  # Save the names of each parameter
  params <- names(sim_results[[1]])
  
  # Initialize a list to contain the averages
  sim_results_avg <- lapply(1:length(params), function(i) list())
  names(sim_results_avg) <- params
  
  for (param in params) {
    if (param != "any_missing") {
      # Save the results across the simulation replications for this parameter
      param_by_rep <- lapply(sim_results, function(sim_iter) sim_iter[[param]])
      
      # Are the results for a matrix or for a vector?
      matrix_or_vec <- dim(param_by_rep[[1]])
      
      # If matrix
      if (!is.null(matrix_or_vec)) {
        # Save the number of sources
        q <- matrix_or_vec[1]
        
        # Initialize a list for the results
        results_for_param <- matrix(list(), nrow = q, ncol = 1)
        
        # Iterate through the results for each source
        for (s in 1:q) {
          # Save the results for source s
          res_by_source <- lapply(param_by_rep, function(sim_iter) sim_iter[[s,1]])
          
          # Calculate the average coverage
          avg_coverage_source <- Reduce("+", lapply(res_by_source, function(sim_iter) sim_iter[[1]]))/observed_counts[[s]]
          
          # Calculate the average MSE
          avg_mse_source <- Reduce("+", lapply(res_by_source, function(sim_iter) sim_iter[[2]]))/observed_counts[[s]]
          
          # Calculate the average CI width
          avg_ci_width_source <- Reduce("+", lapply(res_by_source, function(sim_iter) sim_iter[[3]]))/observed_counts[[s]]
          
          # Save the results
          results_for_param[[s,1]] <- list(avg_coverage = mean(avg_coverage_source),
                                           avg_mse = mean(avg_mse_source),
                                           avg_ci_width = mean(avg_ci_width_source))
        }
        
        # Save the overall results
        sim_results_avg[[param]] <- results_for_param
        
      }
      
      # If not a matrix
      if (is.null(matrix_or_vec)) {
        # Save the results for this parameter
        sim_results_avg[[param]] <- mean(Reduce("+", param_by_rep)/nsim)
      }
    }
  }

  # Return
  sim_results_avg
}