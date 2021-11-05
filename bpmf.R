# -----------------------------------------------------------------------------
# Helper functions for Bayesian PMF
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Bayesian PMF functions
# -----------------------------------------------------------------------------

bpmf <- function(X1, X2, Y = NULL, nuclear_norm_init = TRUE, dims, hyperparameters, nsample = 5000, progress = TRUE) {
  # Gibbs sampling algorithm for sampling the underlying structure and the 
  # regression coefficient vector for a response vector. 
  # 
  # Arguments: 
  # X1, X2 = data matrices 
  # Y = observed response vector. If none provided, set as NULL. Otherwise
  #     binary or normally distributed. 
  # if (nuclear_norm_init) dims = list(p1, p2, n)
  # if (!nuclear_norm_init) dims = list(p1, p2, n, r, r1, r2)
  # dims = all important dimensions: pi x n is the 
  #        dimension of Ri, i=1,2. r is the rank of the joint structure. 
  #        ri is rank of the individual structure of Ri, i=1,2. 
  # hyperparameters = list(sigma21, sigma22, sigma2_joint, sigma2_indiv) list of 
  #             all the variances and other hyperparameters. 
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  p1 <- dims$p1 # number rows X1
  p2 <- dims$p2 # number rows X2
  n <- dims$n # number cols X1 and X2
  
  # ---------------------------------------------------------------------------
  # Extracting the variances and hyperparameters
  # ---------------------------------------------------------------------------
  
  sigma21 <- hyperparameters$sigma21 # error variance for X1
  sigma22 <- hyperparameters$sigma22 # error variance for X2
  sigma2_joint <- hyperparameters$sigma2_joint # variance for joint structure
  sigma2_indiv1 <- hyperparameters$sigma2_indiv1 # variance for individual structure of X1
  sigma2_indiv2 <- hyperparameters$sigma2_indiv2 # variance for individual structure of X2
  Sigma_beta <- hyperparameters$Sigma_beta # variance for coefficients
  shape <- hyperparameters$shape
  rate <- hyperparameters$rate
  
  # ---------------------------------------------------------------------------
  # Is there any missingness in X1 or X2?
  # ---------------------------------------------------------------------------
  
  missingness_in_data <- any(is.na(X1)) | any(is.na(X2))
  
  # If so, which entries are missing?
  if (missingness_in_data) {
    missing_obs_data <- list(X1 = which(is.na(X1)),
                             X2 = which(is.na(X2)))
  }
  
  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------
  
  response_given <- !is.null(Y) 
  
  # If so, what kind of response is it?
  if (response_given) {
    response_type <- if (all(unique(Y) %in% c(0,1,NA))) "binary" else "continuous"
  }
  
  # If there is a response, is there missingness in the outcome?
  if (response_given) {
    missingness_in_response <- any(is.na(Y))
    
    # If there is missingness, which entries are missing?
    if (missingness_in_response) {
      missing_obs_Y <- which(is.na(Y))
    }
  }
  
  # ---------------------------------------------------------------------------
  # Obtaining the ranks 
  # ---------------------------------------------------------------------------
  
  if (nuclear_norm_init) {
    # Preparing the data 
    data <- matrix(list(), nrow = 2, ncol = 1)
    data[1,1][[1]] <- X1
    data[2,1][[1]] <- X2
    nuclear_norm_res <- BIDIFAC(data, rmt = TRUE, pbar = FALSE)
    
    # Print when finished
    print("posterior mode obtained")
    
    # Saving the results
    sigma.mat <- nuclear_norm_res$sigma.mat
    C <- nuclear_norm_res$C
    r <- rankMatrix(C[1,1][[1]]/sigma.mat[1,1]) 
    I <- nuclear_norm_res$I
    r1 <- rankMatrix(I[1,1][[1]]/sigma.mat[1,1])
    r2 <- rankMatrix(I[2,1][[1]]/sigma.mat[2,1])
    
    dims$r <- r
    dims$r1 <- r1
    dims$r2 <- r2
    
    # Scaling X1 and X2 appropriately
    X1 <- X1/sigma.mat[1,1]
    X2 <- X2/sigma.mat[2,1]
  }
  
  if (!nuclear_norm_init) {
    r <- dims$r # number latent components of J1, J2
    r1 <- dims$r1 # number latent components of individual structure A1
    r2 <- dims$r2 # number of latent components of individual structure A2
    n_beta <- 1 + r + r1 + r2 # total number of coefficients
  }
  
  # ---------------------------------------------------------------------------
  # Initialize V, Ui, Vi, Wi, i=1,2.
  # ---------------------------------------------------------------------------
  
  V0 <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
  
  U10 <- matrix(rnorm(p1*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p1, ncol = r)
  U20 <- matrix(rnorm(p2*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = p2, ncol = r)
  
  V10 <- matrix(rnorm(n*r1, mean = 0, sd = sqrt(sigma2_indiv1)), nrow = n, ncol = r1)
  V20 <- matrix(rnorm(n*r2, mean = 0, sd = sqrt(sigma2_indiv2)), nrow = n, ncol = r2)
  
  W10 <- matrix(rnorm(p1*r1, mean = 0, sd = sqrt(sigma2_indiv1)), nrow = p1, ncol = r1)
  W20 <- matrix(rnorm(p2*r2, mean = 0, sd = sqrt(sigma2_indiv2)), nrow = p2, ncol = r2)
  
  VStar0 <- cbind(1, V0, V10, V20)
  
  if (response_given) {
    # Sigma_beta <- diag(c(100, rep(lambda2_joint, r), rep(lambda2_indiv1, r1), rep(lambda2_indiv2, r2)))
    beta0 <- matrix(mvrnorm(1, mu = c(rep(0, n_beta)), Sigma = Sigma_beta))
    
    # If the response is binary, initialize the latent continuous variable
    if (response_type == "binary") {
      Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = 1))
    }
    
    if (response_type == "continuous") {
      tau20 <- 1/rgamma(1, shape = shape, rate = rate)
    }
  }
  
  # If there is missingness in X1 or X2, generate starting values for the missing entries
  if (missingness_in_data) {
    Xm0 <- list(Xm10 = list(U10 %*% t(V0) + W10 %*% t(V10)),
                Xm20 = list(U20 %*% t(V0) + W20 %*% t(V20)))
  }
  
  # If there is missingness in Y, generate starting values for the missing entries
  if (response_given) {
    if (missingness_in_response) {
      # Combine initial values
      VStar0 <- cbind(1, V0, V10, V20)
      
      if (response_type == "continuous") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rnorm(n, mean = VStar0 %*% beta0, sd = sqrt(tau20)))
      }
      
      if (response_type == "binary") {
        # Generate starting values for the missing data
        Ym0 <- matrix(rbinom(n, size = 1, prob = pnorm(VStar0 %*% beta0)))
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------
  
  V.draw <- lapply(1:nsample, function(i) matrix(nrow = n, ncol = r))
  U1.draw <- lapply(1:nsample, function(i) matrix(nrow = p1, ncol = r))
  U2.draw <- lapply(1:nsample, function(i) matrix(nrow = p2, ncol = r))
  V1.draw <- lapply(1:nsample, function(i) matrix(nrow = n, ncol = r1))
  V2.draw <- lapply(1:nsample, function(i) matrix(nrow = n, ncol = r2))
  W1.draw <- lapply(1:nsample, function(i) matrix(nrow = p1, ncol = r1))
  W2.draw <- lapply(1:nsample, function(i) matrix(nrow = p2, ncol = r2))
  
  if (response_given) {
    beta.draw <- lapply(1:nsample, function(i) matrix(nrow = n_beta, ncol = 1))
    
    if (response_type == "binary") {
      Z.draw <- lapply(1:nsample, function(i) matrix(nrow = n, ncol = 1))
    }
    
    if (response_type == "continuous") {
      tau2.draw <- lapply(1:nsample, function(i) matrix(nrow = 1, ncol = 1))
    }
    
    if (missingness_in_response) {
      Ym.draw <- lapply(1:nsample, function(i) matrix(nrow = n, ncol = 1)) 
    }
  }
  
  if (missingness_in_data) {
    Xm.draw <- list(Xm1.draw = lapply(1:nsample, function(i) matrix(nrow = p1, ncol = n)),
                    Xm2.draw = lapply(1:nsample, function(i) matrix(nrow = p2, ncol = n)))
  }
  
  # ---------------------------------------------------------------------------
  # Storing the initial values 
  # ---------------------------------------------------------------------------
  
  V.draw[[1]] <- V0
  U1.draw[[1]] <- U10
  U2.draw[[1]] <- U20
  V1.draw[[1]] <- V10
  V2.draw[[1]] <- V20
  W1.draw[[1]] <- W10
  W2.draw[[1]] <- W20
  
  if (response_given) {
    beta.draw[[1]] <- beta0
    
    if (response_type == "binary") {
      Z.draw[[1]] <- Z0
    }
    
    if (response_type == "continuous") {
      tau2.draw[[1]] <- tau20
    }
    
    if (missingness_in_response) {
      Ym.draw[[1]] <- Ym0
    }
  }
  
  if (missingness_in_data) {
    Xm10 <- U10 %*% t(V0) + W10 %*% t(V10)
    Xm.draw$Xm1.draw[[1]] <- Xm10
    
    Xm20 <- U20 %*% t(V0) + W20 %*% t(V20)
    Xm.draw$Xm2.draw[[1]] <- Xm20
  }
  
  # ---------------------------------------------------------------------------
  # Computing the inverses (all diagonal matrices)
  # ---------------------------------------------------------------------------
  
  if (!response_given) {
    # For V - Combined error variances between X1, X2
    Sigma_V_Inv <- diag(1/c(rep(sigma21, p1), rep(sigma22, p2)))
    
    # For V1 - Combined error variances between X1 
    Sigma_V1_Inv <- diag(1/c(rep(sigma21, p1)))
    
    # For V2 - Combined error variances between X2 
    Sigma_V2_Inv <- diag(1/c(rep(sigma22, p2)))
  }
  
  if (response_given) {
    if (response_type == "binary") {
      # For V - Combined error variances between X1, X2, and Z
      Sigma_V_Inv <- diag(1/c(rep(sigma21, p1), rep(sigma22, p2), 1))
      
      # For V1 - Combined error variances between X1 and Z
      Sigma_V1_Inv <- diag(1/c(rep(sigma21, p1), 1))
      
      # For V2 - Combined error variances between X2 and Z
      Sigma_V2_Inv <- diag(1/c(rep(sigma22, p2), 1))
    } 
    
    # For beta - Combined error variances between intercept and all betas
    Sigma_beta_Inv <- solve(Sigma_beta)
  }
  
  # ---------------------------------------------------------------------------
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------
  
  for (iter in 1:nsample) {
    if (progress) svMisc::progress(iter/(nsample/100))
    
    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------
    
    V.iter <- V.draw[[iter]]
    U1.iter <- U1.draw[[iter]]
    U2.iter <- U2.draw[[iter]]
    V1.iter <- V1.draw[[iter]]
    V2.iter <- V2.draw[[iter]]
    W1.iter <- W1.draw[[iter]]
    W2.iter <- W2.draw[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.draw[[iter]]
      
      # Breaking them down into the intercept, joint effects, individual effects
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE]
      beta_indiv1.iter <- beta.iter[(r+2):(r+1+r1),, drop = FALSE]
      beta_indiv2.iter <- beta.iter[(r+1+r1+1):n_beta,, drop = FALSE]
      
      if (response_type == "binary") {
        Z.iter <- Z.draw[[iter]]
      }
      
      if (response_type == "continuous") {
        tau2.iter <- tau2.draw[[iter]]
      }
      
      if (missingness_in_response) {
        # Save the current imputations for the missing values
        Ym.iter <- Ym.draw[[iter]] 
        
        # Creating the completed outcome vector
        Y_complete <- Y
        
        # Filling in the missing entries for R1 and R2. 
        Y_complete[missing_obs_Y] <- Ym.iter[missing_obs_Y] 
      }
      
      if (!missingness_in_response) {
        Y_complete <- Y
      }
    }
    
    if (missingness_in_data) {
      # Storing the current imputations for the missing values
      Xm.iter <- list(Xm1.iter = Xm.draw$Xm1.draw[[iter]],
                      Xm2.iter = Xm.draw$Xm2.draw[[iter]])
      
      # Creating the completed matrices. 
      X1_complete = X1
      X2_complete = X2
      
      # Fill in the completed matrices with the imputed values
      X1_complete[missing_obs_data$X1] <- Xm.iter$Xm1.iter[missing_obs_data$X1]
      X2_complete[missing_obs_data$X2] <- Xm.iter$Xm2.iter[missing_obs_data$X2]
    }
    
    if (!missingness_in_data) {
      X1_complete <- X1
      X2_complete <- X2
    }
    
    # -------------------------------------------------------------------------
    # Computing the inverse that changes with tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        # For V - Combined error variances between X1, X2, and Y
        Sigma_V_Inv <- diag(1/c(rep(sigma21, p1), rep(sigma22, p2), tau2.iter))
        
        # For V1 - Combined error variances between X1 and Z
        Sigma_V1_Inv <- diag(1/c(rep(sigma21, p1), tau2.iter))
        
        # For V2 - Combined error variances between X2 and Z
        Sigma_V2_Inv <- diag(1/c(rep(sigma22, p2), tau2.iter))
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for U1
    # -------------------------------------------------------------------------
    
    X1.iter <- X1_complete - W1.iter %*% t(V1.iter)
    Bu1 <- solve((1/sigma21) * t(V.iter) %*% V.iter + (1/sigma2_joint) * diag(r))
    U1.draw[[iter+1]] <- t(sapply(1:p1, function(l1) {
      bu1 <- (1/sigma21) * t(V.iter) %*% X1.iter[l1, ]
      
      U1i <- mvrnorm(1, mu = Bu1 %*% bu1, Sigma = Bu1)
      U1i
    }))
    
    # -------------------------------------------------------------------------
    # Posterior sample for U2
    # -------------------------------------------------------------------------
    
    X2.iter <- X2_complete - W2.iter %*% t(V2.iter)
    Bu2 <- solve((1/sigma22) * t(V.iter) %*% V.iter + (1/sigma2_joint) * diag(r))
    U2.draw[[iter+1]] <- t(sapply(1:p2, function(l2) {
      bu2 <- (1/sigma22) * t(V.iter) %*% X2.iter[l2, ]
      
      U2i <- mvrnorm(1, mu = Bu2 %*% bu2, Sigma = Bu2)
      U2i
    }))
    
    # Update the current value of U1 and U2
    U1.iter <- U1.draw[[iter+1]]
    U2.iter <- U2.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for V
    # -------------------------------------------------------------------------
    
    if (!response_given) {
      # Concatenating Ui's together
      U.iter <- rbind(U1.iter, U2.iter)
      
      # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
      tU_Sigma <- crossprod(U.iter, Sigma_V_Inv)
      tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter)
      
      # The combined centered Xis with the latent response vector
      X.iter <- rbind(X1_complete - W1.iter %*% t(V1.iter),
                      X2_complete - W2.iter %*% t(V2.iter))
      
      Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
      
      V.draw[[iter+1]] <- t(sapply(1:n, function(j) {
        bv <-  tU_Sigma %*% X.iter[,j]
        
        Vj <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
        Vj
      }))
    }
    
    if (response_given) {
      # Concatenating Ui's together
      U.iter <- rbind(U1.iter, U2.iter, t(beta_joint.iter))
      
      # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
      tU_Sigma <- crossprod(U.iter, Sigma_V_Inv) 
      tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter)
      
      if (response_type == "binary") {
        Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
        V.draw[[iter+1]] <- t(sapply(1:n, function(i) {
          # The combined centered Xis with the latent response vector
          X.iter <- rbind(X1_complete - W1.iter %*% t(V1.iter),
                          X2_complete - W2.iter %*% t(V2.iter),
                          (Z.iter - c(beta_intercept.iter) - V1.iter %*% beta_indiv1.iter - V2.iter %*% beta_indiv2.iter)[i,])
          
          bv <- tU_Sigma %*% X.iter[,i]
          
          Vj <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
          Vj
        }))
      }
      
      if (response_type == "continuous") {
        Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
        V.draw[[iter+1]] <- t(sapply(1:n, function(i) {
          # The combined centered Xis with the latent response vector
          X.iter <- rbind(X1_complete - W1.iter %*% t(V1.iter),
                          X2_complete - W2.iter %*% t(V2.iter),
                          (Y_complete - c(beta_intercept.iter) - V1.iter %*% beta_indiv1.iter - V2.iter %*% beta_indiv2.iter)[i,])
          
          bv <- tU_Sigma %*% X.iter[,i]
          
          Vj <- mvrnorm(1, mu = Bv %*% bv, Sigma = Bv)
          Vj
        }))
      }
    }
    
    # Updating the value of V
    V.iter <- V.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for V1
    # -------------------------------------------------------------------------
    
    if (!response_given) {
      X1.iter <- X1_complete - U1.iter %*% t(V.iter)
      Bv1 <- solve((1/sigma21) * t(W1.iter) %*% W1.iter + (1/sigma2_indiv1) * diag(r1))
      
      V1.draw[[iter+1]] <- t(sapply(1:n, function(j) {
        bv1 <- (1/sigma21) * t(W1.iter) %*% X1.iter[, j]
        
        V1j <- mvrnorm(1, mu = Bv1 %*% bv1, Sigma = Bv1)
        V1j
      }))
    }
    
    if (response_given) {
      if (response_type == "binary") {
        # Combined centered X1 and Z
        X.iter <- rbind(X1_complete - U1.iter %*% t(V.iter),
                        t(Z.iter - c(beta_intercept.iter) - V.iter %*% beta_joint.iter - V2.iter %*% beta_indiv2.iter))
        
        # Combined Ws and beta
        W.iter <- rbind(W1.iter, t(beta_indiv1.iter))
        
        tW_Sigma <- crossprod(W.iter, Sigma_V1_Inv)
        tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter)
        
        Bv1 <- solve(tW_Sigma_W + (1/sigma2_indiv1) * diag(r1))
        V1.draw[[iter+1]] <- t(sapply(1:n, function(i) {
          bv1 <- tW_Sigma %*% X.iter[, i]
          
          V1j <- mvrnorm(1, mu = Bv1 %*% bv1, Sigma = Bv1)
          V1j
        }))
      }
      
      if (response_type == "continuous") {
        # Combined centered X1 and Y
        Y.iter <- t(Y_complete - c(beta_intercept.iter) - V.iter %*% beta_joint.iter - V2.iter %*% beta_indiv2.iter)
        colnames(Y.iter) <- colnames(X1)
        X.iter <- rbind(X1_complete - U1.iter %*% t(V.iter), Y.iter)
        
        # Combined Ws and beta
        W.iter <- rbind(W1.iter, t(beta_indiv1.iter))
        
        tW_Sigma <- crossprod(W.iter, Sigma_V1_Inv)
        tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter)
        
        Bv1 <- solve(tW_Sigma_W + (1/sigma2_indiv1) * diag(r1))
        V1.draw[[iter+1]] <- t(sapply(1:n, function(i) {
          bv1 <- tW_Sigma %*% X.iter[, i]
          
          V1j <- mvrnorm(1, mu = Bv1 %*% bv1, Sigma = Bv1)
          V1j
        }))
      }
    }
    
    # Update the current value of V1
    V1.iter <- V1.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for V2
    # -------------------------------------------------------------------------
    
    if (!response_given) {
      X2.iter <- X2_complete - U2.iter %*% t(V.iter)
      Bv2 <- solve((1/sigma22) * t(W2.iter) %*% W2.iter + (1/sigma2_indiv2) * diag(r2))
      
      V2.draw[[iter+1]] <- t(sapply(1:n, function(j) {
        bv2 <- (1/sigma22) * t(W2.iter) %*% X2.iter[, j]
        
        V12 <- mvrnorm(1, mu = Bv2 %*% bv2, Sigma = Bv2)
        V12
      }))
    }
    
    if (response_given) {
      if (response_type == "binary") {
        # Combined centered X2 and Z
        X.iter <- rbind(X2_complete - U2.iter %*% t(V.iter),
                        t(Z.iter - c(beta_intercept.iter) - V.iter %*% beta_joint.iter - V1.iter %*% beta_indiv1.iter))
        
        # Combined Ws and beta
        W.iter <- rbind(W2.iter, t(beta_indiv2.iter))
        
        tW_Sigma <- crossprod(W.iter, Sigma_V2_Inv)
        tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter)
        
        Bv2 <- solve(tW_Sigma_W + (1/sigma2_indiv2) * diag(r2))
        V2.draw[[iter+1]] <- t(sapply(1:n, function(i) {
          bv2 <- tW_Sigma %*% X.iter[, i]
          
          V2j <- mvrnorm(1, mu = Bv2 %*% bv2, Sigma = Bv2)
          V2j
        }))
      }
      
      if (response_type == "continuous") {
        # Combined centered X2 and Y
        X.iter <- rbind(X2_complete - U2.iter %*% t(V.iter),
                        t(Y_complete - c(beta_intercept.iter) -  V.iter %*% beta_joint.iter - V1.iter %*% beta_indiv1.iter))
        
        # Combined Ws and beta
        W.iter <- rbind(W2.iter, t(beta_indiv2.iter))
        
        tW_Sigma <- crossprod(W.iter, Sigma_V2_Inv)
        tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter)
        
        Bv2 <- solve(tW_Sigma_W + (1/sigma2_indiv2) * diag(r2))
        V2.draw[[iter+1]] <- t(sapply(1:n, function(i) {
          bv2 <- tW_Sigma %*% X.iter[, i]
          
          V2j <- mvrnorm(1, mu = Bv2 %*% bv2, Sigma = Bv2)
          V2j
        }))
      }
    }
    
    # Update the current value of V2
    V2.iter <- V2.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for W1
    # -------------------------------------------------------------------------
    
    X1.iter <- X1_complete - U1.iter %*% t(V.iter)
    Bw1 <- solve((1/sigma21) * t(V1.iter) %*% V1.iter + (1/sigma2_indiv1) * diag(r1))
    
    W1.draw[[iter+1]] <- t(sapply(1:p1, function(l1) {
      bw1 <- (1/sigma21) * t(V1.iter) %*% X1.iter[l1,] 
      
      W1i <- mvrnorm(1, mu = Bw1 %*% bw1, Sigma = Bw1)
      W1i
    }))
    
    # -------------------------------------------------------------------------
    # Posterior sample for W2
    # -------------------------------------------------------------------------
    
    X2.iter <- X2_complete - U2.iter %*% t(V.iter)
    Bw2 <- solve((1/sigma22) * t(V2.iter) %*% V2.iter + (1/sigma2_indiv2) * diag(r2))
    W2.draw[[iter+1]] <- t(sapply(1:p2, function(l2) {
      bw2 <- (1/sigma22) * t(V2.iter) %*% X2.iter[l2,] 
      
      W2i <- mvrnorm(1, mu = Bw2 %*% bw2, Sigma = Bw2)
      W2i
    }))
    
    # Update the current value of W1 and W2
    W1.iter <- W1.draw[[iter+1]]
    W2.iter <- W2.draw[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      # Combine current values of V, V1, and V2
      VStar.iter <- cbind(1, V.iter, V1.iter, V2.iter)
      
      if (response_type == "continuous") {
        tau2.draw[[iter+1]] <- 1/rgamma(1, shape = shape + (n/2), rate = rate + 0.5 * sum((Y_complete - VStar.iter %*% beta.iter)^2))
        
        # Update the current value of tau2
        tau2.iter <- tau2.draw[[iter+1]]
      }
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Bbeta <- solve(t(VStar.iter) %*% VStar.iter + Sigma_beta_Inv)
        bbeta <- t(VStar.iter) %*% Z.iter
        beta.draw[[iter+1]] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      }
      
      if (response_type == "continuous") {
        Bbeta <- solve((1/tau2.iter) * t(VStar.iter) %*% VStar.iter + Sigma_beta_Inv)
        bbeta <- (1/tau2.iter) * t(VStar.iter) %*% Y_complete
        beta.draw[[iter+1]] <- matrix(mvrnorm(1, mu = Bbeta %*% bbeta, Sigma = Bbeta), ncol = 1)
      }
      
      # Update the current value of beta
      beta.iter <- beta.draw[[iter+1]]
      
      # Breaking them down into the intercept, joint effects, individual effects
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE]
      beta_indiv1.iter <- beta.iter[(r+2):(r+1+r1),, drop = FALSE]
      beta_indiv2.iter <- beta.iter[(r+1+r1+1):n_beta,, drop = FALSE]
    }
    
    # -------------------------------------------------------------------------
    # Posterior sample for latent continuous response Z
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "binary") {
        Z.draw[[iter+1]] <- matrix(sapply(1:n, function(i) {
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
          Ym.draw[[iter+1]] <- matrix(rnorm(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter)), ncol = 1)
        }
        
        if (response_type == "binary") {
          Ym.draw[[iter+1]] <- matrix(rbinom(n, size = 1, prob = pnorm(VStar.iter %*% beta.iter)), ncol = 1)
        }
      }
    }
    
    if (missingness_in_data) {
      E1.iter <- matrix(rnorm(p1*n, 0, sqrt(sigma21)), nrow = p1, ncol = n)
      Xm.draw$Xm1.draw[[iter+1]] <- (U1.iter %*% t(V.iter) + W1.iter %*% t(V1.iter) + E1.iter)[missing_obs_data$X1]
      
      E2.iter <- matrix(rnorm(p2*n, 0, sqrt(sigma22)), nrow = p2, ncol = n)
      Xm.draw$Xm2.draw[[iter+1]] <- (U2.iter %*% t(V.iter) + W2.iter %*% t(V2.iter) + E2.iter)[missing_obs_data$X2]
    }
  }
  
  # Return
  if (!response_given & !missingness_in_data) {
    return(list(X1 = X1, # Returning the scaled version of the data
                X2 = X2, # Returning the scaled version of the data
                sigma.mat = sigma.mat, # Returning the scaling factors
                V.draw = V.draw, 
                U1.draw = U1.draw,
                U2.draw = U2.draw,
                V1.draw = V1.draw,
                V2.draw = V2.draw,
                W1.draw = W1.draw,
                W2.draw = W2.draw))
  }
  
  if (!response_given & missingness_in_data) {
    return(list(X1 = X1, # Returning the scaled version of the data
                X2 = X2, # Returning the scaled version of the data
                sigma.mat = sigma.mat, # Returning the scaling factors
                V.draw = V.draw, 
                U1.draw = U1.draw,
                U2.draw = U2.draw,
                V1.draw = V1.draw,
                V2.draw = V2.draw,
                W1.draw = W1.draw,
                W2.draw = W2.draw,
                Xm.draw = Xm.draw))
  }
  
  if (response_given) {
    if (!missingness_in_response) {
      if (response_type == "continuous" & !missingness_in_data) {
        return(list(X1 = X1, # Returning the scaled version of the data
                    X2 = X2, # Returning the scaled version of the data
                    sigma.mat = sigma.mat, # Returning the scaling factors
                    V.draw = V.draw, 
                    U1.draw = U1.draw,
                    U2.draw = U2.draw,
                    V1.draw = V1.draw,
                    V2.draw = V2.draw,
                    W1.draw = W1.draw,
                    W2.draw = W2.draw,
                    beta.draw = beta.draw,
                    tau2.draw = tau2.draw))
      }
      
      if (response_type == "binary" & !missingness_in_data) {
        return(list(X1 = X1, # Returning the scaled version of the data
                    X2 = X2, # Returning the scaled version of the data
                    sigma.mat = sigma.mat, # Returning the scaling factors
                    V.draw = V.draw, 
                    U1.draw = U1.draw,
                    U2.draw = U2.draw,
                    V1.draw = V1.draw,
                    V2.draw = V2.draw,
                    W1.draw = W1.draw,
                    W2.draw = W2.draw,
                    beta.draw = beta.draw,
                    Z.draw = Z.draw))
      }
      
      if (response_type == "continuous" & missingness_in_data) {
        return(list(X1 = X1, # Returning the scaled version of the data
                    X2 = X2, # Returning the scaled version of the data
                    sigma.mat = sigma.mat, # Returning the scaling factors
                    V.draw = V.draw, 
                    U1.draw = U1.draw,
                    U2.draw = U2.draw,
                    V1.draw = V1.draw,
                    V2.draw = V2.draw,
                    W1.draw = W1.draw,
                    W2.draw = W2.draw,
                    beta.draw = beta.draw,
                    tau2.draw = tau2.draw,
                    Xm.draw = Xm.draw))
      }
      
      if (response_type == "binary" & missingness_in_data) {
        return(list(X1 = X1, # Returning the scaled version of the data
                    X2 = X2, # Returning the scaled version of the data
                    sigma.mat = sigma.mat, # Returning the scaling factors
                    V.draw = V.draw, 
                    U1.draw = U1.draw,
                    U2.draw = U2.draw,
                    V1.draw = V1.draw,
                    V2.draw = V2.draw,
                    W1.draw = W1.draw,
                    W2.draw = W2.draw,
                    beta.draw = beta.draw,
                    Z.draw = Z.draw,
                    Xm.draw = Xm.draw))
      }
    }
    
    if (missingness_in_response) {
      if (response_type == "continuous" & !missingness_in_data) {
        return(list(X1 = X1, # Returning the scaled version of the data
                    X2 = X2, # Returning the scaled version of the data
                    sigma.mat = sigma.mat, # Returning the scaling factors
                    V.draw = V.draw, 
                    U1.draw = U1.draw,
                    U2.draw = U2.draw,
                    V1.draw = V1.draw,
                    V2.draw = V2.draw,
                    W1.draw = W1.draw,
                    W2.draw = W2.draw,
                    beta.draw = beta.draw,
                    tau2.draw = tau2.draw,
                    Ym.draw = Ym.draw))
      }
      
      if (response_type == "binary" & !missingness_in_data) {
        return(list(X1 = X1, # Returning the scaled version of the data
                    X2 = X2, # Returning the scaled version of the data
                    sigma.mat = sigma.mat, # Returning the scaling factors
                    V.draw = V.draw, 
                    U1.draw = U1.draw,
                    U2.draw = U2.draw,
                    V1.draw = V1.draw,
                    V2.draw = V2.draw,
                    W1.draw = W1.draw,
                    W2.draw = W2.draw,
                    beta.draw = beta.draw,
                    Ym.draw = Ym.draw,
                    Z.draw = Z.draw))
      }
      
      if (response_type == "continuous" & missingness_in_data) {
        return(list(X1 = X1, # Returning the scaled version of the data
                    X2 = X2, # Returning the scaled version of the data
                    sigma.mat = sigma.mat, # Returning the scaling factors
                    V.draw = V.draw, 
                    U1.draw = U1.draw,
                    U2.draw = U2.draw,
                    V1.draw = V1.draw,
                    V2.draw = V2.draw,
                    W1.draw = W1.draw,
                    W2.draw = W2.draw,
                    beta.draw = beta.draw,
                    tau2.draw = tau2.draw,
                    Ym.draw = Ym.draw,
                    Xm.draw = Xm.draw))
      }
      
      if (response_type == "binary" & missingness_in_data) {
        return(list(X1 = X1, # Returning the scaled version of the data
                    X2 = X2, # Returning the scaled version of the data
                    sigma.mat = sigma.mat, # Returning the scaling factors
                    V.draw = V.draw, 
                    U1.draw = U1.draw,
                    U2.draw = U2.draw,
                    V1.draw = V1.draw,
                    V2.draw = V2.draw,
                    W1.draw = W1.draw,
                    W2.draw = W2.draw,
                    beta.draw = beta.draw,
                    tau2.draw = tau2.draw,
                    Ym.draw = Ym.draw,
                    Z.draw = Z.draw,
                    Xm.draw = Xm.draw))
      }
    }
  }
}


# -----------------------------------------------------------------------------
# Helper functions for initializing with BIDIFAC
# -----------------------------------------------------------------------------

# packages
library(Matrix)

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
  n.vec=do.call(c, lapply(data[1,], ncol))
  
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
bpmf_data <- function(parameters, hyperparameters, s2n, response, missingness, entrywise, prop_missing) {
  # Generates fake data depending on the dims provided in `parameters`
  # and the model parameters provided in `hyperparameters`
  # s2n = desired signal-to-noise ratio 
  
  # -------------------------------------------------------------------------
  # Setting the dimensions and latent components
  # -------------------------------------------------------------------------

  p1 <- parameters$p1 # number rows X1
  p2 <- parameters$p2 # number rows X2
  n <- parameters$n # number cols X1 and X2
  r <- parameters$r # number latent components of X1, X2
  r1 <- parameters$r1 # number latent components of individual structure X1
  r2 <- parameters$r2 # number of latent components of individual structure X2
  n_beta <- 1 + r + r1 + r2
  
  # Setting the hyperparameters
  sigma21 <- hyperparameters$sigma21 # Error variance for X1
  sigma22 <- hyperparameters$sigma22 # Error variance for X2
  sigma2_joint_model <- hyperparameters$sigma2_joint_model # Model variance of joint structure
  sigma2_indiv1_model <- hyperparameters$sigma2_indiv1_model # Model variance of individual structure for X1
  sigma2_indiv2_model <- hyperparameters$sigma2_indiv2_model # Model variance of individual structure for X2
  sigma2_joint_true <- hyperparameters$sigma2_joint_true # Variance of joint structure for generating data
  sigma2_indiv1_true <- hyperparameters$sigma2_indiv1_true # Variance of individual structure for X1 for generating data
  sigma2_indiv2_true <- hyperparameters$sigma2_indiv2_true # Variance of individual structure for X2 for generating data
  Sigma_beta <- hyperparameters$Sigma_beta # Covariance for betas
  shape <- hyperparameters$shape 
  rate <- hyperparameters$rate
  
  # -------------------------------------------------------------------------
  # Generating the underlying structure
  # -------------------------------------------------------------------------
  
  V <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint_true)), nrow = n, ncol = r)
  U1 <- matrix(rnorm(p1*r, mean = 0, sd = sqrt(sigma2_joint_true)), nrow = p1, ncol = r)
  U2 <- matrix(rnorm(p2*r, mean = 0, sd = sqrt(sigma2_joint_true)), nrow = p2, ncol = r)
  V1 <- matrix(rnorm(n*r1, mean = 0, sd = sqrt(sigma2_indiv1_true)), nrow = n, ncol = r1)
  V2 <- matrix(rnorm(n*r2, mean = 0, sd = sqrt(sigma2_indiv2_true)), nrow = n, ncol = r2)
  W1 <- matrix(rnorm(p1*r1, mean = 0, sd = sqrt(sigma2_indiv1_true)), nrow = p1, ncol = r1)
  W2 <- matrix(rnorm(p2*r2, mean = 0, sd = sqrt(sigma2_indiv2_true)), nrow = p2, ncol = r2)
  
  # Generating the residual error
  E1 <- matrix(rnorm(p1*n, 0, sqrt(sigma21)), nrow = p1, ncol = n)
  E2 <- matrix(rnorm(p2*n, 0, sqrt(sigma22)), nrow = p2, ncol = n)
  
  # Calculating the structure
  X1_joint_structure <- U1 %*% t(V)
  X1_indiv_structure <- W1 %*% t(V1)
  
  X2_joint_structure <- U2 %*% t(V)
  X2_indiv_structure <- W2 %*% t(V2)
  
  # -------------------------------------------------------------------------
  # Standardizing the variance of the signal
  # -------------------------------------------------------------------------
  
  # Calculating the scaling coefficient so that the variance of the underlying structure = s2n * noise variance
  s2n_coef_X1 <- s2n * sd(c(E1))/sd(c(X1_joint_structure + X1_indiv_structure))
  s2n_coef_X2 <- s2n * sd(c(E2))/sd(c(X2_joint_structure + X2_indiv_structure))
  
  # Scaling the underlying structure
  X1_joint_structure_s2n <- s2n_coef_X1 * X1_joint_structure
  X1_indiv_structure_s2n <- s2n_coef_X1 * X1_indiv_structure
  
  X2_joint_structure_s2n <- s2n_coef_X2 * X2_joint_structure
  X2_indiv_structure_s2n <- s2n_coef_X2 * X2_indiv_structure
  
  # -------------------------------------------------------------------------
  # Calculating the observed data
  # -------------------------------------------------------------------------
  
  X1 <- X1_joint_structure_s2n + X1_indiv_structure_s2n + E1
  X2 <- X2_joint_structure_s2n + X2_indiv_structure_s2n + E2
  
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
    # Generate betas
    beta <- matrix(mvrnorm(1, mu = rep(0, n_beta), Sigma = Sigma_beta), ncol = 1)
    
    # Combine the Vs
    VStar <- cbind(1, V, V1, V2)
    
    # True probability of being a case
    Prob.VStar.beta <- pnorm(VStar %*% beta)
    
    if (response == "binary") {
      Y <- matrix(rbinom(n, size = 1, prob = pnorm(VStar %*% beta)), ncol = 1)
      tau2 <- NULL
    }
    
    if (response == "continuous") {
      tau2 <- 1/rgamma(1, shape = shape, rate = rate)
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
      missing_obs_X1 <- NULL
      missing_obs_X2 <- NULL
      
      X1_missing <- NULL
      X2_missing <- NULL
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
      if (entrywise) { # if removing observations entrywise
        # these are counters going down the columns of R. So 9 would be the 9th entry counting down. 
        missing_obs_X1 <- sample(x = 1:length(X1), size = prop_missing*length(X1), replace = FALSE) 
        missing_obs_X2 <- sample(x = 1:length(X2), size = prop_missing*length(X2), replace = FALSE)
        
        # Duplicate X1 and X2 so that I have one with the full data and one with the missing data
        X1_missing <- X1_scaled
        X1_missing[missing_obs_X1] <- NA
        
        X2_missing <- X2_scaled
        X2_missing[missing_obs_X2] <- NA
      } else { # if removing entire columns
        # Gives the column indices to remove
        cols_to_remove_X1 <- sample(x=1:n, size = n*prop_missing, replace = FALSE)
        cols_to_remove_X2 <- sample(x=1:n, size = n*prop_missing, replace = FALSE)
        
        # Enforce different subjects to be missing from either dataset
        any_overlap <- length(intersect(cols_to_remove_X1, cols_to_remove_X2)) != 0
        while (any_overlap) {
          # Sample the columns to remove from X2 again
          cols_to_remove_X2 <- sample(x=1:n, size = n*prop_missing, replace = FALSE)
          
          # Check for overlap
          any_overlap <- length(intersect(cols_to_remove_X1, cols_to_remove_X2)) != 0
        }
        
        X1_missing <- X1_scaled
        X1_missing[,cols_to_remove_X1] <- NA
        
        X2_missing <- X2_scaled
        X2_missing[,cols_to_remove_X2] <- NA
        
        missing_obs_X1 <- which(is.na(X1_missing))
        missing_obs_X2 <- which(is.na(X2_missing))
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # Return
  # -------------------------------------------------------------------------
  
  list(X1 = X1, X2 = X2, # The "observed data"
       Y = Y, # The "observed outcome"
       X1_missing = X1_missing, X2_missing = X2_missing, # Missing data 
       missing_obs_X1 = missing_obs_X1, missing_obs_X2 = missing_obs_X2, # Missing data 
       Y_missing = Y_missing, missing_obs_Y = missing_obs_Y, # Missing data 
       s2n = s2n, s2n_coef_X1 = s2n_coef_X1, s2n_coef_X2 = s2n_coef_X2, # Scaling for s2n
       X1_joint_structure_s2n = X1_joint_structure_s2n, X2_joint_structure_s2n = X2_joint_structure_s2n, # Joint structure
       X1_indiv_structure_s2n = X1_indiv_structure_s2n, X2_indiv_structure_s2n = X2_indiv_structure_s2n, # Individual structure
       beta = beta, tau2 = tau2)
}

# Center the data and the underlying structure 
center_data <- function(data, structure) {
  # Center each entry in data = list(X1, X2)
  data_centered <- lapply(1:length(data), function(i) list())
  means_for_centering <- lapply(1:length(data), function(i) list())
  structure_centered <- lapply(1:length(structure), function(i) list())
  
  for (i in 1:length(data)) {
    # Center the data itself
    data_centered[[i]] <- scale(data[[i]], center = TRUE, scale = FALSE)
    
    # Use the means for centering to center the structure
    means_for_centering[[i]] <- attr(data_centered[[i]], "scaled:center")
    
    # Scale the structure
    structure_centered[[i]] <- list(joint = sweep(structure[[i]]$joint, 2, means_for_centering[[i]]),
                                    indiv = sweep(structure[[i]]$indiv, 2, means_for_centering[[i]]))
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
  sim_results <- lapply(1:n_param, function(i) list(avg_coverage = 0,
                                                    avg_MSE = 0,
                                                    avg_CI_width = 0))
  names(sim_results) <- names(draws)
  
  # Iterate through the parameters, checking the coverage, MSE, and CI width
  for (param in 1:n_param) {
    sim_results[[param]]$avg_coverage <- sim_results[[param]]$avg_coverage + check_coverage(truth[[param]], draws[[param]], burnin = burnin)
    sim_results[[param]]$avg_MSE <- sim_results[[param]]$avg_MSE + mse(truth[[param]], draws[[param]])
    sim_results[[param]]$avg_CI_width <- sim_results[[param]]$avg_CI_width + ci_width(draws[[param]], burnin = burnin)
  }
  
  # Return the results
  sim_results
}

# Count the number of times each observation in each dataset was observed
calculate_observed <- function(sim_results, parameters, response) {
  # sim_results (list) = list of results from the simulation

  # Count how many datasources there were (add 1 for the response)
  n_sources <- sum(substr(names(parameters), start = 1, stop = 1) == "p") + 1
  
  # Create a list of indices for observation
  obs_inds <- lapply(1:n_sources, function(i) {
    if (i != n_sources) 1:(parameters[[i]] * parameters$n) else 1:parameters$n
  })
  
  # Create a list to store the results
  observed_counts <- lapply(1:n_sources, function(i) {
    if (i != n_sources) rep(0, length(obs_inds[[i]])) else rep(0, length(obs_inds[[i]]))
  })
  
  for (i in 1:length(sim_results)) {
    for (j in 1:n_sources) {
      # Save the indices for the missing observations in the jth source (could be the response)
      current_missing_obs <- sim_results[[i]]$missing_obs[[j]]
      observed_counts[[j]] <- observed_counts[[j]] + !(obs_inds[[j]] %in% current_missing_obs)
    }
  }
  
  names(observed_counts) <- c(paste0("X", 1:(n_sources-1)), "Y")
  
  if (is.null(response)) observed_counts <- observed_counts[!(names(observed_counts) %in% "Y")]
  
  # Return
  observed_counts
}

# Average the results before returning from the function
average_results <- function(sim_results, observed_counts, nsim) {
  # results_list (list) = list of results to return that should be averaged
  # observed_counts (list) = list of how many times each observation in each dataset was missing
  # nsim (int) = the number of simulation iterations that were run
  
  n_results <- length(sim_results[[1]])
  sim_results_avg <- lapply(1:n_results, function(i) list())
  names(sim_results_avg) <- names(sim_results[[1]])[1:n_results]
  
  for (i in 1:n_results) {
    # Select the index for the observed counts (change this later)
    k <- if (i == 1 | i == 2) 1 else 2
    
    # Select all results from across sim iters
    all_results_i <- lapply(1:nsim, function(res) {
      sim_results[[res]][[i]]
    })
    
    # Sum across each result type
    sim_results_avg[[i]]$avg_coverage <- mean(Reduce("+", lapply(all_results_i, function(res) {
      res$avg_coverage
    }))/observed_counts[[k]])

    sim_results_avg[[i]]$avg_MSE <- mean(Reduce("+", lapply(all_results_i, function(res) {
      res$avg_MSE
    }))/observed_counts[[k]])
    
    sim_results_avg[[i]]$avg_CI_width <- mean(Reduce("+", lapply(all_results_i, function(res) {
      res$avg_CI_width
    }))/observed_counts[[k]])
  }
  
  # Return
  sim_results_avg
}