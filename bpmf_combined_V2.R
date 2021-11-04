# -----------------------------------------------------------------------------
# Bayesian PMF Gibbs sampling algorithm to sample from the posterior for the
# underlying joint and individual structure for two data sources, X1 and X2. 
# 
# Function allows there to be no response, a binary reseponse, or a normally 
# distributed response. 
# 
# Response may have missingness. X1 and X2 may have missingness. 
# -----------------------------------------------------------------------------

# Editing comments

library(MASS)
library(svMisc)

# -----------------------------------------------------------------------------
# Gibbs sampling algorithms
# -----------------------------------------------------------------------------

bayesian_pmf_2datasources <- function(X1, X2, Y = NULL, nuclear_norm_init = TRUE, dims, hyperparameters, nsample = 5000, progress = TRUE) {
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
# Coverage simulation
# -----------------------------------------------------------------------------

bayesian_pmf_2datasource_coverage_simulation <- function(nsample, 
                                                         parameters,
                                                         hyperparameters,
                                                         nsim = 1000,
                                                         s2n,
                                                         response = NULL, 
                                                         missingness = NULL, 
                                                         prop_missing = NULL, 
                                                         entrywise = NULL,
                                                         nuclear_norm_init = TRUE) {
  # Coverage simulation to check the combined, "all-in-one" model 
  # nsample = number of posterior samples to generate for each component
  # parameters = list(p1, p2, n, r, r1, r2)
  # s2n = signal-to-noise ratio
  # response = "none" or "binary" or "continuous" 
  # missingness = "none" if no imputation, "missingness_in_response", "missingness_in_data", or "both
  # prop_missing = NULL if no missing or in (0,1) if missingness
  # entrywise = NULL if no missing or TRUE if entrywise missingness in X1 and X2 or FALSE if columnwise missingness in X1 and X2
  
  # ---------------------------------------------------------------------------
  # Change the below to be fed in by the parameters and hyperparameters arguments above 
  # ---------------------------------------------------------------------------
  
  # Saving the dimensions and latent components
  p1 <- parameters$p1 # number rows X1
  p2 <- parameters$p2 # number rows X2
  n <- parameters$n # number cols X1 and X2
  
  # Saving the hyperparameters
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

  # Setting the parameters for the Gibbs sampling function arguments
  hyperparameters_for_model <- list(sigma21 = sigma21, sigma22 = sigma22, 
                                    sigma2_joint = sigma2_joint_model, sigma2_indiv1 = sigma2_indiv1_model, 
                                    sigma2_indiv2 = sigma2_indiv2_model, 
                                    Sigma_beta = Sigma_beta, shape = shape, rate = rate)
  
  cl <- makeCluster(n_clust)
  registerDoParallel(cl)
  sim_results <- foreach (sim_iter = 1:nsim, .packages = c("Matrix", "MASS")) %dopar% {
    
    # -------------------------------------------------------------------------
    # Generating the data
    # -------------------------------------------------------------------------
    
    sim_data <- bpmf_data(parameters = parameters,
                          hyperparameters = hyperparameters,
                          s2n = s2n, response = response, missingness = missingness,
                          entrywise = entrywise, prop_missing = prop_missing)
    
    # If there's missingness in the data
    X1_missing <- sim_data$X1_missing
    X2_missing <- sim_data$X2_missing
    
    missing_obs_X1 <- sim_data$missing_obs_X1
    missing_obs_X2 <- sim_data$missing_obs_X2
    
    # The response
    Y <- sim_data$Y
    Y_missing <- sim_data$Y_missing
    missing_obs_Y <- sim_data$missing_obs_Y
    
    # The standardizing coefficients
    s2n_coef_X1 <- sim_data$s2n_coef_X1
    s2n_coef_X2 <- sim_data$s2n_coef_X2
    
    # The response parameters
    beta <- sim_data$beta
    tau2 <- sim_data$tau2
    
    # -------------------------------------------------------------------------
    # Center the data
    # -------------------------------------------------------------------------
    
    centered_data <- center_data(data = list(X1 = sim_data$X1, X2 = sim_data$X2),
                                 structure = list(X1 = list(joint = sim_data$X1_joint_structure_s2n, 
                                                            indiv = sim_data$X1_indiv_structure_s2n),
                                                  X2 = list(joint = sim_data$X2_joint_structure_s2n, 
                                                            indiv = sim_data$X2_indiv_structure_s2n)))
    
    # Saving the centered data
    X1 <- centered_data$data_centered[[1]]
    X2 <- centered_data$data_centered[[2]]
    
    # Saving the centered structure
    X1_joint_structure_s2n_center <- centered_data$structure_centered[[1]]$joint
    X2_joint_structure_s2n_center <- centered_data$structure_centered[[2]]$joint
    
    X1_indiv_structure_s2n_center <- centered_data$structure_centered[[1]]$indiv
    X2_indiv_structure_s2n_center <- centered_data$structure_centered[[2]]$indiv
    
    # -------------------------------------------------------------------------
    # Running the Gibbs sampling algorithm
    # -------------------------------------------------------------------------
    
    # Naming the data going into the Gibbs sampling algorithm
    if (is.null(missingness)) {
      X1_for_model <- X1
      X2_for_model <- X2
      Y_for_model <- Y
    }
    
    if (!is.null(missingness)) {
      if (missingness == "missingness_in_data") {
        X1_for_model <- X1_missing
        X2_for_model <- X2_missing
      }
      
      if (missingness == "missingness_in_response") {
        Y_for_model <- Y_missing
      }
    }
    
    # Gibbs sampling
    res <- bayesian_pmf_2datasources(X1_for_model, X2_for_model, Y = Y_for_model, nuclear_norm_init = nuclear_norm_init, dims = parameters, hyperparameters_for_model, nsample = nsample, progress = FALSE)
    
    # -------------------------------------------------------------------------
    # Extracting the results for each of decomposition matrices
    # -------------------------------------------------------------------------
    U1.draw <- res$U1.draw
    U2.draw <- res$U2.draw
    V.draw <- res$V.draw
    W1.draw <- res$W1.draw
    W2.draw <- res$W2.draw
    V1.draw <- res$V1.draw
    V2.draw <- res$V2.draw
    
    if (!is.null(response)) {
      beta.draw <- res$beta.draw
      
      if (response == "binary") {
        Z.draw <- res$Z.draw
      }
      
      if (response == "continuous") {
        tau2.draw <- res$tau2.draw
      }
    }
    
    if (!is.null(missingness)) {
      if (missingness == "missingness_in_data") {
        Xm.draw <- res$Xm.draw
      }
      
      if (missingness == "missingness_in_response") {
        Ym.draw <- res$Ym.draw
      }
    }
    
    # -------------------------------------------------------------------------
    # Adding a burn-in
    # -------------------------------------------------------------------------
    burnin <- nsample/2
    
    U1.burnin <- U1.draw[burnin:nsample]
    U2.burnin <- U2.draw[burnin:nsample]
    V.burnin <- V.draw[burnin:nsample]
    W1.burnin <- W1.draw[burnin:nsample]
    W2.burnin <- W2.draw[burnin:nsample]
    V1.burnin <- V1.draw[burnin:nsample]
    V2.burnin <- V2.draw[burnin:nsample]
    
    if (!is.null(response)) {
      beta.burnin <- beta.draw[burnin:nsample]
      
      if (response == "binary") {
        Z.burnin <- Z.draw[burnin:nsample]
      }
      
      if (response == "continuous") {
        tau2.burnin <- tau2.draw[burnin:nsample]
      }
    }
    
    if (!is.null(missingness)) {
      if (missingness == "missingness_in_data" | missingness == "both") {
        Xm.burnin <- Xm.draw
        Xm.burnin$Xm1.draw <- Xm.burnin$Xm1.draw[burnin:nsample]
        Xm.burnin$Xm2.draw <- Xm.burnin$Xm2.draw[burnin:nsample]
      }
      
      if (missingness == "missingness_in_response" | missingness == "both") {
        Ym.burnin <- Ym.draw[burnin:nsample]
      }
    }
    
    # Computing the underlying structure for X1 and X2 after burn-in from the sampler
    X1_joint_structure_burnin <- lapply(1:(burnin+1), function(res) {
      U1.burnin[[res]] %*% t(V.burnin[[res]])
    })
    
    X1_indiv_structure_burnin <- lapply(1:(burnin+1), function(res) {
      W1.burnin[[res]] %*% t(V1.burnin[[res]])
    })
    
    X2_joint_structure_burnin <- lapply(1:(burnin+1), function(res) {
      U2.burnin[[res]] %*% t(V.burnin[[res]]) 
    })
    
    X2_indiv_structure_burnin <- lapply(1:(burnin+1), function(res) {
      W2.burnin[[res]] %*% t(V2.burnin[[res]])
    })

    # -------------------------------------------------------------------------
    # Checking the coverage
    # -------------------------------------------------------------------------
    
    # Saving the draws from the sampler together
    draws <- list(X1_joint_structure_burnin = X1_joint_structure_burnin,
                  X1_indiv_structure_burnin = X1_indiv_structure_burnin,
                  X2_joint_structure_burnin = X2_joint_structure_burnin,
                  X2_indiv_structure_burnin = X2_indiv_structure_burnin)
    
    # Saving the truth together for comparison
    truth <- list(X1_joint_structure_s2n_center = X1_joint_structure_s2n_center,
                  X1_indiv_structure_s2n_center = X1_indiv_structure_s2n_center,
                  X2_joint_structure_s2n_center = X2_joint_structure_s2n_center,
                  X2_indiv_structure_s2n_center = X2_indiv_structure_s2n_center)
    
    if (!is.null(missingness)) {
      if (missingness == "missingness_in_data" | missingness == "both") {
        draws$Xm1 = Xm.burnin$Xm1.draw
        draws$Xm2 = Xm.burnin$Xm2.draw
        
        truth$Xm1 = X1
        truth$Xm2 = X2
      }
      
      if (missingness == "missingness_in_response" | missingness == "both") {
        draws$Ym = Ym.burnin
        truth$Ym = Y
      }
    }
    
    if (!is.null(response)) {
      draws$beta = beta.burnin
      truth$beta = beta
      
      if (response == "continuous") {
        draws$tau2 = tau2.burnin
      }
    }
    
    # Adding up the coverage, MSE, and CI width from this iteration
    sim_iter_results <- get_results(truth, draws, burnin)
    
    # ---------------------------------------------------------------------------
    # Returning the results at the end of the loop
    # ---------------------------------------------------------------------------
    
    # Add the indices for the missing values in this sim_iter
    sim_iter_results$missing_obs <- list(X1 = missing_obs_X1,
                                         X2 = missing_obs_X2,
                                         Y = missing_obs_Y)
    
    # Return 
    sim_iter_results
  }
  stopCluster(cl)
  
  # ---------------------------------------------------------------------------
  # Averaging the results
  # ---------------------------------------------------------------------------
  
  # Calculating how many times each entry in the data was observed
  observed_counts <- calculate_observed(sim_results, parameters, response)
  
  # Taking the averages of all the results
  sim_results_avg <- average_results(sim_results, observed_counts, nsim)
  
  # ---------------------------------------------------------------------------
  # Returning the results
  # ---------------------------------------------------------------------------
  
  sim_results_avg
  
}
