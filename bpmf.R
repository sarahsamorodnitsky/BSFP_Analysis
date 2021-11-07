# -----------------------------------------------------------------------------
# Helper functions for Bayesian PMF
# -----------------------------------------------------------------------------

# Packages
library(Matrix)
library(MASS)

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
  missingness_in_data <- any(apply(data, 1, function(source) any(is.na(source[[1]])))) 
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
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
    print("posterior mode obtained")
    
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
    r <- ranks$r 
    r.vec <- unlist(ranks[!(names(ranks) %in% "r")])
  }
  
  r_total <- sum(r.vec)
  n_beta <- 1 + r + r_total
  
  # If a response is given, set up the variance matrix for the prior of the betas using the ranks
  if (response_given) {
    Sigma_beta <- matrix(0, nrow = n_beta, ncol = n_beta)
    diag(Sigma_beta) <- c(beta_vars[1], rep(beta_vars[-1], c(r, r.vec)))
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
  }
  
  # If there is missingness in Y, generate starting values for the missing entries
  if (response_given) {
    if (missingness_in_response) {
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
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  U.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  W.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = q))
  
  if (!response_given) {
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- NULL
  }
  
  if (!missingness_in_data) {
    Xm.draw <- NULL
  }

  if (response_given) {
    beta.draw <- matrix(nrow = nsample, ncol = n_beta)
    
    Z.draw <- matrix(nrow = nsample, ncol = n)
    
    if (response_type == "continuous") {
      tau2.draw <- c()
    }
    
    if (missingness_in_response) {
      Ym.draw <- lapply(1:nsample, function(i) list())
    }
  }
  
  if (missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  # ---------------------------------------------------------------------------
  # Storing the initial values 
  # ---------------------------------------------------------------------------
  
  V.draw[[1]] <- V0
  U.draw[[1]] <- U0
  Vs.draw[[1]] <- Vs0
  W.draw[[1]] <- W0
  
  if (response_given) {
    beta.draw[1,] <- beta0
    Z.draw[1,] <- Z0
    tau2.draw[1] <- tau20

    if (missingness_in_response) {
      Ym.draw[[1]] <- Ym0
    }
  }
  
  if (missingness_in_data) {
    for (s in 1:q) {
      Xm.draw[[1]][[s,1]] <- U.draw[[1]][[s,1]] %*% t(V.draw[[1]][[1,1]]) + W.draw[[1]][[s,s]] %*% t(Vs.draw[[1]][[1,s]])
    }
  }
  
  # ---------------------------------------------------------------------------
  # Computing the inverses (all diagonal matrices)
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
  
  nsample <- nsample-1
  for (iter in 1:nsample) {
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
      beta.indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept, joint, and individual effects
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE]
      beta_indiv.iter.temp <- beta.iter[(r+2):nrow(beta.iter),, drop = FALSE]
      
      for (s in 1:q) {
        if (s == 1) beta.indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
        if (s != 1) beta.indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
      }
      
      if (response_type == "binary") {
        Z.iter <- Z.draw[iter,]
      }
      
      if (response_type == "continuous") {
        tau2.iter <- tau2.draw[iter]
      }
      
      if (missingness_in_response) {
        # Save the current imputations for the missing values
        Ym.iter <- Ym.draw[iter,]
        
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
      U.iter <- rbind(do.call(rbind, U.iter), t(beta_joint.iter))
      
      # Computing the crossprod: t(U.iter) %*% solve(Sigma) %*% U.iter
      tU_Sigma <- crossprod(U.iter, SigmaVInv) 
      tU_Sigma_U <- crossprod(t(tU_Sigma), U.iter)
      
      if (response_type == "binary") {
        Bv <- solve(tU_Sigma_U + (1/sigma2_joint) * diag(r))
        V.draw[[iter+1]][[1]] <- t(sapply(1:n, function(i) {
          # The combined centered Xis with the latent response vector
          X.iter <- rbind(data.rearrange(X_complete)$out - data.rearrange(W.iter)$out %*% do.call(rbind, lapply(Vs.iter, t)),
                          (Z.iter - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% beta_indiv.iter)[i,])
          
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
                          (Y_complete - c(beta_intercept.iter) - do.call(cbind, Vs.iter) %*% beta_indiv.iter)[i,])
          
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
        W.iter <- rbind(W.iter[[s,s]], t(beta.indiv.iter[[s,1]]))
        
        tW_Sigma <- crossprod(W.iter, SigmaVsInv[[s,s]])
        tW_Sigma_W <- crossprod(t(tW_Sigma), W.iter)
        
        if (response_type == "binary") {
          # Combined centered Xs and Z
          Xs.iter <- rbind(X_complete[[s,1]] - U.iter[[s,1]] %*% t(V.iter[[1,1]]),
                          t(Z.iter - c(beta_intercept.iter) - V.iter[[1,1]] %*% beta_joint.iter - 
                              do.call(cbind, Vs.iter[1, !(1:ncol(Vs.iter) %in% s)]) %*% do.call(rbind, beta.indiv.iter[!(c(1:q) %in% s), 1])))
          
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
                               do.call(cbind, Vs.iter[1, !(1:ncol(Vs.iter) %in% s)]) %*% do.call(rbind, beta.indiv.iter[!(c(1:q) %in% s), 1])))
          
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
      beta.indiv.iter <- matrix(list(), ncol = 1, nrow = q)
      
      # Breaking beta down into the intercept, joint, and individual effects
      beta_intercept.iter <- beta.iter[1,, drop = FALSE]
      beta_joint.iter <- beta.iter[2:(r+1),, drop = FALSE]
      beta_indiv.iter.temp <- beta.iter[(r+2):nrow(beta.iter),, drop = FALSE]
      
      for (s in 1:q) {
        if (s == 1) beta.indiv.iter[[s, 1]] <- beta_indiv.iter.temp[1:r.vec[s],, drop = FALSE]
        if (s != 1) beta.indiv.iter[[s, 1]] <- beta_indiv.iter.temp[(r.vec[s-1]+1):(r.vec[s-1] + r.vec[s]),, drop = FALSE]
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
          Ym.draw[[iter+1]] <- matrix(n, mean = VStar.iter %*% beta.iter, sd = sqrt(tau2.iter), ncol = 1)[missing_obs_Y,]
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