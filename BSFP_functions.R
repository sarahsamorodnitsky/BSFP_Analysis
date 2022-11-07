# Bayesian Simultaneous Factorization and Prediction Functions

#' Bayesian Simultaneous Factorization and Prediction (BSFP)
#'
#' Given q sources of data and a continuous or binary outcome measured on n samples,
#' BSFP can decompose variation across the sources and use the estimate factors to 
#' predict the outcome.
#' Use this function to simulate samples from posterior distributions 
#' of joint and individual structures. This function
#' can be used in several ways: 
#' (1) Initialize at the theoretical mode of the decomposition. Priors are fixed
#' at theoretical values. Ranks are determined at initialization. Can include an
#' outcome here or not. If no outcome is included, no burn-in is used. If an outcome
#' is included, a burn-in must be specified or the default is nsample/2. 
#' (2) User specifies ranks to estimate. Model is initialized using priors and 
#' then function samples from posterior distributions. User must specifiy the 
#' hyperparameters for the prior distributions on the structures. 
#' @param data A matrix of lists or a list of matrices that share the same number of
#' columns. The matrices must be oriented in pxn orientation. May contain NAs if 
#' there are missing values in the dataset. 
#' @param Y A matrix of lists or a nx1 matrix of continuous or binary outcome. 
#' May be NULL if no outcome is given. May contain NAs if there are missing outcomes. 
#' @param nninit Boolean determining if nuclear-norm penalized objective is used
#' to initialize the model. If TRUE, ranks = NULL. 
#' @param model_params A list of hyperparameters for the model. 
#' May be left NULL if theoretical defaults are desired. Otherwise, must be in 
#' the following named-list form: model_params = list(error_vars = c(), joint_var = double(),
#' indiv_vars = c(), beta_vars = c(), response_vars = c()). error_vars, indiv_vars are 
#' vectors of length q for each source. response_vars must define the shape and rate for
#' the Inverse-Gamma prior of the response variance. beta_vars must be of length q+1
#' to specify the prior variance on the intercept, the prior variance on each joint factor's
#' contribution to the outcome, and for each individual factor's contribution from each source. 
#' @param ranks A list of length q+1 for the ranks of the joint and individual structures.
#' Leave NULL if nninit=TRUE. 
#' @param scores Matrix with scores estimated by existing factorization method. Use only
#' if desired to run Bayesian linear model with scores estimated separately. Otherwise, 
#' leave NULL. 
#' @param nsample Integer specifying the number of posterior samples to generate. 
#' @param burnin Integer specifying the number of posterior sampling iterations to burn. 
#' Default is nsample/2, which is used if burnin=NULL. 
#' @param progress Should a progress bar be displayed to visualize the progress of the sampler? 
#' @param starting_values List of initial values for Gibbs sampler. If NULL and nninit=TRUE, 
#' fixes at posterior mode. If NULL and nninit=FALSE, simulates from prior distributions. 
#' @export

bsfp <- function(data, Y, nninit = TRUE, model_params = NULL, ranks = NULL, scores = NULL, nsample, burnin = NULL, progress = TRUE, starting_values = NULL) {
  
  # ---------------------------------------------------------------------------
  # Determining type of input data
  # ---------------------------------------------------------------------------
  
  # Was the data input as a list?
  if (!("matrix" %in% class(data))) {
    
    # Save the number of sources
    q <- length(data)
    
    # Initialize new data matrix
    new_data <- matrix(list(), nrow = q, ncol = 1)
    
    # Add in the sources
    for (s in 1:q) {
      new_data[[s,1]] <- data[[s]]
    }
    
    # Rename the new version of the data
    data <- new_data
  }
  
  # If a response was given, was it input as a matrix of lists or as a matrix?
  if (!is.null(Y)) {
    if (class(Y[1,1]) != "list") {
      
      # Create a new version of Y
      new_Y <- matrix(list(), nrow = 1, ncol = 1)
      
      # Input the response
      new_Y[[1,1]] <- Y
      
      # Return to the Y variable
      Y <- new_Y
    }
  }
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  q <- nrow(data) # Number of sources
  p.vec <- apply(data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(data[[1,1]]) # Number of subjects
  
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
  # Extracting the model parameters
  # ---------------------------------------------------------------------------
  
  # If no model parameters are given
  if (is.null(model_params)) {
    
    # Error variances
    error_vars <- sapply(1:q, function(s) 1)
    
    # Variance of joint structure
    lambda_joint <- sqrt(sum(p.vec)) + sqrt(n)
    sigma2_joint <- joint_var <- 1/(lambda_joint)
    
    # Variance of individual structure for Biocrates
    lambda_indiv <- sapply(1:q, function(s) sqrt(p.vec[s]) + sqrt(n))
    sigma2_indiv <- indiv_vars <- 1/(lambda_indiv) 
    
    # For the regression coefficients, beta
    lambda2_intercept <- 1e6
    lambda2_joint <- 1
    lambda2_indiv1 <- 1
    lambda2_indiv2 <- 1
    beta_vars <- c(lambda2_intercept, lambda2_joint, lambda2_indiv1, lambda2_indiv2)
    
    # For the response vector
    shape <- 0.01
    rate <- 0.01
    
    # Putting the model parameters together
    model_params <- list(error_vars = error_vars,
                         joint_var = sigma2_joint,
                         indiv_vars = sigma2_indiv,
                         beta_vars = beta_vars,
                         response_vars = c(shape = shape, rate = rate))
  }
  
  # If model parameters are given
  if (!is.null(model_params)) {
    error_vars <- model_params$error_vars # Error variances
    sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
    sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
    
    # If a prior on the regression betas is given
    if (!is.null(model_params$beta_vars)) {
      beta_vars <- model_params$beta_vars # Variances on betas
    }
    
    # If a prior on the response variance is given
    if (!is.null(model_params$response_vars)) {
      response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
    }
    
    # If there is a response and model_params is not NULL, there must be beta priors and response variance priors
    if ((!is.null(Y[[1,1]]) & is.null(model_params$beta_vars))) {
      stop("A response vector is given so please provide hyperparameters for the regression coefficients in the form
           of model_params$beta_vars = c(intercept_prior_var, prior_var_on_joint_factors, prior_var_on_individual_factors_source_1, ..., prior_var_on_individual_factors_source_q")
    }
    
    # If there is a CONTINUOUS response and model_params is not NULL, there must be beta priors and response variance priors
    if ((!is.null(Y[[1,1]]) & is.null(model_params$response_vars))) {
      if (response_type == "continuous") {
        stop("A response vector is given so please provide hyperparameters for the response variance prior in the form
             of a shape and a rate parameter, i.e. model_params$response_vars = c(shape = a, rate = b).")
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Check for missingness in data
  # ---------------------------------------------------------------------------
  
  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
  # ---------------------------------------------------------------------------
  # Obtaining the ranks 
  # ---------------------------------------------------------------------------
  
  if (nninit) {
    
    # If there is not missing data
    if (!missingness_in_data) {
      rank_init <- BIDIFAC(data, rmt = TRUE, pbar = FALSE, scale_back = FALSE)
    }
    
    # If there is missing data
    if (missingness_in_data) {
      rank_init <- impute.BIDIFAC(data = data, rmt = TRUE, pbar = FALSE, scale_back = FALSE)
    }
    
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
    beta.draw <- tau2.draw <- Z.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (!missingness_in_data) {
    Xm.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  }
  
  if (response_given) {
    beta.draw <- Z.draw <- tau2.draw <- Ym.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
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
    
  }
  
  # If there is missingness in the data, generate starting values for the missing entries
  if (missingness_in_data) {
    Xm0 <- matrix(list(), ncol = 1, nrow = q)
    for (s in 1:q) {
      Xm0[[s,1]] <- rank_init$X[[s,1]][missing_obs[[s]]]
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
  # Scaling the imputed data back to the original data scale
  # ---------------------------------------------------------------------------
  
  # If there is any missingness
  if (missingness_in_data) {
    for (iter in 1:nsample) {
      for (s in 1:q) {
        Xm.draw[[iter]][[s,1]] <- Xm.draw[[iter]][[s,1]] * sigma.mat[s,1]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Aligning the posterior samples
  # ---------------------------------------------------------------------------
  
  # If a burn-in is not given
  if (is.null(burnin)) {
    burnin <- nsample/2
  }
  
  # Save the iterations after burn-in
  iters_burnin <- seq(burnin+1, nsample)
  
  # Combining the results into a list
  BSFP.fit <-   list(data = data, # Returning the scaled version of the data
                     Y = Y, # Return the response vector
                     sigma.mat = sigma.mat, # Scaling factors
                     V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
                     VStar.draw = VStar.draw, # Components that predict Y,
                     Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
                     scores = scores, # Scores if provided by another method 
                     ranks = c(r, r.vec), # Ranks
                     tau2.draw = tau2.draw, beta.draw = beta.draw) # Regression parameters
  
  # Put Y back into the matrix-list format
  new_Y <- matrix(list(), nrow = 1, ncol = 1)
  new_Y[[1,1]] <- Y
  
  # Applying the Match Align algorithm to undo rotational invariance in the results
  aligned_results <- match_align_bpmf(BSFP.fit, y = new_Y,
                                      model_params = model_params, p.vec = p.vec, 
                                      iters_burnin = iters_burnin)
  
  # Return
  list(data = data, # Returning the scaled version of the data
       Y = Y, # Return the response vector
       sigma.mat = sigma.mat, # Scaling factors
       J.draw = J.draw, A.draw = A.draw, S.draw = S.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
       VStar.draw = VStar.draw, # Components that predict Y,
       Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
       scores = scores, # Scores if provided by another method 
       ranks = c(r, r.vec), # Ranks
       tau2.draw = tau2.draw, beta.draw = beta.draw) # Regression parameters
  
}
