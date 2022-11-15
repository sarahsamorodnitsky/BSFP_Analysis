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
  
  # Initialize the indices of features in each source 
  p.ind <- lapply(1:q, function(s) {
    if (s == 1) {
      1:p.vec[s]
    } else {
      (p.vec[s-1] + 1):cumsum(p.vec)[s]
    }
  })
  
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
    
    # Saving the rank vector
    ranks <- c(r, r.vec)
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
  
  # Save the indices for the factors from each source 
  rank.inds <- lapply(1:(q+1), function(s) {

    # For any structure,
    if (ranks[s] == 0) {
      NULL
    }
    
    # For joint structure
    else if (s == 1) {
      if (ranks[s] > 0) {
        1:ranks[s]
      }
    }
    
    # For each individual structure
    else if (s > 1) {
      (sum(ranks[1:(s-1)])+1):sum(ranks[1:s])
    }
    
  })
  
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
    if (nninit) {
      Xm0 <- matrix(list(), ncol = 1, nrow = q)
      for (s in 1:q) {
        Xm0[[s,1]] <- rank_init$X[[s,1]][missing_obs[[s]]]
      }
    }
    
    if (!nninit) {
      Xm0 <- matrix(list(), ncol = 1, nrow = q)
      for (s in 1:q) {
        Xm0[[s,1]] <- rep(0, length(missing_obs[[s]]))
      }
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
  # Calculating the joint and individual structure, scaled to the data
  # ---------------------------------------------------------------------------
  
  # Storing the structures at each Gibbs sampling iteration
  J.draw <- A.draw <- S.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  # Calculating the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  
  if (is.null(scores)) {
    for (iter in 1:nsample) {
      for (s in 1:q) {
        # Calculating the joint structure and scaling by sigma.mat
        J.draw[[iter]][[s,1]] <- (U.draw[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]])) * sigma.mat[s,1]
        
        # Calculating the individual structure and scaling by sigma.mat
        A.draw[[iter]][[s,1]] <- (W.draw[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) * sigma.mat[s,1]
        
        # Calculate the overall structure (joint + individual
        S.draw[[iter]][[s,1]] <- J.draw[[iter]][[s,1]] + A.draw[[iter]][[s,1]]
      }
      
      # Calculate the structure for Y
      if (response_given) {
        if (response_type == "continuous") {
          EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]]
        }
        
        if (response_type == "binary") {
          EY.draw[[iter]][[1,1]] <- pnorm(VStar.draw[[iter]][[1,1]] %*% beta.draw[[iter]][[1,1]])
        }
      }
    }
  }
  
  # Return
  list(data = data, # Returning the scaled version of the data
       Y = Y, # Return the response vector
       sigma.mat = sigma.mat, # Scaling factors
       J.draw = J.draw, A.draw = A.draw, S.draw = S.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.draw = U.draw, W.draw = W.draw, Vs.draw = Vs.draw, # Components of the structure
       Xm.draw = Xm.draw, Ym.draw = Ym.draw, Z.draw = Z.draw, # Missing data imputation
       scores = scores, # Scores if provided by another method 
       ranks = c(r, r.vec), # Ranks
       model_params = model_params, # Parameters used in model
       tau2.draw = tau2.draw, beta.draw = beta.draw) # Regression parameters
}


#' BSFP.Predict
#' 
#' Gibbs sampling algorithm for sampling the underlying structure and the 
#' regression coefficient vector for a response vector for unseen test data.
#' @param bpmf.fit Results from a training bpmf.fit 
#' @param Y_test Column vector with outcome on test samples or NULL
#' @param nninit Should the model be initialized with a nuclear norm penalized 
#' objective? If FALSE, provide ranks.
#' @param model_params (error_vars, joint_vars, indiv_vars, beta_vars = NULL, response_vars). 
#' May be given by BSFP.fit
#' @param nsample Number of Gibbs sampling iterations
#' @param progress Should the progress of the sampler be displayed?
#' @param starting_values list of starting values for V, U, W, Vs. If NULL and
#'  nninit = TRUE, init with UNIFAC. If NULL and nninit = FALSE, init from prior. 
#'  If not NULL, will init with provided starting values unless nninit. 

bsfp.predict <- function(bsfp.fit, test_data, Y_test, nninit = TRUE, model_params = NULL, sparsity = FALSE, nsample, progress = TRUE, starting_values = NULL) {
  
  # ---------------------------------------------------------------------------
  # Extracting the dimensions
  # ---------------------------------------------------------------------------
  
  q <- nrow(test_data) # Number of sources
  p.vec <- apply(test_data, 1, function(source) nrow(source[[1]])) # Number of features per source
  p <- sum(p.vec) # Total number of features
  n <- ncol(test_data[[1,1]]) # Number of subjects
  
  # ---------------------------------------------------------------------------
  # Extracting the model parameters
  # ---------------------------------------------------------------------------
  
  if (is.null(model_params)) {
    model_params <- bsfp.fit$model.params
  }
  
  if (!is.null(model_params)) {
    error_vars <- model_params$error_vars # Error variances
    sigma2_joint <- joint_var <- model_params$joint_var # Variance of joint structure
    sigma2_indiv <- indiv_vars <- model_params$indiv_vars # Variances of individual structure
    beta_vars <- model_params$beta_vars # Variances on betas
    response_vars <- model_params$response_vars; shape <- response_vars[1]; rate <- response_vars[2] # Hyperparameters of variance of response
  }

  # ---------------------------------------------------------------------------
  # Is there a response vector?
  # ---------------------------------------------------------------------------
  
  response_given <- !is.null(Y_test[[1,1]]) 
  
  # If so, what kind of response is it?
  if (response_given) {
    Y_test <- matrix(unlist(Y_test))
    
    response_type <- if (all(unique(Y_test) %in% c(0, 1, NA))) "binary" else "continuous"
  }
  
  # ---------------------------------------------------------------------------
  # Obtaining the ranks 
  # ---------------------------------------------------------------------------
  
  # Saving the ranks from the training data fit
  ranks <- bsfp.fit$ranks
  
  r <- ranks[1]
  r.vec <- ranks[-1]
  sigma.mat <- matrix(nrow = q, 1)
  for (s in 1:q) {
    sigma.mat[s,1] <- 1
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
  # Save the loadings, betas, and outcome variance estimated on the training data
  # ---------------------------------------------------------------------------
  
  # Save the joint loadings
  U.train <- bsfp.fit$U.draw
  
  # Save the individual loadings
  W.train <- bsfp.fit$W.draw
  
  # Save the betas
  beta.train <- bsfp.fit$beta.draw
  
  # Save the response variance
  tau2.train <- bsfp.fit$tau2.draw
  
  # ---------------------------------------------------------------------------
  # Storing the posterior samples
  # ---------------------------------------------------------------------------
  
  V.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  Vs.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = q))
  
  if (!response_given) {
    Z.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  }
  
  if (response_given) {
    Z.draw <- VStar.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1)) 
    
  }
  
  # ---------------------------------------------------------------------------
  # Initialize V, Vs
  # ---------------------------------------------------------------------------
  
  # If no starting values were provided, initialize from priors
  if (is.null(starting_values)) {
    V0 <- matrix(list(), nrow = 1, ncol = 1)
    if (r > 0) {
      V0[[1,1]] <- matrix(rnorm(n*r, mean = 0, sd = sqrt(sigma2_joint)), nrow = n, ncol = r)
    } 
    if (r == 0) {
      V0[[1,1]] <- matrix(0, nrow = n, ncol = 1)
    }
    
    Vs0 <- matrix(list(), nrow = 1, ncol = q)
    
    for (s in 1:q) {
      
      # Initialize W and V
      if (r.vec[s] > 0) {
        Vs0[[1,s]] <- matrix(rnorm(n*r.vec[s], mean = 0, sd = sqrt(sigma2_indiv[s])), nrow = n, ncol = r.vec[s])
        
      } 
      if (r.vec[s] == 0) {
        Vs0[[1,s]] <- matrix(0, nrow = n, ncol = 1)
        
      }
      
    }
  }
  
  # If starting values were provided, use them as initial values
  if (!is.null(starting_values)) {
    V0 <- starting_values$V
    Vs0 <- starting_values$Vs
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
    Z0 <- matrix(rnorm(n, mean = VStar0 %*% beta.train[[1]][[1,1]], sd = 1))
    
  }
  
  # ---------------------------------------------------------------------------
  # Storing the initial values 
  # ---------------------------------------------------------------------------
  
  V.draw[[1]] <- V0
  Vs.draw[[1]] <- Vs0
  
  if (response_given) {
    Z.draw[[1]][[1,1]] <- Z0
    VStar.draw[[1]][[1,1]] <- VStar0
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
  # Start Gibbs sampling!
  # ---------------------------------------------------------------------------
  
  for (iter in 1:(nsample-1)) {
    if (progress) svMisc::progress(iter/((nsample-1)/100))
    
    # ---------------------------------------------------------------------------
    # Storing the current values of the parameters
    # ---------------------------------------------------------------------------
    
    V.iter <- V.draw[[iter]]
    U.iter <- U.train[[iter]]
    Vs.iter <- Vs.draw[[iter]]
    W.iter <- W.train[[iter]]
    
    if (response_given) {
      # The current values of the betas
      beta.iter <- beta.train[[iter]][[1,1]] 
      
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
        tau2.iter <- tau2.train[[iter]][[1,1]]
      }
      
      Y_complete <- Y_test
      
    }
    
    X_complete <- test_data
    
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
    # Save the training U, joint loadings
    # -------------------------------------------------------------------------
    
    U.iter <- U.train[[iter+1]]
    
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
    # Save training sample for W, individual loadings
    # -------------------------------------------------------------------------
    
    # Update the current value of W
    W.iter <- W.train[[iter+1]]
    
    # -------------------------------------------------------------------------
    # Posterior sample for beta
    # -------------------------------------------------------------------------
    
    if (response_given) {
      
      # Update the current value of beta
      beta.iter <- beta.train[[iter+1]][[1,1]]
      
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
    # Save the training tau2
    # -------------------------------------------------------------------------
    
    if (response_given) {
      if (response_type == "continuous") {
        # Update the current value of tau2
        tau2.iter <- tau2.train[[iter+1]][[1,1]]
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
  }
  
  # ---------------------------------------------------------------------------
  # Calculating the joint and individual structure, scaled to the data
  # ---------------------------------------------------------------------------
  
  # Storing the structures at each Gibbs sampling iteration
  J.draw <- A.draw <- S.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = q, ncol = 1))
  
  # Calculating the structure for Y at each Gibbs sampling iteration
  EY.draw <- lapply(1:nsample, function(i) matrix(list(), nrow = 1, ncol = 1))
  
  for (iter in 1:nsample) {
    for (s in 1:q) {
      # Calculating the joint structure and scaling by sigma.mat
      J.draw[[iter]][[s,1]] <- (U.train[[iter]][[s,1]] %*% t(V.draw[[iter]][[1,1]])) * sigma.mat[s,1]
      
      # Calculating the individual structure and scaling by sigma.mat
      A.draw[[iter]][[s,1]] <- (W.train[[iter]][[s,s]] %*% t(Vs.draw[[iter]][[1,s]])) * sigma.mat[s,1]
      
      # Calculate the overall structure (joint + individual
      S.draw[[iter]][[s,1]] <- J.draw[[iter]][[s,1]] + A.draw[[iter]][[s,1]]
    }
    
    # Calculate the structure for Y
    if (response_given) {
      if (response_type == "continuous") {
        EY.draw[[iter]][[1,1]] <- VStar.draw[[iter]][[1,1]] %*% beta.train[[iter]][[1,1]]
      }
      
      if (response_type == "binary") {
        EY.draw[[iter]][[1,1]] <- pnorm(VStar.draw[[iter]][[1,1]] %*% beta.train[[iter]][[1,1]])
      }
    }
  }
  
  # Return
  list(test_data = test_data, # Returning the scaled version of the data
       Y_test = Y_test, # Return the response vector
       sigma.mat = sigma.mat, # Scaling factors
       J.draw = J.draw, A.draw = A.draw, S.draw = S.draw, EY.draw = EY.draw, # Underlying structure
       V.draw = V.draw, U.train = U.train, W.train = W.train, Vs.draw = Vs.draw, # Components of the structure
       VStar.draw = VStar.draw, # Components that predict Y,
       ranks = c(r, r.vec), # Ranks
       tau2.train = tau2.train, beta.train = beta.train) # Regression parameters
  
}


# -----------------------------------------------------------------------------
# Alignment functions
# -----------------------------------------------------------------------------

#' MatchAlign Algorithm (code adapted from infinitefactor package by Evan Poworoznek (2021))
#' 
#' Given results from BSFP, apply Varimax rotation and greedy alignment
#' using the MatchAlign algorithm to adjust for rotational, permutation,
#' and sign invariance. This function is used within the match_align_bsfp
#' function only.
#' @param lambda (list): list of loadings matrices. For BSFP, this will be combined
#'                loadings and betas. 
#' @param eta (list): list of scores matrices
#' @param var_betas (matrix): prior variance matrix for regression coefficients, if
#'                      considering an outcome
#' @param piv (matrix OR int): user-provided pivot. If NULL, use median largest singular value. 
#'                      If matrix, use as pivot. If integer, use corresponding Varimax-rotated posterior sample 
#' @param index (int): allows user to choose a different pivot for sensitivity analysis. 
#'                     0 uses the default pivot. Could be positive or negative. 
#' @export              

jointRot_multi <- function(lambda, eta, piv = NULL, y = NULL, var_betas = NULL, index = 0) {
  
  # Apply the Varimax rotation to the loadings
  vari = lapply(lambda, varimax)
  
  # Save the resulting Varimax-rotated loadings
  loads = lapply(vari, `[[`, 1)
  
  # Save the resulting rotation matrix
  rots = lapply(vari, `[[`, 2)
  
  # Apply the rotation matrix to the corresponding scores
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)
  
  # If no user-provided pivot is given, use default
  if (is.null(piv)) {
    norms = sapply(loads, norm, "2")
    piv = loads[order(norms)][[round(length(lambda)/2) + index]]
    
    # Save the posterior sample after burn-in that was used as pivot
    piv_index_to_return <- c(1:length(loads))[order(norms)][round(length(lambda)/2) + index]
  }
  
  # If user-provided pivot is given
  if (is.matrix(piv)) {
    piv_index_to_return <- NULL
  }
  
  if (!is.matrix(piv)) {
    piv = loads[[piv]]
    
    # Save the posterior sample after burn-in that was used as pivot
    piv_index_to_return <- piv
  }
  
  # Match to the defined pivot
  matches = lapply(loads, msfOUT, piv)
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)
  
  # Return
  return(list(lambda = lamout, eta = etaout, pivot_index = piv_index_to_return))
}

#' match_align_bsfp
#' 
#' Correcting for rotation, permutation, and sign invariance by applying the 
#' MatchAlign algorithm (Poworoznek et al. (2021)) to results from BSFP fit.
#' 
#' @param BSFP.fit (list): a list of results from the BPMF model fit 
#' @param y (matrix of lists): outcome vector if available
#' @param model_params (list): model parameters used in model fitting
#' @param p.vec (vector): vector with number of variables in each source
#' @param iters_burnin (vector): indices for posterior samples after burn-in
#' @param piv.list (list of matrices or ints): list of pivot matrices for joint then 
#' individual structures or indices for pivots from posterior samples
#' @param index (int): for sensitivity of the results, choose a pivot that falls
#' the some indices away from the chosen pivot (index specifies the number of indices)

match_align_bsfp <- function(BSFP.fit, y = NULL, model_params, p.vec, iters_burnin, piv.list = NULL, index = 0) {
  
  # Save the ranks from the model fit
  ranks <- BSFP.fit$ranks
  joint.rank <- ranks[1]
  indiv.ranks <- ranks[-1]
  
  # Save the number of datasets
  q <- nrow(BSFP.fit$data)
  
  # Save the number of posterior samples
  nsample <- length(BSFP.fit$V.draw)
  
  # Save the burnin
  burnin <- length(iters_burnin)
  
  # If no response is given
  if (is.null(y)) {
    joint.betas.final <- individual.betas.final <- joint_var_betas <- NULL
    indiv_var_betas <- lapply(1:q, function(s) NULL)
  }
  
  # Create a vector for the indices of each rank
  if (!is.null(y)) {

    rank.inds <- lapply(1:(q+1), function(s) {
      
      # For any structure,
      if (ranks[s] == 0) {
        NULL
      }
      
      # For joint structure
      else if (s == 1) {
        if (ranks[s] > 0) {
          1:ranks[s]
        }
      }
      
      # For each individual structure
      else if (s > 1) {
        (sum(ranks[1:(s-1)])+1):sum(ranks[1:s])
      }
      
    })
    
    beta.ind <- lapply(1:(q+1), function(s) {
      
      if (is.null(rank.inds[[s]])) {
        NULL
      }
      
      else {
        rank.inds[[s]] + 1
      }
    })
  }
  
  # Checking if a pivot was provided
  if (is.null(piv.list)) {
    piv.list <- lapply(1:(q+1), function(i) NULL)
  }
  
  # Joint structure (plus joint regression coefficients) --
  
  # Save the overall number of features
  p <- sum(p.vec)
  
  # Combine joint loadings and betas (fix betas at 0 if rank is 0)
  if (!is.null(y)) {
    
    if (r > 0) {
      joint.loadings <- lapply(iters_burnin, function(iter) {
        rbind(do.call(rbind, BSFP.fit$U.draw[[iter]]), # Joint loadings  
              t(BSFP.fit$beta.draw[[iter]][[1,1]][beta.ind[[1]],])) # Joint regression coefficients
      })
    }
    
    if (r == 0) {
      joint.loadings <- lapply(iters_burnin, function(iter) {
        rbind(do.call(rbind, BSFP.fit$U.draw[[iter]]), # Joint loadings  
              t(matrix(0, nrow = 1, ncol = 1))) # Joint regression coefficients
      })
    }

  }
  
  if (is.null(y)) {
    joint.loadings <- lapply(iters_burnin, function(iter) {
      rbind(do.call(rbind, BSFP.fit$U.draw[[iter]])) 
    })
  }
  
  # Save joint scores
  joint.scores <- lapply(iters_burnin, function(iter) BSFP.fit$V.draw[[iter]][[1,1]])
  
  # Save the prior variance on the betas for joint factors
  if (!is.null(y)) {
    joint_var_betas <- diag(rep(model_params$beta_vars[2], ranks[1]))
  }
  
  # If we estimated more than 1 factor
  if (joint.rank > 1) {
    # Apply the factor switching method to the joint structure
    joint.results.rotate <- jointRot_multi(joint.loadings, joint.scores, y = y, var_betas = joint_var_betas, index = index, piv = piv.list[[1]])
    
    # Separate the joint loadings from the joint scores
    joint.loadings.final <- lapply(joint.results.rotate$lambda, function(iter) iter[1:p,])
    
    if (!is.null(y)) {
      joint.betas.final <- lapply(joint.results.rotate$lambda, function(iter) t(iter[p+1,,drop=FALSE]))
    }
    
    # Applying the permutation and sign switching to the scores
    joint.scores.final <- joint.results.rotate$eta
    
    # Save the pivot index
    joint_pivot_index <- joint.results.rotate$pivot_index
  }
  
  # If we estimated 1 or fewer factors
  if (joint.rank <= 1) {
    joint.scores.final <- joint.scores
    joint.loadings.final <- lapply(joint.loadings, function(iter) iter[1:p,,drop=FALSE])
    
    if (!is.null(y)) {
      joint.betas.final <- lapply(joint.loadings, function(iter) t(iter[p+1,,drop=FALSE]))
    }
    
    joint_pivot_index <- NULL
  }

  # Individual structure (plus individual regression coefficients) --
  
  # Combine the individual loadings and betas
  if (!is.null(y)) {
    individual.loadings <- lapply(1:q, function(s) {
      
      if (indiv.ranks[s] > 0) {
        lapply(iters_burnin, function(iter) {
          rbind(BSFP.fit$W.draw[[iter]][[s,s]], 
                BSFP.fit$beta.draw[[iter]][[1,1]][beta.ind[[s+1]],])
        })
      } else {
        lapply(iters_burnin, function(iter) {
          rbind(BSFP.fit$W.draw[[iter]][[s,s]], 
                matrix(0, nrow = 1, ncol = 1))
        })
      }

    })
  }
  
  if (is.null(y)) {
    individual.loadings <- lapply(1:q, function(s) lapply(iters_burnin, function(iter) {
      rbind(BSFP.fit$W.draw[[iter]][[s,s]])
    }))
  }
  
  individual.scores <- lapply(1:q, function(s) lapply(iters_burnin, function(iter) BSFP.fit$Vs.draw[[iter]][[1,s]]))
  
  # Save the prior variance on the betas for the individual factors
  if (!is.null(y)) {
    indiv_var_betas <- lapply(1:q, function(s) diag(rep(model_params$beta_vars[s+2], ranks[s+1])))
  }
  
  # Applying the alignment algorithm
  individual.results.rotate <- lapply(1:q, function(s) {
    
    if (indiv.ranks[s] > 1) {
      out <- jointRot_multi(individual.loadings[[s]], individual.scores[[s]], y = y, var_betas = indiv_var_betas[[s]], index = index, piv = piv.list[[s+1]])
      out
    } else {
      list(lambda = individual.loadings[[s]], eta = individual.scores[[s]])
    }
  })
  
  # Save the final individual loadings and betas
  individual.loadings.final <- lapply(1:q, function(s) lapply(individual.results.rotate[[s]]$lambda, function(iter) iter[1:p.vec[s],,drop=FALSE]))
  
  if (!is.null(y)) {
    individual.betas.final <- lapply(1:q, function(s) lapply(individual.results.rotate[[s]]$lambda, function(iter) t(iter[p.vec[s]+1,,drop=FALSE])))
  }
  
  # Applying the permutation and sign switching to the scores
  individual.scores.final <- lapply(1:q, function(s) individual.results.rotate[[s]]$eta)
  
  # Save the pivot index
  individual_pivot_index <- lapply(1:q, function(s) {
    
    if (indiv.ranks[s] > 1) {
      individual.results.rotate[[s]]$pivot_index
    }

  })
  
  # ---------------------------------------------------------------------------
  # Reordering aligned results so that factors are ordered from most-to-least
  # variance explained. 
  # ---------------------------------------------------------------------------
  
  # Order the components by squared Frobenius-norm of the corresponding rank-1 structure after burn-in --
  
  # For each component, calculate the posterior mean of the corresponding rank-1 structure
  if (!is.null(y)) {
    
    if (joint.rank > 0) {
      joint.rank1.structure <- lapply(1:joint.rank, function(r) {
        # Calculate the rank-1 structure for the given component at each iteration
        factor_r_structure <- lapply(1:burnin, function(iter) {
          rbind(joint.loadings.final[[iter]][,r,drop=FALSE], t(joint.betas.final[[iter]][r,,drop=FALSE])) %*% 
            t(joint.scores.final[[iter]][,r,drop=FALSE])
        })
        
        # Calculate the posterior mean
        Reduce("+", factor_r_structure)/length(factor_r_structure)
      })
    }
    
  }
  
  if (is.null(y)) {
    
    if (joint.rank > 0) {
      joint.rank1.structure <- lapply(1:joint.rank, function(r) {
        # Calculate the rank-1 structure for the given component at each iteration
        factor_r_structure <- lapply(1:burnin, function(iter) {
          joint.loadings.final[[iter]][,r,drop=FALSE] %*% 
            t(joint.scores.final[[iter]][,r,drop=FALSE])
        })
        
        # Calculate the posterior mean
        Reduce("+", factor_r_structure)/length(factor_r_structure)
      })
    }
    
  }
  
  # Calculate the norm of each rank-1 structure
  if (joint.rank > 0) {
    joint.structure.norm <- sapply(joint.rank1.structure, function(str) frob(str))
  }
  
  if (joint.rank == 0) {
    joint.structure.norm <- 1
  }
  
  # Order the factors
  joint.factor.order <- order(joint.structure.norm, decreasing = TRUE)
  
  # Reorder the joint scores, loadings, and betas after burn-in
  joint.scores.final.order <- lapply(1:burnin, function(iter) {
    joint.scores.final[[iter]][,joint.factor.order,drop=FALSE]
  })
  joint.loadings.final.order <- lapply(1:burnin, function(iter) {
    joint.loadings.final[[iter]][,joint.factor.order,drop=FALSE]
  })
  
  if (!is.null(y)) {
    joint.betas.final.order <- lapply(1:burnin, function(iter) {
      joint.betas.final[[iter]][joint.factor.order,,drop=FALSE]
    })
  }
  
  # Order the factors in each individual structure by the Frobenius norm of their corresponding rank-1 structure --
  
  # For each source and each component, calculate the posterior mean of the corresponding rank-1 structure
  if (!is.null(y)) {
    
    individual.rank1.structure <- lapply(1:q, function(s) {
      
      if (indiv.ranks[s] > 0) {
        lapply(1:indiv.ranks[s], function(r) {
          # Calculate the rank-1 structure for the given component at each iteration
          factor_r_structure <- lapply(1:burnin, function(iter) {
            rbind(individual.loadings.final[[s]][[iter]][,r,drop=FALSE], t(individual.betas.final[[s]][[iter]][r,,drop=FALSE])) %*% 
              t(individual.scores.final[[s]][[iter]][,r,drop=FALSE])
          })
          
          # Calculate the posterior mean
          Reduce("+", factor_r_structure)/length(factor_r_structure)
        })
      } else {
        NULL
      }

    })
  }
  
  if (is.null(y)) {
    
    individual.rank1.structure <- lapply(1:q, function(s) {
      
      if (indiv.ranks[s] > 0) {
        lapply(1:indiv.ranks[s], function(r) {
          # Calculate the rank-1 structure for the given component at each iteration
          factor_r_structure <- lapply(1:burnin, function(iter) {
            individual.loadings.final[[s]][[iter]][,r,drop=FALSE] %*% 
              t(individual.scores.final[[s]][[iter]][,r,drop=FALSE])
          })
          
          # Calculate the posterior mean
          Reduce("+", factor_r_structure)/length(factor_r_structure)
        })
      } else {
        NULL
      }

    })
  }
  
  # Calculate the norm of each rank-1 structure for each source
  individual.structure.norm <- lapply(1:q, function(s) {
    if (indiv.ranks[s] > 0) {
      sapply(individual.rank1.structure[[s]], function(str) frob(str))
    } else {
      1
    }
  })
  
  # Order the factors
  individual.factor.order <- lapply(1:q, function(s) {
    order(individual.structure.norm[[s]], decreasing = TRUE)
  })
  
  # Reorder the individual scores, loadings, and betas after burn-in
  individual.scores.final.order <- lapply(1:q, function(s) {
    lapply(1:burnin, function(iter) {
      individual.scores.final[[s]][[iter]][,individual.factor.order[[s]],drop=FALSE]
    })
  })
  individual.loadings.final.order <- lapply(1:q, function(s) {
    lapply(1:burnin, function(iter) {
      individual.loadings.final[[s]][[iter]][,individual.factor.order[[s]],drop=FALSE]
    })
  })
  
  if (!is.null(y)) {
    individual.betas.final.order <- lapply(1:q, function(s) {
      lapply(1:burnin, function(iter) {
        individual.betas.final[[s]][[iter]][individual.factor.order[[s]],,drop=FALSE]
      })
    })
  }
  
  # If no response was given
  if (is.null(y)) {
    joint.betas.final.order <- individual.final.burnin <- NULL
  }
  
  # Return the final loadings, scores, and betas
  list(joint.scores.final = joint.scores.final.order, 
       joint.loadings.final = joint.loadings.final.order, 
       joint.betas.final = joint.betas.final.order, 
       joint_pivot_index = joint_pivot_index,
       individual.scores.final = individual.scores.final.order, 
       individual.loadings.final = individual.loadings.final.order, 
       individual.betas.final = individual.betas.final.order,
       individual_pivot_index = individual_pivot_index)
}

# -----------------------------------------------------------------------------
# Summary functions
# -----------------------------------------------------------------------------

#' summarize_factors
#' 
#' Calculate posterior summaries of aligned estimated factors from BSFP. 
#' 
#' @param data A matrix of lists or a list of matrices that share the same number of
#' columns. The matrices must be oriented in pxn orientation. May contain NAs if 
#' there are missing values in the dataset. 
#' @param Y A matrix of lists or a nx1 matrix of continuous or binary outcome. 
#' May be NULL if no outcome is given. May contain NAs if there are missing outcomes. 
#' @param iters_burnin (vector): indices for posterior samples after burn-in
#' @param aligned_results (list): results from match_align_bsfp
#' @param ranks Estimated joint and individual ranks from BSFP
#' @param tau2.draw Posterior samples for error variance in Y if response is given
#' @param Xm.draw Imputed values for missing values in X
#' @param Ym.draw Imputed values for unobserved outcomes in Y
#' @export

summarize_factors <- function(data, Y = NULL, iters_burnin,
                              aligned_results, ranks, tau2.draw = NULL, Xm.draw = NULL, Ym.draw = NULL) {
  
  # ---------------------------------------------------------------------------
  # Save the aligned results
  # ---------------------------------------------------------------------------
  
  joint.scores.final <- aligned_results$joint.scores.final
  joint.loadings.final <- aligned_results$joint.loadings.final
  joint.betas.final <- aligned_results$joint.betas.final
  
  individual.scores.final <- aligned_results$individual.scores.final
  individual.loadings.final <- aligned_results$individual.loadings.final
  individual.betas.final <- aligned_results$individual.betas.final
  
  # ---------------------------------------------------------------------------
  # Save data attributes
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
  
  # Initialize the indices of features in each source 
  p.ind <- lapply(1:q, function(s) {
    if (s == 1) {
      1:p.vec[s]
    } else {
      (p.vec[s-1] + 1):cumsum(p.vec)[s]
    }
  })
  
  # Rank indices
  rank.inds <- lapply(1:(q+1), function(s) {
    
    # For any structure,
    if (ranks[s] == 0) {
      NULL
    }
    
    # For joint structure
    else if (s == 1) {
      if (ranks[s] > 0) {
        1:ranks[s]
      }
    }
    
    # For each individual structure
    else if (s > 1) {
      (sum(ranks[1:(s-1)])+1):sum(ranks[1:s])
    }
    
  })
  
  # Regression coefficient indices
  beta.ind <- lapply(1:(q+1), function(s) {
    
    if (is.null(rank.inds[[s]])) {
      NULL
    }
    
    else {
      rank.inds[[s]] + 1
    }
  })
  
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
  # Check for missingness in data
  # ---------------------------------------------------------------------------
  
  # Check for missingness
  missingness_in_data <- any(sapply(data[,1], function(source) any(is.na(source))))
  
  # Which entries are missing?
  missing_obs <- lapply(data[,1], function(source) which(is.na(source)))
  
  # ---------------------------------------------------------------------------
  # Joint factors
  # ---------------------------------------------------------------------------
  
  if (ranks[1] > 0) {
    # Scores
    joint.scores.summary <- lapply(1:ranks[1], function(r) {
      score.matrix <- do.call(cbind, lapply(joint.scores.final, function(iter) iter[,r,drop=FALSE]))
      summary.mat <- cbind.data.frame(Post.Mean = rowMeans(score.matrix),
                                      Lower.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.025)),
                                      Upper.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.975)))
      rownames(summary.mat) <- paste("Sample", 1:n)
      summary.mat
    })
    names(joint.scores.summary) <- paste0("Joint.Factor.", 1:ranks[1])
    
    # Loadings
    joint.loadings.summary <- lapply(1:ranks[1], function(r) {
      lapply(1:q, function(s) {
        load.matrix <- do.call(cbind, lapply(joint.loadings.final, function(iter) iter[p.ind[[s]],r,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(load.matrix),
                                        Lower.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.975)))
        rownames(summary.mat) <- paste("Source", s, "Feature", 1:p.vec[s])
        summary.mat
      })
    })
    
    # Regression coefficients
    if (response_given) {
      joint.betas.summary <- lapply(1:ranks[1], function(r) {
        beta.mat <- do.call(cbind, lapply(joint.betas.final, function(iter) iter[r,,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(beta.mat),
                                        Lower.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.975)))
        summary.mat
      })
      names(joint.betas.summary) <- paste0("Joint.Factor.Regression.Coefficient.", 1:ranks[1])
    }
    
    if (!response_given) {
      joint.betas.summary <- NULL
    }
  }
  
  if (ranks[1] == 0) {
    # Scores
    score.matrix <- do.call(cbind, lapply(joint.scores.final, function(iter) iter[,1,drop=FALSE]))
    joint.scores.summary <- cbind.data.frame(Post.Mean = rowMeans(score.matrix),
                                             Lower.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.025)),
                                             Upper.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.975)))
    rownames(joint.scores.summary) <- paste("Sample", 1:n)
    names(joint.scores.summary) <- paste0("Joint.Factor.", 1:ranks[1])
    
    # Loadings
    joint.loadings.summary <- lapply(1:q, function(s) {
      load.matrix <- do.call(cbind, lapply(joint.loadings.final, function(iter) iter[p.ind[[s]],1,drop=FALSE]))
      summary.mat <- cbind.data.frame(Post.Mean = rowMeans(load.matrix),
                                      Lower.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.025)),
                                      Upper.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.975)))
      rownames(summary.mat) <- paste("Source", s, "Feature", 1:p.vec[s])
      summary.mat
    })
    
    # Regression coefficients
    if (response_given) {
      beta.mat <- do.call(cbind, lapply(joint.betas.final, function(iter) iter[1,,drop=FALSE]))
      joint.betas.summary <- cbind.data.frame(Post.Mean = rowMeans(beta.mat),
                                              Lower.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.025)),
                                              Upper.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.975)))
      joint.betas.summary
    }
    
    if (!response_given) {
      joint.betas.summary <- NULL
    }
    
  }
  
  # ---------------------------------------------------------------------------
  # Individual factors
  # ---------------------------------------------------------------------------
  
  # Scores
  individual.scores.summary <- lapply(1:q, function(s) {
    
    if (ranks[s+1] > 0) {
      source.scores.list <- lapply(1:ranks[s+1], function(rs) {
        score.matrix <- do.call(cbind, lapply(individual.scores.final[[s]], function(iter) iter[,rs,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(score.matrix),
                                        Lower.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.975)))
        rownames(summary.mat) <- paste("Sample", 1:n)
        summary.mat
      })
      names(source.scores.list) <- paste0("Source.", s, ".Individual.Factor.", 1:ranks[s+1])
      source.scores.list
    } else {
      score.matrix <- do.call(cbind, lapply(individual.scores.final[[s]], function(iter) iter[,1,drop=FALSE]))
      summary.mat <- cbind.data.frame(Post.Mean = rowMeans(score.matrix),
                                      Lower.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.025)),
                                      Upper.95.CI = apply(score.matrix, 1, function(row) quantile(row, 0.975)))
      rownames(summary.mat) <- paste("Sample", 1:n)
      summary.mat
    }
    
  })
  names(individual.scores.summary) <- paste0("Source.", 1:q)
  
  # Loadings
  individual.loadings.summary <- lapply(1:q, function(s) {
    
    if (ranks[s+1] > 0) {
      source.load.list <- lapply(1:ranks[s+1], function(rs) {
        load.matrix <- do.call(cbind, lapply(individual.loadings.final[[s]], function(iter) iter[,rs,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(load.matrix),
                                        Lower.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.975)))
        rownames(summary.mat) <- paste("Source", s, "Feature", 1:p.vec[s])
        summary.mat
      })
      names(source.load.list) <- paste0("Source.", s, ".Individual.Factor.", 1:ranks[s+1])
      source.load.list
    } else {
      load.matrix <- do.call(cbind, lapply(individual.loadings.final[[s]], function(iter) iter[,1,drop=FALSE]))
      summary.mat <- cbind.data.frame(Post.Mean = rowMeans(load.matrix),
                                      Lower.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.025)),
                                      Upper.95.CI = apply(load.matrix, 1, function(row) quantile(row, 0.975)))
      rownames(summary.mat) <- paste("Source", s, "Feature", 1:p.vec[s])
      summary.mat
    }
    
  })
  names(individual.loadings.summary) <- paste0("Source.", 1:q)
  
  # Regression coefficients
  if (response_given) {
    individual.betas.summary <- lapply(1:q, function(s) {
      
      if (ranks[s+1] > 0) {
        source.beta.list <- lapply(1:ranks[s+1], function(rs) {
          beta.mat <- do.call(cbind, lapply(individual.betas.final[[s]], function(iter) iter[rs,,drop=FALSE]))
          summary.mat <- cbind.data.frame(Post.Mean = rowMeans(beta.mat),
                                          Lower.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.025)),
                                          Upper.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.975)))
          summary.mat
        })
        names(source.beta.list) <- paste0("Source.", s, ".Individual.Factor.Regression.Coefficient.", 1:ranks[s+1])
        source.beta.list
      } else {
        beta.mat <- do.call(cbind, lapply(individual.betas.final[[s]], function(iter) iter[1,,drop=FALSE]))
        summary.mat <- cbind.data.frame(Post.Mean = rowMeans(beta.mat),
                                        Lower.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.025)),
                                        Upper.95.CI = apply(beta.mat, 1, function(row) quantile(row, 0.975)))
        summary.mat
      }
      
    })
    names(individual.betas.summary) <- paste0("Source.", 1:q)
  }
  
  if (!response_given) {
    individual.betas.summary <- NULL 
  }
  
  # ---------------------------------------------------------------------------
  # Summarizing the outcome variance
  # ---------------------------------------------------------------------------
  
  if (response_given) {
    tau2.summary <- cbind.data.frame(Post.Mean = mean(unlist(tau2.draw)[iters_burnin]),
                                     Lower.95.CI = quantile(unlist(tau2.draw)[iters_burnin], 0.025),
                                     Upper.95.CI = quantile(unlist(tau2.draw)[iters_burnin], 0.975))
  }
  
  if (!response_given) {
    tau2.summary <- NULL
  }
  
  # ---------------------------------------------------------------------------
  # Summarizing the imputed values
  # ---------------------------------------------------------------------------
  
  # For each source
  Xm.summary <- lapply(1:q, function(s) {
    
    if (length(missing_obs[[s]]) > 0) {
      # Create matrix of imputed values for each missing sample
      Xm.s <- do.call(cbind, lapply(iters_burnin, function(iter) {
        Xm.draw[[iter]][[s,1]]
      }))
      
      # Summarize the posterior mean and 95% CI
      Xm.s.summary <- cbind.data.frame(Post.Mean = rowMeans(Xm.s),
                                       Lower.95.CI = apply(Xm.s, 1, function(row) quantile(row, 0.025)),
                                       Upper.95.CI = apply(Xm.s, 1, function(row) quantile(row, 0.975)))
      
      # Add rownames
      rownames(Xm.s.summary) <- paste0("Source ", s, " Sample ", missing_obs[[s]])
      Xm.s.summary
    } else {
      NULL
    }
  })
  
  # For the response
  if (response_given) {
    if (missingness_in_response) {
      Ym.mat <- do.call(cbind, lapply(iters_burnin, function(iter) Ym.draw[[iter]][[1,1]]))
      Ym.summary <- cbind.data.frame(Post.Mean = rowMeans(Ym.mat),
                                     Lower.95.CI = apply(Ym.mat, 1, function(row) quantile(row, 0.025)),
                                     Upper.95.CI = apply(Ym.mat, 1, function(row) quantile(row, 0.975)))
      rownames(Ym.summary) <- paste0("Sample ", missing_obs_Y)
    }
    
    if (!missingness_in_response) {
      Ym.summary <- NULL
    }
  }
  
  if (!response_given) {
    Ym.summary <- NULL
  }
  
  
  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  
  list(
    # Summaries of the joint factors
    joint.scores.summary = joint.scores.summary,
    joint.loadings.summary = joint.loadings.summary,
    
    # Summaries of the individual factors
    individual.scores.summary = individual.scores.summary,
    individual.loadings.summary = individual.loadings.summary,
    
    # Summaries of the factor contributions in predictive model
    joint.betas.summary = joint.betas.summary,
    individual.betas.summary = individual.betas.summary,
    
    # Summaries of the imputed values
    Xm.summary = Xm.summary, Ym.summary = Ym.summary, # Missing data imputation
    
    # Miscellaneous
    ranks = c(r, r.vec), # Ranks
    tau2.summary = tau2.summary) # Regression parameters
  
}

#' log_joint_density: 
#' 
#' Calculate the log-joint density of the estimated model at a given posterior
#' sampling iteration to facilitate checking convergence. 
#' 
#' @param data A matrix of lists or a list of matrices that share the same number of
#' columns. The matrices must be oriented in pxn orientation. May contain NAs if 
#' there are missing values in the dataset. 
#' @param Y A matrix of lists or a nx1 matrix of continuous or binary outcome. 
#' May be NULL if no outcome is given. May contain NAs if there are missing outcomes. 
#' @param U.iter Posterior draw for joint loadings, U, at a given iteration
#' @param V.iter Posterior draw for joint scores, V, at a given iteration
#' @param W.iter Posterior draw for individual loadings, W, at a given iteration
#' @param Vs.iter Posterior draw for individual scores, Vs, at a given iteration
#' @param ranks Estimated joint and individual ranks from BSFP
#' @param beta.iter Posterior draw for regression coefficients at a given iteration
#' @param tau2.iter Posterior draw for estimated error variance in Y at a given iteration
#' @param Xm.iter Imputed values for unobserved X values at a given iteration
#' @param Ym.iter Imputed values for unobserved Y values at a given iteration
#' @export

log_joint_density <- function(data, Y = NULL, U.iter, V.iter, W.iter, Vs.iter, model_params, ranks, beta.iter = NULL, tau2.iter = NULL, Xm.iter = NULL, Ym.iter = NULL) {
  
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
      like <- like + sum(sapply(1:n_beta, function(rs) {
        dnorm(beta.iter[[1,1]][rs,], mean = 0, sd = sqrt(Sigma_beta[rs,rs]), log = TRUE)
      }))

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
    }
  }
  
  # Return
  like
}

# -----------------------------------------------------------------------------
# Miscellaneous
# -----------------------------------------------------------------------------

#' bsfp_data
#' 
#' Generate fake data according to BSF/BSFP model
#' 
#' @param p.vec (vector) vector of with number of features per source
#' @param n (int) number of samples
#' @param ranks (vector) vector with joint and individual ranks (first index 
#' is joint, remaining are individual for each source)
#' @param s2nX (dbl) signal-to-noise ratio in X (ratio of variance in X to variance
#' in error). >1 means higher signal in X, <1 means lower signal in X
#' @param s2nY (dbl) signal-to-noise ratio in Y
#' @param response NULL if no response is desired, "continuous", or "binary" 
#' @param missingness NULL if no missingness, or "missingness_in_data" or "missingness_in_response"
#' @param missing_data_type NULL if missingness = NULL; otherwise, = entrywise 
#' if randomly sample across sources; columnwise if randomly set samples to missing;
#'  MNAR is randomly set features below a threshold to missing
#' @param prop_missing NULL if missingness = NULL; otherwise the proportion of 
#' entries or columns set to missing
#' @param sparsity (Boolean) TRUE if generate regression coefficients under 
#' spike-and-slab prior, FALSE otherwise
#' @param identically_zero (Boolean) TRUE if generating response with sparsity 
#' and want spike coefficients to be exactly 0
#' @param num_in_spike (vector) number of coefficients to be generated from spike 
#' if generating response with sparsity. NULL if sparsity=FALSE. Should have an 
#' integer for joint and each individual structure. May be set to NULL and 
#' the number of factors in the spike will be determined randomly. 
#' @export

bsfp_data <- function(p.vec, n, ranks, true_params, s2nX = NULL, s2nY = NULL, response = NULL, missingness = NULL, missing_data_type = NULL, prop_missing = NULL, sparsity, identically_zero = FALSE, num_in_spike = NULL) {
  
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
  response_vars <- true_params$response_vars
  
  if (!is.null(response_vars)) {
    shape <- response_vars[1]; rate <- response_vars[2]
  }
  
  # -------------------------------------------------------------------------
  # Generating the underlying structure
  # -------------------------------------------------------------------------
  
  joint.structure <- indiv.structure <- overall.structure <- matrix(list(), ncol = 1, nrow = q)
  
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
    overall.structure[[s,1]] <- joint.structure[[s,1]] + indiv.structure[[s,1]]
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
      overall.structure[[s,1]] <- joint.structure[[s,1]] + indiv.structure[[s,1]]
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
      
      # If a true error variance is provided, set tau2 to this value
      if (length(error_vars) > q) {
        tau2 <- matrix(list(), nrow = 1, ncol = 1)
        tau2[[1,1]] <- error_vars[q+1]
      }
      
      # If a true error variance is not provided, use prior 
      if (length(error_vars) == q) {
        tau2 <- matrix(list(), nrow = 1, ncol = 1)
        tau2[[1,1]] <- matrix(1/rgamma(1, shape = shape, rate = rate)) 
      }
      
      Y <- EY <- matrix(list(), nrow = 1, ncol = 1)
      EY[[1,1]] <- VStar %*% beta[[1,1]]
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
      
      missing_obs <- missing_data <- missing_cols <- missing_rows <- matrix(list(), nrow = q, ncol = 1)
      
      if (missing_data_type == "entrywise") { # if removing observations entrywise
        for (s in 1:q) {
          # these are counters going down the columns of R. So 9 would be the 9th entry counting down. 
          missing_obs[[s,1]] <- sort(sample(x = 1:length(data[[s,1]]), size = prop_missing*length(data[[s,1]]), replace = FALSE))
          
          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][missing_obs[[s,1]]] <- NA
        }
      } 
      
      if (missing_data_type == "columnwise") { # if removing entire columns
        # Gives the column indices to remove
        
        for (s in 1:q) {
          # These are counters going down the COLUMNS of X. So 9 would be the 9th column.
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
      
      if (missing_data_type == "MNAR") { # Restrict missing to within features
        for (s in 1:q) {
          # Sort the entries in source s
          sorted_entries_s <- sort(c(data[[s,1]]), decreasing = FALSE)
          
          # Calculate the threshold at which the bottom prop_missing are set to NA
          number_below_lod <- prop_missing * length(sorted_entries_s)
          LOD <- ceiling(max(sorted_entries_s[1:number_below_lod]))
          
          # Duplicate Xs so that I have one with the full data and one with the missing data
          missing_data[[s,1]] <- data[[s,1]]
          missing_data[[s,1]][missing_data[[s,1]] < LOD] <- NA
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
       overall.structure = overall.structure, # Joint + Individual structure
       V = V, U = U, Vs = Vs, W = W, # Components of the structure
       beta = beta, EY = EY, tau2 = tau2, gamma = gamma, p.prior = p.prior)
}

# -----------------------------------------------------------------------------
# UNIFAC/BIDIFAC functions (Credit to Jun Young Park (2020))
# -----------------------------------------------------------------------------

BIDIFAC=function(data,rmt=T, sigma=NULL,
                 start=NULL, out=FALSE,
                 eps=1e-3, max.iter=1000, pbar=TRUE, seed=NULL, scale_back = TRUE, ...){
  if (!is.null(seed)){set.seed(seed)}
  # if (!rmt & class(sigma)!="matrix") stop("sigma must be a matrix.")
  if (!rmt & !("matrix" %in% class(sigma))) stop("sigma must be a matrix.")
  
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

impute.BIDIFAC=function(data, 
                        rmt=T, sigma=NULL,
                        pbar=TRUE,
                        start=NULL, max.iter.impute=100, 
                        eps.impute=1e-3, scale_back=TRUE,...){
  dim.data=dim(data)
  p=dim.data[1]; q=dim.data[2]
  
  dim.list=do.call(cbind,lapply(data, dim))
  mvec=apply(matrix(dim.list[1,],p),1,unique)
  nvec=apply(matrix(dim.list[2,],p),2,unique)
  if (class(mvec)=="list" ) stop("the number of rows do not match")
  if (class(nvec)=="list" ) stop("the number of columns do not match")
  
  impute.index=matrix(list(), nrow = p, ncol=q)
  if (is.null(sigma)) sigma=matrix(1,p,q)
  
  for (i in 1:p){
    for (j in 1:q){
      
      if (any(is.na(data[[i,j]]))) {
        fillmat=fill.matrix(data[[i,j]])
        impute.index[[i,j]]=fillmat$na.ind
        if (rmt) sigma[i,j]=sigma.rmt(fillmat$X.fill)
        data[[i,j]]=fillmat$X.fill/sigma[i,j]
      } else {
        data[[i,j]] <- data[[i,j]]/sigma[i,j]
      }

    }
  }
  
  bool2=TRUE
  count2=1; impute.vec=0
  if (pbar) pb = txtProgressBar(min = 0, max=max.iter.impute, initial=0, char="-", style = 3)
  while (bool2){
    if (pbar){  setTxtProgressBar(pb, count2)  }
    impute.vec.old=impute.vec
    fit=BIDIFAC(data, rmt=F, sigma=matrix(1,p,q),start=start, pbar = F)
    
    start=list(
      data.rearrange(fit$G)$out, data.rearrange(fit$R)$out,
      data.rearrange(fit$C)$out, data.rearrange(fit$I)$out)
    
    impute.vec=NULL
    for (i in 1:p){
      for (j in 1:q){
        imp=fit$S[[i,j]][impute.index[[i,j]]]
        data[[i,j]][impute.index[[i,j]]]=imp
        impute.vec=c(impute.vec,imp)
      }
    }
    
    crit2=frob(impute.vec.old-impute.vec)/length(impute.vec)
    if (crit2<eps.impute){ bool2=FALSE }
    else if (count2==max.iter.impute){ 
      bool2=FALSE 
      cat("\n The algorithm did not converge. \n")
      cat("Try to increase max.iter.impute or relax eps/eps.impute.")
    }
    else{ count2=count2+1 }
  }
  
  if (scale_back) {
    for (i in 1:p){
      for (j in 1:q){
        fit$X[[i,j]]=fit$X[[i,j]]*sigma[i,j]
        fit$S[[i,j]]=fit$S[[i,j]]*sigma[i,j]
        fit$G[[i,j]]=fit$G[[i,j]]*sigma[i,j]
        fit$R[[i,j]]=fit$R[[i,j]]*sigma[i,j]
        fit$C[[i,j]]=fit$C[[i,j]]*sigma[i,j]
        fit$I[[i,j]]=fit$I[[i,j]]*sigma[i,j]
      }
    }
  }
  
  if (!scale_back) {
    for (i in 1:p){
      for (j in 1:q){
        fit$X[[i,j]]=fit$X[[i,j]]#*sigma[i,j]
        fit$S[[i,j]]=fit$S[[i,j]]#*sigma[i,j]
        fit$G[[i,j]]=fit$G[[i,j]]#*sigma[i,j]
        fit$R[[i,j]]=fit$R[[i,j]]#*sigma[i,j]
        fit$C[[i,j]]=fit$C[[i,j]]#*sigma[i,j]
        fit$I[[i,j]]=fit$I[[i,j]]#*sigma[i,j]
      }
    }
  }
  
  fit$sigma.mat=sigma
  
  return(fit)
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

fill.matrix=function(X){
  na.ind.original=na.ind=which(is.na(X),arr.ind = T)
  bool=T
  while (bool){
    impute.X=rep(NA, nrow(na.ind))
    for (j in 1:nrow(na.ind)){
      colmean=mean(X[,na.ind[j,2]], na.rm=T)
      rowmean=mean(X[na.ind[j,1],], na.rm=T)
      impute.X[j]=mean(c(colmean,rowmean), na.rm=T)
    }
    X[na.ind]=impute.X
    
    na.ind=which(is.na(X),arr.ind = T)
    if (length(na.ind)==0) bool=F
  }
  
  return(list(X.fill=X, na.ind=na.ind.original))
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
