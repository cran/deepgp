
# Function Contents -----------------------------------------------------------
# Internal:
#   krig: prediction/posterior sampling using typical kriging equations
#   krig_vec: prediction/posterior sampling using vecchia approximation

# krig ------------------------------------------------------------------------

krig <- function(y, 
                 xdmat = NULL, xdmat_new = NULL, xdmat_cross = NULL, # used for sep = FALSE
                 x = NULL, x_new = NULL, # used for sep = TRUE or grad (predict or enhance)
                 tau2 = 1, theta, g = 0, v = 2.5, 
                 sep = FALSE, # use separable lengthscales?  If so, theta must match dimension
                 s2 = FALSE, # calculate and return point-wise variances?
                 sigma = FALSE, # calculate and return full covariance?
                 grad = FALSE, # include predictions of gradients?  Inherits mean/s2, requires sigma = FALSE
                 grad_enhance = FALSE,
                 nsamples = 0, # draw full samples from posterior?  If so, how many?
                 prior_mean = 0,
                 prior_mean_new = 0) {
  
  out <- list()
  
  # Check settings
  if (s2 & sigma) stop("only s2 or sigma should be true, not both")
  if ((grad_enhance | grad) & v != 999) stop("grad options require 'exp2' kernel")
  if (sep | grad_enhance | grad)
    if (is.null(x) | is.null(x_new))
      stop("x and x_new arguments are required")
  if (nsamples > 0) {
    if (s2) stop("posterior sampling uses s2 = FALSE and sigma = TRUE")
    s2 <- FALSE # NOT required for posterior sampling
    sigma <- TRUE # required for posterior sampling
  }
  if (grad_enhance | grad) { # need to get grad_indices
    n <- nrow(x)
    d <- ncol(x)
    n_new <- nrow(x_new)
    if (grad_enhance) {
      x <- bind(x, d)
      grad_indx <- rep(0:d, each = n)
    } else grad_indx <- rep(0, times = n)
    if (grad) {
      yi <- 1:n_new
      dyi <- (n_new + 1):(n_new*(d+1))
      x_new <- bind(x_new, d)
      grad_indx_new <- rep(0:d, each = n_new)
    } else grad_indx_new <- rep(0, times = n_new)
  } 

  # Pre-calculate C and C_cross matrices
  if (v == 999) { 
    if (sep) {
      if (grad_enhance | grad) {
        C <- Exp2SepGrad(x, x, grad_indx, grad_indx, 1, theta, g)
        C_cross <- Exp2SepGrad(x_new, x, grad_indx_new, grad_indx, 1, theta, 0) # no g in rectangular matrix
      } else {
        C <- Exp2Sep(x, x, 1, theta, g)
        C_cross <- Exp2Sep(x_new, x, 1, theta, 0) # no g in rectangular matrix
      }
    } else {
      if (grad_enhance | grad) {
        C <- Exp2Grad(x, x, grad_indx, grad_indx, 1, theta, g)
        C_cross <- Exp2Grad(x_new, x, grad_indx_new, grad_indx, 1, theta, 0) # no g in rectangular matrix
      } else {
        C <- Exp2(xdmat, 1, theta, g)
        C_cross <- Exp2(xdmat_cross, 1, theta, 0) # no g in rectangular matrix
      }
    }
  } else { # no grad or grad_enhance options
    if (sep) {
      C <- MaternSep(x, x, 1, theta, g, v) 
      C_cross <- MaternSep(x_new, x, 1, theta, 0, v)
    } else {
      C <- Matern(xdmat, 1, theta, g, v) 
      C_cross <- Matern(xdmat_cross, 1, theta, 0, v)
    }
  }

  # Posterior mean
  C_inv <- invdet(C)$Mi
  C_inv_y <- C_inv %*% (y - prior_mean)
  mu <- prior_mean_new + C_cross %*% C_inv_y 
  if (grad) {
    out$mean <- mu[yi]
    out$grad_mean <- matrix(mu[dyi], ncol = d, byrow = FALSE)
  } else out$mean <- mu
  
  # Point-wise variances
  if (s2) {
    if (grad) {
      if (length(theta) == 1) theta_vec <- rep(theta, d) else theta_vec <- theta
      diag_term <- c(rep(1, n_new), rep(2/theta_vec, each = n_new))
      s2_all <- tau2*(diag_term + g - diag_quad_mat(C_cross, C_inv))
      out$s2 <- s2_all[yi]
      out$grad_s2 <- matrix(s2_all[dyi], ncol = d, byrow = FALSE)
    } else {
      out$s2 <- tau2*(1 + g - diag_quad_mat(C_cross, C_inv))
    }
  } 
  
  # Full covariance
  if (sigma) {
    if (v == 999) {
      if (sep) {
        if (grad_enhance | grad) {
          C_new <- Exp2SepGrad(x_new, x_new, grad_indx_new, grad_indx_new, 1, theta, g)
        } else {
          C_new <- Exp2Sep(x_new, x_new, 1, theta, g)
        }
      } else {
        if (grad_enhance | grad) {
          C_new <- Exp2Grad(x_new, x_new, grad_indx_new, grad_indx_new, 1, theta, g)        
        } else {
          C_new <- Exp2(xdmat_new, 1, theta, g)
        }
      }
    } else { # no grad options
      if (sep) {
        C_new <- MaternSep(x_new, x_new, 1, theta, g, v) 
      } else {
        C_new <- Matern(xdmat_new, 1, theta, g, v) 
      }
    }
    quad_term <- C_cross %*% C_inv %*% t(C_cross)
    out$sigma <- tau2*(C_new - quad_term) # possibly includes gradients, not separated
  }
  
  # Posterior samples?  mean and sigma have already been calculated
  if (nsamples > 0) {
    smooth_sigma <- out$sigma - diag(tau2*g - eps, nrow(out$sigma)) # leave jitter on diagonal for numerical stability
    samples <- mvtnorm::rmvnorm(nsamples, mean = mu, sigma = smooth_sigma, checkSymmetry = FALSE)
    if (grad) {
      out$samples <- samples[, yi]
      out$grad_samples <- array(samples[, dyi], dim = c(nsamples, n_new, d))
    } else out$samples <- samples
  }
  
  return(out)
}

# krig_vecchia ----------------------------------------------------------------
# grad_enhance information is contained within approx

krig_vec <- function(y, approx, 
                     tau2 = 1, theta, g = 0, v = 2.5, 
                     sep = FALSE, # use separable lengthscales?  If so, theta must match dimension
                     s2 = FALSE, # calculate and return point-wise variances?
                     sigma = FALSE, # calculate and return full covariance?
                     grad = FALSE, # include predictions of gradients?
                     nsamples = 0, # draw full samples from posterior?  If so, how many?
                     prior_mean = 0,
                     prior_mean_new = 0) {
  
  out <- list()
  grad_enhance <- any(approx$grad_indx > 0)

  # Check settings
  if (nsamples > 0) {
    if (is.null(approx$ord_new)) stop("posterior samples require approx with lite = FALSE")
  } else {
    lite <- is.null(approx$ord_new)
    if (lite & sigma) stop("mismatch: approx is lite, but sigma is FALSE")
    if (!lite & s2) stop("mismatch: approx is not lite, but s2 is TRUE")
    if (s2 & sigma) stop("only s2 or sigma should be true, not both")
    if ((grad_enhance | grad) & v != 999) stop("grad options require 'exp2' kernel")
    if (grad & !lite) stop("no grad option for lite = FALSE")
  }

  # Get ordered and prior-mean-adjusted y
  yo <- y[approx$ord] - ifel(length(prior_mean) == 1, prior_mean, prior_mean[approx$ord])

  if (nsamples > 0) {

    n <- length(y)
    n_new <- nrow(approx$x_ord) - n
    d <- ncol(approx$x_ord)
    samples_all <- matrix(nrow = nsamples, ncol = n_new)
    
    for (i in 1:n_new) {
      
      NN <- rev(approx$NNarray[n + i, ]) # now last index is n + i (the prediction location)
      NN <- NN[!is.na(NN)]
      ncond <- length(NN) - 1
      x_combined <- approx$x_ord[NN, , drop = FALSE] # last entry is index of predictive location
      if (grad_enhance | grad) { # only offered for v = 999
        grad_indx_combined <- approx$grad_indx[NN]
        if (sep) {
          K <- Exp2SepGrad(x_combined, x_combined, grad_indx_combined, grad_indx_combined, 
                           1, theta, g)
        } else {
          K <- Exp2Grad(x_combined, x_combined, grad_indx_combined, grad_indx_combined, 
                             1, theta, g)
        }
      } else {
        if (v == 999) {
          if (sep) {
            K <- Exp2Sep(x_combined, x_combined, 1, theta, g)
          } else K <- Exp2(sq_dist(x_combined), 1, theta, g)
        } else {
          if (sep) {
            K <- MaternSep(x_combined, x_combined, 1, theta, g, v)
          } else K <- Matern(sq_dist(x_combined), 1, theta, g, v)
        }
      }
      L <- t(chol(K))
      s2_store <- tau2*(L[ncond + 1, ncond + 1]^2) # does not depend on sampled y
      for (k in 1:nsamples) {
        if (i == 1) {
          yo_updated <- yo # no sampled points yet
        } else yo_updated <- c(yo, samples_all[k, 1:i]) # TODO: subtract prior mean here?
        mean_store <- L[ncond + 1, 1:ncond] %*% forwardsolve(L[1:ncond, 1:ncond], yo_updated[NN[1:ncond]])
        samples_all[k, i] <- rnorm(1, mean_store, sqrt(s2_store))
      }
    } # end of i for loop
    samples_all <- samples_all[, approx$rev_ord_new, drop = FALSE]
    samples_all <- samples_all + prior_mean_new
    if (grad) { # re-arrange samples from their long stacked vector form
      n_new <- sum(approx$grad_indx[!approx$observed] == 0) # number of unique locations
      yi <- 1:n_new
      dyi <- (n_new + 1):(n_new*(d+1))
      out$samples <- samples_all[, yi]
      out$grad_samples <- array(samples_all[, dyi], dim = c(nsamples, n_new, d))
    } else out$samples <- samples_all
  
  } else { # Only get posterior moments directly if samples are not wanted
  
    if (lite) { # Independent predictions
      
      n_new <- nrow(approx$x_new) 
      d <- ncol(approx$x_new)
      m <- approx$m_new
      mu <- vector(length = n_new)
      if (s2) s2_all <- vector(length = n_new)
      
      for (i in 1:n_new) {
        NN <- approx$NNarray_new[i, ]
        x_combined <- rbind(approx$x_ord[NN, , drop = FALSE], approx$x_new[i, , drop = FALSE])
        if (grad_enhance | grad) { # only offered for v = 999
          grad_indx_combined <- c(approx$grad_indx[NN], approx$grad_indx_new[i])
          if (sep) {
            K <- Exp2SepGrad(x_combined, x_combined, grad_indx_combined, grad_indx_combined, 
                             1, theta, g)
          } else {
            K <- Exp2Grad(x_combined, x_combined, grad_indx_combined, grad_indx_combined, 
                               1, theta, g)
          }
        } else {
          if (v == 999) {
            if (sep) {
              K <- Exp2Sep(x_combined, x_combined, 1, theta, g)
            } else K <- Exp2(sq_dist(x_combined), 1, theta, g)
          } else {
            if (sep) {
              K <- MaternSep(x_combined, x_combined, 1, theta, g, v)
            } else K <- Matern(sq_dist(x_combined), 1, theta, g, v)
          }
        }
        L <- t(chol(K))
        mu[i] <- L[m + 1, 1:m] %*% forwardsolve(L[1:m, 1:m], yo[NN])
        if (s2) s2_all[i] <- tau2*(L[m + 1, m + 1]^2)
      } # end of i for loop
      mu <- mu + prior_mean_new
      if (grad) { # re-arrange mean and s2 from their long stacked vector form
        n_new <- sum(approx$grad_indx_new == 0) # number of unique locations
        yi <- 1:n_new
        dyi <- (n_new + 1):(n_new*(d+1))
        out$mean <- mu[yi]
        out$grad_mean <- matrix(mu[dyi], ncol = d, byrow = FALSE)
        if (s2) {
          out$s2 <- s2_all[yi]
          out$grad_s2 <- matrix(s2_all[dyi], ncol = d, byrow = FALSE)
        }
      } else {
        out$mean <- mu
        if (s2) out$s2 <- s2_all
      }

    } else { # Joint predictions (no gradient option)
        
      U_mat <- create_U(approx, theta = theta, g = g, v = v, sep = sep)
      Upp <- U_mat[!approx$observed, !approx$observed]
      Uppinv <- Matrix::solve(Upp, sparse = TRUE)
      Winv <- Matrix::crossprod(Uppinv)
      Uop <- U_mat[approx$observed, !approx$observed]
      UopTy <- Matrix::crossprod(Uop, yo)
      mu_ordered <- -Matrix::crossprod(Uppinv, UopTy)
      out$mean <- prior_mean_new + mu_ordered[approx$rev_ord_new] 
      if (sigma) out$sigma <- as.matrix(tau2*Winv[approx$rev_ord_new, approx$rev_ord_new])

    } 
  }
    
  return(out)
}