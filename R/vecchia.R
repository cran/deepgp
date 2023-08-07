
# Function Contents -----------------------------------------------------------
# Internal:
#   rand_mvn_vec: samples from MVN Gaussian using vecchia approximation
#   fill_final_row_vec: uses kriging to fill final row of w_0/z_0
#   create_U: uses C++ to create sparse Matrix U
#   create_approx: creates vecchia approximation object
#   update_obs_in_approx: updates x_ord inside approx object
#   add_pred_to_approx: incorporates predictive locations inside vecchia approx
#   krig_vec: prediction using vecchia

# Random MVN ------------------------------------------------------------------

rand_mvn_vec <- function(approx, theta, v, mean = rep(0, nrow(approx$x_ord)),
                         g = eps, scale = 1) {
  
  z <- rnorm(nrow(approx$x_ord))
  Uraw <- create_U(approx, g, theta, v, raw_form = TRUE)
  sample <- forward_solve_raw(Uraw, z, approx$NNarray)
  sample <- sample[approx$rev_ord_obs] + mean
  return(sample)
}

# Fill Final Row --------------------------------------------------------------

fill_final_row_vec <- function(x, w_0, D, theta_w_0, v, m) { 
  
  n <- nrow(x)
  new_w <- vector(length = D)
  old_x <- as.matrix(x[1:(n - 1), ])
  new_x <- matrix(x[n, ], nrow = 1)
  for (i in 1:D) 
    new_w[i] <- krig_vec(w_0[, i], theta_w_0[i], eps, v = v, m = m,
                         x = old_x, x_new = new_x)$mean
  
  return(rbind(w_0, new_w))
}

# Create U --------------------------------------------------------------------

create_U <- function(approx, g, theta, v, sep = FALSE,
                     raw_form = FALSE) {
  
  n <- nrow(approx$x_ord)
  if (sep) {
    U <- U_entries_sep(approx$n_cores, approx$x_ord, approx$revNN, 
                       1, theta, g, v)
  } else U <- U_entries(approx$n_cores, approx$x_ord, approx$revNN, 
                        1, theta, g, v)
  
  if (raw_form) {
    m <- ncol(approx$NNarray) - 1
    U <- rev_matrix(U)
    for (i in 1:m) {
      zeros <- (U[i, ] == 0)
      U[i, ] <- c(U[i, !zeros], rep(0, times = sum(zeros))) 
    }
    return(U)
  } else {
    U <- c(t(U))[approx$notNA]
    U <- Matrix::sparseMatrix(i = approx$pointers[, 2], 
                              j = approx$pointers[, 1], x = U, 
                              dims = c(n, n))
    return(U)
  } 
}

# Create Approximation --------------------------------------------------------

create_approx <- function(x, m) {
  
  n <- nrow(x)
  rev_ord_obs <- sample(1:n, n, replace = FALSE) 
  ord <- order(rev_ord_obs)
  x_ord <- x[ord, , drop = FALSE]
 
  NNarray <- GpGp::find_ordered_nn(x_ord, m)
  revNN <- rev_matrix(NNarray)
  notNA <- as.vector(t(!is.na(NNarray)))
  pointers <- row_col_pointers(NNarray)
  n_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
  
  out <- list(m = m, ord = ord, NNarray = NNarray, revNN = revNN, 
              notNA = notNA, pointers = pointers,
              n_cores = n_cores, rev_ord_obs = rev_ord_obs, x_ord = x_ord)
  return(out)
}

# Update Observation in Approx ------------------------------------------------

update_obs_in_approx <- function(approx, x_new, col_index = NULL) {
  
  if (is.null(col_index)) {
    approx$x_ord <- x_new[approx$ord, , drop = FALSE]
  } else {
    approx$x_ord[, col_index] <- x_new[approx$ord]
  }
  
  return(approx)
}

# Add Pred to Approx ----------------------------------------------------------

add_pred_to_approx <- function(approx, x_pred, m) {
  
  n <- nrow(approx$x_ord)
  n_pred <- nrow(x_pred)
  ord_pred <- sample(1:n_pred, n_pred, replace = FALSE)
  rev_ord_pred <- order(ord_pred)
  
  ord <- c(approx$ord, ord_pred + n) # observed data FIRST
  x_ord <- rbind(approx$x_ord, x_pred[ord_pred, , drop = FALSE])
  observed <- c(rep(TRUE, n), rep(FALSE, n_pred))
  
  NNarray <- GpGp::find_ordered_nn(x_ord, m)
  revNN <- rev_matrix(NNarray)
  notNA <- as.vector(t(!is.na(NNarray)))
  pointers <- row_col_pointers(NNarray)

  out <- list(m = m, ord = ord, NNarray = NNarray, revNN = revNN, 
              notNA = notNA, pointers = pointers,
              n_cores = approx$n_cores, rev_ord_obs = approx$rev_ord_obs,
              rev_ord_pred = rev_ord_pred, observed = observed, x_ord = x_ord)
  return(out)
}

# Krig Vecchia ----------------------------------------------------------------

krig_vec <- function(y, theta, g, tau2 = 1, 
                     s2 = FALSE, sigma = FALSE, v, m = NULL,
                     x = NULL, x_new =  NULL, # inputs required for sigma = FALSE 
                     NNarray_pred = NULL, # optional input for sigma = FALSE
                     approx = NULL, # inputs required for sigma = TRUE
                     sep = FALSE,
                     prior_mean = rep(0, length(y)), 
                     prior_mean_new = 0) {
  
  out <- list()
  
  if (!sigma) { # lite = TRUE 
    
    n_new <- nrow(x_new)
    if (is.null(NNarray_pred)) # calculate NN if not provided
      NNarray_pred <- FNN::get.knnx(x, x_new, m)$nn.index # NN DO NOT USE SCALED LENGTHSCALES IN SEPARABLE GP
    out$mean <- vector(length = n_new)
    if (s2) out$s2 <- vector(length = n_new)
    for (i in 1:n_new) {
      NN <- NNarray_pred[i, ]
      x_combined <- rbind(x[NN, , drop = FALSE], x_new[i, , drop = FALSE])
      if (v == 999) {
        if (sep) {
          K <- Exp2Sep(x_combined, x_combined, 1, theta, g)
        } else K <- Exp2(sq_dist(x_combined), 1, theta, g)
      } else {
        if (sep) {
          K <- MaternSep(x_combined, x_combined, 1, theta, g, v)
        } else K <- Matern(sq_dist(x_combined), 1, theta, g, v)
      }
      L <- t(chol(K))
      out$mean[i] <- L[m + 1, 1:m] %*% forwardsolve(L[1:m, 1:m], y[NN] - prior_mean[NN])
      if (s2) out$s2[i] <- tau2 * (L[m + 1, m + 1] ^ 2)
    }
    out$mean <- out$mean + prior_mean_new
    
  } else { # lite = FALSE
    
    if (is.null(approx$observed)) # add pred to approximation if not provided
      approx <- add_pred_to_approx(approx, x_new, m)
    prior_mean_ordered <- prior_mean[approx$ord[approx$observed]]
    yo <- y[approx$ord[approx$observed]] - prior_mean_ordered
    U_mat <- create_U(approx, g, theta, v, sep = sep)
    Upp <- U_mat[!approx$observed, !approx$observed]
    Uppinv <- Matrix::solve(Upp, sparse = TRUE)
    Winv <- Matrix::crossprod(Uppinv)
    Uop <- U_mat[approx$observed, !approx$observed]
    UopTy <- Matrix::crossprod(Uop, yo)
    mu_ordered <- -Matrix::crossprod(Uppinv, UopTy)
    out$mean <- prior_mean_new + mu_ordered[approx$rev_ord_pred]
    out$sigma <- as.matrix(tau2 * Winv[approx$rev_ord_pred, approx$rev_ord_pred])
    
  }
  
  return(out)
}
