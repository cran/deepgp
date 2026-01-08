
# Function Contents -----------------------------------------------------------
# External:
#   to_vec: converts a "gp", "dgp2", or "dgp3" class object to its Vecchia equivalent
# Internal:
#   rand_mvn_vec: samples from MVN Gaussian using Vecchia approximation
#   create_U: uses C++ to create sparse U matrix
#   create_approx: creates Vecchia approximation object
#   add_pred_to_approx: incorporates predictive locations inside Vecchia approx
#   clean_pred_from_approx: remove predictive information from an approx object

# to_vec ----------------------------------------------------------------------
#' @title Converts non-Vecchia object to its Vecchia version
#' @description Converts an object of class "gp", "dgp2", or "dgp3" to its 
#'              Vecchia equivalent ("gpvec", "dgp2vec", or "dgp3vec").
#' 
#' @details Creates and appends Vecchia-approximation objects to the 
#'     provided fit.  Defaults to random ordering and nearest neighbors
#'     conditioning.  Useful for speeding up predictions when the testing
#'     size is large but the training size is small.
#'
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer} with \code{vecchia = FALSE}
#' @param m size of Vecchia conditioning sets, defaults to the lower of 25 or 
#'        maximum available
#' @param ord optional ordering for Vecchia approximation, defaults to random
#' @param cores number of cores to use for OpenMP parallelization.  Defaults 
#'        to \code{min(4, maxcores - 1)} where \code{maxcores} is the number 
#'        of detectable available cores.
#'
#' @return The same object, with \code{x_approx}, \code{w_approx}, and \code{z_approx}
#'         appended (depending on the number of layers).
#'
#' @examples
#' # See ?fit_one_layer for an example
#' 
#' @export

to_vec <- function(object, m = NULL, ord = NULL, cores = NULL) {
  
  n <- nrow(object$x)
  d <- ncol(object$x)
  grad_enhance <- !is.null(object$dydx)
  if (is.null(m)) m <- min(25, ifel(grad_enhance, n*(d + 1) - 1, n - 1))
  test <- check_vecchia(n, d, m, ord, grad_enhance)
  if (inherits(object, "gp")) {
    object$x_approx <- create_approx(object$x, m, ord, grad_enhance, cores)
    class(object) <- "gpvec"
  } else if (inherits(object, "dgp2")) {
    object$x_approx <- create_approx(object$x, m, ord, grad_enhance, cores)
    object$w_approx <- create_approx(object$w[object$nmcmc, 1:n, ], m, ord, 
                                     grad_enhance, cores)
    class(object) <- "dgp2vec"
  } else if (inherits(object, "dgp3")) { # no grad_enhance option
    object$x_approx <- create_approx(object$x, m, ord, FALSE, cores)
    object$z_approx <- create_approx(object$z[object$nmcmc, , ], m, ord, FALSE, cores)
    object$w_approx <- create_approx(object$w[object$nmcmc, , ], m, ord, FALSE, cores)
    class(object) <- "dgp3vec"
  }
  return(object)
}

# rand_mvn_vec ----------------------------------------------------------------

rand_mvn_vec <- function(approx, tau2 = 1, theta, g = eps, v, grad = FALSE,
                         prior_mean = 0) {
  
  z <- rnorm(nrow(approx$x_ord))
  if (grad) {
    if (v != 999) stop("grad requires 'exp2' kernel") # Matern might be implemented later

    n <- sum(approx$grad_indx == 0)
    d <- ncol(approx$x_ord)
    
    if (length(theta) > 1) {
      U <- U_entries_sep_grad(approx$cores, approx$x_ord, approx$NNarray,
                              approx$grad_indx, tau2, theta, g, v)
    } else {
      U <- U_entries_grad(approx$cores, approx$x_ord, approx$NNarray,
                          approx$grad_indx, tau2, theta, g, v)
    }
    sample <- forward_solve_raw(U, z, approx$NNarray)
    sample <- sample[approx$rev_ord_obs] + prior_mean # prior_mean is NOT ordered
    return(sample)
  } else {
    if (length(theta) > 1) {
      U <- U_entries_sep(approx$cores, approx$x_ord, approx$NNarray, 
                         tau2, theta, g, v)
    } else {
      U <- U_entries(approx$cores, approx$x_ord, approx$NNarray, 
                          tau2, theta, g, v)
    }
    sample <- forward_solve_raw(U, z, approx$NNarray)
    sample <- sample[approx$rev_ord_obs] + prior_mean # prior_mean is NOT ordered
    return(sample)
  }
}

# create_U --------------------------------------------------------------------

create_U <- function(approx, tau2 = 1, theta, g, v, sep = FALSE) {
  
  n <- nrow(approx$x_ord)
  
  if (is.null(approx$grad_indx)) { # no gradients
    if (sep) {
      U <- U_entries_sep(approx$cores, approx$x_ord, approx$NNarray, 
                         tau2, theta, g, v)
    } else {
      U <- U_entries(approx$cores, approx$x_ord, approx$NNarray, 
                     tau2, theta, g, v)
    }
  } else { # gradients
    if (sep) {
      U <- U_entries_sep_grad(approx$cores, approx$x_ord, approx$NNarray, 
                              approx$grad_indx, tau2, theta, g, v)
    } else {
      U <- U_entries_grad(approx$cores, approx$x_ord, approx$NNarray, 
                          approx$grad_indx, tau2, theta, g, v)
    }
  }
  
  U <- U[-approx$na_indx]
  U <- Matrix::sparseMatrix(i = approx$pointers[, 1], j = approx$pointers[, 2], 
                            x = U, dims = c(n, n))
  return(U)
}

# create_approx ---------------------------------------------------------------

create_approx <- function(x, m, ord = NULL, grad_enhance = FALSE, cores = NULL) {
  
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  if (is.null(cores)) cores <- check_cores(cores)

  if (is.null(ord)) ord <- sample(1:n, n, replace = FALSE) # random

  if (grad_enhance) { 
    grad_indx <- rep(0:d, each = n)
    ord <- rep(ord, times = d + 1) + n*grad_indx
    x <- bind(x, d)
  } else grad_indx <- rep(0, times = n)

  x_ord <- x[ord, , drop = FALSE]
  NNarray <- GpGp::find_ordered_nn(x_ord, m) # TODO: avoid duplicate NN calculations for gradient locations?
  rev_ord_obs <- order(ord) # includes partial derivative indices
  na_indx <- which(is.na(t(NNarray))) # vector indices of NA values (for converting U to sparse matrix)
  pointers <- row_col_pointers(NNarray)
  
  out <- list(m = m, ord = ord, NNarray = NNarray, na_indx = na_indx,
              pointers = pointers, rev_ord_obs = rev_ord_obs, 
              x_ord = x_ord, cores = cores, grad_indx = grad_indx)
  return(out)
}

# add_pred_to_approx ----------------------------------------------------------

add_pred_to_approx <- function(approx, x_new, m, lite, grad = FALSE, ord_new = NULL) {
  
  approx$m_new <- m
  if (!is.matrix(x_new)) x_new <- as.matrix(x_new)
  n_new <- nrow(x_new)
  d <- ncol(x_new)
  
  if (lite) {
    if (grad) { # do we want predictions of the gradient?
      x_new <- bind(x_new, d)
      grad_indx_new <- rep(0:d, each = n_new)
      approx$grad_indx_new <- grad_indx_new
    } else approx$grad_indx_new <- rep(0, times = n_new)
    approx$x_new <- x_new
    approx$NNarray_new <- FNN::get.knnx(approx$x_ord, x_new, m)$nn.index
  } else { 
    if (is.null(ord_new)) { # TODO: don't update for every t in T?
      ord_new <- sample(1:n_new, n_new, replace = FALSE)
    }
    if (grad) {
      x_new <- bind(x_new, d)
      grad_indx_new <- rep(0:d, each = n_new)
      ord_new <- rep(ord_new, times = d + 1) + n_new*grad_indx_new
      approx$grad_indx <- c(approx$grad_indx, grad_indx_new)
    } else approx$grad_indx <- c(approx$grad_indx, rep(0, times = n_new))
    approx$ord_new <- ord_new
    approx$rev_ord_new <- order(ord_new)
    approx$observed <- c(rep(TRUE, times = length(approx$ord)), 
                         rep(FALSE, times = length(approx$ord_new)))
   
    # OVERWRITE THE FOLLOWING
    approx$x_ord <- rbind(approx$x_ord, x_new[ord_new, , drop = FALSE])
    approx$NNarray <- GpGp::find_ordered_nn(approx$x_ord, m) # TODO: don't recalculate NN?
    approx$na_indx <- which(is.na(t(approx$NNarray)))
    approx$pointers <- row_col_pointers(approx$NNarray)

  }
  
  return(approx)
}

# clean_pred_from_approx ------------------------------------------------------

clean_pred_from_approx <- function(approx) {

  if (!is.null(approx$m_new)) approx["m_new"] <- NULL
  lite <- is.null(approx$observed)
  if (lite) {
    if (!is.null(approx$x_new)) approx["x_new"] <- NULL
    if (!is.null(approx$NNarray_new)) approx["NNarray_new"] <- NULL
    if (!is.null(approx$grad_indx_new)) approx["grad_indx_new"] <- NULL
  } else {
    n <- sum(approx$observed)
    if (!is.null(approx$ord_new)) approx["ord_new"] <- NULL
    if (!is.null(approx$rev_ord_new)) approx["rev_ord_new"] <- NULL
    approx$grad_indx <- approx$grad_indx[approx$observed]
    approx$x_ord <- approx$x_ord[approx$observed, , drop = FALSE]
    approx$NNarray <- approx$NNarray[approx$observed, 1:(approx$m+1), drop = FALSE]
    approx$na_indx <- which(is.na(t(approx$NNarray)))
    approx$pointers <- row_col_pointers(approx$NNarray)
  }
  return(approx)
}
