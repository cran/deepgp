
# Function Contents -----------------------------------------------------------
# Internal:
#   eps: minimum nugget value (1.5e-8)
#   invdet: calculates inverse and log determinant of a matrix using C
#           (credit given to the "laGP package, R.B. Gramacy & Furong Sun)
#   fill_final_rows: uses kriging to fill-in w_0/z_0
#   exp_improv: calculates expected improvement
#   calc_entropy: calculates entropy for two classes (pass/fail)
#   ifel: returns second or third element depending on first argument
#   monotransform: uses cumulative sum on a reference grid to make warping monotonic
#   bind: appends copies of a matrix to itself (for gradient-enhancement)
#   get_dydw: solves for dydw given dwdx and dydx (for gradient-enhancement)
#   get_prior_mean: provides identity warping with appropriate gradients
# External (see documentation below):
#   sq_dist
#   score
#   rmse
#   crps

# eps -------------------------------------------------------------------------

eps <- sqrt(.Machine$double.eps)

# invdet ----------------------------------------------------------------------

invdet <- function(M) {

  n <- nrow(M)
  out <- .C("inv_det_R",
            n = as.integer(n),
            M = as.double(M),
            Mi = as.double(diag(n)),
            ldet = double(1),
            PACKAGE = "deepgp")

  return(list(Mi = matrix(out$Mi, ncol=n), ldet = out$ldet))
}

# sq_dist ---------------------------------------------------------------------
#' @title Calculates squared pairwise distances
#' @description Calculates squared pairwise euclidean distances using C.
#' 
#' @details C code derived from the "laGP" package (Robert B Gramacy and 
#'     Furong Sun).
#'     
#' @param X1 matrix of input locations
#' @param X2 matrix of second input locations (if \code{NULL}, distance is 
#'        calculated between \code{X1} and itself)
#'        
#' @return symmetric matrix of squared euclidean distances
#' 
#' @references 
#' Gramacy, RB and F Sun. (2016). laGP: Large-Scale Spatial Modeling via Local 
#'     Approximate Gaussian Processes in R. \emph{Journal of Statistical 
#'     Software 72} (1), 1-46. doi:10.18637/jss.v072.i01
#' 
#' @examples 
#' x <- seq(0, 1, length = 10)
#' d2 <- sq_dist(x)
#' 
#' @export

sq_dist <- function(X1, X2 = NULL) {

  X1 <- as.matrix(X1)
  n1 <- nrow(X1)
  m <- ncol(X1)

  if(is.null(X2)) {
    outD <- .C("distance_symm_R",
               X = as.double(t(X1)),
               n = as.integer(n1),
               m = as.integer(m),
               D = double(n1 * n1),
               PACKAGE = "deepgp")
    return(matrix(outD$D, ncol = n1, byrow = TRUE))
  } else {
    X2 <- as.matrix(X2)
    n2 <- nrow(X2)
    if(ncol(X1) != ncol(X2)) stop("dimension of X1 & X2 do not match")
    outD <- .C("distance_R",
               X1 = as.double(t(X1)),
               n1 = as.integer(n1),
               X2 = as.double(t(X2)),
               n2 = as.integer(n2),
               m = as.integer(m),
               D = double(n1 * n2),
               PACKAGE = "deepgp")
    return(matrix(outD$D, ncol = n2, byrow = TRUE))
  }
}

# score -----------------------------------------------------------------------
#' @title Calculates score
#' @description Calculates score, proportional to the multivariate normal log
#'     likelihood.  Higher scores indicate better fits.  Only 
#'     applicable to noisy data.  Requires full covariance matrix 
#'     (e.g. \code{predict} with \code{lite = FALSE}).
#'     
#' @param y response vector
#' @param mu predicted mean
#' @param sigma predicted covariance
#' 
#' @references 
#' Gneiting, T, and AE Raftery. 2007. Strictly Proper Scoring Rules, Prediction, 
#'     and Estimation. \emph{Journal of the American Statistical Association 102} 
#'     (477), 359-378.
#'     
#' @examples
#' # Additional examples including real-world computer experiments are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' 
#' @export

score <- function(y, mu, sigma) {
  
  id <- invdet(sigma)
  score <- (- id$ldet - t(y - mu) %*% id$Mi %*% (y - mu)) / length(y)
  return(c(score))
}

# rmse ------------------------------------------------------------------------
#' @title Calculates RMSE
#' @description Calculates root mean square error (lower RMSE indicate better 
#'     fits).
#' 
#' @param y response vector
#' @param mu predicted mean
#' 
#' @examples
#' # Additional examples including real-world computer experiments are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' 
#' @export

rmse <- function(y, mu) {
  return(sqrt(mean((y - mu)^2)))
}

# crps ------------------------------------------------------------------------
#' @title Calculates CRPS
#' @description Calculates continuous ranked probability score (lower CRPS indicate
#' better fits, better uncertainty quantification).
#' 
#' @param y response vector
#' @param mu predicted mean
#' @param s2 predicted point-wise variances
#' 
#' @references 
#' Gneiting, T, and AE Raftery. 2007. Strictly Proper Scoring Rules, Prediction, 
#'     and Estimation. \emph{Journal of the American Statistical Association 102} 
#'     (477), 359-378.
#'     
#' @examples
#' # Additional examples including real-world computer experiments are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' 
#' @export

crps <- function(y, mu, s2) {
  sigma <- sqrt(s2)
  z <- (y - mu)/sigma
  return(mean(sigma*(-(1/sqrt(pi)) + 2*dnorm(z) + z*(2*pnorm(z) - 1))))
}

# fill_final_rows -------------------------------------------------------------

fill_final_rows <- function(w, x, tau2, theta, v, pmx = FALSE, vecchia = FALSE,
                            m = nrow(x) - 1) {

  n <- nrow(w)
  n_new <- nrow(x)
  if (n_new <= n) stop("no additional rows to fill")
  D <- ncol(w)
  if (pmx & (ncol(x) != D)) stop("pmx requires ncol(x) = D")
  old_x <- x[1:n, , drop = FALSE]
  new_x <- x[(n + 1):n_new, , drop = FALSE]
  new_w <- matrix(nrow = n_new - n, ncol = D)
  
  if (vecchia) {
    temp_approx <- create_approx(old_x, m)
    temp_approx <- add_pred_to_approx(temp_approx, new_x, m, lite = TRUE)
    for (i in 1:D) {
      new_w[i] <- krig_vec(w[, i], 
                           approx = temp_approx,
                           tau2 = tau2,
                           theta = ifel(length(theta) == 1, theta, theta[i]), 
                           g = eps,
                           v = v,
                           prior_mean = ifel(pmx, old_x[, i], 0),
                           prior_mean_new = ifel(pmx, new_x[, i], 0))$mean
    }
  } else {
    xdmat <- sq_dist(old_x)
    xdmat_cross <- sq_dist(new_x, old_x)
    for (i in 1:D) {
      new_w[, i] <- krig(w[, i], 
                         xdmat = xdmat, 
                         xdmat_cross = xdmat_cross, 
                         tau2 = tau2, 
                         theta = ifel(length(theta) == 1, theta, theta[i]), 
                         g = eps, 
                         v = v,
                         prior_mean = ifel(pmx, old_x[, i], 0),
                         prior_mean_new = ifel(pmx, new_x[, i], 0))$mean 
    }
  }
  return(rbind(w, new_w))
}

# exp_improv ------------------------------------------------------------------

exp_improv <- function(mu, sig2, f_min) {
  
  s <- sqrt(sig2)
  z <- (f_min - mu)/s
  ei <- (f_min - mu)*pnorm(z) + s*dnorm(z)
  
  return(ei)
}

# calc_entropy ----------------------------------------------------------------

calc_entropy <- function(mu, sig2, limit) {
  
  fail_prob <- pnorm((mu - limit) / sqrt(sig2))
  ent <- -(1 - fail_prob)*log(1 - fail_prob) - fail_prob*log(fail_prob)
  ent[which(is.nan(ent))] <- 0
  
  return(ent)
}

# ifel ------------------------------------------------------------------------

ifel <- function(logical, yes, no) {
  if (logical) {
    return(yes)
  } else return(no)
}


# monotransform ---------------------------------------------------------------

monotransform <- function(x, x_grid, w_grid) {
  if (!all(diff(x_grid) >= 0)) stop("x_grid must be increasing")
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(w_grid)) w_grid <- as.matrix(w_grid)
  d <- ncol(x)
  wwarp <- matrix(nrow = nrow(x), ncol = d)
  for (i in 1:d) {
    wwarp_grid <- c(0, cumsum(abs(diff(w_grid[, i]))))
    wwarp[, i] <- suppressWarnings(approx(x_grid, wwarp_grid, x[, i])$y)
  }
  return(wwarp)
}

# bind ------------------------------------------------------------------------

bind <- function(x, d) {
  # Appends x to itself d-many times
  if (!is.matrix(x)) x <- as.matrix(x)
  xbind <- x
  for (i in 1:d)
    xbind <- rbind(xbind, x)
  return(xbind)
}

# get_dydw --------------------------------------------------------------------

get_dydw <- function(w, dydx) {
  # Solves for dydw given dwdx (latter rows of w) and dydx
  if (!is.matrix(w)) w <- as.matrix(w)
  n <- nrow(dydx)
  d <- ncol(dydx)
  dydw <- matrix(nrow = n, ncol = d)
  for (i in 1:n) {
    irows <- seq(n+i, n*(d+1), by = n)
    dydw[i, ] <- solve(w[irows, , drop = FALSE], dydx[i, ])
  }
  return(dydw)
}

# get_prior_mean --------------------------------------------------------------

get_prior_mean <- function(x) {
  # Calculates prior mean for latent layer W when grad_enhance and pmx are TRUE
  d <- ncol(x)
  n <- nrow(x)
  prior_mean <- x
  for (i in 1:d) {
    zeromat <- matrix(0, nrow = n, ncol = d)
    zeromat[, i] <- rep(1, n)
    prior_mean <- rbind(prior_mean, zeromat)
  }
  return(prior_mean)
}
