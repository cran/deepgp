
# Function Contents -----------------------------------------------------------
# Internal:
#   eps: minimum nugget value (1.5e-8)
#   krig: calculates posterior mean, sigma, and tau2 (optionally f_min)
#   krig_sep: same as above, but with a separable kernel
#   invdet: calculates inverse and log determinant of a matrix using C
#           (credit given to the "laGP package, R.B. Gramacy & Furong Sun)
#   fill_final_row: uses kriging to fill final row of w_0/z_0
#   exp_improv: calculates expected improvement
#   calc_entropy: calculates entropy for two classes (pass/fail)
#   ifel: returns second or third element depending on first argument
#   monowarp_ref: uses cumulative sum on a reference grid to make warping monotonic
#   fo_approx: linearly approximates along a reference grid
# External (see documentation below):
#   sq_dist
#   score
#   rmse
#   crps

eps <- sqrt(.Machine$double.eps)

# Krig ------------------------------------------------------------------------

krig <- function(y, dx, d_new = NULL, d_cross = NULL, theta, g, tau2 = 1,
                 s2 = FALSE, sigma = FALSE, f_min = FALSE, v = 2.5,
                 prior_mean = 0, prior_mean_new = 0) {
  
  out <- list()
  if (v == 999) {
    C <- Exp2(dx, 1, theta, g)
    C_cross <- Exp2(d_cross, 1, theta, 0) # no g in rectangular matrix
  } else {
    C <- Matern(dx, 1, theta, g, v) 
    C_cross <- Matern(d_cross, 1, theta, 0, v)
  }
  C_inv <- invdet(C)$Mi
  out$mean <- prior_mean_new + C_cross %*% C_inv %*% (y - prior_mean)
  
  if (f_min) { # predict at observed locations, return min expected value
    if (v == 999) {
      C_cross_observed_only <- Exp2(dx, 1, theta, 0)
    } else C_cross_observed_only <- Matern(dx, 1, theta, 0, v)
    out$f_min <- min(C_cross_observed_only %*% C_inv %*% y)
  }
  
  if (s2) {
    C_new <- rep(1 + g, times = nrow(d_cross))
    out$s2 <- tau2 * (C_new - diag_quad_mat(C_cross, C_inv))
  }
  
  if (sigma) {
    quad_term <- C_cross %*% C_inv %*% t(C_cross)
    if (v == 999) {
      C_new <- Exp2(d_new, 1, theta, g)
    } else  C_new <- Matern(d_new, 1, theta, g, v) 
    out$sigma <- tau2 * (C_new - quad_term)
  }
  return(out)
}

# Krig SEPARABLE --------------------------------------------------------------

krig_sep <- function(y, x, x_new, theta, g, tau2 = 1,
                     s2 = FALSE, sigma = FALSE, f_min = FALSE, v = 2.5) {
  
  out <- list()
  if (v == 999) {
    C <- Exp2Sep(x, x, 1, theta, g)
    C_cross <- Exp2Sep(x_new, x, 1, theta, 0) # no g in rectangular matrix
  } else {
    C <- MaternSep(x, x, 1, theta, g, v) 
    C_cross <- MaternSep(x_new, x, 1, theta, 0, v)
  }
  C_inv <- invdet(C)$Mi
  out$mean <- C_cross %*% C_inv %*% y
  
  if (f_min) { # predict at observed locations, return min expected value
    if (v == 999) {
      C_cross_observed_only <- Exp2Sep(x, x, 1, theta, 0)
    } else C_cross_observed_only <- MaternSep(x, x, 1, theta, 0, v)
    out$f_min <- min(C_cross_observed_only %*% C_inv %*% y)
  }
  
  if (s2) {
    C_new <- rep(1 + g, times = nrow(x_new))
    out$s2 <- tau2 * (C_new - diag_quad_mat(C_cross, C_inv))
  }
  
  if (sigma) {
    quad_term <- C_cross %*% C_inv %*% t(C_cross)
    if (v == 999) {
      C_new <- Exp2Sep(x_new, x_new, 1, theta, g)
    } else  C_new <- MaternSep(x_new, x_new, 1, theta, g, v) 
    out$sigma <- tau2 * (C_new - quad_term)
  }
  return(out)
}

# Matrix Inverse --------------------------------------------------------------

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

# Squared Distance-------------------------------------------------------------
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

# Score -----------------------------------------------------------------------
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

# RMSE ------------------------------------------------------------------------
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
  return(sqrt(mean((y - mu) ^ 2)))
}

# CRPS ------------------------------------------------------------------------
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
  z <- (y - mu) / sigma
  return(mean(sigma * (-(1 / sqrt(pi)) + 2 * dnorm(z) + z * (2 * pnorm(z) - 1))))
}

# Fill Final Row --------------------------------------------------------------

fill_final_row <- function(x, w_0, D, theta_w_0, v) {
  n <- nrow(x)
  dx <- sq_dist(x)
  new_w <- vector(length = D)
  old_x <- x[1:(n - 1), ]
  new_x <- matrix(x[n, ], nrow = 1)
  for (i in 1:D) {
    new_w[i] <- krig(w_0[, i], dx[1:(n - 1), 1:(n - 1)], 
                     d_cross = sq_dist(new_x, old_x), theta = theta_w_0[i], 
                     g = eps, v = v)$mean
  }
  return(rbind(w_0, new_w))
}

# EI --------------------------------------------------------------------------

exp_improv <- function(mu, sig2, f_min) {
  
  s <- sqrt(sig2)
  z <- (f_min - mu) / s
  ei <- (f_min - mu)*pnorm(z) + s*dnorm(z)
  
  return(ei)
}

# Entropy ---------------------------------------------------------------------

calc_entropy <- function(mu, sig2, limit) {
  
  fail_prob <- pnorm((mu - limit) / sqrt(sig2))
  ent <- -(1 - fail_prob) * log(1 - fail_prob) - fail_prob * log(fail_prob)
  ent[which(is.nan(ent))] <- 0
  
  return(ent)
}

# If else ---------------------------------------------------------------------

ifel <- function(logical, yes, no) {
  if (logical) {
    return(yes)
  } else return(no)
}

# Monowarp --------------------------------------------------------------------

monowarp_ref <- function(x, xg, wg, index) {
  # x: matrix of input locations for returned w
  # xg: matrix of grid locations
  # wg: matrix of w values at xg locations to be warped
  if (!is.matrix(x)) x <- matrix(x, ncol = 1)
  if (!is.matrix(xg)) xg <- matrix(xg, ncol = 1)
  if (!is.matrix(wg)) wg <- matrix(wg, ncol = 1)
  if (!is.matrix(index)) index <- matrix(index, ncol = 1)
  
  w <- matrix(nrow = nrow(x), ncol = ncol(wg))
  for (i in 1:ncol(wg)) {
    wg[, i] <- exp(wg[, i])
    wg[, i] <- cumsum(wg[, i])
    r <- range(wg[, i])
    wg[, i] <- (wg[, i] - r[1]) / (r[2] - r[1])
    w[, i] <- fo_approx(xg[, i], wg[, i], x[, i], index[, i])
  }
  return(w)
}

# Fixed order approx ----------------------------------------------------------
# Implements linear approximation assuming that input locations are unchanged,
# so their ordering may be pre-calculated.  Uses linear extrapolation outside
# the range of the data.

fo_approx <- function(xg, wg, x, index = NULL) {
  # calculate slopes and intercepts between every adjacent pair
  n <- length(wg)
  slopes <- icepts <- rep(NA, n + 1)
  slopes[-c(1, n + 1)] <- (wg[2:n] - wg[1:(n - 1)]) / (xg[2:n] - xg[1:(n - 1)])
  icepts[-c(1, n + 1)] <- -slopes[-c(1, n + 1)] * xg[1:(n-1)] + wg[1:(n-1)] 
    
  # accommodating linear extrapolation outside the range
  slopes[1] <- slopes[2]
  slopes[n + 1] <- slopes[n]
  icepts[1] <- icepts[2]
  icepts[n + 1] <- icepts[n]

  if (is.null(index)) index <- fo_approx_init(as.matrix(xg), as.matrix(x))
      
  # linear interpolation
  return(x * slopes[index] + icepts[index])
}
