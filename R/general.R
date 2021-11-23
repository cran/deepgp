
# Function Contents -----------------------------------------------------------
# Internal:
#   krig: calculates posterior mean, sigma, and tau2 (optionally f_min)
#   invdet: calculates inverse and log determinant of a matrix using C
#   fill_final_row: uses kriging to fill final row of w_0/z_0
# External (see documentation below):
#   sq_dist
#   score
#   rmse

eps <- sqrt(.Machine$double.eps)

# Kriging Function ------------------------------------------------------------
# Calculates posterior mean, sigma/s2 and estimate of tau2 using kriging equations
# Optionally calculates f_min (min predicted mean at observed data locations)

krig <- function(y, dx, d_new = NULL, d_cross = NULL, theta, g, mean = TRUE, 
                 s2 = FALSE, sigma = FALSE, tau2 = FALSE, f_min = FALSE, 
                 cov = "matern", v = 2.5) {
  
  if (s2 & sigma) s2 <- FALSE # don't calculate diagonal separately
  
  out <- list()
  if (cov == "matern") {
    C <- MaternFun(dx, c(1, theta, g, v)) 
  } else C <- ExpFun(dx, c(1, theta, g))
  C_inv <- invdet(C)$Mi
  
  if (mean) {
    if (cov == "matern") {
      C_cross <- MaternFun(d_cross, c(1, theta, 0, v))
    } else C_cross <- ExpFun(d_cross, c(1, theta, 0)) # no g in rectangular matrix
    out$mean <- C_cross %*% C_inv %*% y
  }
  
  if (f_min) { # predict at observed locations, return min expected value
    if (cov == "matern") {
      C_cross_observed_only <- MaternFun(dx, c(1, theta, 0, v))
    } else C_cross_observed_only <- ExpFun(dx, c(1, theta, 0))
    out$f_min <- min(C_cross_observed_only %*% C_inv %*% y)
  }
  
  if (tau2) {
    scale <- c(t(y) %*% C_inv %*% y)/length(y)
    out$tau2 <- scale
  } else scale <- 1
  
  if (s2) {
    if (!mean) {
      if (cov == "matern") {
        C_cross <- MaternFun(d_cross, c(1, theta, 0, v))
      } else C_cross <- ExpFun(d_cross, c(1, theta, 0))
    }
    C_new <- rep(1 + g, times = nrow(d_new))
    out$s2 <- scale * (C_new - diag(C_cross %*% C_inv %*% t(C_cross)))
  }
  
  if (sigma) {
    if (!mean) {
      if (cov == "matern") {
        C_cross <- MaternFun(d_cross, c(1, theta, 0, v))
      } else C_cross <- ExpFun(d_cross, c(1, theta, 0))
    }
    if (cov == "matern") {
      C_new <- MaternFun(d_new, c(1, theta, g, v)) 
    } else C_new <- ExpFun(d_new, c(1, theta, g))
    out$sigma <- scale * (C_new - (C_cross %*% C_inv %*% t(C_cross)))
  }
  
  return(out)
}

# Matrix Inverse C Function ---------------------------------------------------
# Calculates matrix inverse and log determinant using C
# Credit given to the "laGP" package (Robert B Gramacy & Furong Sun)

invdet <- function(M) {
  
  n <- nrow(M)
  out <- .C("inv_det_R",
            n = as.integer(n),
            M = as.double(M),
            Mi = as.double(diag(n)),
            ldet = double(1))
  
  return(list(Mi = matrix(out$Mi, ncol=n), ldet = out$ldet))
}

# Distance C Function ---------------------------------------------------------
#' @title Calculates squared pairwise distances
#' @description Calculates squared pairwise euclidean distances using C.
#' 
#' @details C code derived from the "laGP" package (Robert B Gramacy and 
#'     Furong Sun).
#' @param X1 matrix of input locations
#' @param X2 matrix of second input locations (if \code{NULL}, distance is 
#'        calculated between \code{X1} and itself)
#' @return symmetric matrix of squared euclidean distances
#' 
#' @references 
#' Gramacy, RB and F Sun. (2016). laGP: Large-Scale Spatial Modeling via Local 
#'     Approximate Gaussian Processes in R. \emph{Journal of Statistical Software 
#'     72} (1), 1-46. doi:10.18637/jss.v072.i01
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
               D = double(n1 * n1))
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
               D = double(n1 * n2))
    return(matrix(outD$D, ncol = n2, byrow = TRUE))
  }
}

# Calculate Score Function ----------------------------------------------------
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
#' @export

score <- function(y, mu, sigma) {
  
  id <- invdet(sigma)
  score <- (- id$ldet - t(y - mu) %*% id$Mi %*% (y - mu)) / length(y)
  return(c(score))
}

# Calculate RMSE Function -----------------------------------------------------
#' @title Calculates RMSE
#' @description Calculates root mean square error (lower RMSE indicate better fits).
#' @param y response vector
#' @param mu predicted mean
#' 
#' @examples
#' # See "deepgp-package", "fit_one_layer", "fit_two_layer", or "fit_three_layer"
#' # for an example
#' 
#' @export

rmse <- function(y, mu) {
  return(sqrt(mean((y - mu) ^ 2)))
}

# Fill final row function -----------------------------------------------------
# Uses kriging prediction to fill the final row of w_0 or z_0
# Used in sequential design

fill_final_row <- function(x, w_0, D, theta_w_0, cov, v) {
  n <- nrow(x)
  dx <- sq_dist(x)
  new_w <- vector(length = D)
  old_x <- x[1:(n - 1), ]
  new_x <- matrix(x[n, ], nrow = 1)
  for (i in 1:D) {
    new_w[i] <- krig(w_0[, i], dx[1:(n - 1), 1:(n - 1)], 
                     d_cross = sq_dist(new_x, old_x), theta = theta_w_0[i], 
                     g = eps, cov = cov, v = v)$mean
  }
  return(rbind(w_0, new_w))
}
