
# Function Contents -----------------------------------------------------------
# Internal:
#   rand_mvn: generates random draw from MVN distribution
#   krig: calculates posterior mean, sigma, and tau2
#   invdet: calculates inverse and log determinant of a matrix using C
#   calc_tau2: calculates MLE for tau2
# External (see documentation below):
#   calc_K
#   sq_dist
#   score
#   rmse

# Covariance Function ---------------------------------------------------------
#' @title Calculates covariance matrix
#' @description Calculates covariance matrix based on inverse exponentiated
#'     squared euclidean distance with specified hyperparameters.
#' @param d2 matrix of squared distances among input locations
#' @param theta length scale parameter determining strength of correlation
#' @param g nugget parameter determining noise level (only used if \code{d2} 
#'        is square)
#' @return symmetric covariance matrix
#' 
#' @examples 
#' x <- seq(0, 1, length = 10)
#' K <- calc_K(sq_dist(x), theta = 0.1, g = 0.01)
#' 
#' @export

calc_K <- function(d2, theta, g = NULL) {

  if (is.null(g)) g <- sqrt(.Machine$double.eps)
  if (nrow(d2) == ncol(d2)) return(exp(- d2 / theta) + diag(g, ncol(d2)))
  return(exp(- d2 / theta))
}

# Random MVN Sample Function --------------------------------------------------
# Draws a random sample from multivariate normal distribution
# Credit given to the "mvtnorm" package (Alan Genz et al)

rand_mvn <- function(n, mean = rep(0, nrow(sigma)), sigma) {
  
  suppressWarnings({ L <- chol(sigma, pivot = TRUE) })
  L <- L[, order(attr(L, "pivot"))]
  sample <- matrix(rnorm(n * nrow(sigma), 0, 1), nrow = n) %*% L
  sample <- sweep(sample, 2, mean, "+")
  return(t(sample)) # output is an nrow(sigma) by n matrix
}

# Kriging Function ------------------------------------------------------------
# Calculates posterior mean, sigma, and estimate of tau2 using kriging equations

krig <- function(y, dx, d_new = NULL, d_cross = NULL, theta, g, mean = TRUE,
                 sigma = TRUE, tau2 = TRUE) {

  out <- list()
  C <- calc_K(dx, theta, g)
  C_inv <- invdet(C)$Mi

  if (tau2) {
    scale <- c(t(y) %*% C_inv %*% y)/length(y)
    out$tau2 <- scale
  } else {
    scale <- 1
  }

  if (mean) {
    C_cross <- calc_K(d_cross, theta) # nugget not included in rectangular matrix
    out$mean <- C_cross %*% C_inv %*% y
  }

  if (sigma) {
    if (!mean) C_cross <- calc_K(d_cross, theta) # nugget not included in rectangular matrix
    C_new <- calc_K(d_new, theta, g)
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
            ldet = double(1),
            PACKAGE = "deepgp")

  return(list(Mi = matrix(out$Mi, ncol=n), ldet = out$ldet))
}

# Distance C Function ---------------------------------------------------------
#' @title Calculates squared pairwise distances
#' @description Calculates squared pairwise euclidean distances using C.
#' 
#' @details C code derived from the "laGP" package (Robert B Gramacy and Furong Sun).
#' @param X1 matrix of input locations
#' @param X2 matrix of second input locations (if \code{NULL}, distance is calculated
#'        between \code{X1} and itself)
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

  # coerse arguments and extract dimensions
  X1 <- as.matrix(X1)
  n1 <- nrow(X1)
  m <- ncol(X1)

  if(is.null(X2)) {
    # calculate D
    outD <- .C("distance_symm_R",
               X = as.double(t(X1)),
               n = as.integer(n1),
               m = as.integer(m),
               D = double(n1 * n1),
               PACKAGE = "deepgp")

    # return the distance matrix
    return(matrix(outD$D, ncol = n1, byrow = TRUE))

  } else {
    # coerse arguments and extract dimensions
    X2 <- as.matrix(X2)
    n2 <- nrow(X2)

    # check inputs
    if(ncol(X1) != ncol(X2)) stop("dimension of X1 & X2 do not match")

    # calculate D
    outD <- .C("distance_R",
               X1 = as.double(t(X1)),
               n1 = as.integer(n1),
               X2 = as.double(t(X2)),
               n2 = as.integer(n2),
               m = as.integer(m),
               D = double(n1 * n2),
               PACKAGE = "deepgp")

    # return the distance matrix
    return(matrix(outD$D, ncol = n2, byrow = TRUE))
  }
}

# Calculate Tau Hat Function --------------------------------------------------
# Calculates MLE estimate for tau2

calc_tau2 <- function(y, x, theta, g) {
  
  n <- length(y) # sample size
  
  K <- calc_K(sq_dist(x), theta, g)
  id <- invdet(K)
  tau2 <- (t(y) %*% id$Mi %*% y)/n
  
  return(c(tau2))
}

# Calculate Score Function ----------------------------------------------------
#' @title Calculates score
#' @description Calculates score (higher scores indicate better fits).  Only 
#'     applicable to noisy data.
#' @param y response vector
#' @param mu predicted mean
#' @param sigma predicted covariance
#' 
#' @references 
#' Gneiting, T, and AE Raftery. 2007. “Strictly Proper Scoring Rules, Prediction, 
#'     and Estimation.” \emph{Journal of the American Statistical Association 102} 
#'     (477), 359-378.
#' 
#' @examples
#' # See "deepgp-package", "fit_one_layer", "fit_two_layer", or "fit_three_layer"
#' # for an example
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

