
# Function Contents -----------------------------------------------------------
# Internal:
#   eps: minimum nugget value (1.5e-8)
#   krig: calculates posterior mean, sigma, and tau2 (optionally f_min)
#   invdet: calculates inverse and log determinant of a matrix using C
#           (credit given to the "laGP package, R.B. Gramacy & Furong Sun)
#   fill_final_row: uses kriging to fill final row of w_0/z_0
#   clean_prediction: removes prediction elements from object
# External (see documentation below):
#   sq_dist
#   score
#   rmse

eps <- sqrt(.Machine$double.eps)

# Krig ------------------------------------------------------------------------

krig <- function(y, dx, d_new = NULL, d_cross = NULL, theta, g, tau2 = 1,
                 s2 = FALSE, sigma = FALSE, f_min = FALSE, v = 2.5) {
  
  out <- list()
  if (v == 999) {
    C <- Exp2Fun(dx, c(1, theta, g))
    C_cross <- Exp2Fun(d_cross, c(1, theta, 0)) # no g in rectangular matrix
  } else {
    C <- MaternFun(dx, c(1, theta, g, v)) 
    C_cross <- MaternFun(d_cross, c(1, theta, 0, v))
  }
  C_inv <- invdet(C)$Mi
  out$mean <- C_cross %*% C_inv %*% y
  
  if (f_min) { # predict at observed locations, return min expected value
    if (v == 999) {
      C_cross_observed_only <- Exp2Fun(dx, c(1, theta, 0))
    } else C_cross_observed_only <- MaternFun(dx, c(1, theta, 0, v))
    out$f_min <- min(C_cross_observed_only %*% C_inv %*% y)
  }
  
  if (s2) {
    quadterm <- C_cross %*% C_inv %*% t(C_cross)
    C_new <- rep(1 + g, times = nrow(d_new))
    out$s2 <- tau2 * (C_new - diag(quadterm))
  }
  
  if (sigma) {
    quadterm <- C_cross %*% C_inv %*% t(C_cross)
    if (v == 999) {
      C_new <- Exp2Fun(d_new, c(1, theta, g))
    } else  C_new <- MaternFun(d_new, c(1, theta, g, v)) 
    out$sigma <- tau2 * (C_new - quadterm)
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
            ldet = double(1))

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
#' @export

rmse <- function(y, mu) {
  return(sqrt(mean((y - mu) ^ 2)))
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

# Clean Prediction ------------------------------------------------------------

clean_prediction <- function(object) {
  
  if (!is.null(object$x_new)) 
    object <- object[-which(names(object) == "x_new")]
  if (!is.null(object$mean)) 
    object <- object[-which(names(object) == "mean")]
  if (!is.null(object$s2)) 
    object <- object[-which(names(object) == "s2")]
  if (!is.null(object$s2_smooth)) 
    object <- object[-which(names(object) == "s2_smooth")]
  if (!is.null(object$Sigma)) 
    object <- object[-which(names(object) == "Sigma")]
  if (!is.null(object$Sigma_smooth)) 
    object[-which(names(object) == "Sigma_smooth")]
  if (!is.null(object$EI)) 
    object[-which(names(object) == "EI")]
  
  return(object)
}
