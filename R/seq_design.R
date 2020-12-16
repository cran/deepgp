
# Function Contents -----------------------------------------------------------
# Internal:
#   alc.C: calculates ALC using C
#   Wij.C: calculates Wij matrix using C
# External (see documentation below):
#   ALC (S3 method for gp, dgp2, dgp3 classes)
#   IMSE (S3 method for gp, dgp2, dgp3 classes)
#   EI (S3 method for gp, dgp2, dgp3 classes)

# Imported Functions ----------------------------------------------------------
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach

# ALC C Function --------------------------------------------------------------
# Calculates vector of ALC values in C

alc.C <- function(X, Ki, theta, g, Xcand, Xref, tau2, verb = 0) {
  result <- .C("alcGP_R",
               X = as.double(t(X)),
               n = as.integer(nrow(X)),
               col = as.integer(ncol(X)),
               Ki = as.double(t(Ki)),
               d = as.double(theta),
               g = as.double(g),
               ncand = as.integer(nrow(Xcand)),
               Xcand = as.double(t(Xcand)),
               nref = as.integer(nrow(Xref)),
               Xref = as.double(t(Xref)),
               phi = as.double(tau2),
               verb = as.integer(verb),
               alc = double(nrow(Xcand)),
               PACKAGE = "deepgp")
  return(result$alc)
}

# Define ALC for S3 Objects ---------------------------------------------------
#' @title Active Learning Cohn for Sequential Design
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object. 
#'    Calculates ALC over the input locations \code{x_new} using specified
#'    reference grid.  If no reference grid is specified, \code{x_new} is used as
#'    the reference.  Optionally utilizes SNOW parallelization.  User should 
#'    select the point with the highest ALC to add to the design.   
#'    
#' @details All iterations in the object are used in the calculation, so samples 
#'     should be burned-in.  Thinning the samples using \code{trim} will
#'     speed up computation.  This function may be used in two ways:
#'     \itemize{
#'         \item called on an object with only MCMC iterations, in which case 
#'               \code{x_new} must be specified
#'         \item called on an object that has been predicted over, in which 
#'               case the \code{x_new} from \code{predict} is used
#'     }
#'     In \code{dgp2} and \code{dgp3} objects that have been run through \code{predict},
#'     the stored \code{w_new} mappings are used.  Through \code{predict}, the user may
#'     specify a mean mapping (\code{mean_map = TRUE}) or a full sample from the MVN distribution
#'     over \code{w_new} (\code{mean_map = FALSE}).  When the object has not yet been predicted
#'     over, the mean mapping is used.
#'     
#'     SNOW parallelization reduces computation time but requires more memory
#'     storage.  C code derived from the "laGP" package (Robert B Gramacy and Furong Sun).
#' 
#' @param object object of class \code{gp}, \code{dgp2}, or \code{dgp3}
#' @param x_new matrix of possible input locations, if object has been run 
#'        through \code{predict} the previously stored \code{x_new} is used
#' @param ref optional reference grid for ALC approximation, if \code{ref = NULL} 
#'        then \code{x_new} is used
#' @param cores number of cores to utilize in parallel, by default no parallelization
#'        is used
#' @return list with elements:
#' \itemize{
#'   \item \code{value}: vector of ALC values, indices correspond to \code{x_new}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." arXiv:2012.08015. \cr\cr
#' Seo, S, M Wallat, T Graepel, and K Obermayer. 2000. “Gaussian Process Regression:
#'     Active Data Selection and Test Point Rejection.” In Mustererkennung 2000, 
#'     27–34. New York, NY: Springer–Verlag.\cr\cr
#' Gramacy, RB and F Sun. (2016). laGP: Large-Scale Spatial Modeling via Local 
#'     Approximate Gaussian Processes in R. \emph{Journal of Statistical Software 
#'     72} (1), 1-46. doi:10.18637/jss.v072.i01
#' 
#' @examples
#' # See "deepgp-package" or "fit_two_layer" for an example
#' 
#' @rdname ALC
#' @export

ALC <- function(object, x_new, ref, cores)
  UseMethod("ALC", object)

# ALC One Layer Function ------------------------------------------------------
#' @rdname ALC
#' @export

ALC.gp <- function(object, x_new = NULL, ref = NULL, cores = 1) {

  tic <- proc.time()[3]
  
  if (is.null(object$x_new)) predicted <- FALSE else predicted <- TRUE
  
  # use previously stored x_new or newly specified x_new
  if (predicted) {
    if (!is.null(x_new)) warning('using x_new that was previously stored')
    x_new <- object$x_new
  } else {
    if (is.null(x_new)) stop('x_new has not been specified')
    if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  }
  
  # check that ref is a matrix
  if (!is.null(ref) & !is.matrix(ref)) stop('ref must be a matrix')
  if (is.null(ref)) ref <- x_new
  
  dx <- sq_dist(object$x)
  
  if (cores == 1) {
    alc <- rep(0, times = nrow(x_new))
    
    for(t in 1:object$nmcmc) {
      
      K <- calc_K(dx, theta = object$theta[t], g = object$g[t])
      Ki <- invdet(K)$Mi
      if (predicted) { tau2 <- object$tau2[t] 
      } else tau2 <- krig(object$y, dx, theta = object$theta[t], g = object$g[t], 
                          mean = FALSE, sigma = FALSE)$tau2
      
      alc <- alc + alc.C(object$x, Ki, object$theta[t], object$g[t], x_new, ref, tau2)
    } # end of t for loop
  } else {
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  
    alc <- foreach(t = 1:object$nmcmc, .combine = '+') %dopar% {
    
      K <- calc_K(dx, theta = object$theta[t], g = object$g[t])
      Ki <- invdet(K)$Mi
      if (predicted) { tau2 <- object$tau2[t] 
      } else tau2 <- krig(object$y, dx, theta = object$theta[t], g = object$g[t], 
                        mean = FALSE, sigma = FALSE)$tau2
    
      return(alc.C(object$x, Ki, object$theta[t], object$g[t], x_new,
                  ref, tau2))
    } # end of foreach statement
  
    stopCluster(cl)
  } # end of else statement

  toc <- proc.time()[3]
  return(list(value = alc / object$nmcmc, time = toc - tic))
}

# ALC Two Layer Function ------------------------------------------------------
#' @rdname ALC
#' @export

ALC.dgp2 <- function(object, x_new = NULL, ref = NULL, cores = 1) {

  tic <- proc.time()[3]
  
  if (is.null(object$x_new)) predicted <- FALSE else predicted <- TRUE
  
  # use previously stored x_new or newly specified x_new
  if (predicted) {
    if (!is.null(x_new)) warning('using x_new that was previously stored')
    x_new <- object$x_new
  } else {
    if (is.null(x_new)) stop('x_new has not been specified')
    if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  }
  
  # check that ref is a matrix
  if (!is.null(ref) & !is.matrix(ref)) stop('ref must be a matrix')
  
  # specify pre-calculations if predicted is FALSE
  if (!predicted) {
    m <- nrow(x_new)
    D <- ncol(object$w[[1]])
    dx <- sq_dist(object$x)
    d_cross <- sq_dist(x_new, object$x)
  }
  
  if (cores == 1) {
    alc <- rep(0, times = nrow(x_new))
    
    for(t in 1:object$nmcmc) {
      w <- object$w[[t]]
      
      if (predicted) {
        w_new <- object$w_new[[t]]
        tau2 <- object$tau2[t]
      } else {
        # calculate w_new using conditional uncertainty
        w_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D) {
          w_new[, i] <- krig(w[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                             g = NULL, sigma = FALSE, tau2 = FALSE)$mean
        }
        # calculate tau2
        tau2 <- krig(object$y, sq_dist(w), theta = object$theta_y[t], 
                     g = object$g[t], mean = FALSE, sigma = FALSE)$tau2
      } 
      
      if (is.null(ref)) ref <- w_new
      
      K <- calc_K(sq_dist(w), theta = object$theta_y[t], g = object$g[t])
      Ki <- invdet(K)$Mi
      
      alc <- alc + alc.C(w, Ki, object$theta_y[t], object$g[t], w_new, ref, tau2)
    } # end of t for loop
  } else {
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)

    alc <- foreach(t = 1:object$nmcmc, .combine = '+') %dopar% {
      w <- object$w[[t]]

      if (predicted) {
        w_new <- object$w_new[[t]]
        tau2 <- object$tau2[t]
      } else {
        # calculate w_new using conditional uncertainty
        w_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D) {
          w_new[, i] <- krig(w[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                             g = NULL, sigma = FALSE, tau2 = FALSE)$mean
        }
        # calculate tau2
        tau2 <- krig(object$y, sq_dist(w), theta = object$theta_y[t], 
                     g = object$g[t], mean = FALSE, sigma = FALSE)$tau2
      } 
  
      if (is.null(ref)) ref <- w_new
    
      K <- calc_K(sq_dist(w), theta = object$theta_y[t], g = object$g[t])
      Ki <- invdet(K)$Mi
    
      return(alc.C(w, Ki, object$theta_y[t], object$g[t], w_new,
                   ref, tau2))
    } # end of foreach statement
  
    stopCluster(cl)
  } # end of else statement

  toc <- proc.time()[3]
  return(list(value = alc / object$nmcmc, time = toc - tic))
}

# ALC Three Layer Function ----------------------------------------------------
#' @rdname ALC
#' @export

ALC.dgp3 <- function(object, x_new = NULL, ref = NULL, cores = 1) {
  
  tic <- proc.time()[3]
  
  if (is.null(object$x_new)) predicted <- FALSE else predicted <- TRUE
  
  # use previously stored x_new or newly specified x_new
  if (predicted) {
    if (!is.null(x_new)) warning('using x_new that was previously stored')
    x_new <- object$x_new
  } else {
    if (is.null(x_new)) stop('x_new has not been specified')
    if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  }
  
  # check that w_ref is a matrix
  if (!is.null(ref) & !is.matrix(ref)) stop('ref must be a matrix')
  
  # specify pre-calculations if predicted is FALSE
  if (!predicted) {
    m <- nrow(x_new)
    D <- ncol(object$w[[1]])
    dx <- sq_dist(object$x)
    d_cross <- sq_dist(x_new, object$x)
  }
  
  if (cores == 1) {
    alc <- rep(0, times = nrow(x_new))
    
    for(t in 1:object$nmcmc) {
      w <- object$w[[t]]
      
      if (predicted) {
        w_new <- object$w_new[[t]]
        tau2 <- object$tau2[t]
      } else {
        # calculate z_new using conditional uncertainty
        z <- object$z[[t]]
        z_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D) {
          z_new[, i] <- krig(z[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                             g = NULL, sigma = FALSE, tau2 = FALSE)$mean
        }
        
        # calculate w_new using conditional uncertainty
        w_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D) {
          w_new[, i] <- krig(w[, i], sq_dist(z), NULL, sq_dist(z_new, z), 
                             object$theta_w[t, i], g = NULL, sigma = FALSE, 
                             tau2 = FALSE)$mean
        }
        # calculate tau2
        tau2 <- krig(object$y, sq_dist(w), theta = object$theta_y[t], 
                     g = object$g[t], mean = FALSE, sigma = FALSE)$tau2
      } 
      
      if (is.null(ref)) ref <- w_new
      
      K <- calc_K(sq_dist(w), theta = object$theta_y[t], g = object$g[t])
      Ki <- invdet(K)$Mi
      
      alc <- alc + alc.C(w, Ki, object$theta_y[t], object$g[t], w_new, ref, tau2)
    } # end of t for loop
  } else {
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  
    alc <- foreach(t = 1:object$nmcmc, .combine = '+') %dopar% {
      w <- object$w[[t]]
    
      if (predicted) {
        w_new <- object$w_new[[t]]
        tau2 <- object$tau2[t]
      } else {
        # calculate z_new using conditional uncertainty
        z <- object$z[[t]]
        z_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D) {
          z_new[, i] <- krig(z[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                             g = NULL, sigma = FALSE, tau2 = FALSE)$mean
        }
      
        # calculate w_new using conditional uncertainty
        w_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D) {
          w_new[, i] <- krig(w[, i], sq_dist(z), NULL, sq_dist(z_new, z), 
                             object$theta_w[t, i], g = NULL, sigma = FALSE, 
                             tau2 = FALSE)$mean
        }
        # calculate tau2
        tau2 <- krig(object$y, sq_dist(w), theta = object$theta_y[t], 
                     g = object$g[t], mean = FALSE, sigma = FALSE)$tau2
      } 
    
      if (is.null(ref)) ref <- w_new
    
      K <- calc_K(sq_dist(w), theta = object$theta_y[t], g = object$g[t])
      Ki <- invdet(K)$Mi
    
      return(alc.C(w, Ki, object$theta_y[t], object$g[t], w_new,
                   ref, tau2))
    } # end of foreach statement
  
    stopCluster(cl)
  } # end of else statement
  
  toc <- proc.time()[3]
  return(list(value = alc / object$nmcmc, time = toc - tic))
}


# Wij C Function --------------------------------------------------------------
# Calculates the Wij matrix in C

Wij.C <- function(x1, x2, theta, a, b){
  W <- matrix(1, nrow = nrow(x1), ncol = nrow(x2))
  result <- .C("Wij_R",
               X1 = as.double(t(x1)),
               n1 = as.integer(nrow(x1)),
               X2 = as.double(t(x2)),
               n2 = as.integer(nrow(x2)),
               col = as.integer(ncol(x1)),
               theta = as.double(theta),
               a = as.double(a),
               b = as.double(b),
               W = as.double(t(W)),
               PACKAGE = "deepgp")
  return(matrix(result$W, nrow = nrow(x1), ncol = nrow(x2)))
}

# Define IMSE for S3 Objects -------------------------------------------------
#' @title Integrated Mean-Squared (prediction) Error for Sequential Design
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.
#'     Calculates IMSE over the input locations \code{x_new}.  Optionally 
#'     utilizes SNOW parallelization.  User should select the point with the 
#'     lowest IMSE to add to the design.
#'     
#' @details All iterations in the object are used in the calculation, so samples
#'     should be burned-in.  Thinning the samples using \code{trim} will speed 
#'     up computation.  This function may be used in two ways:
#'     \itemize{
#'         \item called on an object with only MCMC iterations, in which case 
#'               \code{x_new} must be specified
#'         \item called on an object that has been predicted over, in which case
#'         the \code{x_new} from \code{predict} is used
#'     }
#'     In \code{dgp2} and \code{dgp3} objects that have been run through \code{predict},
#'     the stored \code{w_new} mappings are used.  Through \code{predict}, the user may
#'     specify a mean mapping (\code{mean_map = TRUE}) or a full sample from the MVN distribution
#'     over \code{w_new} (\code{mean_map = FALSE}).  When the object has not yet been predicted
#'     over, the mean mapping is used.
#'     
#'     SNOW parallelization reduces computation time but requires more memory storage.
#' 
#' @param object object of class \code{gp}, \code{dgp2}, or \code{dgp3}
#' @param x_new matrix of possible input locations, if object has been run 
#'        through \code{predict} the previously stored \code{x_new} is used
#' @param cores number of cores to utilize in parallel, by default no 
#'        parallelization is used
#' @return list with elements:
#' \itemize{
#'   \item \code{value}: vector of IMSE values, indices correspond to \code{x_new}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." arXiv:2012.08015. \cr\cr
#' Binois, M, J Huang, RB Gramacy, and M Ludkovski. 2019. “Replication or Exploration? 
#'     Sequential Design for Stochastic Simulation Experiments.” \emph{Technometrics 
#'     61}, 7-23. Taylor & Francis. doi:10.1080/00401706.2018.1469433.
#' 
#' @examples
#' # See "deepgp-package" or "fit_three_layer" for an example
#' 
#' @rdname IMSE
#' @export

IMSE <- function(object, x_new, cores)
  UseMethod("IMSE", object)

# IMSE One Layer Function ----------------------------------------------------
#' @rdname IMSE
#' @export

IMSE.gp <- function(object, x_new = NULL, cores = 1) {

  tic <- proc.time()[3]

  if (is.null(object$x_new)) predicted <- FALSE else predicted <- TRUE
  
  # use previously stored x_new or newly specified x_new
  if (predicted) {
    if (!is.null(x_new)) warning('using x_new that was previously stored')
    x_new <- object$x_new
  } else {
    if (is.null(x_new)) stop('x_new has not been specified')
    if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  }

  n <- nrow(object$x)
  m <- nrow(x_new)
  dx <- sq_dist(object$x)
  
  # define bounds
  a <- apply(x_new, 2, min)
  b <- apply(x_new, 2, max)

  Knew_inv <- matrix(nrow = n + 1, ncol = n + 1)
  Wijs <- matrix(nrow = n + 1, ncol = n + 1)
  
  if (cores == 1) {
    imse <- rep(0, times = nrow(x_new))
    
    for (t in 1:object$nmcmc) {
      Kn <- calc_K(dx, theta = object$theta[t], g = object$g[t])
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
      
      # precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(object$x, object$x, object$theta[t], a, b)
      
      if (predicted) { tau2 <- object$tau2[t] 
      } else tau2 <- krig(object$y, dx, theta = object$theta[t], g = object$g[t], 
                          mean = FALSE, sigma = FALSE)$tau2
      
      imse_store <- vector(length = m)
      
      for (i in 1:m) {
        # specify new design point
        x_star <- matrix(x_new[i, ], nrow = 1)
        
        # calculate new Ki matrix
        k <- calc_K(sq_dist(object$x, x_star), theta = object$theta[t])
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v
        
        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(object$x, x_star, 
                                                  object$theta[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(x_star, x_star, object$theta[t], a, b)
        imse_store[i] <- tau2 * prod(b - a) * (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
      imse <- imse + imse_store
    } # end of t for loop
  } else {
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)

    imse <- foreach(t = 1:object$nmcmc, .combine = '+') %dopar% {
    
      Kn <- calc_K(dx, theta = object$theta[t], g = object$g[t])
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
  
      # precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(object$x, object$x, object$theta[t], a, b)
      
      if (predicted) { tau2 <- object$tau2[t] 
      } else tau2 <- krig(object$y, dx, theta = object$theta[t], g = object$g[t], 
                          mean = FALSE, sigma = FALSE)$tau2
  
      imse_store <- vector(length = m)
  
      for (i in 1:m) {
        # specify new design point
        x_star <- matrix(x_new[i, ], nrow = 1)
    
        # calculate new Ki matrix
        k <- calc_K(sq_dist(object$x, x_star), theta = object$theta[t])
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v
    
        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(object$x, x_star, 
                                                  object$theta[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(x_star, x_star, object$theta[t], a, b)
        imse_store[i] <- tau2 * prod(b - a) * (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
      return(imse_store)
    } # end of foreach statement
  
    stopCluster(cl)
  } # end of else statement
  
  toc <- proc.time()[3]

  return(list(value = imse / object$nmcmc, time = toc - tic))
}

# IMSE Two Layer Function ----------------------------------------------------
#' @rdname IMSE
#' @export

IMSE.dgp2 <- function(object, x_new = NULL, cores = 1) {

  tic <- proc.time()[3]

  if (is.null(object$x_new)) predicted <- FALSE else predicted <- TRUE
  
  # use previously stored x_new or newly specified x_new
  if (predicted) {
    if (!is.null(x_new)) warning('using x_new that was previously stored')
    x_new <- object$x_new
  } else {
    if (is.null(x_new)) stop('x_new has not been specified')
    if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  }

  n <- nrow(object$x)
  m <- nrow(x_new)
  if (!predicted) {
    D <- ncol(object$w[[1]])
    dx <- sq_dist(object$x)
    d_cross <- sq_dist(x_new, object$x)
  }

  Knew_inv <- matrix(nrow = n + 1, ncol = n + 1)
  Wijs <- matrix(nrow = n + 1, ncol = n + 1)
  
  if (cores == 1) {
    imse <- rep(0, times = nrow(x_new))
    
    for (t in 1:object$nmcmc) {
      w <- object$w[[t]]
      Kn <- calc_K(sq_dist(w), theta = object$theta_y[t], g = object$g[t])
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
      
      if (predicted) {
        w_new <- object$w_new[[t]]
        tau2 <- object$tau2[t]
      } else {
        w_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D)
          w_new[, i] <- krig(w[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                             g = NULL, sigma = FALSE, tau2 = FALSE)$mean
        tau2 <- krig(object$y, sq_dist(w), theta = object$theta_y[t], 
                     g = object$g[t], mean = FALSE, sigma = FALSE)$tau2
      }
      
      # define bounds
      a <- apply(w_new, 2, min)
      b <- apply(w_new, 2, max)
      
      # precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(w, w, object$theta_y[t], a, b)
      
      imse_store <- vector(length = m)
      
      for (i in 1:m) {
        # specify new design point
        w_star <- matrix(w_new[i, ], nrow = 1)
        
        # calculate new Ki matrix
        k <- calc_K(sq_dist(w, w_star), theta = object$theta_y[t])
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v
        
        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(w, w_star, 
                                                  object$theta_y[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(w_star, w_star, object$theta_y[t], a, b)
        imse_store[i] <- tau2 * prod(b - a) * (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
      imse <- imse + imse_store
    } # end of t for loop
  } else {
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  
    imse <- foreach(t = 1:object$nmcmc, .combine = '+') %dopar% {
      w <- object$w[[t]]
      Kn <- calc_K(sq_dist(w), theta = object$theta_y[t], g = object$g[t])
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
      
      if (predicted) {
        w_new <- object$w_new[[t]]
        tau2 <- object$tau2[t]
      } else {
        w_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D)
          w_new[, i] <- krig(w[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                             g = NULL, sigma = FALSE, tau2 = FALSE)$mean
        tau2 <- krig(object$y, sq_dist(w), theta = object$theta_y[t],
                     g = object$g[t], mean = FALSE, sigma = FALSE)$tau2
      }
      
      # define bounds
      a <- apply(w_new, 2, min)
      b <- apply(w_new, 2, max)

      # precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(w, w, object$theta_y[t], a, b)

      imse_store <- vector(length = m)

      for (i in 1:m) {
        # specify new design point
        w_star <- matrix(w_new[i, ], nrow = 1)

        # calculate new Ki matrix
        k <- calc_K(sq_dist(w, w_star), theta = object$theta_y[t])
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v

        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(w, w_star, object$theta_y[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(w_star, w_star, object$theta_y[t], a, b)
        imse_store[i] <- tau2 * prod(b - a) * (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
      return(imse_store)
    } # end of foreach statement
  
    stopCluster(cl)
  } # end of else statement
  
  toc <- proc.time()[3]

  return(list(value = imse / object$nmcmc, time = toc - tic))
}

# IMSE Three Layer Function ---------------------------------------------------
#' @rdname IMSE
#' @export

IMSE.dgp3 <- function(object, x_new = NULL, cores = 1) {
  
  tic <- proc.time()[3]
  
  if (is.null(object$x_new)) predicted <- FALSE else predicted <- TRUE
  
  # use previously stored x_new or newly specified x_new
  if (predicted) {
    if (!is.null(x_new)) warning('using x_new that was previously stored')
    x_new <- object$x_new
  } else {
    if (is.null(x_new)) stop('x_new has not been specified')
    if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  }
  
  n <- nrow(object$x)
  m <- nrow(x_new)
  if (!predicted) {
    D <- ncol(object$w[[1]])
    dx <- sq_dist(object$x)
    d_cross <- sq_dist(x_new, object$x)
  }
  
  Knew_inv <- matrix(nrow = n + 1, ncol = n + 1)
  Wijs <- matrix(nrow = n + 1, ncol = n + 1)
  
  if (cores == 1) {
    imse <- rep(0, times = nrow(x_new))
    
    for (t in 1:object$nmcmc) {
      w <- object$w[[t]]
      Kn <- calc_K(sq_dist(w), theta = object$theta_y[t], g = object$g[t])
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
      
      if (predicted) {
        w_new <- object$w_new[[t]]
        tau2 <- object$tau2[t]
      } else {
        # calculate z_new using conditional uncertainty
        z <- object$z[[t]]
        z_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D)
          z_new[, i] <- krig(z[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                             g = NULL, sigma = FALSE, tau2 = FALSE)$mean
        w_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D)
          w_new[, i] <- krig(w[, i], sq_dist(z), NULL, sq_dist(z_new, z), 
                             object$theta_w[t, i], g = NULL, sigma = FALSE, 
                             tau2 = FALSE)$mean
        tau2 <- krig(object$y, sq_dist(w), theta = object$theta_y[t], 
                     g = object$g[t], mean = FALSE, sigma = FALSE)$tau2
      }
      
      # define bounds
      a <- apply(w_new, 2, min)
      b <- apply(w_new, 2, max)
      
      # precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(w, w, object$theta_y[t], a, b)
      
      imse_store <- vector(length = m)
      
      for (i in 1:m) {
        # specify new design point
        w_star <- matrix(w_new[i, ], nrow = 1)
        
        # calculate new Ki matrix
        k <- calc_K(sq_dist(w, w_star), theta = object$theta_y[t])
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v
        
        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(w, w_star, 
                                                  object$theta_y[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(w_star, w_star, object$theta_y[t], a, b)
        imse_store[i] <- tau2 * prod(b - a) * (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
      imse <- imse + imse_store
    } # end of t for loop
  } else {
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    imse <- foreach(t = 1:object$nmcmc, .combine = '+') %dopar% {
      w <- object$w[[t]]
      Kn <- calc_K(sq_dist(w), theta = object$theta_y[t], g = object$g[t])
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
      
      if (predicted) {
        w_new <- object$w_new[[t]]
        tau2 <- object$tau2[t]
      } else {
        # calculate z_new using conditional uncertainty
        z <- object$z[[t]]
        z_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D)
          z_new[, i] <- krig(z[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                             g = NULL, sigma = FALSE, tau2 = FALSE)$mean
        w_new <- matrix(nrow = m, ncol = D)
        for (i in 1:D)
          w_new[, i] <- krig(w[, i], sq_dist(z), NULL, sq_dist(z_new, z), 
                             object$theta_w[t, i], g = NULL, sigma = FALSE, 
                             tau2 = FALSE)$mean
        tau2 <- krig(object$y, sq_dist(w), theta = object$theta_y[t], 
                     g = object$g[t], mean = FALSE, sigma = FALSE)$tau2
      }
      
      # define bounds
      a <- apply(w_new, 2, min)
      b <- apply(w_new, 2, max)
      
      # precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(w, w, object$theta_y[t], a, b)
      
      imse_store <- vector(length = m)
      
      for (i in 1:m) {
        # specify new design point
        w_star <- matrix(w_new[i, ], nrow = 1)
        
        # calculate new Ki matrix
        k <- calc_K(sq_dist(w, w_star), theta = object$theta_y[t])
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v
        
        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(w, w_star, 
                                                  object$theta_y[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(w_star, w_star, object$theta_y[t], a, b)
        imse_store[i] <- tau2 * prod(b - a) * (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
      return(imse_store)
    } # end of foreach statement
    
    stopCluster(cl)
  } # end of else statement
  
  toc <- proc.time()[3]
  
  return(list(value = imse / object$nmcmc, time = toc - tic))
}


# Define EI for S3 Objects ----------------------------------------------------
#' @title Expected Improvement for Sequential Design
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.  
#'     Calculates expected improvement over input locations \code{x_new} with 
#'     the goal of MINIMIZING the function.  Optionally utilizes SNOW 
#'     parallelization.  User should select the point with the highest EI to 
#'     add to the design.
#'     
#' @details The object must be an output of \code{predict} with 
#'     \code{lite = TRUE} and \code{store_all = TRUE}.  This will store the posterior 
#'     mean and point-wise variance for every iteration.  Once prediction is done, computation is 
#'     relatively quick, so SNOW parallelization is only recommended when 
#'     \code{nmcmc} is large.
#' 
#' @param object object of class \code{gp}, \code{dgp2}, or \code{dgp3} that has been 
#'        predicted over \code{x_new} with \code{lite = TRUE} and \code{store_all = TRUE}
#' @param cores number of cores to utilize in parallel, by default no 
#'        parallelization is used
#' @return list with elements:
#' \itemize{
#'   \item \code{value}: vector of EI values, indices correspond to \code{x_new}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Jones, DR, M Schonlau, and WJ Welch. 1998. "Efficient Global Optimization of 
#'     Expensive Black-Box Functions." \emph{Journal of Global Optimization 13}, 
#'     455-492. doi:10.1023/A:1008306431147.
#' 
#' @examples
#' # See "deepgp-package" or "fit_one_layer" for an example
#' 
#' @rdname EI
#' @export

EI <- function(object, cores)
  UseMethod("EI", object)


# EI One Layer Function -------------------------------------------------------
#' @rdname EI
#' @export

EI.gp <- function(object, cores = 1) {

  tic <- proc.time()[3]

  if (is.null(object$mu_t)) stop('object must be predicted with store_all = TRUE')
  
  if (cores == 1) {
    ei <- rep(0, times = nrow(object$x_new))
    for (t in 1:object$nmcmc) {
      mu <- object$mu_t[,t]
      sig <- sqrt(object$s2_t[,t])
      f_min <- min(mu)
      
      ei <- ei + (f_min - mu) * pnorm(f_min, mean = mu, sd = sig) +
                  sig * dnorm(f_min, mean = mu, sd = sig)
    }
  } else {
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
  
    ei <- foreach(t = 1:object$nmcmc, .combine = '+') %dopar% {

      mu <- object$mu_t[, t]
      sig <- sqrt(object$s2_t[, t])
      f_min <- min(mu) # use min of predicted mean

      return((f_min - mu) * pnorm(f_min, mean = mu, sd = sig) +
              sig * dnorm(f_min, mean = mu, sd = sig))
    }
  
    stopCluster(cl)
  } # end of else statement
  
  toc <- proc.time()[3]

  return(list(value = ei / object$nmcmc, time = toc - tic))
}

# EI Two Layer Function (same as one layer) -----------------------------------
#' @rdname EI
#' @export

EI.dgp2 <- EI.gp

# EI Three Layer Function (same as one layer) ---------------------------------
#' @rdname EI
#' @export

EI.dgp3 <- EI.gp


