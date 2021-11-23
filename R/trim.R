
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   trim.gp
#   trim.dgp2
#   trim.dgp3

# Define Trim for S3 Objects --------------------------------------------------
#' @title Trim/Thin MCMC iterations
#' @description Acts on a \code{gp}, \code{dgp2}, or \code{dgp3} object.
#'    Removes the specified number of MCMC iterations (starting at the first 
#'    iteration).  After these samples are removed, the remaining samples may 
#'    be thinned.
#' 
#' @details The resulting object will have \code{nmcmc} equal to the previous 
#'     \code{nmcmc} minus \code{burn} divided by \code{thin}.  It is recommended 
#'     to start an MCMC fit then investigate trace plots to assess burn-in.  Once
#'     burn-in has been achieved, use this function to remove the starting 
#'     iterations.  Thinning reduces the size of the resulting object while 
#'     accounting for the high correlation between consecutive iterations.
#'     
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer}
#' @param burn integer specifying number of iterations to cut off as burn-in
#' @param thin integer specifying amount of thinning (\code{thin = 1} keeps all 
#'        iterations, \code{thin = 2} keeps every other iteration, 
#'        \code{thin = 10} keeps every tenth iteration, etc.)
#' @return object of the same class with the selected iterations removed
#' 
#' @examples 
#' # See "deepgp-package", "fit_one_layer", "fit_two_layer", or "fit_three_layer"
#' # for an example
#' 
#' @rdname trim
#' @export

trim <- function(object, burn, thin)
  UseMethod("trim", object)


# Trim One Layer Function -----------------------------------------------------
#' @rdname trim
#' @export

trim.gp <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[3]
  
  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  object$g <- object$g[indx]
  object$theta <- object$theta[indx]
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

# Trim Two Layer Function -----------------------------------------------------
#' @rdname trim
#' @export

trim.dgp2 <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[3]
  
  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  object$g <- object$g[indx]
  object$theta_y <- object$theta_y[indx]
  object$theta_w <- as.matrix(object$theta_w[indx,])
  object$w <- object$w[indx]
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

# Trim Three Layer Function ---------------------------------------------------
#' @rdname trim
#' @export

trim.dgp3 <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[3]
  
  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  object$g <- object$g[indx]
  object$theta_y <- object$theta_y[indx]
  object$theta_w <- as.matrix(object$theta_w[indx,])
  object$theta_z <- as.matrix(object$theta_z[indx,])
  object$w <- object$w[indx]
  object$z <- object$z[indx]
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

