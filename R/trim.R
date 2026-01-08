
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   trim.gp == trim.gpvec
#   trim.dgp2 == trim.dgp2vec
#   trim.dgp3 == trim.dgp3vec

# trim S3 class ---------------------------------------------------------------
#' @title Trim/Thin MCMC iterations
#' @description Acts on a \code{gp}, \code{gpvec}, \code{dgp2}, \code{dgp2vec},
#'    \code{dgp3vec}, or \code{dgp3} object.
#'    Removes the specified number of MCMC iterations (starting at the first 
#'    iteration).  After these samples are removed, the remaining samples are
#'    optionally thinned.
#' 
#' @details The resulting object will have \code{nmcmc} equal to the previous 
#'     \code{nmcmc} minus \code{burn} divided by \code{thin}.  It is 
#'     recommended to start an MCMC fit then investigate trace plots to assess 
#'     burn-in.  Once burn-in has been achieved, use this function to remove 
#'     the starting iterations.  Thinning reduces the size of the resulting 
#'     object while accounting for the high correlation between consecutive 
#'     iterations.
#'     
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer}
#' @param burn integer specifying number of iterations to cut off as burn-in
#' @param thin integer specifying amount of thinning (\code{thin = 1} keeps all 
#'        iterations, \code{thin = 2} keeps every other iteration, 
#'        \code{thin = 10} keeps every tenth iteration, etc.)
#'        
#' @return object of the same class with the selected iterations removed
#' 
#' @examples 
#' # See ?fit_one_layer, ?fit_two_layer, or ?fit_three_layer
#' # for examples
#' 
#' @rdname trim
#' @export

trim <- function(object, burn, thin)
  UseMethod("trim", object)

# trim.gp ---------------------------------------------------------------------
#' @rdname trim
#' @export

trim.gp <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[3]
  
  if (burn >= object$nmcmc) stop("burn must be less than nmcmc")
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  if (length(object$g) > 1) object$g <- object$g[indx, drop = FALSE]
  if (object$settings$sep) {
    object$theta <- object$theta[indx, , drop = FALSE]
  } else object$theta <- object$theta[indx, drop = FALSE]
  object$tau2 <- object$tau2[indx, drop = FALSE]
  object$ll <- object$ll[indx, drop = FALSE]
  
  toc <- proc.time()[3]
  object$time <- object$time + unname(toc - tic)
  
  return(object)
}

# trim.gpvec ------------------------------------------------------------------
#' @rdname trim
#' @export

trim.gpvec <- trim.gp

# Trim Two Layer --------------------------------------------------------------
#' @rdname trim
#' @export

trim.dgp2 <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[3]
  
  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  if (length(object$g) > 1) object$g <- object$g[indx, drop = FALSE]
  object$tau2_y <- object$tau2_y[indx, drop = FALSE]
  object$theta_y <- object$theta_y[indx, drop = FALSE]
  if (object$settings$monowarp) {
    if (!is.null(object$tau2_w)) {
      object$tau2_w <- object$tau2_w[indx, , drop = FALSE]
    }
    object$w_grid <- object$w_grid[indx, , , drop = FALSE]
  }
  object$theta_w <- object$theta_w[indx, , drop = FALSE]
  object$w <- object$w[indx, , , drop = FALSE]
  object$ll <- object$ll[indx, drop = FALSE]
  
  toc <- proc.time()[3]
  object$time <- object$time + unname(toc - tic)
  
  return(object)
}

# trim.dgp2vec ----------------------------------------------------------------
#' @rdname trim
#' @export

trim.dgp2vec <- trim.dgp2

# trim.dgp3 -------------------------------------------------------------------
#' @rdname trim
#' @export

trim.dgp3 <- function(object, burn, thin = 1) {
  
  tic <- proc.time()[3]
  
  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')
  
  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]
  
  object$nmcmc <- length(indx)
  if (length(object$g) > 1) object$g <- object$g[indx, drop = FALSE]
  object$tau2_y <- object$tau2_y[indx, drop = FALSE]
  object$theta_y <- object$theta_y[indx, drop = FALSE]
  object$theta_w <- object$theta_w[indx, , drop = FALSE]
  object$theta_z <- object$theta_z[indx, , drop = FALSE]
  object$w <- object$w[indx, , , drop = FALSE]
  object$z <- object$z[indx, , , drop = FALSE]
  object$ll <- object$ll[indx, drop = FALSE]
  
  toc <- proc.time()[3]
  object$time <- object$time + unname(toc - tic)
  
  return(object)
}

# trim.dgp3vec ----------------------------------------------------------------
#' @rdname trim
#' @export

trim.dgp3vec <- trim.dgp3
