
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   continue.gp
#   continue.dgp2
#   continue.dgp3

# Define Continue for S3 Objects ----------------------------------------------
#' @title Continues MCMC sampling
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.  
#'     Continues MCMC sampling of hyperparameters and hidden layers and appends
#'     results to existing object.
#'     
#' @details See "\code{fit_one_layer}", "\code{fit_two_layer}", or 
#'     "\code{fit_three_layer}" for details on MCMC.  The resulting object 
#'     will have \code{nmcmc} equal to the previous \code{nmcmc} plus 
#'     \code{new_mcmc}.  It is recommended to start an MCMC fit then 
#'     investigate trace plots to assess burn-in.  The primary use of this 
#'     function is to gather more MCMC iterations in order to obtain burned-in 
#'     samples.
#' 
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer}
#' @param new_mcmc number of MCMC iterations to conduct and append
#' @param trace logical indicating whether to print iteration progress
#' @return object of the same class with the new iterations appended
#' 
#' @examples 
#' # See "deepgp-package" or "fit_two_layer" for an example
#' 
#' @rdname continue
#' @export

continue <- function(object, new_mcmc, trace)
  UseMethod("continue", object)

# Continue One Layer Function -------------------------------------------------
#' @rdname continue
#' @export

continue.gp <- function(object, new_mcmc = 1000, trace = TRUE) {
  
  tic <- proc.time()[3]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) {
    true_g = object$g[1] 
  } else true_g = NULL
  
  # Run continuing MCMC iterations
  new_fit <- fit_one_layer(object$x, object$y, nmcmc = new_mcmc, trace = trace,
                           g_0 = object$g[object$nmcmc], 
                           theta_0 = object$theta[object$nmcmc],
                           true_g = true_g, settings = object$settings)
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, new_fit$g)
  object$theta <- c(object$theta, new_fit$theta)
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}

# Continue Two Layer Function -------------------------------------------------
#' @rdname continue
#' @export

continue.dgp2 <- function(object, new_mcmc = 1000, trace = TRUE) {
  
  tic <- proc.time()[3]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) {
    true_g = object$g[1] 
  } else true_g = NULL
  
  # Run continuing MCMC iterations
  new_fit <- fit_two_layer(object$x, object$y, D = ncol(object$w[[1]]),
                           nmcmc = new_mcmc, trace = trace,
                           w_0 = object$w[[object$nmcmc]], 
                           g_0 = object$g[object$nmcmc],
                           theta_y_0 = object$theta_y[object$nmcmc],
                           theta_w_0 = object$theta_w[object$nmcmc, ], 
                           true_g = true_g,
                           settings = object$settings)
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, new_fit$g)
  object$theta_y <- c(object$theta_y, new_fit$theta_y)
  object$theta_w <- rbind(object$theta_w, new_fit$theta_w)
  object$w <- c(object$w, new_fit$w)
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}

# Continue Three Layer Function -----------------------------------------------
#' @rdname continue
#' @export

continue.dgp3 <- function(object, new_mcmc = 1000, trace = TRUE) {
  
  tic <- proc.time()[3]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) {
    true_g = object$g[1] 
  } else true_g = NULL
  
  # Run continuing MCMC iterations
  new_fit <- fit_three_layer(object$x, object$y, D = ncol(object$w[[1]]),
                             nmcmc = new_mcmc, trace = trace,
                             w_0 = object$w[[object$nmcmc]], 
                             z_0 = object$z[[object$nmcmc]],
                             g_0 = object$g[object$nmcmc], 
                             theta_y_0 = object$theta_y[object$nmcmc],
                             theta_w_0 = object$theta_w[object$nmcmc, ],
                             theta_z_0 = object$theta_z[object$nmcmc, ], 
                             true_g = true_g,
                             settings = object$settings)
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, new_fit$g)
  object$theta_y <- c(object$theta_y, new_fit$theta_y)
  object$theta_w <- rbind(object$theta_w, new_fit$theta_w)
  object$theta_z <- rbind(object$theta_z, new_fit$theta_z)
  object$w <- c(object$w, new_fit$w)
  object$z <- c(object$z, new_fit$z)
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}
