
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   continue.gp
#   continue.dgp2
#   continue.dgp3
#   continue.gpvec
#   continue.dgp2vec
#   continue.dgp3vec

# Define Continue for S3 Objects ----------------------------------------------
#' @title Continues MCMC sampling
#' @description Acts on a \code{gp}, \code{gpvec}, \code{dgp2}, 
#'     \code{dgp2vec}, \code{dgp3}, or \code{dgp3vec} object.  
#'     Continues MCMC sampling of hyperparameters and hidden layers using 
#'     settings from the original object.  Appends new samples to existing
#'     samples.  When \code{vecchia = TRUE}, this function provides the option
#'     to update Vecchia ordering/conditioning sets based on latent layer
#'     warpings through the specification of \code{re_approx = TRUE}.
#'     
#' @details See \code{fit_one_layer}, \code{fit_two_layer}, or 
#'     \code{fit_three_layer} for details on MCMC.  The resulting 
#'     object will have \code{nmcmc} equal to the previous \code{nmcmc} plus 
#'     \code{new_mcmc}.  It is recommended to start an MCMC fit then 
#'     investigate trace plots to assess burn-in.  The primary use of this 
#'     function is to gather more MCMC iterations in order to obtain burned-in 
#'     samples.
#'     
#'     Specifying \code{re_approx = TRUE} updates random orderings and 
#'     nearest-neighbor conditioning sets (only for \code{vecchia = TRUE} 
#'     fits).  In one-layer, there is no latent warping but the Vecchia 
#'     approximation is still re-randomized and nearest-neighbors are adjusted 
#'     accordingly.  In two- and three-layers, the latest samples of hidden 
#'     layers are used to update nearest-neighbors.  If you update the 
#'     Vecchia approximation, you should later remove previous samples 
#'     (updating the approximation effectively starts a new chain).  When 
#'     \code{re_approx = FALSE} the previous orderings and conditioning sets 
#'     are used (maintaining the continuity of the previous chain).
#' 
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer}
#' @param new_mcmc number of new MCMC iterations to conduct and append
#' @param verb logical indicating whether to print iteration progress
#' @param re_approx logical indicating whether to re-randomize the ordering 
#'        and update Vecchia nearest-neighbor conditioning sets (only for fits 
#'        with \code{vecchia = TRUE})
#' @param ... N/A
#'        
#' @return object of the same class with the new iterations appended
#' 
#' @examples 
#' # See ?fit_two_layer for an example
#' 
#' @rdname continue
#' @export

continue <- function(object, new_mcmc, verb, re_approx, ...)
  UseMethod("continue", object)

# Continue One Layer ----------------------------------------------------------
#' @rdname continue
#' @export

continue.gp <- function(object, new_mcmc = 1000, verb = TRUE, ...) {
  
  tic <- proc.time()[[3]]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) true_g <- object$g[1] else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(theta = object$theta[object$nmcmc],
                  g = object$g[object$nmcmc],
                  tau2 = object$tau2[object$nmcmc])
  
  # Run continuing MCMC iterations
  samples <- gibbs_one_layer(x = object$x, 
                             y = object$y, 
                             nmcmc = new_mcmc + 1, 
                             verb = verb, 
                             initial = initial, 
                             true_g = true_g, 
                             settings = object$settings, 
                             v = object$v)

  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, samples$g[-1])
  object$theta <- c(object$theta, samples$theta[-1])
  object$tau2 <- c(object$tau2, samples$tau2[-1])
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# Continue Two Layer ----------------------------------------------------------
#' @rdname continue
#' @export

continue.dgp2 <- function(object, new_mcmc = 1000, verb = TRUE, ...) {
  
  tic <- proc.time()[[3]]

  monowarp <- (!is.null(object$x_grid))
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) true_g <- object$g[1] else true_g <- NULL
  
  if (monowarp) {
    initial <- list(w = monowarp_ref(object$x, object$x_grid, object$w_grid[[object$nmcmc]]),
                    theta_y = object$theta_y[object$nmcmc, ],
                    theta_w = object$theta_w[object$nmcmc, ],
                    g = object$g[object$nmcmc],
                    tau2 = object$tau2[object$nmcmc])
    samples <- gibbs_two_layer_mono(x = object$x,
                                    y = object$y,
                                    x_grid = object$x_grid,
                                    nmcmc = new_mcmc + 1,
                                    D = ncol(object$theta_y),
                                    verb = verb,
                                    initial = initial,
                                    true_g = true_g,
                                    settings = object$settings,
                                    v = object$v)
  } else {
    initial <- list(w = object$w[[object$nmcmc]],
                    theta_y = object$theta_y[object$nmcmc],
                    theta_w = object$theta_w[object$nmcmc, ],
                    g = object$g[object$nmcmc],
                    tau2 = object$tau2[object$nmcmc])
    # Run continuing MCMC iterations
    samples <- gibbs_two_layer(x = object$x, 
                               y = object$y, 
                               nmcmc = new_mcmc + 1, 
                               D = ncol(object$w[[1]]), 
                               verb = verb, 
                               initial = initial, 
                               true_g = true_g, 
                               settings = object$settings, 
                               v = object$v)
  }
  
  

  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, samples$g[-1])
  if (monowarp) {
    object$theta_y <- rbind(object$theta_y, samples$theta_y[-1, , drop = FALSE])
    object$w_grid <- c(object$w_grid, samples$w_grid[-1])
  } else {
    object$theta_y <- c(object$theta_y, samples$theta_y[-1])
    object$w <- c(object$w, samples$w[-1])
  }
  object$theta_w <- rbind(object$theta_w, samples$theta_w[-1, , drop = FALSE])
  object$tau2 <- c(object$tau2, samples$tau2[-1])
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# Continue Three Layer --------------------------------------------------------
#' @rdname continue
#' @export

continue.dgp3 <- function(object, new_mcmc = 1000, verb = TRUE, ...) {
  
  tic <- proc.time()[[3]]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) true_g <- object$g[1] else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(w = object$w[[object$nmcmc]],
                  z = object$z[[object$nmcmc]],
                  theta_y = object$theta_y[object$nmcmc],
                  theta_w = object$theta_w[object$nmcmc, ],
                  theta_z = object$theta_z[object$nmcmc, ],
                  g = object$g[object$nmcmc],
                  tau2 = object$tau2[object$nmcmc])
  
  # Run continuing MCMC iterations
  samples <- gibbs_three_layer(x = object$x, 
                               y = object$y, 
                               nmcmc = new_mcmc + 1, 
                               D = ncol(object$w[[1]]), 
                               verb = verb, 
                               initial = initial, 
                               true_g = true_g, 
                               settings = object$settings, 
                               v = object$v)
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, samples$g[-1])
  object$theta_y <- c(object$theta_y, samples$theta_y[-1])
  object$theta_w <- rbind(object$theta_w, samples$theta_w[-1, , drop = FALSE])
  object$theta_z <- rbind(object$theta_z, samples$theta_z[-1, , drop = FALSE])
  object$w <- c(object$w, samples$w[-1])
  object$z <- c(object$z, samples$z[-1])
  object$tau2 <- c(object$tau2, samples$tau2[-1])
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# Continue One Layer Vecchia --------------------------------------------------
#' @rdname continue
#' @export

continue.gpvec <- function(object, new_mcmc = 1000, verb = TRUE, 
                           re_approx = FALSE, ...) {
  
  tic <- proc.time()[[3]]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) true_g <- object$g[1] else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(theta = object$theta[object$nmcmc],
                  g = object$g[object$nmcmc],
                  tau2 = object$tau2[object$nmcmc])
  
  # Optionally redo approximation
  if (re_approx) {
    x_approx <- NULL
  } else x_approx <- object$x_approx
  
  # Run continuing MCMC iterations
  samples <- gibbs_one_layer_vec(x = object$x, 
                                 y = object$y, 
                                 nmcmc = new_mcmc + 1, 
                                 verb = verb, 
                                 initial = initial, 
                                 true_g = true_g, 
                                 settings = object$settings, 
                                 v = object$v, 
                                 m = object$m, 
                                 ordering = object$ordering,
                                 x_approx = x_approx)
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, samples$g[-1])
  object$theta <- c(object$theta, samples$theta[-1])
  object$tau2 <- c(object$tau2, samples$tau2[-1])
  object$x_approx <- samples$x_approx
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# Continue Two Layer Vecchia --------------------------------------------------
#' @rdname continue
#' @export

continue.dgp2vec <- function(object, new_mcmc = 1000, verb = TRUE, 
                             re_approx = FALSE, ...) {
  
  tic <- proc.time()[[3]]

  monowarp <- (!is.null(object$x_grid))
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) true_g <- object$g[1] else true_g <- NULL
  
  # Optionally redo approximation
  if (re_approx) {
    if (!monowarp) x_approx <- NULL
    w_approx <- NULL
  } else {
    if (!monowarp) x_approx <- object$x_approx
    w_approx <- object$w_approx
  }

  # Start MCMC at last samples
  if (monowarp) {
    inital <- list(w = monowarp_ref(object$x, object$xgrid, object$w_grid[[object$nmcmc]]),
                   theta_y = object$theta_y[object$nmcmc, ],
                   theta_w = object$theta_w[object$nmcmc, ],
                   g = object$g[object$nmcmc],
                   tau2 = object$tau2[object$nmcmc])

    samples <- gibbs_two_layer_vec_mono(x = object$x,
                                        y = object$y,
                                        x_grid = object$x_grid,
                                        nmcmc = new_mcmc + 1,
                                        D = ncol(object$theta_y),
                                        verb = verb,
                                        initial = initial,
                                        true_g = true_g,
                                        settings = object$settings,
                                        v = object$v,
                                        m = object$m,
                                        ordering = object$ordering,
                                        w_approx = w_approx)
  } else {
    initial <- list(w = object$w[[object$nmcmc]],
                    theta_y = object$theta_y[object$nmcmc],
                    theta_w = object$theta_w[object$nmcmc, ],
                    g = object$g[object$nmcmc],
                    tau2 = object$tau2[object$nmcmc])

    samples <- gibbs_two_layer_vec(x = object$x, 
                                   y = object$y, 
                                   nmcmc = new_mcmc + 1, 
                                   D = ncol(object$w[[1]]), 
                                   verb = verb, 
                                   initial = initial, 
                                   true_g = true_g, 
                                   settings = object$settings, 
                                   v = object$v, 
                                   m = object$m, 
                                   ordering = object$ordering,
                                   x_approx = x_approx, 
                                   w_approx = w_approx)
  }
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, samples$g[-1])
  if (monowarp) {
    object$theta_y <- rbind(object$theta_y, samples$theta_y[-1, , drop = FALSE])
    object$w_grid <- c(object$w_grid, samples$w_grid[-1])
  } else {
    object$theta_y <- c(object$theta_y, samples$theta_y[-1])
    object$w <- c(object$w, samples$w[-1])
  }
  object$theta_w <- rbind(object$theta_w, samples$theta_w[-1, , drop = FALSE])
  object$tau2 <- c(object$tau2, samples$tau2[-1])
  if (monowarp) object$x_approx <- samples$x_approx
  object$w_approx <- samples$w_approx
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# Continue Three Layer Vecchia ------------------------------------------------
#' @rdname continue
#' @export

continue.dgp3vec <- function(object, new_mcmc = 1000, verb = TRUE, 
                             re_approx = FALSE, ...) {
  
  tic <- proc.time()[[3]]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) true_g <- object$g[1] else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(w = object$w[[object$nmcmc]],
                  z = object$z[[object$nmcmc]],
                  theta_y = object$theta_y[object$nmcmc],
                  theta_w = object$theta_w[object$nmcmc, ],
                  theta_z = object$theta_z[object$nmcmc, ],
                  g = object$g[object$nmcmc],
                  tau2 = object$tau2[object$nmcmc])
  
  # Optionally redo approximation
  if (re_approx) {
    x_approx <- NULL
    z_approx <- NULL
    w_approx <- NULL
  } else {
    x_approx <- object$x_approx
    z_approx <- object$z_approx
    w_approx <- object$w_approx
  }
  
  # Run continuing MCMC iterations
  samples <- gibbs_three_layer_vec(x = object$x, 
                                   y = object$y, 
                                   nmcmc = new_mcmc + 1, 
                                   D = ncol(object$w[[1]]),
                                   verb = verb, 
                                   initial = initial, 
                                   true_g = true_g, 
                                   settings = object$settings, 
                                   v = object$v, 
                                   m = object$m, 
                                   ordering = object$ordering,
                                   x_approx = x_approx, 
                                   z_approx = z_approx, 
                                   w_approx = w_approx)
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, samples$g[-1])
  object$theta_y <- c(object$theta_y, samples$theta_y[-1])
  object$theta_w <- rbind(object$theta_w, samples$theta_w[-1, , drop = FALSE])
  object$theta_z <- rbind(object$theta_z, samples$theta_z[-1, , drop = FALSE])
  object$w <- c(object$w, samples$w[-1])
  object$z <- c(object$z, samples$z[-1])
  object$tau2 <- c(object$tau2, samples$tau2[-1])
  object$x_approx <- samples$x_approx
  object$z_approx <- samples$z_approx
  object$w_approx <- samples$w_approx
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}