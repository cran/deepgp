
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   continue.gp
#   continue.dgp2
#   continue.dgp3
#   continue.gpvec
#   continue.dgp2vec
#   continue.dgp3vec

# continue S3 class -----------------------------------------------------------
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
#' @param ord optional ordering to be used in Vecchia re-approximation (only
#'        for fits with \code{vecchia = TRUE} when \code{re_approx = TRUE}')
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

# continue.gp -----------------------------------------------------------------
#' @rdname continue
#' @export

continue.gp <- function(object, new_mcmc = 1000, verb = TRUE, ...) {
  
  tic <- proc.time()[[3]]
  
  # Use true nugget if it was specified
  if (length(object$g) == 1) true_g <- object$g else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(theta = ifel(object$settings$sep, object$theta[object$nmcmc, ], 
                               object$theta[object$nmcmc]),
                  g = ifel(is.null(true_g), object$g[object$nmcmc], NULL))
  
  # Run continuing MCMC iterations
  samples <- gibbs_one_layer(x = object$x,
                             y = object$y,
                             dydx = object$dydx,
                             nmcmc = new_mcmc + 1,
                             verb = verb,
                             initial = initial,
                             true_g = true_g,
                             settings = object$settings,
                             v = object$v)

  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$tau2 <- c(object$tau2, samples$tau2[-1])
  if (object$settings$sep) {
    object$theta <- rbind(object$theta, samples$theta[-1, , drop = FALSE])
  } else object$theta <- c(object$theta, samples$theta[-1])
  if (is.null(true_g)) object$g <- c(object$g, samples$g[-1])
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# continue.dgp2 ---------------------------------------------------------------
#' @rdname continue
#' @export

continue.dgp2 <- function(object, new_mcmc = 1000, verb = TRUE, ...) {
  
  tic <- proc.time()[[3]]
  grad_enhance <- (!is.null(object$dydx))

  # Use true nugget if it was specified
  if (length(object$g) == 1) true_g <- object$g else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(theta_y = object$theta_y[object$nmcmc],
                  theta_w = object$theta_w[object$nmcmc, ],
                  g = ifel(is.null(true_g), object$g[object$nmcmc], NULL))
  if (object$settings$monowarp) {
    initial$w <- NULL
  } else {
    initial$w <- object$w[object$nmcmc, , ]
  }
  if (grad_enhance) initial$dwdx <- object$dwdx


  if (object$settings$monowarp) {
    samples <- gibbs_two_layer_mono(x = object$x,
                                    y = object$y,
                                    x_grid = object$x_grid,
                                    nmcmc = new_mcmc + 1,
                                    verb = verb,
                                    initial = initial,
                                    true_g = true_g,
                                    settings = object$settings,
                                    v = object$v)
  } else if (grad_enhance) {
    samples <- gibbs_two_layer_grad(x = object$x,
                                    y = object$y,
                                    dydx = object$dydx,
                                    nmcmc = new_mcmc + 1,
                                    verb = verb,
                                    initial = initial,
                                    true_g = true_g,
                                    settings = object$settings,
                                    v = object$v)
  } else {
    samples <- gibbs_two_layer(x = object$x, 
                               y = object$y, 
                               nmcmc = new_mcmc + 1, 
                               D = dim(object$w)[3], 
                               verb = verb, 
                               initial = initial, 
                               true_g = true_g, 
                               settings = object$settings, 
                               v = object$v)
  }
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$tau2_y <- c(object$tau2_y, samples$tau2_y[-1])
  object$theta_y <- c(object$theta_y, samples$theta_y[-1])
  if (object$settings$monowarp) {
    if (!is.null(samples$tau2_w)) {
      object$tau2_w <- rbind(object$tau2_w, samples$tau2_w[-1, , drop = FALSE])
    }
    object$w_grid <- abind::abind(object$w_grid, samples$w_grid[-1, , , drop = FALSE], along = 1)
  }
  object$theta_w <- rbind(object$theta_w, samples$theta_w[-1, , drop = FALSE])
  object$w <- abind::abind(object$w, samples$w[-1, , , drop = FALSE], along = 1)
  if (is.null(true_g)) object$g <- c(object$g, samples$g[-1])
  object$ll <- c(object$ll, samples$ll[-1])

  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# continue.dgp3 ---------------------------------------------------------------
#' @rdname continue
#' @export

continue.dgp3 <- function(object, new_mcmc = 1000, verb = TRUE, ...) {
  
  tic <- proc.time()[[3]]
  
  # Use true nugget if it was specified
  if (length(object$g) == 1) true_g <- object$g else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(w = object$w[object$nmcmc, , ],
                  z = object$z[object$nmcmc, , ],
                  theta_y = object$theta_y[object$nmcmc],
                  theta_w = object$theta_w[object$nmcmc, ],
                  theta_z = object$theta_z[object$nmcmc, ],
                  g = ifel(is.null(true_g), object$g[object$nmcmc], NULL))
  
  # Run continuing MCMC iterations
  samples <- gibbs_three_layer(x = object$x, 
                               y = object$y, 
                               nmcmc = new_mcmc + 1, 
                               D = dim(object$w)[3], 
                               verb = verb, 
                               initial = initial, 
                               true_g = true_g, 
                               settings = object$settings, 
                               v = object$v)
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$tau2_y <- c(object$tau2_y, samples$tau2_y[-1])
  object$theta_y <- c(object$theta_y, samples$theta_y[-1])
  object$theta_w <- rbind(object$theta_w, samples$theta_w[-1, , drop = FALSE])
  object$theta_z <- rbind(object$theta_z, samples$theta_z[-1, , drop = FALSE])
  object$w <- abind::abind(object$w, samples$w[-1, , , drop = FALSE], along = 1)
  object$z <- abind::abind(object$z, samples$z[-1, , , drop = FALSE], along = 1)
  if (is.null(true_g)) object$g <- c(object$g, samples$g[-1])
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# continue.gpvec --------------------------------------------------------------
#' @rdname continue
#' @export

continue.gpvec <- function(object, new_mcmc = 1000, verb = TRUE, 
                           re_approx = FALSE, ord = NULL, ...) {
  
  tic <- proc.time()[[3]]

  # Use true nugget if it was specified
  if (length(object$g) == 1) true_g <- object$g else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(theta = ifel(object$settings$sep, object$theta[object$nmcmc, ], 
                               object$theta[object$nmcmc]),
                  g = ifel(is.null(true_g), object$g[object$nmcmc], NULL))
  
  # Run continuing MCMC iterations
  samples <- gibbs_one_layer_vec(x = object$x, 
                                 y = object$y, 
                                 dydx = object$dydx,
                                 nmcmc = new_mcmc + 1, 
                                 verb = verb, 
                                 initial = initial, 
                                 true_g = true_g, 
                                 settings = object$settings, 
                                 v = object$v, 
                                 m = object$x_approx$m, 
                                 ord = ifel(re_approx, ord, NULL), # only used if re_approx = TRUE
                                 cores = object$x_approx$cores,
                                 x_approx = ifel(re_approx, NULL, object$x_approx))
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$tau2 <- c(object$tau2, samples$tau2[-1])
  if (object$settings$sep) {
    object$theta <- rbind(object$theta, samples$theta[-1, , drop = FALSE])
  } else object$theta <- c(object$theta, samples$theta[-1])
  if (is.null(true_g)) object$g <- c(object$g, samples$g[-1])
  object$x_approx <- samples$x_approx
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# continue.dgp2vec ------------------------------------------------------------
#' @rdname continue
#' @export

continue.dgp2vec <- function(object, new_mcmc = 1000, verb = TRUE, 
                             re_approx = FALSE, ord = NULL, ...) {
  
  tic <- proc.time()[[3]]
  grad_enhance <- (!is.null(object$dydx))

  # Use true nugget if it was specified
  if (length(object$g) == 1) true_g <- object$g else true_g <- NULL

  # Start MCMC at last samples
  initial <- list(theta_y = object$theta_y[object$nmcmc],
                  theta_w = object$theta_w[object$nmcmc, ],
                  g = ifel(is.null(true_g), object$g[object$nmcmc], NULL))
  if (object$settings$monowarp) {
    initial$w <- NULL
  } else {
    initial$w <- object$w[object$nmcmc, , ]
  }
  if (grad_enhance) initial$dwdx <- object$dwdx
  
  if (object$settings$monowarp) {
    samples <- gibbs_two_layer_vec_mono(x = object$x,
                                        y = object$y,
                                        x_grid = object$x_grid,
                                        nmcmc = new_mcmc + 1,
                                        verb = verb,
                                        initial = initial,
                                        true_g = true_g,
                                        settings = object$settings,
                                        v = object$v,
                                        m = object$w_approx$m,
                                        ord = ifel(re_approx, ord, NULL), # only used if re_approx = TRUE
                                        cores = object$w_approx$cores,
                                        w_approx = ifel(re_approx, NULL, object$w_approx))
  } else if (grad_enhance) {
    samples <- gibbs_two_layer_vec_grad(x = object$x,
                                        y = object$y,
                                        dydx = object$dydx,
                                        nmcmc = new_mcmc + 1,
                                        verb = verb,
                                        initial = initial,
                                        true_g = true_g,
                                        settings = object$settings,
                                        v = object$v,
                                        m = object$w_approx$m,
                                        ord = ifel(re_approx, ord, NULL), # only used if re_approx = TRUE
                                        cores = object$w_approx$cores,
                                        x_approx = ifel(re_approx, NULL, object$x_approx),
                                        w_approx = ifel(re_approx, NULL, object$w_approx))
  } else {
    samples <- gibbs_two_layer_vec(x = object$x, 
                                   y = object$y, 
                                   nmcmc = new_mcmc + 1, 
                                   D = dim(object$w)[3], 
                                   verb = verb, 
                                   initial = initial, 
                                   true_g = true_g, 
                                   settings = object$settings, 
                                   v = object$v, 
                                   m = object$x_approx$m, 
                                   ord = ifel(re_approx, ord, NULL), # only used if re_approx = TRUE
                                   cores = object$x_approx$cores,
                                   x_approx = ifel(re_approx, NULL, object$x_approx), 
                                   w_approx = ifel(re_approx, NULL, object$w_approx))
  }
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$tau2_y <- c(object$tau2_y, samples$tau2[-1])
  object$theta_y <- c(object$theta_y, samples$theta_y[-1])
  object$theta_w <- rbind(object$theta_w, samples$theta_w[-1, , drop = FALSE])
  if (object$settings$monowarp) {
    if (!is.null(samples$tau2_w)) {
      object$tau2_w <- rbind(object$tau2_w, samples$tau2_w[-1, , drop = FALSE])
    }
    object$w_grid <- abind::abind(object$w_grid, samples$w_grid[-1, , , drop = FALSE], along = 1)
  }
  object$w <- abind::abind(object$w, samples$w[-1, , , drop = FALSE], along = 1)
  if (is.null(true_g)) object$g <- c(object$g, samples$g[-1])
  if (!object$settings$monowarp) object$x_approx <- samples$x_approx
  object$w_approx <- samples$w_approx
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}

# continue.dgp3vec ------------------------------------------------------------
#' @rdname continue
#' @export

continue.dgp3vec <- function(object, new_mcmc = 1000, verb = TRUE, 
                             re_approx = FALSE, ord = NULL, ...) {
  
  tic <- proc.time()[[3]]
  
  # Use true nugget if it was specified
  if (length(object$g) == 1) true_g <- object$g else true_g <- NULL
  
  # Start MCMC at last samples
  initial <- list(w = object$w[object$nmcmc, , ],
                  z = object$z[object$nmcmc, , ],
                  theta_y = object$theta_y[object$nmcmc],
                  theta_w = object$theta_w[object$nmcmc, ],
                  theta_z = object$theta_z[object$nmcmc, ],
                  g = ifel(is.null(true_g), object$g[object$nmcmc], NULL))
  
  # Run continuing MCMC iterations
  samples <- gibbs_three_layer_vec(x = object$x, 
                                   y = object$y, 
                                   nmcmc = new_mcmc + 1, 
                                   D = dim(object$w)[3],
                                   verb = verb, 
                                   initial = initial, 
                                   true_g = true_g, 
                                   settings = object$settings, 
                                   v = object$v, 
                                   m = object$x_approx$m, 
                                   ord = ifel(re_approx, ord, NULL), # only used if re_approx = TRUE
                                   cores = object$x_approx$cores,
                                   x_approx = ifel(re_approx, NULL, object$x_approx), 
                                   z_approx = ifel(re_approx, NULL, object$z_approx), 
                                   w_approx = ifel(re_approx, NULL, object$w_approx))
  
  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$tau2_y <- c(object$tau2_y, samples$tau2_y[-1])
  object$theta_y <- c(object$theta_y, samples$theta_y[-1])
  object$theta_w <- rbind(object$theta_w, samples$theta_w[-1, , drop = FALSE])
  object$theta_z <- rbind(object$theta_z, samples$theta_z[-1, , drop = FALSE])
  if (is.null(true_g)) object$g <- c(object$g, samples$g[-1])
  object$w <- abind::abind(object$w, samples$w[-1, , , drop = FALSE], along = 1)
  object$z <- abind::abind(object$z, samples$z[-1, , , drop = FALSE], along = 1)
  object$x_approx <- samples$x_approx
  object$z_approx <- samples$z_approx
  object$w_approx <- samples$w_approx
  object$ll <- c(object$ll, samples$ll[-1])
  
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  return(object)
}