
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   plot.gp == plot.gpvec
#   plot.dgp2 == plot.dgp2vec
#   plot.dgp3 == plot.dgp3vec

# Define Plot for S3 Objects --------------------------------------------------
#' @name plot
#' @title Plots object from \code{deepgp} package
#' 
#' @description Acts on a \code{gp}, \code{gpvec}, \code{dgp2}, \code{dgp2vec},
#'     \code{dgp3}, or \code{dgp3vec} object.  
#'     Generates trace plots for length scale and nugget hyperparameters.
#'     Generates plots of hidden layers for one-dimensional inputs.  Generates
#'     plots of the posterior mean and estimated 90\% prediction intervals for 
#'     one-dimensional inputs; generates heat maps of the posterior mean and 
#'     point-wise variance for two-dimensional inputs.
#'     
#' @details Trace plots are useful in assessing burn-in.  Hidden layer plots 
#'     are colored on a gradient - red lines represent earlier iterations and 
#'     yellow lines represent later iterations - to help assess burn-in of the 
#'     hidden layers.  These plots are meant to help in model fitting and 
#'     visualization.
#' 
#' @param x object of class \code{gp}, \code{gpvec}, \code{dgp2}, 
#'        \code{dgp2vec}, \code{dgp3}, or \code{dgp3vec}
#' @param trace logical indicating whether to generate trace plots (default is
#'        TRUE if the object has not been through \code{predict})
#' @param hidden logical indicating whether to generate plots of hidden layers
#'        (two or three layer only, default is FALSE)
#' @param predict logical indicating whether to generate posterior predictive 
#'        plot (default is TRUE if the object has been through \code{predict})
#' @param ... N/A
#' 
#' @examples 
#' # See "fit_one_layer", "fit_two_layer", or "fit_three_layer"
#' # for an example
#' 
#' @rdname plot
NULL

# Plot One Layer --------------------------------------------------------------
#' @rdname plot
#' @export

plot.gp <- function(x, trace = NULL, predict = NULL, ...) {
  
  # save and restore par settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  Dx <- ncol(x$x)
  if (is.null(x$mean)) {
    if (is.null(trace)) trace <- TRUE
    if (is.null(predict)) predict <- FALSE
  } else {
    if (is.null(trace)) trace <- FALSE
    if (is.null(predict)) predict <- TRUE
  }
  
  if (trace) {
    par(mfrow = c(1, 2), mar = c(5, 4, 2, 2))
    plot(x$g, type = 'l', ylab = 'g', xlab = 'Iteration',
         main = 'Trace Plot of g')
    plot(x$theta, type = 'l', ylab = 'theta_y', xlab = 'Iteration',
         main = 'Trace Plot of theta')
  }
  
  if (predict) {
    if (Dx == 1) {
      par(mfrow = c(1, 1))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, x$Sigma_smooth))
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = 'l', xlab = 'X', ylab = 'Y', 
           ylim = c(min(q1), max(q3)),
           col = 'blue', ...)
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = 'grey', lty = 1)
        lines(x$x_new[o], x$mean[o], col = 'blue')
      }
      lines(x$x_new[o], q1[o], col = 'blue', lty = 2)
      lines(x$x_new[o], q3[o], col = 'blue', lty = 2)
      points(x$x, x$y, pch = 20)
    } else if (Dx == 2) {
      if (!requireNamespace("interp", quietly = TRUE)) {
        stop("Package \"interp\" needed for this plot. Please install it.",
             call. = FALSE)
      }
      cols <- heat.colors(128)
      i1 <- interp::interp(x$x_new[, 1], x$x_new[, 2], x$mean)
      if (is.null(x$Sigma)) {
        i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(x$s2))
      } else i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(diag(x$Sigma)))
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
      image(i1, col = cols, main = 'Posterior Mean', xlab = 'X1', ylab = 'X2')
      contour(i1, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = 'Posterior Variance', xlab = 'X1', 
            ylab = 'X2')
      contour(i2, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat('Dimension of X too large for default plotting')
  }
}

# Plot One Layer Vecchia ------------------------------------------------------
#' @rdname plot
#' @export

plot.gpvec <- plot.gp

# Plot Two Layer --------------------------------------------------------------
#' @rdname plot
#' @export

plot.dgp2 <- function(x, trace = NULL, hidden = NULL, predict = NULL, ...) {
  
  # save and restore par settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  Dx <- ncol(x$x)
  D <- ncol(x$w[[1]])
  if (is.null(x$mean)) {
    if (is.null(trace)) trace <- TRUE
    if (is.null(predict)) predict <- FALSE
  } else {
    if (is.null(trace)) trace <- FALSE
    if (is.null(predict)) predict <- TRUE
  }
  if (is.null(hidden)) hidden <- FALSE
  
  if (trace) {
    par(mfrow = c(1, D + 2), mar = c(5, 4, 2, 2))
    plot(x$g, type = 'l', ylab = 'g', xlab = 'Iteration',
         main = 'Trace Plot of g')
    plot(x$theta_y, type = 'l', ylab = 'theta_y', xlab = 'Iteration',
         main = 'Trace Plot of theta_y')
    for (i in 1:D)
      plot(x$theta_w[, i], type = 'l', ylab = 'theta_w', xlab = 'Iteration',
           main = paste0('Trace Plot of theta_w [', i, ']'))
  }
  
  if (hidden) {
    if (Dx == 1 & D == 1) {
      # specify the hidden layers to plot
      indx <- floor(seq(from = 1, to = x$nmcmc, length = 100))
      if (indx[1] == 0) indx[1] = 1
      col <- heat.colors(100 + 10) # add ten to avoid using colors that are too light
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
      
      # plot x to w
      o <- order(x$x)
      plot(x$x[o], x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), type = 'l', xlab = 'X', 
           ylab = 'W', col = col[1], main = paste0('MCMC samples of X to W'), 
           ylim = c(min(unlist(x$w[indx])), max(unlist(x$w[indx]))))
      for (j in 2:length(indx)) {
        lines(x$x[o], x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), col = col[j])
      }
    } else cat('Default plotting not prepared for these dimensions')
  }
  
  if (predict) {
    if (Dx == 1){
      par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, x$Sigma_smooth))
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = 'l', xlab = 'X', ylab = 'Y', 
           ylim = c(min(q1), max(q3)),
           col = 'blue', ...)
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = 'grey', lty = 1)
        lines(x$x_new[o], x$mean[o], col = 'blue')
      }
      lines(x$x_new[o], q1[o], col = 'blue', lty = 2)
      lines(x$x_new[o], q3[o], col = 'blue', lty = 2)
      points(x$x, x$y, pch = 20)
      
    } else if (Dx == 2) {
      if (!requireNamespace("interp", quietly = TRUE)) {
        stop("Package \"interp\" needed for this plot. Please install it.",
             call. = FALSE)
      }
      cols <- heat.colors(128)
      i1 <- interp::interp(x$x_new[, 1], x$x_new[, 2], x$mean)
      if (is.null(x$Sigma)) {
        i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(x$s2))
      } else i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(diag(x$Sigma)))
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
      image(i1, col = cols, main = 'Posterior Mean', xlab = 'X1', ylab = 'X2')
      contour(i1, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = 'Posterior Variance', xlab = 'X1', ylab = 'X2')
      contour(i2, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat('Dimension of X too large for default plotting')
  }
}

# Plot Two Layer Vecchia ------------------------------------------------------
#' @rdname plot
#' @export

plot.dgp2vec <- plot.dgp2

# Plot Three Layer ------------------------------------------------------------
#' @rdname plot
#' @export

plot.dgp3 <- function(x, trace = NULL, hidden = NULL, predict = NULL, ...) {
  
  # save and restore par settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  Dx <- ncol(x$x)
  D <- ncol(x$w[[1]])
  if (is.null(x$mean)) {
    if (is.null(trace)) trace <- TRUE
    if (is.null(predict)) predict <- FALSE
  } else {
    if (is.null(trace)) trace <- FALSE
    if (is.null(predict)) predict <- TRUE
  }
  if (is.null(hidden)) hidden <- FALSE
  
  if (trace) {
    par(mfrow = c(2, D + 1), mar = c(5, 4, 2, 2))
    plot(x$g, type = 'l', ylab = 'g', xlab = 'Iteration',
         main = 'Trace Plot of g')
    plot(x$theta_y, type = 'l', ylab = 'theta_y', xlab = 'Iteration',
         main = 'Trace Plot of theta_y')
    for (i in 1:D)
      plot(x$theta_w[,i], type = 'l', ylab = 'theta_w', xlab = 'Iteration',
           main = paste0('Trace Plot of theta_w [', i, ']'))
    for (i in 1:D)
      plot(x$theta_z[,i], type = 'l', ylab = 'theta_z', xlab = 'Iteration',
           main = paste0('Trace Plot of theta_z [', i, ']'))
  }
  
  if (hidden) {
    if (Dx == 1 & D == 1) {
      # specify the hidden layers to plot
      indx <- floor(seq(from = 1, to = x$nmcmc, length = 100))
      if (indx[1] == 0) indx[1] <- 1
      col <- heat.colors(100 + 10)
      par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
      
      # plot z to w
      o <- order(x$z[[indx[1]]])
      plot(x$z[[indx[1]]][o] - mean(x$z[[indx[1]]]), 
           x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), type = 'l', xlab = 'Z', 
           ylab = 'W', col = col[1], main = paste0('MCMC samples of Z to W'), 
           xlim = c(min(unlist(x$z[indx])), max(unlist(x$z[indx]))),
           ylim = c(min(unlist(x$w[indx])), max(unlist(x$w[indx]))))
      for (j in 2:length(indx)) {
        o <- order(x$z[[indx[j]]])
        lines(x$z[[indx[j]]][o] - mean(x$z[[indx[j]]]), 
              x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), col = col[j])
      }
      
      # plot x to z
      o <- order(x$x)
      plot(x$x[o], x$z[[indx[1]]][o] - mean(x$z[[indx[1]]]), type = 'l', 
           xlab = 'X', ylab = 'Z', col = col[1], 
           main = paste0('MCMC samples of X to Z'), 
           ylim = c(min(unlist(x$z[indx])), max(unlist(x$z[indx]))))
      for (j in 2:length(indx)) {
        lines(x$x[o], x$z[[indx[j]]][o] - mean(x$z[[indx[j]]]), col = col[j])
      }
    } else cat('Default plotting not prepared for these dimensions')
  }
  
  if (predict) {
    if (Dx == 1) {
      par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, x$Sigma_smooth))
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = 'l', xlab = 'X', ylab = 'Y', 
           ylim = c(min(q1), max(q3)),
           col = 'blue', ...)
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = 'grey', lty = 1)
        lines(x$x_new[o], x$mean[o], col = 'blue')
      }
      lines(x$x_new[o], q1[o], col = 'blue', lty = 2)
      lines(x$x_new[o], q3[o], col = 'blue', lty = 2)
      points(x$x, x$y, pch = 20)
    } else if (Dx == 2) {
      if (!requireNamespace("interp", quietly = TRUE)) {
        stop("Package \"interp\" needed for this function to work. Please install it.",
             call. = FALSE)
      }
      cols <- heat.colors(128)
      i1 <- interp::interp(x$x_new[, 1], x$x_new[, 2], x$mean)
      if (is.null(x$Sigma)) {
        i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(x$s2))
      } else i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(diag(x$Sigma)))
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
      image(i1, col = cols, main = 'Posterior Mean', xlab = 'X1', ylab = 'X2')
      contour(i1, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = 'Posterior Variance', xlab = 'X1', ylab = 'X2')
      contour(i2, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat('Dimension of X too large for default plotting')
  }
}

# Plot Three Layer Vecchia ----------------------------------------------------
#' @rdname plot
#' @export

plot.dgp3vec <- plot.dgp3
