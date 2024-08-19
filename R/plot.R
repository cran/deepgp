
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
#'     Generates trace plots for outer log likelihood, length scale,
#'     and nugget hyperparameters.
#'     Generates plots of hidden layers for one-dimensional inputs or monotonic
#'     warpings.  Generates
#'     plots of the posterior mean and estimated 90\% prediction intervals for 
#'     one-dimensional inputs; generates heat maps of the posterior mean and 
#'     point-wise variance for two-dimensional inputs.
#'     
#' @details Trace plots are useful in assessing burn-in.  If there are too
#'     many hyperparameters to plot them all, then it is most useful to 
#'     visualize the log likelihood (e.g., \code{plot(fit$ll, type = "l")}).
#'
#'     Hidden layer plots are colored on a gradient - red lines represent 
#'     earlier iterations and yellow lines represent later iterations - to 
#'     help assess burn-in of the hidden layers.  Only every 100th sample
#'     is plotted.
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
#' # See ?fit_one_layer, ?fit_two_layer, or ?fit_three_layer
#' # for examples
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
    fixed_g <- (length(unique(x$g)) == 1)
    if (is.matrix(x$theta)) { # separable lengthscale
      if (fixed_g) {
        par(mfrow = c(1, ncol(x$theta) + 1), mar = c(5, 4, 2, 2))
        plot(x$ll, type = "l", ylab = "logl", xlab = "Iteration",
             main = "logl")
        for (i in 1:ncol(x$theta))
          plot(x$theta[, i], type = "l", ylab = paste0("theta[", i, "]"), 
               xlab = "Iteration", main = paste0("theta[", i, "]"))
      } else {
        par(mfrow = c(1, ncol(x$theta) + 2), mar = c(5, 4, 2, 2))
        plot(x$ll, type = "l", ylab = "logl", xlab = "Iteration",
             main = "logl")
        plot(x$g, type = "l", ylab = "g", xlab = "Iteration",
             main = "g")
        for (i in 1:ncol(x$theta))
          plot(x$theta[, i], type = "l", ylab = paste0("theta[", i, "]"), 
               xlab = "Iteration", main = paste0("theta[", i, "]"))
      }
    } else { # isotropic lengthscale
      if (fixed_g) {
        par(mfrow = c(1, 2), mar = c(5, 4, 2, 2))
        plot(x$ll, type = "l", ylab = "logl", xlab = "Iteration",
             main = "logl")
        plot(x$theta, type = "l", ylab = "theta", xlab = "Iteration",
             main = "theta")
      } else {
        par(mfrow = c(1, 3), mar = c(5, 4, 2, 2))
        plot(x$ll, type = "l", ylab = "logl", xlab = "Iteration",
             main = "logl")
        plot(x$g, type = "l", ylab = "g", xlab = "Iteration",
             main = "g")
        plot(x$theta, type = "l", ylab = "theta", xlab = "Iteration",
             main = "theta")
      }
    }
  }
  
  if (predict) {
    if (Dx == 1) {
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        Sigma_smooth <- x$Sigma - diag(mean(x$g * x$tau2), nrow(x$x_new))
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, Sigma_smooth))
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = "l", xlab = "x", ylab = "y", 
           ylim = c(min(q1), max(q3)),
           col = "blue", ...)
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = "grey", lty = 1)
        lines(x$x_new[o], x$mean[o], col = "blue")
      }
      lines(x$x_new[o], q1[o], col = "blue", lty = 2)
      lines(x$x_new[o], q3[o], col = "blue", lty = 2)
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
      image(i1, col = cols, main = "Posterior Mean", xlab = "x1", ylab = "x2")
      contour(i1, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = "Posterior Std. Dev.", xlab = "x1", 
            ylab = "x2")
      contour(i2, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat("Dimension of x too large for default plotting")
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

  monowarp <- (!is.null(x$x_grid))
  
  Dx <- ncol(x$x)
  if (monowarp) D <- ncol(x$w_grid[[1]]) else D <- ncol(x$w[[1]])
  if (is.null(x$mean)) {
    if (is.null(trace)) trace <- TRUE
    if (is.null(predict)) predict <- FALSE
  } else {
    if (is.null(trace)) trace <- FALSE
    if (is.null(predict)) predict <- TRUE
  }
  if (is.null(hidden)) hidden <- FALSE
  
  if (trace) {
    fixed_g <- (length(unique(x$g)) == 1)
    if (fixed_g) {
      if (monowarp) {
        par(mfrow = c(1, 2*D + 1), mar = c(5, 4, 2, 2))
      } else par(mfrow = c(1, D + 2), mar = c(5, 4, 2, 2))
    } else {
      if (monowarp) {
        par(mfrow = c(1, 2*D + 2), mar = c(5, 4, 2, 2))
      } else par(mfrow = c(1, D + 3), mar = c(5, 4, 2, 2))
    }
    plot(x$ll, type = "l", ylab = "outer logl", xlab = "Iteration",
          main = "outer logl")
    if (!fixed_g) 
      plot(x$g, type = "l", ylab = "g", xlab = "Iteration",
           main = "g")
    if (monowarp) {
      for (i in 1:D)
        plot(x$theta_y[, i], type = "l", ylab = paste0("theta_y[", i, "]"), 
             xlab = "Iteration", main = paste0("theta_y[", i, "]"))
    } else {
      plot(x$theta_y, type = "l", ylab = "theta_y", xlab = "Iteration",
             main = "theta_y")
    }
    for (i in 1:D)
      plot(x$theta_w[, i], type = "l", ylab = paste0("theta_w[", i, "]"), 
           xlab = "Iteration", main = paste0("theta_w[", i, "]"))
  }
  
  if (hidden) {
    # specify the hidden layers to plot
    indx <- floor(seq(from = 1, to = x$nmcmc, length = 100))
    if (indx[1] == 0) indx[1] = 1
    col <- heat.colors(100 + 10) # add ten to avoid using colors that are too light
    if (monowarp) { 
      grid_index <- fo_approx_init(x$x_grid, x$x)
      par(mfrow = c(1, D), mar = c(4, 4, 2, 2))
      for (i in 1:D) {
        o <- order(x$x[, i])
        plot(x$x[o, i], monowarp_ref(x$x[, i], x$x_grid[, i], x$w_grid[[indx[1]]][, i], 
                                     grid_index[, i])[o], 
             type = "l", xlab = paste0("x", i), ylab = paste0("w", i), col = col[1], 
             main = "monowarped ESS samples", ylim = c(0, 1))
        for (j in 2:length(indx)) 
          lines(x$x[o, i], monowarp_ref(x$x[, i], x$x_grid[, i], x$w_grid[[indx[j]]][, i], 
                                        grid_index[, i])[o], col = col[j])
        abline(0, 1, lwd = 2, lty = 2)
      }
    } else if (Dx == 1 & D == 1) {
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
      o <- order(x$x)
      plot(x$x[o], x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), type = "l", xlab = "x", 
           ylab = "w", col = col[1], main = "ESS samples", 
           ylim = c(min(unlist(x$w[indx])), max(unlist(x$w[indx]))))
      for (j in 2:length(indx)) {
        lines(x$x[o], x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), col = col[j])
      }
    } else cat("Default plotting not prepared for these dimensions")
  } 

  if (predict) {
    if (Dx == 1){
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        Sigma_smooth <- x$Sigma - diag(mean(x$g * x$tau2), nrow(x$x_new))
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, Sigma_smooth))
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = "l", xlab = "x", ylab = "y", 
           ylim = c(min(q1), max(q3)),
           col = "blue", ...)
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = "grey", lty = 1)
        lines(x$x_new[o], x$mean[o], col = "blue")
      }
      lines(x$x_new[o], q1[o], col = "blue", lty = 2)
      lines(x$x_new[o], q3[o], col = "blue", lty = 2)
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
      image(i1, col = cols, main = "Posterior Mean", xlab = "x1", ylab = "x2")
      contour(i1, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = "Posterior Std. Dev.", xlab = "x1", ylab = "x2")
      contour(i2, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat("Dimension of x too large for default plotting")
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
    fixed_g <- (length(unique(x$g)) == 1)
    if (fixed_g) {
      par(mfrow = c(2, D + 1), mar = c(5, 4, 2, 2))
      plot(x$ll, type = "l", ylab = "outer logl", xlab = "Iteration",
           main = "outer logl")
      plot(x$theta_y, type = "l", ylab = "theta_y", xlab = "Iteration",
           main = "theta_y")
      for (i in 1:D)
        plot(x$theta_w[,i], type = "l", ylab = paste0("theta_w[", i, "]"), 
             xlab = "Iteration", main = paste0("theta_w[", i, "]"))
      for (i in 1:D)
        plot(x$theta_z[,i], type = "l", ylab = paste0("theta_z[", i, "]"), 
             xlab = "Iteration", main = paste0("theta_z[", i, "]"))
    } else {
      par(mfrow = c(2, D + 2), mar = c(5, 4, 2, 2))
      plot(x$ll, type = "l", ylab = "outer logl", xlab = "Iteration",
           main = "outer logl")
      plot(x$g, type = "l", ylab = "g", xlab = "Iteration",
           main = "g")
      plot(x$theta_y, type = "l", ylab = "theta_y", xlab = "Iteration",
           main = "theta_y")
      for (i in 1:D)
        plot(x$theta_w[,i], type = "l", ylab = paste0("theta_w[", i, "]"), 
             xlab = "Iteration", main = paste0("theta_w[", i, "]"))
      for (i in 1:D)
        plot(x$theta_z[,i], type = "l", ylab = paste0("theta_z[", i, "]"), 
             xlab = "Iteration", main = paste0("theta_z[", i, "]"))
    }
  }
  
  if (hidden) {
    if (Dx == 1 & D == 1) {
      # specify the hidden layers to plot
      indx <- floor(seq(from = 1, to = x$nmcmc, length = 100))
      if (indx[1] == 0) indx[1] <- 1
      col <- heat.colors(100 + 10)
      par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
      # plot x to z
      o <- order(x$x)
      plot(x$x[o], x$z[[indx[1]]][o] - mean(x$z[[indx[1]]]), type = "l", 
           xlab = "x", ylab = "z", col = col[1], 
           main = paste0("ESS samples (inner layer)"), 
           ylim = c(min(unlist(x$z[indx])), max(unlist(x$z[indx]))))
      for (j in 2:length(indx)) {
        lines(x$x[o], x$z[[indx[j]]][o] - mean(x$z[[indx[j]]]), col = col[j])
      }
      # plot z to w
      zmin <- 100; zmax <- -100 # arbitrary extreme values
      wmin <- 100; wmax <- -100
      for (i in 1:length(indx)) {
        zmin <- min(zmin, x$z[[indx[i]]] - mean(x$z[[indx[i]]]))
        zmax <- max(zmax, x$z[[indx[i]]] - mean(x$z[[indx[i]]]))
        wmin <- min(wmin, x$w[[indx[i]]] - mean(x$w[[indx[i]]]))
        wmax <- max(wmax, x$w[[indx[i]]] - mean(x$w[[indx[i]]]))
      }
      o <- order(x$z[[indx[1]]])
      plot(x$z[[indx[1]]][o] - mean(x$z[[indx[1]]]), 
           x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), 
           type = "l", xlab = "z (centered at zero)", 
           ylab = "w (centered at zero)", col = col[1], 
           main = paste0("ESS samples (middle layer)"), xlim = c(zmin, zmax),
           ylim = c(wmin, wmax))
      for (j in 2:length(indx)) {
        o <- order(x$z[[indx[j]]])
        lines(x$z[[indx[j]]][o] - mean(x$z[[indx[j]]]), 
              x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), col = col[j])
      }
    } else cat("Default plotting not prepared for these dimensions")
  }
  
  if (predict) {
    if (Dx == 1) {
      par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        Sigma_smooth <- x$Sigma - diag(mean(x$g * x$tau2), nrow(x$x_new))
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, Sigma_smooth))
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = "l", xlab = "x", ylab = "y", 
           ylim = c(min(q1), max(q3)),
           col = "blue", ...)
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = "grey", lty = 1)
        lines(x$x_new[o], x$mean[o], col = "blue")
      }
      lines(x$x_new[o], q1[o], col = "blue", lty = 2)
      lines(x$x_new[o], q3[o], col = "blue", lty = 2)
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
      image(i1, col = cols, main = "Posterior Mean", xlab = "x1", ylab = "x2")
      contour(i1, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = "Posterior Std. Dev.", xlab = "x1", ylab = "x2")
      contour(i2, add = TRUE)
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat("Dimension of X too large for default plotting")
  }
}

# Plot Three Layer Vecchia ----------------------------------------------------
#' @rdname plot
#' @export

plot.dgp3vec <- plot.dgp3

