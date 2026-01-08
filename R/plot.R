
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   plot.gp == plot.gpvec
#   plot.dgp2 == plot.dgp2vec
#   plot.dgp3 == plot.dgp3vec

# plot S3 class ---------------------------------------------------------------
#' @name plot
#' @title Plots object from \code{deepgp} package
#' 
#' @description Acts on a \code{gp}, \code{gpvec}, \code{dgp2}, \code{dgp2vec},
#'     \code{dgp3}, or \code{dgp3vec} object.  
#'     Generates trace plots for outer log likelihood, length scale,
#'     and nugget hyperparameters.
#'     Generates plots of hidden layers for low dimensions or 
#'     monotonic warpings.  Generates
#'     plots of the posterior mean and estimated 90\% prediction intervals for 
#'     one-dimensional inputs; generates heat maps of the posterior mean and 
#'     point-wise variance for two-dimensional inputs.
#'     
#' @details Trace plots are useful in assessing burn-in.  If there are too
#'     many hyperparameters to plot them all, then it is most useful to 
#'     visualize the log likelihood (e.g., \code{plot(fit$ll, type = "l")}).
#'
#'     In one dimension, hidden layer plots show 100 evenly distributed samples.
#'     In two dimensions (two layer only), hidden layer plots show 3 samples.
#'     Hidden layer plots are colored on a gradient - red lines represent 
#'     earlier iterations and yellow lines represent later iterations - to 
#'     help assess burn-in of the hidden layers.  
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

# plot.gp ---------------------------------------------------------------------
#' @rdname plot
#' @export

plot.gp <- function(x, trace = NULL, predict = NULL, ...) {
  
  # save and restore par settings 
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  d <- ncol(x$x)
  if (is.null(x$mean)) {
    if (is.null(trace)) trace <- TRUE
    if (is.null(predict)) predict <- FALSE
  } else {
    if (is.null(trace)) trace <- FALSE
    if (is.null(predict)) predict <- TRUE
  }
  
  if (trace) {
    fixed_g <- (length(x$g) == 1)
    if (is.matrix(x$theta)) {
      nplots <- ncol(x$theta) + 2 + as.numeric(!fixed_g)
    } else nplots <- 3 + as.numeric(!fixed_g)
    if (nplots > 4) {
      par(mfrow = c(2, ceiling(nplots/2)), mar = c(5, 4, 2, 2))
    } else par(mfrow = c(1, nplots), mar = c(5, 4, 2, 2))
    plot(x$ll, type = "l", ylab = "logl", xlab = "Iteration", main = "logl")
    plot(x$tau2, type = "l", ylab = "tau2", xlab = "Iteration", main = "tau2")
    if (!fixed_g) plot(x$g, type = "l", ylab = "g", xlab = "Iteration", main = "g")
    if (is.matrix(x$theta)) { # separable lengthscale
      for (i in 1:ncol(x$theta)) {
        plot(x$theta[, i], type = "l", ylab = paste0("theta[", i, "]"), 
             xlab = "Iteration", main = paste0("theta[", i, "]"))
      }
    } else { # isotropic lengthscale
      plot(x$theta, type = "l", ylab = "theta", xlab = "Iteration", main = "theta")
    }
  }
  
  if (predict) {
    if (d == 1) {
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        Sigma_smooth <- x$Sigma - diag(mean(x$g * x$tau2), nrow(x$x_new))
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, Sigma_smooth, checkSymmetry = FALSE))
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = "l", xlab = "x", ylab = "y", 
           ylim = c(min(q1), max(q3)), col = "blue", ...)
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = "grey", lty = 1)
        lines(x$x_new[o], x$mean[o], col = "blue")
      }
      lines(x$x_new[o], q1[o], col = "blue", lty = 2)
      lines(x$x_new[o], q3[o], col = "blue", lty = 2)
      points(x$x, x$y, pch = 20)
    } else if (d == 2) {
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

# plot.gpvec ------------------------------------------------------------------
#' @rdname plot
#' @export

plot.gpvec <- plot.gp

# plot.dgp2 -------------------------------------------------------------------
#' @rdname plot
#' @export

plot.dgp2 <- function(x, trace = NULL, hidden = NULL, predict = NULL, ...) {
  
  # save and restore par settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  n <- nrow(x$x)
  d <- ncol(x$x)
  D <- dim(x$w)[3]
  if (is.null(x$mean)) {
    if (is.null(trace)) trace <- TRUE
    if (is.null(predict)) predict <- FALSE
  } else {
    if (is.null(trace)) trace <- FALSE
    if (is.null(predict)) predict <- TRUE
  }
  if (is.null(hidden)) hidden <- FALSE
  
  if (trace) {
    fixed_g <- (length(x$g) == 1)
    nplots <- D + 3 + as.numeric(!fixed_g) + ifel(x$settings$monowarp & d > 1, D, 0)
    if (nplots > 4) {
      par(mfrow = c(2, ceiling(nplots/2)), mar = c(5, 4, 2, 2))
    } else par(mfrow = c(1, nplots), mar = c(5, 4, 2, 2))
    plot(x$ll, type = "l", ylab = "outer logl", xlab = "Iteration", main = "outer logl")
    plot(x$tau2_y, type = "l", ylab = "tau2_y", xlab = "Iteration", main = "tau2_y")
    if (!fixed_g) plot(x$g, type = "l", ylab = "g", xlab = "Iteration", main = "g")
    plot(x$theta_y, type = "l", ylab = "theta_y", xlab = "Iteration", main = "theta_y")
    if (x$settings$monowarp & d > 1) {
      for (i in 1:D) {
        plot(x$tau2_w[, i], type = "l", ylab = paste0("tau2_w[", i, "]"), 
           xlab = "Iteration", main = paste0("tau2_w[", i, "]"))
      }
    }
    for (i in 1:D) {
      plot(x$theta_w[, i], type = "l", ylab = paste0("theta_w[", i, "]"), 
           xlab = "Iteration", main = paste0("theta_w[", i, "]"))
    }
  }
  
  if (hidden) { 
    col <- heat.colors(100 + 30)

    if (x$settings$monowarp) { 
      # select 100 lines to plot (don't want things to be too cluttered)
      indx <- floor(seq(from = 1, to = x$nmcmc, length = 100))
      w <- x$w[indx, , , drop = FALSE]
      par(mfrow = c(1, D), mar = c(4, 4, 2, 2))
      for (i in 1:D) {
        o <- order(x$x[, i])
        matplot(x$x[o, i], t(w[, o, i] - apply(w[, , i], 1, mean)), 
                type = "l", xlab = paste0("x", i),
                ylab = paste0("w", i), col = col, lty = 1,
                main = "monowarped ESS samples")
      }
    } else if (d == 1 & D == 1) {
      # select 100 lines to plot (don't want things to be too cluttered)
      indx <- floor(seq(from = 1, to = x$nmcmc, length = 100))
      o <- order(x$x)
      w <- x$w[indx, , 1]
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
      matplot(x$x[o], t(w[, o] - apply(w, 1, mean)), type = "l", xlab = "x", 
           ylab = "w", col = col, lty = 1, main = "ESS samples")
    } else if (d == 2 & D == 2) {
      if (!requireNamespace("interp", quietly = TRUE)) {
        stop("Package \"interp\" needed for this plot. Please install it.",
             call. = FALSE)
      }
      par(mfcol = c(2, 3), mar = c(4, 4, 2, 2))
      # select 3 samples to plot
      for (i in floor(seq(from = 1, to = x$nmcmc, length = 3))) {
        ii <- interp::interp(x$x[, 1], x$x[, 2], x$w[i, 1:n, 1] - mean(x$w[i, 1:n, 1]))
        image(ii, col = heat.colors(128), xlab = "x1", ylab = "x2", 
            main = paste0("w1, iteration ", i))
        points(x$x, pch = 20, cex = 0.5)
        ii <- interp::interp(x$x[, 1], x$x[, 2], x$w[i, 1:n, 2] - mean(x$w[i, 1:n, 1]))
        image(ii, col = heat.colors(128), xlab = "x1", ylab = "x2", 
            main = paste0("w2, iteration ", i))
        points(x$x, pch = 20, cex = 0.5)
      }
    } else cat("Default plotting not prepared for these dimensions")
  } 

  if (predict) {
    n <- length(x$y)
    if (d == 1) {
      par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        Sigma_smooth <- x$Sigma - diag(mean(x$g * x$tau2_y), nrow(x$x_new))
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, Sigma_smooth, checkSymmetry = FALSE))
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
    } else if (d == 2) {
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

# plot.dgp2vec ----------------------------------------------------------------
#' @rdname plot
#' @export

plot.dgp2vec <- plot.dgp2

# plot.dgp3 -------------------------------------------------------------------
#' @rdname plot
#' @export

plot.dgp3 <- function(x, trace = NULL, hidden = NULL, predict = NULL, ...) {
  
  # save and restore par settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  d <- ncol(x$x)
  D <- dim(x$w)[3]
  if (is.null(x$mean)) {
    if (is.null(trace)) trace <- TRUE
    if (is.null(predict)) predict <- FALSE
  } else {
    if (is.null(trace)) trace <- FALSE
    if (is.null(predict)) predict <- TRUE
  }
  if (is.null(hidden)) hidden <- FALSE
  
  if (trace) {
    fixed_g <- (length(x$g) == 1)
    nplots <- 2*D + 3 + as.numeric(!fixed_g)
    if (nplots > 4) {
      par(mfrow = c(2, ceiling(nplots/2)), mar = c(5, 4, 2, 2))
    } else par(mfrow = c(1, nplots), mar = c(5, 4, 2, 2))
    plot(x$ll, type = "l", ylab = "outer logl", xlab = "Iteration", main = "outer logl")
    plot(x$tau2_y, type = "l", ylab = "tau2_y", xlab = "Iteration", main = "tau2_y")
    plot(x$theta_y, type = "l", ylab = "theta_y", xlab = "Iteration", main = "theta_y")
    for (i in 1:D) {
      plot(x$theta_w[,i], type = "l", ylab = paste0("theta_w[", i, "]"), 
           xlab = "Iteration", main = paste0("theta_w[", i, "]"))
    }
    for (i in 1:D) {
      plot(x$theta_z[,i], type = "l", ylab = paste0("theta_z[", i, "]"), 
           xlab = "Iteration", main = paste0("theta_z[", i, "]"))
    }
  }
  
  if (hidden) {
    if (d == 1 & D == 1) {
      # randomly select 100 lines to plot (don't want things to be too cluttered)
      indx <- floor(seq(from = 1, to = x$nmcmc, length = 100))
      col <- heat.colors(100 + 30)
      z <- x$z[indx, , 1]
      par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
      o <- order(x$x)
      matplot(x$x[o], t(z[, o] - apply(z, 1, mean)), type = "l", xlab = "x", 
              ylab = "z", col = col, lty = 1, main = "ESS samples (inner layer)")
      w <- x$w[indx, , 1]
      o <- order(z[1, ])
      plot(z[1, o], w[1, o] - mean(w[1, ]), type = "l", xlab = "z (centered at zero)", 
           ylab = "w (centered at zero)", col = col[1], 
           main = "ESS samples (middle layer)", xlim = c(min(z), max(z)),
           ylim = c(min(w), max(w)))
      for (j in 2:length(indx)) {
        o <- order(z[j, ])
        lines(z[j, o], w[j, o] - mean(w[j, ]), col = col[j])
      }
    } else cat("Default plotting not prepared for these dimensions")
  }
  
  if (predict) {
    if (d == 1) {
      par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        Sigma_smooth <- x$Sigma - diag(mean(x$g * x$tau2_y), nrow(x$x_new))
        y_samples <- t(mvtnorm::rmvnorm(50, x$mean, Sigma_smooth, checkSymmetry = FALSE))
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
    } else if (d == 2) {
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

# plot.dgp3vec ----------------------------------------------------------------
#' @rdname plot
#' @export

plot.dgp3vec <- plot.dgp3

