
# Function Contents -----------------------------------------------------------
# Internal:
#   check_settings: checks hyperparameter proposal/prior specifications
#   check_inputs: errors/warnings for input dimensions/class/etc

check_settings <- function(settings, layers = 1) {
  if (is.null(settings$l)) settings$l <- 1
  if (is.null(settings$u)) settings$u <- 2
  if (is.null(settings$alpha$g)) settings$alpha$g <- 1.5
  if (is.null(settings$beta$g)) settings$beta$g <- 3.9
  if (layers == 1) {
    if (is.null(settings$alpha$theta)) settings$alpha$theta <- 1.5
    if (is.null(settings$beta$theta)) settings$beta$theta <- 3.9/1.5
  } else if (layers == 2) {
    if (is.null(settings$alpha$theta_w)) settings$alpha$theta_w <- 1.5
    if (is.null(settings$alpha$theta_y)) settings$alpha$theta_y <- 1.5
    if (is.null(settings$beta$theta_w)) settings$beta$theta_w <- 3.9/4
    if (is.null(settings$beta$theta_y)) settings$beta$theta_y <- 3.9/6
  } else if (layers == 3) {
    if (is.null(settings$alpha$theta_z)) settings$alpha$theta_z <- 1.5
    if (is.null(settings$alpha$theta_w)) settings$alpha$theta_w <- 1.5
    if (is.null(settings$alpha$theta_y)) settings$alpha$theta_y <- 1.5
    if (is.null(settings$beta$theta_z)) settings$beta$theta_z <- 3.9/4
    if (is.null(settings$beta$theta_w)) settings$beta$theta_w <- 3.9/12
    if (is.null(settings$beta$theta_y)) settings$beta$theta_y <- 3.9/6
  }
  return(settings)
}

check_inputs <- function(x, y, true_g, w_0 = NULL, D = NULL, z_0 = NULL) {
  if (!is.vector(y)) stop("y must be a vector")
  if (nrow(x) != length(y)) stop("dimensions of x and y do not match")
  if (min(x) < -1 | min(x) > 1 | max(x) < 0 | max(x) > 2) 
    warning("this function is designed for x over the range [0, 1]")
  if (is.null(true_g) & (mean(y) < -5 | mean(y) > 5 | var(y) < 0.1 | var(y) > 10))
    warning("this function is designed for y scaled to mean zero and variance 1")
  if (!is.null(w_0))
    if (ncol(w_0) != D) stop("dimension of w_0 does not match D")
  if (!is.null(z_0))
    if (ncol(z_0) != D) stop("dimension of z_0 does not match D")
  return(NULL)
}