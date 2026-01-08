
# Function Contents -----------------------------------------------------------
# Internal:
#   check_settings: checks hyperparameter proposal/prior specifications
#   check_initialization: sets default initial values for hyperparameters and 
#                         hidden layers
#   check_inputs: errors/warnings for input dimensions/scaling/etc
#   check_vecchia: errors/warnings for Vecchia settings
#   check_gradients: errors/warnings for gradient enhancement
#   check_cores: checks if requested cores is too high, sets default values

# check_settings --------------------------------------------------------------

check_settings <- function(settings, layers = 1, noisy = FALSE, monowarp = FALSE) {
  
  # Metropolis Hastings proposal sliding windows
  if (is.null(settings$l)) settings$l <- 1
  if (is.null(settings$u)) settings$u <- 2
  
  # Prior distributions for g and theta
  if (noisy) {
    if (is.null(settings$g$alpha)) settings$g$alpha <- 1.01
    if (is.null(settings$g$beta)) settings$g$beta <- 10
      # mode = 0.001, 95% quantile = 0.3
    
    if (layers == 1) {
      if (is.null(settings$theta$alpha)) settings$theta$alpha <- 1.2
      if (is.null(settings$theta$beta)) settings$theta$beta <- 4
        # mode = 0.05, 95% quantile = 0.84
    } else if (layers == 2) {
      if (monowarp) {
        if (is.null(settings$theta_w$alpha)) settings$theta_w$alpha <- 1.2
        if (is.null(settings$theta_w$beta)) settings$theta_w$beta <- 0.2
      } else {
        if (is.null(settings$theta_w$alpha)) settings$theta_w$alpha <- 1.2
        if (is.null(settings$theta_w$beta)) settings$theta_w$beta <- 2
      }
        # mode = 0.1, 95% quantile = 1.7  
      if (is.null(settings$theta_y$alpha)) settings$theta_y$alpha <- 1.2
      if (is.null(settings$theta_y$beta)) settings$theta_y$beta <- 1
        # mode = 0.2, 95% quantile = 3.4
    } else if (layers == 3) {
      if (is.null(settings$theta_z$alpha)) settings$theta_z$alpha <- 1.2
      if (is.null(settings$theta_z$beta)) settings$theta_z$beta <- 2
      if (is.null(settings$theta_w$alpha)) settings$theta_w$alpha <- 1.2
      if (is.null(settings$theta_w$beta)) settings$theta_w$beta <- 0.8
      if (is.null(settings$theta_y$beta)) settings$theta_y$beta <- 1.2
      if (is.null(settings$theta_y$alpha)) settings$theta_y$alpha <- 1
    }
  } else { # Original default values 
    if (layers == 1) {
      if (is.null(settings$theta$alpha)) settings$theta$alpha <- 1.5
      if (is.null(settings$theta$beta)) settings$theta$beta <- 3.9/1.5
    } else if (layers == 2) {
      if (monowarp) {
        if (is.null(settings$theta_w$alpha)) settings$theta_w$alpha <- 1.5
        if (is.null(settings$theta_w$beta)) settings$theta_w$beta <- 0.39/4
      } else {
        if (is.null(settings$theta_w$alpha)) settings$theta_w$alpha <- 1.5
        if (is.null(settings$theta_w$beta)) settings$theta_w$beta <- 3.9/4
      }  
      if (is.null(settings$theta_y$alpha)) settings$theta_y$alpha <- 1.5
      if (is.null(settings$theta_y$beta)) settings$theta_y$beta <- 3.9/6
    } else if (layers == 3) {
      if (is.null(settings$theta_z$alpha)) settings$theta_z$alpha <- 1.5
      if (is.null(settings$theta_z$beta)) settings$theta_z$beta <- 3.9/4
      if (is.null(settings$theta_w$alpha)) settings$theta_w$alpha <- 1.5
      if (is.null(settings$theta_w$beta)) settings$theta_w$beta <- 3.9/12
      if (is.null(settings$theta_y$alpha)) settings$theta_y$alpha <- 1.5
      if (is.null(settings$theta_y$beta)) settings$theta_y$beta <- 3.9/6
    }
  }

  # Fixed values of tau2 for latent layers
  if (layers == 2) {
    if (is.null(settings$tau2_w)) settings$tau2_w <- 1
  } else if (layers == 3) {
    if (is.null(settings$tau2_w)) settings$tau2_w <- 1
    if (is.null(settings$tau2_z)) settings$tau2_z <- 1
  }

  return(settings)
}

# check_initialization --------------------------------------------------------

check_initialization <- function(initial, layers, D, grad_enhance = FALSE, 
                                 x = NULL, v = v, pmx = FALSE, 
                                 vecchia = FALSE, m = NULL) {
  
  n <- nrow(x)
  d <- ncol(x)

  if (layers == 2) {
    
    if (!is.null(initial$theta_w))
      if (length(initial$theta_w) == 1)
        initial$theta_w <- rep(initial$theta_w, times = D)
    
    if (grad_enhance) { # D == d already checked
      if (is.null(initial$w)) { # start at identity warping (with appropriate derivatives)
        initial$w <- get_prior_mean(x)
      } else if (nrow(initial$w) != (n*(d+1)))
        stop("dimension of w_0 must be n*(d+1)")
    } else {
      if (is.null(initial$w))
        initial$w <- suppressWarnings(matrix(x, nrow = n, ncol = D))
      if (nrow(initial$w) != n)
        initial$w <- fill_final_rows(initial$w, x, tau2 = 1, theta = initial$theta_w, 
                                     v = v, pmx = pmx, vecchia = vecchia, m = m)
    }
    
    if (!is.matrix(initial$w)) 
      initial$w <- as.matrix(initial$w)
    if (ncol(initial$w) != D) 
      stop("dimension of initial$w does not match D")
    
  } else if (layers == 3) {
    
    if (is.null(initial$z)) 
      initial$z <- suppressWarnings(matrix(x, nrow = n, ncol = D))
    if (!is.matrix(initial$z)) 
      initial$z <- as.matrix(initial$z)
    if (ncol(initial$z) != D) 
      stop("dimension of initial$z does not match D")
    if (is.null(initial$w)) 
      initial$w <- suppressWarnings(matrix(x, nrow = n, ncol = D))
    if (!is.matrix(initial$w)) 
      initial$w <- as.matrix(initial$w)
    if (ncol(initial$w) != D) 
      stop("dimension of initial$w does not match D")
    if (nrow(initial$z) != n)
      initial$z <- fill_final_rows(initial$z, x, tau2 = 1, theta = initial$theta_z, 
                                   v = v, vecchia = vecchia, m = m)
    if (nrow(initial$w) != n)
      initial$w <- fill_final_rows(initial$w, initial$z, tau2 = 1, theta = initial$theta_w, 
                                   v = v, vecchia = vecchia, m = m)
  }
  return(initial)
}

# check_inputs ----------------------------------------------------------------

check_inputs <- function(x, y, true_g, nmcmc) {
  
  if (!is.vector(y)) 
    stop("y must be a vector")
  if (nrow(x) != length(y)) 
    stop("dimensions of x and y do not match")
  if (min(x) < -5 | min(x) > 5 | max(x) < -4 | max(x) > 6) 
    message("Warning: this function is designed for x over the range [0, 1]")
  if (is.null(true_g) & (mean(y) < -10 | mean(y) > 10 | var(y) < 0.1 | var(y) > 10))
    message("Warning: this function is designed for y scaled to mean zero and variance 1")
  if (nmcmc <= 1) 
    stop("nmcmc must be greater than 1")

  return(NULL)
}

# check_vecchia ---------------------------------------------------------------

check_vecchia <- function(n, d, m, ord, grad_enhance = FALSE) {
  
  test <- check_omp() # returns NULL
  if (grad_enhance) {
    if (m >= n*(d+1)) stop("m must be less than n*(d+1)")
  } else if (m >= n) stop("m must be less than n")
  if (!is.null(ord)) {
    if (min(ord) != 1)
      stop("ordering must begin at index 1")
    if (sum(duplicated(ord)) > 0)
      stop("ordering must not have any duplicates")
    if (grad_enhance) {
      if (length(ord) != n | max(ord) != n)
        stop(paste0("ord must be integers 1:nrow(x), gradients ",
              "are appended in the same order"))
    } else {
      if (length(ord) != n | max(ord) != n) 
        stop("ord must be integers 1:nrow(x)")
    }
  }
  
  return(NULL)
}

# check_gradients -------------------------------------------------------------

check_gradients <- function(n, d, dydx, cov, true_g, D = NULL, vecchia = FALSE, 
                            monowarp = FALSE) {
  
  if (n != nrow(dydx))
    stop("nrow(x) must match nrow(dydx)")
  if (d != ncol(dydx))
    stop("ncol(x) must match ncol(dydx)")
  if (!is.null(D))
    if (D != d) 
      stop("grad_enhance requires D = ncol(x)")
  if (monowarp | is.numeric(monowarp)) 
    stop("monowarp not implemented for grad_enhance")
  if (is.null(true_g)) {
    message("Warning: gradient-enhancement is designed for deterministic data.")
    message("Setting true_g to a small value is recommended.")
    message("Gradients are assumed to be noise-free.")
  }
  if (!vecchia & (d+1)*n > 200)
    message("vecchia = TRUE recommended for faster computation")
  if (cov == "matern")
    stop("gradient-enhancement requires cov = 'exp2")
    
  return(NULL)
}

# check_cores -----------------------------------------------------------------

check_cores <- function(cores, nmcmc = NULL) {
  
  maxcores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
  if (is.null(cores)) {
    cores <- min(4, maxcores - 1)
  } else {
    if (cores > maxcores) {
      message("cores is greater than available nodes, overwriting with maxcores-1")
      cores <- maxcores - 1
    }
  }
  if (!is.null(nmcmc)) {
    if (cores > nmcmc) {
      message("cores is greater than nmcmc, overwriting with nmcmc")
      cores <- nmcmc
    }
  }
    
  return(cores)
}

