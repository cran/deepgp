
# Function Contents -----------------------------------------------------------
# Internal:
#   check_settings: checks hyperparameter proposal/prior specifications
#   check_initialization: sets default initial values for hyperparameters and 
#                         hidden layers
#   check_inputs: errors/warnings for input dimensions/class/etc
#   check_ordering: errors/warnings for user-specifed Vecchia orderings

# Check Settings --------------------------------------------------------------

check_settings <- function(settings, layers = 1, D = NULL) {
  
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

# Check Initialization --------------------------------------------------------

check_initialization <- function(initial, layers = 2, x = NULL, D = NULL,
                                 vecchia = NULL, v = NULL, m = NULL) {
  
  if (is.null(initial$tau2)) initial$tau2 <- 1
  
  if (layers == 2) {
    
    if (is.null(initial$w)) 
      initial$w <- suppressWarnings(matrix(x, nrow = nrow(x), ncol = D))
    if (!is.matrix(initial$w)) 
      initial$w <- as.matrix(initial$w)
    if (ncol(initial$w) != D) 
      stop("dimension of initial$w does not match D")
    if (length(initial$theta_w) == 1) 
      initial$theta_w <- rep(initial$theta_w, D)
    
    # If initial$w is from previous sequential step, predict at new point
    if (nrow(initial$w) == nrow(x) - 1) {
      if (vecchia) {
        initial$w <- fill_final_row_vec(x, initial$w, D, initial$theta_w, v, m)
      } else {
        initial$w <- fill_final_row(x, initial$w, D, initial$theta_w, v)
      }
    }
    
  } else if (layers == 3) {
    
    if (is.null(initial$w)) 
      initial$w <- suppressWarnings(matrix(x, nrow = nrow(x), ncol = D))
    if (!is.matrix(initial$w)) 
      initial$w <- as.matrix(initial$w)
    if (ncol(initial$w) != D) 
      stop("dimension of initial$w does not match D")
    if (is.null(initial$z)) 
      initial$z <- suppressWarnings(matrix(x, nrow = nrow(x), ncol = D))
    if (!is.matrix(initial$z)) 
      initial$z <- as.matrix(initial$z)
    if (ncol(initial$z) != D) 
      stop("dimension of initial$z does not match D")
    if (length(initial$theta_w) == 1) 
      initial$theta_w <- rep(initial$theta_w, D)
    if (length(initial$theta_z) == 1) 
      initial$theta_z <- rep(initial$theta_z, D)
    
    # If initial$z/initial$w are from previous sequential step, predict at new point
    if (nrow(initial$z) == nrow(x) - 1) {
      if (vecchia) {
        initial$z <- fill_final_row_vec(x, initial$z, D, initial$theta_z, v, m)
      } else {
        initial$z <- fill_final_row(x, initial$z, D, initial$theta_z, v)
      }
    }
    if (nrow(initial$w) == nrow(x) - 1) {
      if (vecchia) {
        initial$w <- fill_final_row_vec(initial$z, initial$w, D, initial$theta_w, v, m)
      } else {
        initial$w <- fill_final_row(initial$z, initial$w, D, initial$theta_w, v)
      }
    }
  }
  return(initial)
}

# Check Inputs ----------------------------------------------------------------

check_inputs <- function(x, y, true_g) {
  
  if (!is.vector(y)) 
    stop("y must be a vector")
  if (nrow(x) != length(y)) 
    stop("dimensions of x and y do not match")
  if (min(x) < -5 | min(x) > 5 | max(x) < -4 | max(x) > 6) 
    warning("this function is designed for x over the range [0, 1]")
  if (is.null(true_g) & (mean(y) < -10 | mean(y) > 10 | var(y) < 0.1 | var(y) > 10))
    warning("this function is designed for y scaled to mean zero and variance 1")
  
  return(NULL)
}

# Check Ordering --------------------------------------------------------------

check_ordering <- function(ordering, n) {
  
  if (length(ordering) != n) 
    stop("length(ordering) must match dimension of x for fitting or x_new for predicting")
  if (min(ordering) != 1)
    stop("ordering must begin at index 1")
  if (max(ordering) != n)
    stop("ordering must end at index nrow(x) for fitting or nrow(x_new) for predicting")
  if (sum(duplicated(ordering)) > 0)
    stop("ordering must not have any duplicates")
  
  return(NULL)
}

