
# Imported Functions ----------------------------------------------------------
#' @importFrom grDevices heat.colors
#' @importFrom graphics image lines matlines par plot points contour
#' @importFrom stats cov dgamma dnorm pnorm qnorm rnorm runif var
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm

# Package Documentation -------------------------------------------------------
#' @useDynLib deepgp, .registration = TRUE
#' @title Package deepgp
#' @author Annie Sauer \email{anniees@vt.edu}
#' @docType package
#' @name deepgp-package
#'
#' @description Performs model fitting and sequential design for deep Gaussian
#'     processes following Sauer, Gramacy, and Higdon (2020) <arXiv:2012.08015>.  
#'     Models extend up to three layers deep; a one layer model is equivalent 
#'     to typical Gaussian process regression.  Both Matern and squared exponential
#'     kernels are implemented.  Sequential design criteria 
#'     include integrated mean-squared error (IMSE), active learning Cohn (ALC), 
#'     and expected improvement (EI).  Applicable to both noisy and 
#'     deterministic functions.  Incorporates SNOW parallelization and 
#'     utilizes C and C++ under the hood.
#' 
#' @section Important Functions:
#' \itemize{
#'   \item \code{\link[deepgp]{fit_one_layer}}: conducts MCMC sampling of 
#'   hyperparameters for a one layer GP
#'   \item \code{\link[deepgp]{fit_two_layer}}: conducts MCMC sampling of 
#'   hyperparameters and hidden layer for a two layer deep GP
#'   \item \code{\link[deepgp]{fit_three_layer}}: conducts MCMC sampling of 
#'   hyperparameters and hidden layers for a three layer deep GP
#'   \item \code{\link[deepgp]{continue}}: collects additional MCMC samples
#'   \item \code{\link[deepgp]{trim}}: cuts off burn-in and optionally thins 
#'   samples
#'   \item \code{\link[deepgp]{predict}}: calculates posterior mean and 
#'   variance over a set of input locations (optionally calculates EI)
#'   \item \code{\link[deepgp]{plot}}: produces trace plots, hidden layer 
#'   plots, and posterior plots
#'   \item \code{\link[deepgp]{ALC}}: calculates active learning Cohn over 
#'   set of input locations using reference grid
#'   \item \code{\link[deepgp]{IMSE}}: calculates integrated mean-squared error
#'    over set of input locations
#' }
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." \emph{Technometrics, to appear;} arXiv:2012.08015. \cr\cr
#' Binois, M, J Huang, RB Gramacy, and M Ludkovski. 2019. Replication or 
#'     Exploration? Sequential Design for Stochastic Simulation Experiments. 
#'     \emph{Technometrics 61}, 7-23. Taylor & Francis. 
#'     doi:10.1080/00401706.2018.1469433.\cr\cr
#' Gramacy, RB. \emph{Surrogates: Gaussian Process Modeling, Design, and 
#'     Optimization for the Applied Sciences}. Chapman Hall, 2020.\cr\cr
#' Jones, DR, M Schonlau, and WJ Welch. 1998. "Efficient Global Optimization 
#'     of Expensive Black-Box Functions." \emph{Journal of Global Optimization 
#'     13}, 455-492. doi:10.1023/A:1008306431147.\cr\cr
#' Murray, I, RP Adams, and D MacKay. 2010. "Elliptical slice sampling." 
#'     \emph{Journal of Machine Learning Research 9}, 541-548.\cr\cr
#' Seo, S, M Wallat, T Graepel, and K Obermayer. 2000. Gaussian Process 
#'     Regression: Active Data Selection and Test Point Rejection. In 
#'     Mustererkennung 2000, 27-34. New York, NY: Springer Verlag.
#' 
#' @examples \donttest{
#' # 1. One Layer and EI ---------------------------------------------------------
#' 
#' f <- function(x) {
#'   sin(5 * pi * x) / (2 * x) + (x - 1) ^ 4
#' }
#' 
#' # Training data
#' x <- seq(0.5, 2, length = 30)
#' y <- f(x) + rnorm(30, 0, 0.01)
#' 
#' # Testing data
#' xx <- seq(0.5, 2, length = 100)
#' yy <- f(xx)
#' 
#' # Standardize inputs and outputs
#' xx <- (xx - min(x)) / (max(x) - min(x))
#' x <- (x - min(x)) / (max(x) - min(x))
#' yy <- (yy - mean(y)) / sd(y)
#' y <- (y - mean(y)) / sd(y)
#' 
#' # Conduct MCMC
#' fit <- fit_one_layer(x, y, nmcmc = 10000)
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#' 
#' # Predict and calculate EI
#' fit <- predict(fit, xx, EI = TRUE)
#' 
#' # Visualize Fit
#' plot(fit)
#' par(new = TRUE) # overlay EI
#' plot(xx, fit$EI, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.max(fit$EI)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' 
#' # 2. Two Layer and ALC ------------------------------------------------------
#' 
#' f <- function(x) {
#'   exp(-10 * x) * (cos(10 * pi * x - 1) + sin(10 * pi * x - 1)) * 5 - 0.2
#' }
#' 
#' # Training data
#' x <- seq(0, 1, length = 30)
#' y <- f(x) + rnorm(30, 0, 0.05)
#' 
#' # Testing data
#' xx <- seq(0, 1, length = 100)
#' yy <- f(xx)
#' 
#' # Conduct MCMC
#' fit <- fit_two_layer(x, y, D = 1, nmcmc = 9000, cov = "exp2")
#' fit <- continue(fit, 1000)
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#' 
#' # Option 1 - calculate ALC from MCMC iterations
#' alc <- ALC(fit, xx)
#' 
#' # Option 2 - calculate ALC after predictions
#' fit <- predict(fit, xx, store_latent = TRUE)
#' alc <- ALC(fit)
#' 
#' # Visualize fit
#' plot(fit)
#' par(new = TRUE) # overlay ALC
#' plot(xx, alc$value, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.max(alc$value)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' 
#' # 3. Three Layer and IMSE ------------------------------------------------------
#' 
#' f <- function(x) {
#'   i <- which(x <= 0.48)
#'   x[i] <- 2 * sin(pi * x[i] * 4) + 0.4 * cos(pi * x[i] * 16)
#'   x[-i] <- 2 * x[-i] - 1
#'   return(x)
#' }
#' 
#' # Training data
#' x <- seq(0, 1, length = 30)
#' y <- f(x) + rnorm(30, 0, 0.05)
#' 
#' # Testing data
#' xx <- seq(0, 1, length = 100)
#' yy <- f(xx)
#' 
#' # Conduct MCMC
#' fit <- fit_three_layer(x, y, D = 1, nmcmc = 10000, cov = "exp2")
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#' 
#' # Option 1 - calculate IMSE from only MCMC iterations
#' imse <- IMSE(fit, xx)
#' 
#' # Option 2 - calculate IMSE after predictions
#' fit <- predict(fit, xx, store_latent = TRUE)
#' imse <- IMSE(fit)
#' 
#' # Visualize fit
#' plot(fit)
#' par(new = TRUE) # overlay IMSE
#' plot(xx, imse$value, col = 2, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.min(imse$value)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' }

NULL