
# Package Documentation -------------------------------------------------------
#' @useDynLib deepgp, .registration = TRUE
#' @title Package deepgp
#' @author Annie Sauer \email{anniees@vt.edu}
#' @docType package
#' @name deepgp-package
#'
#' @description Performs model fitting and sequential design for deep Gaussian
#'     processes using MCMC and elliptical slice sampling.  Models extend up to 
#'     three layers deep; a one layer model is equivalent to typical Gaussian 
#'     process regression.  Sequential design criteria include integrated mean 
#'     square prediction error (IMSPE), active learning Cohn (ALC), and expected 
#'     improvement (EI).  Covariance structure is based on inverse exponentiated 
#'     squared euclidean distance.  Applicable to noisy and deterministic functions.  
#'     Incorporates SNOW parallelization and utilizes C under the hood.  Manuscript 
#'     forthcoming; see Damianou and Lawrence (2013) <arXiv:1211.0358> for deep 
#'     Gaussian process models and Murray, Adams, and MacKay (2010) <arXiv:1001.0175> 
#'     for elliptical slice sampling.
#' 
#' @section Important Functions:
#' \itemize{
#'   \item \code{\link[deepgp]{fit_one_layer}}: conducts MCMC sampling of 
#'   hyperparameters for a one layer GP
#'   \item \code{\link[deepgp]{fit_two_layer}}: conducts MCMC sampling of 
#'   hyperparameters for a two layer deep GP
#'   \item \code{\link[deepgp]{fit_three_layer}}: conducts MCMC sampling of 
#'   hyperparameters for a three layer deep GP
#'   \item \code{\link[deepgp]{continue}}: collects additional MCMC samples
#'   \item \code{\link[deepgp]{trim}}: cuts off burn-in and optionally thins 
#'   samples
#'   \item \code{\link[deepgp]{predict}}: calculates posterior mean and 
#'   covariance over a set of input locations
#'   \item \code{\link[deepgp]{plot}}: produces trace plots, hidden layer 
#'   plots, and posterior plots
#'   \item \code{\link[deepgp]{ALC}}: calculates active learning Cohn over 
#'   set of input locations using reference grid
#'   \item \code{\link[deepgp]{IMSPE}}: calculates integrated mean square 
#'   prediction error over set of input locations
#'   \item \code{\link[deepgp]{EI}}: calculates expected improvement over set 
#'   of input locations
#' }
#' 
#' @references 
#' Binois, M, J Huang, RB Gramacy, and M Ludkovski. 2019. “Replication or 
#'     Exploration? Sequential Design for Stochastic Simulation Experiments.” 
#'     \emph{Technometrics 61}, 7-23. Taylor & Francis. 
#'     doi:10.1080/00401706.2018.1469433.\cr\cr
#' Damianou, A and N Lawrence. (2013). "Deep gaussian processes." 
#'     \emph{Artificial Intelligence and Statistics}, 207-215.\cr\cr
#' Gramacy, RB. \emph{Surrogates: Gaussian Process Modeling, Design, and 
#'     Optimization for the Applied Sciences}. Chapman Hall, 2020.\cr\cr
#' Jones, DR, M Schonlau, and WJ Welch. 1998. "Efficient Global Optimization 
#'     of Expensive Black-Box Functions." \emph{Journal of Global Optimization 
#'     13}, 455-492. doi:10.1023/A:1008306431147.\cr\cr
#' Murray, I, RP Adams, and D MacKay. 2010. "Elliptical slice sampling." 
#'     \emph{Journal of Machine Learning Research 9}, 541-548.\cr\cr
#' Seo, S, M Wallat, T Graepel, and K Obermayer. 2000. “Gaussian Process 
#'     Regression: Active Data Selection and Test Point Rejection.” In 
#'     Mustererkennung 2000, 27–34. New York, NY: Springer–Verlag.
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
#' fit <- predict(fit, xx, lite = FALSE)
#' ei <- EI(fit)
#'   
#' # Visualize Fit
#' plot(fit)
#' par(new = TRUE) # overlay EI
#' plot(xx, ei$value, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.max(ei$value)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' score(yy, fit$mean, fit$Sigma) # higher is better
#' 
#' # 2. Two Layer and ALC --------------------------------------------------------
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
#' fit <- fit_two_layer(x, y, D = 1, nmcmc = 9000)
#' fit <- continue(fit, 1000)
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#' 
#' # Option 1 - calculate ALC from MCMC iterations
#' alc <- ALC(fit, xx)
#' 
#' # Option 2 - calculate ALC after predictions
#' fit <- predict(fit, xx)
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
#' score(yy, fit$mean, fit$Sigma) # higher is better
#' 
#' # 3. Three Layer and IMSPE ----------------------------------------------------
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
#' fit <- fit_three_layer(x, y, D = 1, nmcmc = 10000)
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#' 
#' # Option 1 - calculate IMSPE from only MCMC iterations
#' imspe <- IMSPE(fit, xx)
#' 
#' # Option 2 - calculate IMSPE after predictions
#' fit <- predict(fit, xx)
#' imspe <- IMSPE(fit)
#' 
#' # Visualize fit
#' plot(fit)
#' par(new = TRUE) # overlay IMSPE
#' plot(xx, imspe$value, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.min(imspe$value)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' score(yy, fit$mean, fit$Sigma) # higher is better
#' }

NULL