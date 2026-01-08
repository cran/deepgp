
# Imported Functions ----------------------------------------------------------
#' @importFrom Matrix t solve
#' @importFrom grDevices heat.colors
#' @importFrom graphics image lines par plot points contour abline matplot matlines
#' @importFrom stats cov dgamma dnorm pnorm qnorm rnorm runif var approx
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm
#' @importFrom FNN get.knnx
#' @importFrom GpGp find_ordered_nn
#' @importFrom fields rdist
#' @importFrom abind abind

# Package Documentation -------------------------------------------------------
#' @useDynLib deepgp, .registration = TRUE
#' @title Package deepgp
#' @author Annie S. Booth \email{annie_booth@ncsu.edu}
#' @name deepgp-package
#'
#' @description Performs Bayesian posterior inference for deep Gaussian 
#' processes following Sauer, Gramacy, and Higdon (2023).  
#' See Sauer (2023) for comprehensive 
#' methodological details and \url{https://bitbucket.org/gramacylab/deepgp-ex/} for 
#' a variety of coding examples. Models are trained through MCMC including 
#' elliptical slice sampling of latent Gaussian layers and Metropolis-Hastings 
#' sampling of kernel hyperparameters.  Gradient-enhancement and gradient
#' predictions are offered following Booth (2025).  Vecchia approximation for faster 
#' computation is implemented following Sauer, Cooper, and Gramacy 
#' (2023).  Optional monotonic warpings are implemented following 
#' Barnett et al. (2025).  Downstream tasks include sequential design 
#' through active learning Cohn/integrated mean squared error (ALC/IMSE; Sauer, 
#' Gramacy, and Higdon, 2023), optimization through expected improvement 
#' (EI; Gramacy, Sauer, and Wycoff, 2022), and contour 
#' location through entropy (Booth, Renganathan, and Gramacy, 
#' 2025).  Models extend up to three layers deep; a one 
#' layer model is equivalent to typical Gaussian process regression.  
#' Incorporates OpenMP and SNOW parallelization and utilizes C/C++ under 
#' the hood.
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
#'   variance over a set of input locations (optionally calculates EI or entropy)
#'   \item \code{\link[deepgp]{plot}}: produces trace plots, hidden layer 
#'   plots, and posterior predictive plots
#'   \item \code{\link[deepgp]{ALC}}: calculates active learning Cohn over 
#'   set of input locations using reference grid
#'   \item \code{\link[deepgp]{IMSE}}: calculates integrated mean-squared error
#'    over set of input locations
#' }
#' 
#' @references 
#' Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
#'      *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
#'      \url{http://hdl.handle.net/10919/114845}
#'      \cr\cr
#' Booth, A. S. (2025). Deep Gaussian processes with gradients. arXiv:2512.18066
#'      \cr\cr
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'      Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'      \cr\cr
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
#'      processes for computer experiments. 
#'      *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
#'      \cr\cr
#' Gramacy, R. B., Sauer, A. & Wycoff, N. (2022). Triangulation candidates for Bayesian 
#'     optimization.  *Advances in Neural Information Processing Systems (NeurIPS), 35,* 
#'     35933-35945.  arXiv:2112.07457
#'     \cr\cr
#' Booth, A., Renganathan, S. A. & Gramacy, R. B. (2025). Contour location for 
#'     reliability in airfoil simulation experiments using deep Gaussian 
#'     processes. *Annals of Applied Statistics, 19*(1), 191-211. arXiv:2308.04420
#'	   \cr\cr
#' Barnett, S., Beesley, L. J., Booth, A. S., Gramacy, R. B., & Osthus D. (2025). 
#'     Monotonic warpings for additive and deep Gaussian processes. 
#'     *Statistics and Computing, 35*(3), 65. arXiv:2408.01540
#'     
#' @examples 
#' # See vignette, ?fit_one_layer, ?fit_two_layer, ?fit_three_layer, 
#' # ?ALC, or ?IMSE for examples
#' # Many more examples including real-world computer experiments are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' 
"_PACKAGE"