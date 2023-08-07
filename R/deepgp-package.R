
# Imported Functions ----------------------------------------------------------
#' @importFrom Matrix t solve
#' @importFrom grDevices heat.colors
#' @importFrom graphics image lines matlines par plot points contour
#' @importFrom stats cov dgamma dnorm pnorm qnorm rnorm runif var
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm
#' @importFrom FNN get.knnx

# Package Documentation -------------------------------------------------------
#' @useDynLib deepgp, .registration = TRUE
#' @title Package deepgp
#' @author Annie S. Booth \email{annie_booth@ncsu.edu}
#' @docType package
#' @name deepgp-package
#'
#' @description Performs Bayesian posterior inference for deep Gaussian processes following 
#' Sauer, Gramacy, and Higdon (2023, <arXiv:2012.08015>).  See Sauer (2023, 
#' <http://hdl.handle.net/10919/114845>) for comprehensive methodological details and 
#' <https://bitbucket.org/gramacylab/deepgp-ex/> for a variety of coding examples. 
#' Models are trained through MCMC including elliptical slice sampling of latent 
#' Gaussian layers and Metropolis-Hastings sampling of kernel hyperparameters.  
#' Vecchia-approximation for faster computation is implemented following 
#' Sauer, Cooper, and Gramacy (2022, <arXiv:2204.02904>).  Downstream tasks 
#' include sequential design through active learning Cohn/integrated mean squared 
#' error (ALC/IMSE; Sauer, Gramacy, and Higdon, 2023), optimization through 
#' expected improvement (EI; Gramacy, Sauer, and Wycoff, 2021 <arXiv:2112.07457>), 
#' and contour location through entrop(Sauer, 2023).  Models extend up to three 
#' layers deep; a one layer model is equivalent to typical Gaussian process 
#' regression.  Incorporates OpenMP and SNOW parallelization and utilizes 
#' C/C++ under the hood.
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
#'      \cr\cr
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'      Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'      \cr\cr
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2022). Vecchia-approximated deep Gaussian 
#'      processes for computer experiments. 
#'      *Journal of Computational and Graphical Statistics,* 1-14.  arXiv:2204.02904
#'      \cr\cr
#' Gramacy, R. B., Sauer, A. & Wycoff, N. (2022). Triangulation candidates for Bayesian 
#'     optimization.  *Advances in Neural Information Processing Systems (NeurIPS), 35,* 
#'     35933-35945.  arXiv:2112.07457
#'     
#' @examples 
#' # See "fit_one_layer", "fit_two_layer", "fit_three_layer", 
#' # "ALC", or "IMSE" for examples
#' # Examples of real-world implementations are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' 
NULL