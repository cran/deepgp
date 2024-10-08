% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deepgp-package.R
\docType{package}
\name{deepgp-package}
\alias{deepgp-package}
\title{Package deepgp}
\description{
Performs Bayesian posterior inference for deep Gaussian 
processes following Sauer, Gramacy, and Higdon (2023).  
See Sauer (2023) for comprehensive 
methodological details and \url{https://bitbucket.org/gramacylab/deepgp-ex/} for 
a variety of coding examples. Models are trained through MCMC including 
elliptical slice sampling of latent Gaussian layers and Metropolis-Hastings 
sampling of kernel hyperparameters.  Vecchia-approximation for faster 
computation is implemented following Sauer, Cooper, and Gramacy 
(2023).  Optional monotonic warpings are implemented following 
Barnett et al. (2024).  Downstream tasks include sequential design 
through active learning Cohn/integrated mean squared error (ALC/IMSE; Sauer, 
Gramacy, and Higdon, 2023), optimization through expected improvement 
(EI; Gramacy, Sauer, and Wycoff, 2022), and contour 
location through entropy (Booth, Renganathan, and Gramacy, 
2024).  Models extend up to three layers deep; a one 
layer model is equivalent to typical Gaussian process regression.  
Incorporates OpenMP and SNOW parallelization and utilizes C/C++ under 
the hood.
}
\section{Important Functions}{

\itemize{
  \item \code{\link[deepgp]{fit_one_layer}}: conducts MCMC sampling of 
  hyperparameters for a one layer GP
  \item \code{\link[deepgp]{fit_two_layer}}: conducts MCMC sampling of 
  hyperparameters and hidden layer for a two layer deep GP
  \item \code{\link[deepgp]{fit_three_layer}}: conducts MCMC sampling of 
  hyperparameters and hidden layers for a three layer deep GP
  \item \code{\link[deepgp]{continue}}: collects additional MCMC samples
  \item \code{\link[deepgp]{trim}}: cuts off burn-in and optionally thins 
  samples
  \item \code{\link[deepgp]{predict}}: calculates posterior mean and 
  variance over a set of input locations (optionally calculates EI or entropy)
  \item \code{\link[deepgp]{plot}}: produces trace plots, hidden layer 
  plots, and posterior predictive plots
  \item \code{\link[deepgp]{ALC}}: calculates active learning Cohn over 
  set of input locations using reference grid
  \item \code{\link[deepgp]{IMSE}}: calculates integrated mean-squared error
   over set of input locations
}
}

\examples{
# See vignette, ?fit_one_layer, ?fit_two_layer, ?fit_three_layer, 
# ?ALC, or ?IMSE for examples
# Many more examples including real-world computer experiments are available at: 
# https://bitbucket.org/gramacylab/deepgp-ex/

}
\references{
Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
     *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
     \url{http://hdl.handle.net/10919/114845}
     \cr\cr
Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
     Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
     \cr\cr
Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
     processes for computer experiments. 
     *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
     \cr\cr
Gramacy, R. B., Sauer, A. & Wycoff, N. (2022). Triangulation candidates for Bayesian 
    optimization.  *Advances in Neural Information Processing Systems (NeurIPS), 35,* 
    35933-35945.  arXiv:2112.07457
    \cr\cr
Booth, A., Renganathan, S. A. & Gramacy, R. B. (2024). Contour location for 
    reliability in airfoil simulation experiments using deep Gaussian 
    processes. *In Review.* arXiv:2308.04420

Barnett, S., Beesley, L. J., Booth, A. S., Gramacy, R. B., & Osthus D. (2024). 
    Monotonic warpings for additive and deep Gaussian processes. *In Review.* arXiv:2408.01540
}
\author{
Annie S. Booth \email{annie_booth@ncsu.edu}
}
