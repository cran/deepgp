% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R, R/predict_vecchia.R
\name{predict}
\alias{predict}
\alias{predict.gp}
\alias{predict.dgp2}
\alias{predict.dgp3}
\alias{predict.gpvec}
\alias{predict.dgp2vec}
\alias{predict.dgp3vec}
\title{Predict posterior mean and variance/covariance}
\usage{
\method{predict}{gp}(
  object,
  x_new,
  lite = TRUE,
  return_all = FALSE,
  EI = FALSE,
  entropy_limit = NULL,
  cores = 1,
  ...
)

\method{predict}{dgp2}(
  object,
  x_new,
  lite = TRUE,
  store_latent = FALSE,
  mean_map = TRUE,
  return_all = FALSE,
  EI = FALSE,
  entropy_limit = NULL,
  cores = 1,
  ...
)

\method{predict}{dgp3}(
  object,
  x_new,
  lite = TRUE,
  store_latent = FALSE,
  mean_map = TRUE,
  return_all = FALSE,
  EI = FALSE,
  entropy_limit = NULL,
  cores = 1,
  ...
)

\method{predict}{gpvec}(
  object,
  x_new,
  m = object$m,
  ordering_new = NULL,
  lite = TRUE,
  return_all = FALSE,
  EI = FALSE,
  entropy_limit = NULL,
  cores = 1,
  ...
)

\method{predict}{dgp2vec}(
  object,
  x_new,
  m = object$m,
  ordering_new = NULL,
  lite = TRUE,
  store_latent = FALSE,
  mean_map = TRUE,
  return_all = FALSE,
  EI = FALSE,
  entropy_limit = NULL,
  cores = 1,
  ...
)

\method{predict}{dgp3vec}(
  object,
  x_new,
  m = object$m,
  ordering_new = NULL,
  lite = TRUE,
  store_latent = FALSE,
  mean_map = TRUE,
  return_all = FALSE,
  EI = FALSE,
  entropy_limit = NULL,
  cores = 1,
  ...
)
}
\arguments{
\item{object}{object from \code{fit_one_layer}, \code{fit_two_layer}, or 
\code{fit_three_layer} with burn-in already removed}

\item{x_new}{matrix of predictive input locations}

\item{lite}{logical indicating whether to calculate only point-wise 
variances (\code{lite = TRUE}) or full covariance 
(\code{lite = FALSE})}

\item{return_all}{logical indicating whether to return mean and point-wise
variance prediction for ALL samples (only available for \code{lite = TRUE})}

\item{EI}{logical indicating whether to calculate expected improvement 
(for minimizing the response)}

\item{entropy_limit}{optional limit state for entropy calculations (separating
passes and failures), default value of \code{NULL} bypasses entropy
calculations}

\item{cores}{number of cores to utilize in parallel}

\item{...}{N/A}

\item{store_latent}{logical indicating whether to store and return mapped 
values of latent layers (two or three layer models only)}

\item{mean_map}{logical indicating whether to map hidden layers using 
conditional mean (\code{mean_map = TRUE}) or using a random sample
from the full MVN distribution (two or three layer models only),
\code{mean_map = FALSE} is not yet implemented for fits with 
\code{vecchia = TRUE}}

\item{m}{size of Vecchia conditioning sets (only for fits with 
\code{vecchia = TRUE}), defaults to the \code{m} used for MCMC}

\item{ordering_new}{optional ordering for Vecchia approximation, must correspond
to rows of \code{x_new}, defaults to random, is applied to all layers
in deeper models}
}
\value{
object of the same class with the following additional elements:
\itemize{
  \item \code{x_new}: copy of predictive input locations
  \item \code{mean}: predicted posterior mean, indices correspond to 
        \code{x_new} locations
  \item \code{s2}: predicted point-wise variances, indices correspond to 
        \code{x_new} locations (only returned when \code{lite = TRUE})
  \item \code{mean_all}: predicted posterior mean for each sample (column
        indices), only returned when \code{return_all = TRUE}
  \item \code{s2_all} predicted point-wise variances for each sample (column
        indices), only returned when \code{return-all = TRUE}
  \item \code{Sigma}: predicted posterior covariance, indices correspond to 
        \code{x_new} locations (only returned when \code{lite = FALSE})
  \item \code{EI}: vector of expected improvement values, indices correspond 
        to \code{x_new} locations (only returned when \code{EI = TRUE})
  \item \code{entropy}: vector of entropy values, indices correspond to 
        \code{x_new} locations (only returned when \code{entropy_limit} is
        numeric)
  \item \code{w_new}: list of hidden layer mappings (only returned when 
        \code{store_latent = TRUE}), list index corresponds to iteration and 
        row index corresponds to \code{x_new} location (two or three layer 
        models only)
  \item \code{z_new}: list of hidden layer mappings (only returned when 
        \code{store_latent = TRUE}), list index corresponds to iteration and 
        row index corresponds to \code{x_new} location (three layer models only) 
}
Computation time is added to the computation time of the existing object.
}
\description{
Acts on a \code{gp}, \code{dgp2}, or \code{dgp3} object.
    Calculates posterior mean and variance/covariance over specified input 
    locations.  Optionally calculates expected improvement (EI) or entropy 
    over candidate inputs.  Optionally utilizes SNOW parallelization.
}
\details{
All iterations in the object are used for prediction, so samples 
    should be burned-in.  Thinning the samples using \code{trim} will speed 
    up computation.  Posterior moments are calculated using conditional 
    expectation and variance.  As a default, only point-wise variance is 
    calculated.  Full covariance may be calculated using \code{lite = FALSE}. 
    
    Expected improvement is calculated with the goal of minimizing the 
    response.  See Chapter 7 of Gramacy (2020) for details.  Entropy is 
    calculated based on two classes separated by the specified limit.  
    See Sauer (2023, Chapter 3) for details.
    
    SNOW parallelization reduces computation time but requires 
    more memory storage.
}
\examples{
# See ?fit_one_layer, ?fit_two_layer, or ?fit_three_layer
# for examples

}
\references{
Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
     *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
     \cr\cr
Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
     Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
     \cr\cr
Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
     processes for computer experiments. 
     *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
     \cr\cr
Barnett, S., Beesley, L. J., Booth, A. S., Gramacy, R. B., & Osthus D. (2024). Monotonic 
     warpings for additive and deep Gaussian processes. *In Review.* arXiv:2408.01540
}
