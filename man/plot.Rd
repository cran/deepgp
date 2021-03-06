% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot}
\alias{plot}
\alias{plot.gp}
\alias{plot.dgp2}
\alias{plot.dgp3}
\title{Plots object from "\code{deepgp}" package}
\usage{
\method{plot}{gp}(x, trace = TRUE, predict = TRUE, ...)

\method{plot}{dgp2}(x, trace = TRUE, hidden = FALSE, predict = TRUE, ...)

\method{plot}{dgp3}(x, trace = TRUE, hidden = FALSE, predict = TRUE, ...)
}
\arguments{
\item{x}{object of class \code{gp}, \code{dgp2}, or \code{dgp3}}

\item{trace}{logical indicating whether to generate trace plots}

\item{predict}{logical indicating whether to generate posterior predictive 
plot}

\item{...}{N/A}

\item{hidden}{logical indicating whether to generate plots of hidden layers
("\code{dgp2}" or "\code{dgp3}" only)}
}
\description{
Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.  
    Generates trace plots for length scale and nugget hyperparameters.
    Generates plots of hidden layers for one-dimensional inputs.  Generates
    plots of the posterior mean and estimated 95\% prediction intervals for 
    one-dimensional inputs; generates heat maps of the posterior mean and 
    point-wise variance for two-dimensional inputs.
}
\details{
Trace plots are useful in assessing burn-in.  Hidden layer plots 
    are colored on a gradient - red lines represent earlier iterations and 
    yellow lines represent later iterations - to help assess burn-in of the 
    hidden layers.  These plots are meant to help in model fitting and 
    visualization.
}
\examples{
# See "deepgp-package", "fit_one_layer", "fit_two_layer", or 
# "fit_three_layer" for an example

}
