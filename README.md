# deepgp Package

Maintainer: Annie Sauer <anniees@vt.edu>

Performs posterior inference for deep Gaussian processes following Sauer, Gramacy, and Higdon (2020) <arXiv:2012.08015>.  Models are trained through MCMC including elliptical slice sampling of latent Gaussian layers and Metropolis-Hastings sampling of kernel hyperparameters.  Vecchia-approximation for faster computation is implemented following Sauer, Cooper, and Gramacy (2022) <arXiv:2204.02904>.  Downstream tasks include sequential design through active learning Cohn/integrated mean squared error (ALC/IMSE; Sauer, Gramacy, and Higdon, 2020) and optimization through expected improvement (EI; Gramacy, Sauer, and Wycoff, 2021 <arXiv:2112.07457>).  Models extend up to three layers deep; a one layer model is equivalent to typical Gaussian process regression.  Covariance kernel options are matern (default) and squared exponential.  Applicable to both noisy and deterministic functions.  Incorporates SNOW parallelization; utilizes C and C++ with OpenMP parallelization under the hood.

Run `help("deepgp-package")` or `help(package = "deepgp")` for more information.

## References

Sauer, A., Gramacy, R.B., & Higdon, D. (2020). Active learning for deep Gaussian process surrogates. *Technometrics*, (just-accepted), 1-39.

Sauer, A., Cooper, A., & Gramacy, R. B. (2022). Vecchia-approximated deep Gaussian processes for computer experiments. *pre-print on arXiv:2204.02904*

## Version History

What's new in version 1.0.1?

* Minor bug fixes/improvements.
* New warning message when OpenMP parallelization is not utilized for the Vecchia approximation.  This happens when the package is downloaded from CRAN on a Mac.  To set up OpenMP, download package source and compile with gcc/g++ instead of clang.

What's new in version 1.0.0?

* Models may now leverage the Vecchia approximation (through the specification of `vecchia = TRUE` in fit functions) for faster computation.  The speed of this implementation relies on OpenMP parallelization (make sure the `-fopenmp` flag is present with package installation).
* SNOW parallelization now uses less memory/storage.
* `tau2` is now calculated at the time of MCMC, not at the time of prediction.  This avoids some extra calculations.

What's new in version 0.3.0?

* The Matern kernel is now the default covariance. The smoothness parameter is user-adjustable but must be either `v = 0.5`, `v = 1.5`, or `v = 2.5` (default). The squared exponential kernel is still required for use with ALC and IMSE (set `cov = "exp2"` in fit functions).
* Expected improvement (EI) may now be computed concurrently with predictions. Set `EI = TRUE` inside `predict` calls. EI calculations are nugget-free and are for *minimizing* the response (negate `y` if maximization is desired).
* To save memory, hidden layer mappings used in predictions are no longer stored and returned by default. To store them, set `store_latent = TRUE` inside predict.


