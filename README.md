# deepgp Package

Maintainer: Annie S. Booth <annie_booth@ncsu.edu>

Performs Bayesian posterior inference for deep Gaussian processes following Sauer, Gramacy, and Higdon (2023, <arXiv:2012.08015>).  See Sauer (2023, <http://hdl.handle.net/10919/114845>) for comprehensive methodological details and <https://bitbucket.org/gramacylab/deepgp-ex/> for a variety of coding examples. Models are trained through MCMC including elliptical slice sampling of latent Gaussian layers and Metropolis-Hastings sampling of kernel hyperparameters.  Vecchia-approximation for faster computation is implemented following Sauer, Cooper, and Gramacy (2022, <arXiv:2204.02904>).  Downstream tasks include sequential design through active learning Cohn/integrated mean squared error (ALC/IMSE; Sauer, Gramacy, and Higdon, 2023), optimization through expected improvement (EI; Gramacy, Sauer, and Wycoff, 2021 <arXiv:2112.07457>), and contour location through entropy (Sauer, 2023).  Models extend up to three layers deep; a one layer model is equivalent to typical Gaussian process regression.  Incorporates OpenMP and SNOW parallelization and utilizes C/C++ under the hood.

Run `help("deepgp-package")` or `help(package = "deepgp")` for more information.

## References

Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*

Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015

Sauer, A., Cooper, A., & Gramacy, R. B. (2022). Vecchia-approximated deep Gaussian processes for computer experiments. *Journal of Computational and Graphical Statistics,* 1-14.  arXiv:2204.02904

Gramacy, R. B., Sauer, A. & Wycoff, N. (2022). Triangulation candidates for Bayesian optimization.  *Advances in Neural Information Processing Systems (NeurIPS), 35,* 35933-35945.  arXiv:2112.07457

## Version History

What's new in version 1.1.1?

* Entropy calculations for contour locating sequential designs are offered through the specification of an  `entropy_limit` in any of the `predict` functions.
* In posterior predictions, there is now an option to return point-wise mean and variance estimates for all MCMC samples through the specification of `return_all = TRUE`.
* To save on memory and storage, `predict` functions no longer return `s2_smooth` or `Sigma_smooth`.  If desired, these quantities may be calculated by subtracting `tau2*g` from the diagonal.
* The `vecchia = TRUE` option may now utilize either the Matern (`cov = "matern"`)  or squared exponential kernel (`cov = "exp2"`").
* Performance improvements for `cores = 1` in `predict`, `ALC`, and `IMSE` functions (helps to avoid a SNOW conflict when running multiple instances on the same machine).
* Fit functions now return the outer log likelihood value along with MCMC samples.  Used in trace plots to assess burn-in.
* In `fit_two_layer`, the intermediate latent layer may now have either a prior mean of zero (default) or a prior mean equal to `x` (`pmx = TRUE`).  If `pmx` is set to a constant, this will be the scale parameter on the inner Gaussian layer.

What's new in version 1.1.0?

* Package vignette
* Option to specify `sep = TRUE` in `fit_one_layer` to fit a GP with separable/anisotropic lengthscales.
* Default cores in predict are now 1 (this avoids a conflict when running multiple sessions simultaneously on a single machine).

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


