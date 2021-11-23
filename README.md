# deepgp Package

Maintainer: Annie Sauer <anniees@vt.edu>

Performs model fitting and sequential design for deep Gaussian
processes following Sauer, Gramacy, and Higdon (2020) <arXiv:2012.08015>.  
Models extend up to three layers deep; a one layer model is equivalent 
to typical Gaussian process regression.  Covariance kernel options are 
Matern (default) and squared exponential.  Sequential design criteria 
include integrated mean-squared error (IMSE), active learning Cohn (ALC), 
and expected improvement (EI).  Applicable to both noisy and 
deterministic functions.  Incorporates SNOW parallelization and 
utilizes C and C++ under the hood.

Run `help("deepgp-package")` or `help(package = "deepgp")` for more information.

What's new in version 0.3.0?

+ The Matern kernel is now the default covariance.  The smoothness parameter 
is user-adjustable but must be either `v = 0.5`, `v = 1.5`, or `v = 2.5` (default).
The squared exponential kernel is still required for use with `ALC` and `IMSE` 
(set `cov = "exp2"` in fit functions).
+ Expected improvement (EI) may now be computed concurrently with 
predictions.  Set `EI = TRUE` inside `predict` calls.  EI calculations are
nugget-free and are for *minimizing* the response (negate `y` if maximization
is desired).
+ To save memory, hidden layer mappings used in predictions are no longer
stored and returned by default.  To store them, set `store_latent = TRUE` inside
`predict`.

## Reference

Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
Process Surrogates." *Technometrics, to appear;* arXiv:2012.08015.


    
