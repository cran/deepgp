# deepgp Package

Maintainer: Annie Sauer <anniees@vt.edu>

Performs model fitting and sequential design for deep Gaussian processes 
following Sauer, Gramacy, and Higdon (2020).  Models extend up to three layers 
deep; a one layer model is equivalent to typical Gaussian process regression.  
Sequential design criteria include integrated mean-squared error (IMSE), active 
learning Cohn (ALC), and expected improvement (EI).  Covariance structure is 
based on inverse exponentiated squared euclidean distance.  Applicable to 
noisy and deterministic functions.  Incorporates SNOW parallelization and 
utilizes C under the hood.

View `deepgp-package` help file for more information.

## Reference

Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
Process Surrogates." arXiv:2012.08015.


    
