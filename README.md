# deepgp Package

Performs model fitting and sequential design for deep Gaussian
processes using MCMC and elliptical slice sampling.  Models extend up to 
three layers deep; a one layer model is equivalent to typical Gaussian 
process regression.  Sequential design criteria include integrated mean 
square prediction error (IMSPE), active learning Cohn (ALC), and expected 
improvement (EI).  Covariance structure is based on inverse exponentiated 
squared euclidean distance.  Applicable to noisy and deterministic functions.  
Incorporates SNOW parallelization and utilizes C under the hood.  Manuscript 
forthcoming.

View `deepgp-package` help file for more information.

## CRAN submission notes:

This is a re-submission.  The following changes were made:
* Added references to description field of DESCRIPTION file
* Added toy examples that run in less than 5 seconds to ?fit_one_layer, ?fit_two_layer, and ?fit_three_layer documentation
* Implemented `on.exit` to restore par options after plotting

### Test environments
* local OS X, R 3.6.3
* ubuntu 18.04.4, R 3.6.3
* win-builder (devel)

### R CMD check results
There were no ERRORs or WARNINGs or NOTEs.

    
