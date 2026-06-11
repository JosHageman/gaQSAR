# gaQSAR 1.2.3

## Improvements

* Improved Williams plot diagnostics and documentation for double cross-validation.
* Improved handling of single `gaQSAR` objects versus lists of objects.
* Updated permutation-test p-values with a plus-one correction.
* Improved double cross-validation predictor selection frequencies.
* Updated examples, vignettes and documentation.

# gaQSAR 1.2.1

## Improvements

* Fix: permutation-test validation for `gaQSAR_dcv` now reconstructs outer predictions from saved fold models (removes spurious validation warnings).
* Improved reproducibility: `bestSeed` tracking and outer fold `seed` are stored and used consistently.
* Clarified documentation and messages for validation and CV bookkeeping.

# gaQSAR 1.2.0

## New Features

* Added nested cross-validation with outer CV (LOO/k-fold), inner LOOCV fitness, robust error handling, AD thresholds, verbose output option, and predictor stability metrics.

# gaQSAR 1.0.1

## New Features

* Added VIF (Variable Inflation Factor) scores 
* Added counting of Williams plot objects with high residuals and leverages for both training and test sets
* GA-selected models now report adjusted R-squared alongside training R-squared

# gaQSAR 1.0.0

* Initial release
