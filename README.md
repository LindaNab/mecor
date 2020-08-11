<!-- README.md is generated from README.Rmd. Please edit that file -->

The mecor Package
=================

This package for R implements measurement error correction methods for
measurement error in a continuous covariate or outcome in a linear
regression with a continuous outcome.

Installation
============

It can be installed via

``` r
devtools::install_github("LindaNab/mecor")
```

Quick demo
==========

``` r
library(mecor)
data(icvs)
mecor(Y ~ MeasError(X_star, reference = X) + Z, data = icvs, method = "rc")
```

References
==========

Nab L, van Smeden M, Keogh RH, Groenwold RHH. mecor: an R package for
measurement error correction

Bartlett JW, Stavola DBL, Frost C. Linear mixed models for replication
data to efficiently allow for covariate measurement error. Statistics in
Medicine. 2009:28(25):3158–3178. Carroll RJ, Ruppert D, Stefanski LA,
Crainiceanu CM. Measurement error in non-linear models: A modern
perspective. 2006, 2nd edition. Chapman & Hall/CRC, Boca Raton. Keogh
RH, Carroll RJ, Tooze JA, Kirkpatrick SI, Freedman LS. Statistical
issues related to dietary intake as the response variable in
intervention trials. Statistics in Medicine. 2016:35(25):4493–4508.
Keogh RH, White IR. A toolkit for measurement error correction, with a
focus on nutritional epidemiology. Statistics in Medicine
2014:33(12):2137–2155. Nab L, Groenwold RHH, Welsing PMJ, van Smeden M.
Measurement error in continuous endpoints in randomised trials: Problems
and solutions. Statistics in Medicine. 2019:38(27):5182-5196.
