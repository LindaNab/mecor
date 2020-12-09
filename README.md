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
# load the internal covariate validation study
data("icvs", package = "mecor")
# correct the biased exposure-outcome association
mecor(Y ~ MeasError(X_star, reference = X) + Z, data = icvs, method = "standard")
```

References
==========

Key reference
-------------

-   Nab L, van Smeden M, Keogh RH, Groenwold RHH. mecor: an R package
    for measurement error correction in linear models with a continuous
    outcome.

References to methods implemented in the package
------------------------------------------------

-   Bartlett JW, Stavola DBL, Frost C. Linear mixed models for
    replication data to efficiently allow for covariate measurement
    error. Statistics in Medicine. 2009:28(25):3158–3178.
    \[<a href="doi:10.1002/sim.3713" class="uri">doi:10.1002/sim.3713</a>\]

-   Buonaccorsi JP. Measurement error: Models, methods, and
    applications. 2010. Chapman & Hall/CRC, Boca Raton.

-   Carroll RJ, Ruppert D, Stefanski LA, Crainiceanu CM. Measurement
    error in non-linear models: A modern perspective. 2006, 2nd edition.
    Chapman & Hall/CRC, Boca Raton.

-   Keogh RH, Carroll RJ, Tooze JA, Kirkpatrick SI, Freedman LS.
    Statistical issues related to dietary intake as the response
    variable in intervention trials. Statistics in Medicine.
    2016:35(25):4493–4508.
    \[<a href="doi:10.1002/sim.7011" class="uri">doi:10.1002/sim.7011</a>\]

-   Keogh RH, White IR. A toolkit for measurement error correction, with
    a focus on nutritional epidemiology. Statistics in Medicine
    2014:33(12):2137–2155.
    \[<a href="doi:10.1002/sim.6095" class="uri">doi:10.1002/sim.6095</a>\]

-   Nab L, Groenwold RHH, Welsing PMJ, van Smeden M. Measurement error
    in continuous endpoints in randomised trials: Problems and
    solutions. Statistics in Medicine. 2019:38(27):5182-5196.
    \[<a href="doi:10.1002/sim.8359" class="uri">doi:10.1002/sim.8359</a>\]

-   Rosner B, Spiegelman D, Willett WC. Correction of logistic
    regression relative risk estimates and confidence intervals for
    measurement error: The case of multiple covariates measured with
    error. 1990:132(4):734-745.
    \[<a href="doi:10.1093/oxfordjournals.aje.a115715" class="uri">doi:10.1093/oxfordjournals.aje.a115715</a>\]

-   Rosner B, Spiegelman D, Willett WC. Correction of logistic
    regression relative risk estimates and confidence intervals for
    random within-person measurement error. American Journal of
    Epidemiology. 1992:136(11):1400-1413.
    \[<a href="doi:10.1093/oxfordjournals.aje.a116453" class="uri">doi:10.1093/oxfordjournals.aje.a116453</a>\]

-   Spiegelman D, Carroll RJ, Kipnis V. Efficient regression calibration
    for logistic regression in main study/internal validation study
    designs with an imperfect reference intstrument. Statistics in
    Medicine. 2001:20(1):139-160.
    \[<a href="doi:10.1002/1097-0258(20010115)20:1" class="uri">doi:10.1002/1097-0258(20010115)20:1</a>\<139::AID-SIM644\>3.0.CO;2-K\]
