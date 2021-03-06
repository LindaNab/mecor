% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mecor.R
\name{mecor}
\alias{mecor}
\title{mecor: a Measurement Error Correction Package}
\usage{
mecor(formula, data, method = "standard", B = 0)
}
\arguments{
\item{formula}{an object of class \link[stats]{formula} (or one that is
coerced to that class): a symbolic description of the regression model
containing a \link[mecor]{MeasError}, \link[mecor]{MeasErrorExt} or
\link[mecor]{MeasErrorRandom} object in one of the covariates or the outcome.}

\item{data}{a data.frame, list or environment (or object coercible by
as.data.frame to a data frame) containing the variables in the model
specified in \code{formula}.}

\item{method}{a character string indicating the method used to correct for
the measurement error, either "standard" (regression calibration for
covariate measurement error and method of moments for outcome measurement
error), "efficient" (efficient regression calibration for covariate
measurement error and efficient method of moments for outcome measurement
error), "valregcal" (validation regression calibration) or "mle" (maximum
likelihood estimation). Defaults to "standard".}

\item{B}{number of bootstrap samples, defaults to 0.}
}
\value{
\code{mecor} returns an object of \link[base]{class} "mecor".

An object of class \code{mecor} is a list containing the following components:

\item{corfit}{a list containing the corrected fit, including the coefficients
of the corrected fit (\code{coef}) and the variance--covariance matrix of the
coefficients of the corrected fit obtained by the delta method (\code{vcov}),
and more depending on the method used.}
\item{uncorfit}{an \link[stats]{lm.fit} object of the uncorrected fit.}
}
\description{
mecor provides correction methods for measurement error in a continuous
covariate or outcome in linear regression models with a continuous outcome
}
\examples{
## measurement error in a covariate:
# internal covariate-validation study
data(icvs)
out <-
mecor(Y ~ MeasError(X_star, reference = X) + Z,
      data = icvs,
      method = "standard",
      B = 999)
# replicates study
data(rs)
mecor(Y ~ MeasError(X_star_1, replicate = cbind(X_star_2, X_star_3)) + Z1 + Z2,
      data = rs,
      method = "mle")
# covariate-calibration study
data(ccs)
mecor(Y ~ MeasError(X_star, replicate = cbind(X_star_1, X_star_2)) + Z,
      data = ccs,
      method = "efficient")
# external covariate-validation study
data(ecvs)
calmod_fit <- lm(X ~ X_star + Z, data = ecvs)
data(icvs) # suppose reference X is not available
mecor(Y ~ MeasErrorExt(X_star, model = calmod_fit) + Z,
      data = icvs)
# sensitivity analyses
data(icvs) # suppose reference X is not available
# guesstimate the coefficients of the calibration model:
mecor(Y ~ MeasErrorExt(X_star, model = list(coef = c(0, 0.9, 0.2))) + Z,
      data = icvs)
# assume random measurement error in X_star of magnitude 0.25:
mecor(Y ~ MeasErrorRandom(X_star, variance = 0.25) + Z,
      data = icvs)
data(rs) # suppose replicates X_star_2 and X_star_2 are not available
mecor(Y ~ MeasErrorRandom(X_star_1, variance = 0.25) + Z1 + Z2,
      data = rs)

## measurement error in the outcome:
# internal outcome-validation study
data(iovs)
mecor(MeasError(Y_star, reference = Y) ~ X + Z,
      data = iovs,
      method = "standard")
# external outcome-validation study
data(eovs)
memod_fit <- lm(Y_star ~ Y, data = eovs)
data(iovs) # suppose reference Y is not available
mecor(MeasErrorExt(Y_star, model = memod_fit) ~ X + Z,
      data = iovs,
      method = "standard")
# sensitivity analyses
data(iovs) # suppose reference Y is not available
# guesstimate the coefficients of the measurement error model:
mecor(MeasErrorExt(Y_star, model = list(coef = c(0, 0.5))) ~ X + Z,
      data = iovs,
      method = "standard")

## differential measurement error in the outcome:
# internal outcome-validation study
data(iovs_diff)
mecor(MeasError(Y_star, reference = Y, differential = X) ~ X,
      data = iovs_diff,
      method = "standard")
# sensitivity analysis
data(iovs_diff) # suppose reference Y is not available
# guesstimate the coefficients of the measurement error model:
mecor(MeasErrorExt(Y_star, model = list(coef = c(0, 0.5, 1, 1))) ~ X,
      data = iovs_diff,
      method = "standard")
}
\references{
L. Nab, R.H.H. Groenwold, P.M.J. Welsing, and  M. van Smeden.
Measurement error in continuous endpoints in randomised trials: problems and
solutions

L. Nab, M. van Smeden, R.H. Keogh, and R.H.H. Groenwold.
mecor: an R package for measurement error correction in linear models with
continuous outcomes
}
\author{
Linda Nab, \email{l.nab@lumc.nl}
}
