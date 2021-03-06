% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeasErrorExt.R
\name{MeasErrorExt}
\alias{MeasErrorExt}
\title{Create an External Measurement Error Object}
\usage{
MeasErrorExt(substitute, model)
}
\arguments{
\item{substitute}{a vector containing the error-prone measure}

\item{model}{a fitted linear model of class \link[stats]{lm} or a named
\link[base]{list}. The \link[base]{list} contains a vector named \code{coef}:
the coefficients of the calibration model or measurement error model and an
optional matrix named \code{vcov}: the variance--covariance matrix of the
coefficients}
}
\value{
\code{MeasErrorExt} returns an object of \link[base]{class}
"MeasErrorExt".

An object of class \code{MeasErrorExt} is a list containing the substitute
variable and the fitted calibration model or measurement error model and has
attributes input (the name of the substitute variable) and call (the matched
call).
}
\description{
This function creates an external measurement error object, usually used as
a covariate or the outcome in the \code{formula} argument of
\link[mecor]{mecor} if one wants to correct for the measurement error in that
variable using external data or externally estimated coefficients of the
calibration model (covariate-measurement error) or measurement error model
(outcome-measurement error)
}
\examples{
## measurement error in a covariate:
# external covariate-validation study
data(ecvs)
# calibration model
calmod_fit <- lm(X ~ X_star + Z, data = ecvs)
# the external covariate-validation study can be used to correct for the
# measurement error in X_star in the dataset 'icvs', using the fitted
# calibration model
data(icvs)
with (icvs, MeasErrorExt(substitute = X_star,
                         model = calmod_fit))
# identical to:
calmod_coef <- coefficients(calmod_fit)
calmod_vcov <- vcov(calmod_fit)
with (icvs, MeasErrorExt(substitute = X_star,
                         model = list(coef = calmod_coef,
                                     vcov = calmod_vcov)))
# when no external data is available, guesstimations of the coefficients of
# the calibration model can be used instead:
with (icvs, MeasErrorExt(substitute = X_star,
                         model = list(coef = c('(Intercept)' = 0,
                                               'X_star' = 0.8,
                                               'Z' = 0.2))))
## measurement error in the outcome:
# external outcome-validation study
data(eovs)
memod_fit <- lm(Y_star ~ Y, data = eovs)
# the external outcome-validation study can be used to correct for the
# measurement error in Y_star in the dataset 'iovs', using the fitted
# measurement error model
with (iovs, MeasErrorExt(substitute = Y_star,
                         model = memod_fit))
# identical to:
memod_coef <- coefficients(memod_fit)
memod_vcov <- vcov(memod_fit)
with (iovs, MeasErrorExt(substitute = Y_star,
                         model = list(coef = memod_coef,
                                      vcov = memod_vcov)))
# when no external data is available, guesstimations of the coefficients of
# the measurement error model can be used instead:
with (iovs, MeasErrorExt(substitute = Y_star,
                         model = list(coef = c('(Intercept)' = 0,
                                               'Y' = 0.5))))
}
\author{
Linda Nab, \email{l.nab@lumc.nl}
}
