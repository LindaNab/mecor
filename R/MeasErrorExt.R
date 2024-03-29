#' @title Create an External Measurement Error Object
#'
#' @description
#' This function creates an external measurement error object, usually used as
#' a covariate or the outcome in the \code{formula} argument of
#' \link[mecor]{mecor} if one wants to correct for the measurement error in that
#' variable using external data or externally estimated coefficients of the
#' calibration model (covariate-measurement error) or measurement error model
#' (outcome-measurement error)
#'
#' @param substitute a vector containing the error-prone measure
#' @param model a fitted linear model of class \link[stats]{lm} or a named
#' \link[base]{list}. The \link[base]{list} contains a vector named \code{coef}:
#' the coefficients of the calibration model or measurement error model and an
#' optional matrix named \code{vcov}: the variance--covariance matrix of the
#' coefficients
#'
#' @return \code{MeasErrorExt} returns an object of \link[base]{class}
#' "MeasErrorExt".
#'
#' An object of class \code{MeasErrorExt} is a list containing the substitute
#' variable and the fitted calibration model or measurement error model and has
#' attributes input (the name of the substitute variable) and call (the matched
#' call).
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ## measurement error in a outcome:
#' # external outcome-validation study
#' data(haemoglobin_ext)
#' # calibration model
#' calmod_fit <- lm(capillary ~ venous, data = haemoglobin)
#' # the external covariate-validation study can be used to correct for the
#' # measurement error in X_star in the dataset 'icvs', using the fitted
#' # calibration model
#' data(haemoglobin)
#' with (haemoglobin, MeasErrorExt(substitute = capillary,
#'                                 model = calmod_fit))
#' # identical to:
#' calmod_coef <- coefficients(calmod_fit)
#' calmod_vcov <- vcov(calmod_fit)
#' with (haemoglobin, MeasErrorExt(substitute = capillary,
#'                                 model = list(coef = calmod_coef,
#'                                              vcov = calmod_vcov)))
#' # when no external data is available, guesstimations of the coefficients of
#' # the calibration model can be used instead:
#' with (haemoglobin, MeasErrorExt(substitute = capillary,
#'                                 model = list(coef = c('(Intercept)' = -7,
#'                                                       'venous' = 1.1))))
#' @export
MeasErrorExt <- function(substitute,
                         model) {
  UseMethod("MeasErrorExt", model)
}
#' @export
MeasErrorExt.lm <- function(substitute,
                            model) {
  if (missing(substitute))
    stop("'substitute' in the MeasErrorMatrix object is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the MeasErrorMatrix object is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the MeasErrorMatrix object cannot contain missing values")
  out <- list(substitute = substitute,
              model = model)
  input <- c(substitute = as.list(match.call())$substitute)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  class(out) <- c("MeasErrorExt", "list")
  out
}
#' @export
MeasErrorExt.list <- function(substitute,
                              model) {
  if (missing(substitute))
    stop("'substitute' in the MeasErrorExt object is missing")
  if (is.null(model$coef))
    stop("'coef' in the list 'model' is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the MeasErrorExt object is not a vector")
  if (!is.vector(model$coef))
    stop("'coef' in the list 'model' is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the MeasErrorExt object cannot contain missing values")
  if (length(model$coef) == 1)
    stop("length of 'coef' should be 2 or greater than 2")
  if (!is.null(model$vcov)) {
    if (!is.matrix(model$vcov)) {
      stop("'vcov' in the list 'model' should be a matrix")
    } else if (length(model$coef) != nrow(model$vcov) &
               length(model$coef) != ncol(model$vcov))
      stop(
        "dimensions of 'vcov' in the list 'model' should be equal to length of 'coef' in 'calfit'"
      )
  }
  out <- list(substitute = substitute,
              model = model)
  input <- c(substitute = as.list(match.call())$substitute)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  class(out) <- c("MeasErrorExt", "list")
  out
}
