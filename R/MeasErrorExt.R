#' @title Create External Measurement Error Object
#'
#' @description
#' This function creates an external measurement error object, usually used as
#' the independent or dependent variable in a regression if one wants to correct
#' for the measurement error in that variable using external data or externally
#' estimated coefficients of the calibration model (independent variable) or
#' measurement error model (depedent variable)
#'
#' @param substitute a vector containing the error-prone measure
#' @param model a fitted linear model of class \link[stats]{lm} or a list
#' containing the coefficients of the calibration model or measurement error
#' model
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
#' ## measurement error in independent variable
#' data(ivs)
#' calmod_fit <- lm(X ~ X_star + Z, data = ivs)
#' calmod_coef <- coefficients(calmod_fit)
#' calmod_vcov <- vcov(calmod_fit)
#' me_ext <- with (ivs, MeasErrorExt(substitute = X_star,
#'                                   model = calmod_fit))
#' # identical to
#' me_ext <- with (ivs, MeasErrorExt(substitute = X_star,
#'                                   model = list(coef = calmod_coef,
#'                                                vcov = calmod_vcov)))
#' ## measurement error in dependent variable
#' data(ivs_o)
#' measerror_fit <- lm(Y_star ~ Y, data = ivs_o)
#' measerror_coef <- coefficients(measerror_fit)
#' measerror_vcov <- vcov(measerror_fit)
#' me_ext <- with (ivs_o, MeasErrorExt(substitute = Y_star,
#'                                     model = measerror_fit))
#' # identical to:
#' me_ext <- with (ivs_o, MeasErrorExt(substitute = Y_star,
#'                                     model = list(coef = measerror_coef))
#' me_ext <- with (ivs_o, MeasErrorExt(substitute = Y_star,
#'                                     model = list(coef = c(0, 0.5))))
#' @export
MeasErrorExt <- function(substitute,
                         model){
  UseMethod("MeasErrorExt", model)
}
#' @export
MeasErrorExt.lm <- function(substitute,
                            model){
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
                              model){
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
  if (!is.null(model$vcov)){
    if (!is.matrix(model$vcov)){
      stop("'vcov' in the list 'model' should be a matrix")
    } else if (length(model$coef) != nrow(model$vcov) &
               length(model$coef) != ncol(model$vcov))
      stop("dimensions of 'vcov' in the list 'model' should be equal to length of 'coef' in 'calfit'")
  }
  out <- list(substitute = substitute,
              model = model)
  input <- c(substitute = as.list(match.call())$substitute)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  class(out) <- c("MeasErrorExt", "list")
  out
}
