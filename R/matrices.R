#' @title Create a Measurement Error Matrix Object
#'
#' @description
#' This function creates a measurement error matrix object, usually used as the
#' dependent variable in a regression if one wants to correct for the
#' measurement error in that variable, by defining the coefficients of the
#' measurement error matrix (dependent variable)
#'
#' @param substitute a vector containing the error-prone measure
#' @param mefit a fitted linear model object from \link[stats]{lm} or a named list
#' containing 'coef': the coefficients of the measurement error
#' model (c(theta0, theta1) if there is non-differential measurement error
#' and c(theta00, theta01, theta10, theta11) if there is differential
#' measurement error) and a non-mandatory 'vcov': the variance--covariance matrix
#' of the coefficients
#' @param n number of independent variables in the regression model. Defaults to
#' one. If MeasErrorMatrix is used inside \link[mecor]{mecor}, n is neglected
#' and is calculated automatically.
#'
#' @return \code{MeasErrorMatrix} returns an object of \link[base]{class}
#' "MeasErrorMatrix".
#'
#' An object of class \code{MeasErrorMatrix} is a list containing the substitute
#' variable and the coefficients of the measurement error model and the
#' variance--covariance matrix of the coefficients of the measurement error
#' model and the measurement error matrix. A \code{MeasErrorMatrix} is typically
#' used in \link[mecor]{mecor} in case of measurement error in the dependent
#' variable. The \code{MeasErrorMatrix} is implemented in \code{mecor} to
#' provide measurement error correction using external data and for conducting
#' sensitivity analyses. An MeasErrorMatrix object has attributes input (the
#' name of the substitute variable) and n (number of dependent variables) and
#' call (the matched call).
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ## measurement error in dependent variable
#' data(ivs_o)
#' measerror_fit <- lm(Y_star ~ Y, data = ivs_o)
#' measerror_coef <- coefficients(measerror_fit)
#' measerror_vcov <- vcov(measerror_fit)
#' measerror_matrix <- with (ivs_o, MeasErrorMatrix(substitute = Y_star,
#'                                                  mefit = measerror_fit,
#'                                                  n = 2))
#' # identical to:
#' measerror_matrix <- with (ivs_o, MeasErrorMatrix(substitute = Y_star,
#'                                                  mefit = list(coef = measerror_coef),
#'                                                  n = 2))
#' measerror_matrix <- with (ivs_o, MeasErrorMatrix(substitute = Y_star,
#'                                                  mefit = list(coef = c(0, 0.5)),
#'                                                  n = 2))
#' @export
MeasErrorMatrix <- function(mefit, ...){
  UseMethod("MeasErrorMatrix", mefit)
}
#' @export
MeasErrorMatrix.list <- function(substitute,
                                 mefit,
                                 n = 1){
  if (missing(substitute))
    stop("'substitute' in the MeasErrorMatrix object is missing")
  if (is.null(mefit$coef))
    stop("'coef' in 'mefit' object is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the MeasErrorMatrix object is not a vector")
  if (!is.vector(mefit$coef))
    stop("'coef' in 'mefit' is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the MeasErrorMatrix object cannot contain missing values")
  if (n < 1)
    stop("'n' should be equal or greater than one")
  if (length(mefit$coef) == 2){
    type <- "dep"
    if (!is.null(mefit$vcov)){
      if (!is.matrix(mefit$vcov)){
        stop("'vcov' should be a matrix")
      } else if (nrow(mefit$vcov) != 2 | ncol(mefit$vcov) != 2)
      stop("'vcov' in 'mefit' should be a 2x2 matrix")
    }
  } else if (length(mefit$coef) == 4){
    type <- "dep_diff"
    if (!is.null(mefit$vcov)){
      if (!is.matrix(mefit$vcov)){
        stop("'vcov' should be a matrix")
      } else if (nrow(mefit$vcov) != 4 | ncol(mefit$vcov) != 4)
        stop("'vcov' in 'mefit' should be a 4x4 matrix")
    }
  } else stop("length of coef should be 2 (non-differential outcome measurement error) or 4 (differential outcome measurement error)")
  matrix <- mecor:::regcal_get_calmod_matrix(coef_calmod = mefit$coef,
                                             n = n + 1,
                                             type = type)
  out <- list(substitute = substitute,
              coef = mefit$coef,
              matrix = matrix)
  if (!is.null(mefit$vcov)) out$vcov <- mefit$vcov
  input <- c(substitute = as.list(match.call())$substitute)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  attr(out, "fit_class") <- "list"
  attr(out, "n") <- n
  class(out) <- c("MeasErrorMatrix", "list")
  out
}
#' @export
MeasErrorMatrix.lm <- function(substitute,
                               mefit,
                               n = 1){
  if (missing(substitute))
    stop("'substitute' in the MeasErrorMatrix object is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the MeasErrorMatrix object is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the MeasErrorMatrix object cannot contain missing values")
  if (n < 1)
    stop("'n' should be equal or greater than one")
  coef <- coefficients(mefit)
  vcov <- vcov(mefit)
  if (length(coef) == 2){
    type <- "dep"
    if (!is.null(vcov)){
      if (!is.matrix(vcov)){
        stop("'vcov' should be a matrix")
      } else if (nrow(vcov) != 2 | ncol(vcov) != 2)
        stop("'vcov' in 'mefit' should be a 2x2 matrix")
    }
  } else if (length(coef) == 4){
    type <- "dep_diff"
    if (!is.null(vcov)){
      if (!is.matrix(vcov)){
        stop("'vcov' should be a matrix")
      } else if (nrow(vcov) != 4 | ncol(vcov) != 4)
        stop("'vcov' in 'mefit' should be a 4x4 matrix")
    }
  } else stop("length of coef should be 2 (non-differential outcome measurement error) or 4 (differential outcome measurement error)")
  matrix <- mecor:::regcal_get_calmod_matrix(coef_calmod = coef,
                                             n = n + 1,
                                             type = type)
  out <- list(substitute = substitute,
              coef = mefit$coef,
              matrix = matrix)
  out$vcov <- vcov
  input <- c(substitute = as.list(match.call())$substitute)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  attr(out, "fit_class") <- "lm"
  attr(out, "n") <- n
  class(out) <- c("MeasErrorMatrix", "list")
  out
}
#' @title Create a Calibration Model Matrix Object
#'
#' @description
#' This function creates a calibration model matrix object, usually used as the
#' independent variable in a regression if one wants to correct for
#' the measurement error in that variable, by defining the coefficients of the
#' calibration matrix.
#'
#' @param substitute a vector containing the error-prone measure
#' @param calfit a fitted linear model object from \link[stats]{lm} or a
#' named list containing 'coef': the coefficients of the calibration model
#' c(lambda0, lambdaxstar, lambdaz) and a non-mandatory 'vcov': the
#' variance--covariance matrix of the coefficients
#'
#' @return \code{CalModMatrix} returns an object of \link[base]{class}
#' "CalModMatrix".
#'
#' An object of class \code{CalModMatrix} is a list containing the substitute
#' variable and the coefficients of the calibration model and the
#' variance--covariance matrix of the coefficients of calibration model and the
#' calibration model matrix. A \code{CalModMatrix} is typically used in
#' \link[mecor]{mecor} in case of measurement error in the independent variable.
#' The CalModMatrix is implemented in \code{mecor} to provide measurement error
#' correction using external data and for conducting sensitivity analyses. A
#' \code{CalModMatrix} object has attributes input (the name of the substitute
#' variable) and call (the matched call).
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ## measurement error in dependent variable
#' data(ivs)
#' calmod_fit <- lm(X ~ X_star + Z, data = ivs)
#' calmod_coef <- coefficients(calmod_fit)
#' calmod_vcov <- vcov(calmod_fit)
#' calmod_matrix <- with (ivs, CalModMatrix(substitute = X_star,
#'                                          calfit = calmod_fit))
#' # identical to
#' calmod_matrix <- with (ivs, CalModMatrix(substitute = X_star,
#'                                          calfit = list(coef = calmod_coef,
#'                                                         vcov = calmod_vcov)))
#' @export
CalModMatrix <- function(calfit, ...){
  UseMethod("CalModMatrix", calfit)
}
#' @export
CalModMatrix.list <- function(substitute,
                              calfit){
   if (missing(substitute))
    stop("'substitute' in the CalModMatrix object is missing")
  if (is.null(calfit$coef))
    stop("'coef' in 'calfit' object is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the CalModMatrix object is not a vector")
  if (!is.vector(calfit$coef))
    stop("'coef' in 'calfit' is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the CalModMatrix object cannot contain missing values")
  if (length(calfit$coef) == 1)
    stop("length of 'coef' should be 2 or greater than 2")
  if (!is.null(calfit$vcov)){
    if (!is.matrix(calfit$vcov)){
      stop("'vcov' in 'calfit' should be a matrix")
    } else if (length(calfit$coef) != nrow(calfit$vcov) &
               length(calfit$coef) != ncol(calfit$vcov))
      stop("dimension of 'vcov' in 'calfit' should be equal to length of 'coef' in 'calfit'")
  }
  n <- length(calfit$coef)
  coef <- mecor:::change_order_coefs(calfit$coef)
  matrix <- mecor:::regcal_get_calmod_matrix(coef_calmod = coef,
                                             n = n,
                                             type = "indep")
  out <- list(substitute = substitute,
              coef = calfit$coef,
              matrix = matrix)
  if (!is.null(calfit$vcov)) out$vcov <- calfit$vcov
  input <- c(substitute = as.list(match.call())$substitute)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  attr(out, "fit_class") <- "list"
  attr(out, "n") <- n
  class(out) <- c("CalModMatrix", "list")
  out
}
#' @export
CalModMatrix.lm <- function(substitute,
                            calfit){
  if (missing(substitute))
    stop("'substitute' in the CalModMatrix object is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the CalModMatrix object is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the CalModMatrix object cannot contain missing values")
  coef <- mecor:::change_order_coefs(coefficients(calfit))
  vcov <- vcov(calfit)
  n <- length(coef)
  matrix <- mecor:::regcal_get_calmod_matrix(coef_calmod = coef,
                                             n = n,
                                             type = "indep")
  out <- list(substitute = substitute,
              coef = calfit$coef,
              matrix = matrix)
  out$vcov <- vcov
  input <- c(substitute = as.list(match.call())$substitute)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  attr(out, "fit_class") <- "lm"
  attr(out, "n") <- n
  class(out) <- c("CalModMatrix", "list")
  out
}
