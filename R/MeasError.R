#' @title Create a Measurement Error Object
#'
#' @description
#' This function creates a measurement error object, usually used in a regression where one
#' wants to correct for the measurement error in the measured variable using a
#' reference measure.
#'
#' @param test a vector containing the error-prone measurement
#' @param reference a vector containing the reference measurement
#'
#' @return \code{MeasError} returns an object of \link[base]{class} "MeasError".
#'
#' An object of class \code{mefit} is a formula and additionally has attributes called input
#' (vectors containing the test and reference variable) and call (the matched call).
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ##generation of external validation set
#' Xval <- c(rep(0, 5), rep(1, 5))
#' Yval <- Xval + rnorm(10, 0, 1)
#' Vval_cme <- Yval + rnorm(10, 0, 3) #classical measurement error (cme)
#' Vval_sme <- 1 + 2 * Yval + rnorm(10, 0, 3) #systematic measurement error (sme)
#' Vval_dme <- 2 + 2 * Xval + 3 * Yval + 2 * Xval * Yval + rnorm(10, 0, 3 * (1 - Xval) + 2 * Xval) #systematic differential measurement error (dme)
#'
#' cme <- MeasError(test = Vval_cme, reference = Yval)
#' sme <- MeasError(test = Vval_sme, reference = Yval)
#' dme <- MeasError(test = Vval_dme, reference = Yval)
#' @export
MeasError <- function(test,
                      reference){
  if(missing(test)) stop("'test' is missing in the MeasError object")
  if(!is.vector(test)) stop("'test' is not a vector")
  if(missing(reference)) stop("'reference' is missing in the MeasError object")
  if(!is.vector(reference)) stop("'reference' is not a vector")
  out <- test ~ reference
  input <- c(test = as.list(match.call())$test, reference = as.list(match.call())$reference)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  class(out) <- c("MeasError", "formula")
  out
}
#' @export
print.MeasError <- function(x){
  cat("\nCall:\n", paste(deparse(attributes(x)$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nThe error prone variable", attributes(x)$input$test,
      "is correctly measured by", attributes(x)$input$reference, "\n")
  if(!is.null(attributes(x)$model)){
    cat("The assumed measurement error model is:", attributes(x)$model)}
  invisible(x)
}

update.MeasError <- function(old, new, ...){
  temp <- attributes(old)
  out <- update.formula(old, new, ...)
  attributes(out) <- temp
  class(out) <- c("MeasError", "formula")
  out
}
