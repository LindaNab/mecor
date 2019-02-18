#' @title Create a Measurement Error Object
#'
#' @description
#' This function creates a measurement error object, usually used in a regression where one
#' wants to correct for the measurement error in the measured variable using a
#' reference measure.
#'
#' @param test a vector containing error-prone measurements or a matrix
#' with replicate error-prone measurements
#' @param reference a vector containing reference measurements
#'
#' @return \code{MeasError} returns an object of \link[base]{class} "MeasError".
#'
#' An object of class \code{MeasError} is a data.frame and additionally has
#' attributes input (the name of the test and reference variables) and call (the matched call).
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ##measurement error in endpoint
#' X <- c(rep(0, 50), rep(1, 50))
#' Y <- X + rnorm(100, 0, 1)
#' Vcme <- Y + rnorm(100, 0, 3) #classical measurement error (cme)
#' Vsme <- 1 + 2*Y + rnorm(100, 0, 3) #systematic measurement error (sme)
#' Vdme <- 2 + 2*X + 3*Y + 2*X*Y + rnorm(10, 0, 3*(1-X) + 2*X) #systematic differential measurement error (dme)
#'
#' ##measurement error in exposure
#' X <- rnorm(100, 0, 1)
#' Y <- 1.5*X + rnorm(100, 0, 0.5)
#' W1 <- X + rnorm(100, 0, 0.5)
#' W2 <- X + rnorm(100, 0, 0.5)
#'
#' cme <- MeasError(Vcme, Y)
#' sme <- MeasError(Vsme, Y)
#' dme <- MeasError(Vdme, Y)
#' cme2 <- (MeasError(W1, X))
#' cme3 <- MeasError(cbind(W1, W2), NA)
#'
#' @export
MeasError <- function(test,
                      reference){
  if(missing(test)) stop("'test' is missing in the MeasError object")
  if(!is.vector(test) & !is.matrix(test)) stop("'test' is not a vector or matrix")
  if(missing(reference)) reference = NA
  if(all(is.na(reference)) & !is.matrix(test))
    stop("if there is no reference, multiple test measurements need to be provided")
  if(!all(is.na(reference)) & !is.vector(reference)) stop("'reference' is not a vector")
  out <- data.frame(test = test, reference = reference)
  input <- c(test = as.list(match.call())$test, reference = as.list(match.call())$reference)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  class(out) <- c("MeasError", "data.frame")
  out
}
#' @export
print.MeasError <- function(x){
  cat("\nCall:\n", deparse(attributes(x)$call), "\n", sep = "")
  cat("\nThe error-prone variable", deparse((attr(x, 'input')$test)),
      "is correctly measured by",   deparse((attr(x, 'input')$reference)))
  invisible(x)
}
