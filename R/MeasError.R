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
#' ##measurement error in exposure
#' nobs <- 1e3
#' Z <- rnorm(nobs, 0, 1)
#' X <- Z + rnorm(nobs, 0, 1)
#' Y <- 0.5 * X + 2 * Z + rnorm(nobs, 0, 1)
#' W <- X + rnorm(nobs, 0, 0.5)
#' X <- ifelse(rbinom(nobs, 0, 0.9) == 1, NA, X)
#' data <- data.frame(Z, X, W, Y)
#' W2 <- X + rnorm(nobs, 0, 0.5)
#' data2 <- data.frame(Z, W, W2, Y)
#'
#' MeasError(W, X)
#' MeasError(cbind(W, W2), NA)
#' @export
MeasError <- function(test,
                      reference){
  if(missing(test)) stop("'test' is missing in the MeasError object")
  if(!is.vector(test) & !is.matrix(test)) stop("'test' is not a vector or matrix")
  if(missing(reference)) reference = NA
  if(all(is.na(reference)) & !is.matrix(test))
    stop("if there is no reference, replicate test measurements are needed")
  if(!all(is.na(reference)) & !is.vector(reference)) stop("'reference' is not a vector")
  if(NCOL(test) > 2){
    stop("it is currently not possible to have more than two test measures")}
  if(all(is.na(reference)) & is.matrix(test)){
    test <- data.frame(test)
    colnames(test) <- c("test1", "test2")
    if(any(is.na(test$test1)) == TRUE){
      stop("the first replicate measure cannot contain missing values")
    }
    out <- list(test = test)
    input <- as.list(match.call()$test)[2:3]
    input <- list(test = input)
    names(input$test) <- c("test1", "test2")
    type <- "replicate"}
  else {
    out <- data.frame(test = test, reference = reference)
    input <- c(test = as.list(match.call())$test, reference = as.list(match.call())$reference)
    type <- "internal"}
  attr(out, "input") <- input
  attr(out, "type") <- type
  attr(out, "call") <- match.call()
  class(out) <- c("MeasError", "data.frame")
  out
}
#' @export
print.MeasError <- function(x){
  cat("\nCall:\n", deparse(attributes(x)$call), "\n", sep = "")
  if(is.vector(x$test) & is.vector(x$reference)){
  cat("\nThe error-prone variable", deparse((attr(x, 'input')$test)),
      "is correctly measured by",   deparse((attr(x, 'input')$reference)))}
  if(is.data.frame(x$test) & is.null(x$reference)){
    cat("\n", deparse((attr(x, 'input')$test$test1)), " and ",
              deparse((attr(x, 'input')$test$test2)),
              " are replicate measures with classical measurement error.", sep = "")
  }
  invisible(x)
}
