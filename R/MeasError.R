#' @title Create a Measurement Error Object
#'
#' @description
#' This function creates a measurement error object, usually used as the
#' independent or dependent variable in a regression if one wants to correct for
#' the measurement error in that variable using a reference variable.
#'
#' @param test a vector containing error-prone measurements or a matrix with
#' replicate error-prone measurements
#' @param reference a vector containing reference measurements
#'
#' @return \code{MeasError} returns an object of \link[base]{class} "MeasError".
#'
#' An object of class \code{MeasError} is a list containing the test (and
#' reference) variables and additionally has attributes input (the name of the
#' test and reference variables), type (the study type that provides
#' information to correct for the measurement error) and call (the matched
#' call).
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ## measurement error in exposure
#' # internal validation study
#' data(ivs)
#' with (ivs, MeasError(W, X))
#' # replicates study
#' data(rs)
#' with (rs, MeasError(test = cbind(W, W2), reference = NA))
#' @export
MeasError <- function(test,
                      reference){
  print(match.call())
  if(missing(test))
    stop("'test' is missing in the MeasError object")
  if(!is.vector(test) & !is.matrix(test))
    stop("'test' is not a vector or matrix")
  if(missing(reference))
    reference = NA
  if(all(is.na(reference)) & !is.matrix(test))
    stop("if there is no reference, replicate test measurements are needed")
  if(!all(is.na(reference)) & !is.vector(reference))
    stop("'reference' is not a vector")
  if(all(is.na(reference)) & is.matrix(test)){
    test <- data.frame(test)
    nrep <- ncol(test)
    nrep_range <- 1:nrep
    colnames(test) <- sapply(nrep_range, FUN = function(x) paste0("test", x))
    if(any(is.na(test$test1)) == TRUE){
      stop("the first replicate measure cannot contain missing values")}
    out <- list(test = test)
    input <- as.list(match.call()$test)[-1]
    input <- list(test = input)
    names(input$test) <- colnames(test)
    type <- "replicate"}
  else {
    out <- data.frame(test = test, reference = reference)
    attr(out, "row.names") <- NULL
    input <- c(test = as.list(match.call())$test,
               reference = as.list(match.call())$reference)
    type <- "internal"
  }
  attr(out, "input") <- input
  attr(out, "type") <- type
  attr(out, "call") <- match.call()
  class(out) <- c("MeasError", "data.frame")
  out
}
#' @export
print.MeasError <- function(x){
  cat("\nCall:\n", deparse(attributes(x)$call), "\n", sep = "")
  if(attr(x, "type") == "internal"){
    cat("\nThe error-prone variable", deparse((attr(x, 'input')$test)),
        "is correctly measured by",   deparse((attr(x, 'input')$reference)))
  }
  if(attr(x, "type") == "replicate"){
    cat("\n", deparse((attr(x, 'input')$test$test1)), " and ",
              deparse((attr(x, 'input')$test$test2)),
              " are replicate measures with classical measurement error", sep = "")
  }
  invisible(x)
}
