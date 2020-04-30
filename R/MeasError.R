#' @title Create a Measurement Error Object
#'
#' @description
#' This function creates a measurement error object, usually used as the
#' independent or dependent variable in a regression if one wants to correct for
#' the measurement error in that variable using a reference variable.
#'
#' @param substitute a vector containing error-prone measurements or a matrix with
#' replicate error-prone measurements
#' @param reference a vector containing reference measurements
#'
#' @return \code{MeasError} returns an object of \link[base]{class} "MeasError".
#'
#' An object of class \code{MeasError} is a list containing the substitute (and
#' reference) variables and additionally has attributes input (the name of the
#' substitute and reference variables), type (the study type that provides
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
#' with (rs, MeasError(substitute = cbind(W, W2), reference = NA))
#' @export
MeasError <- function(substitute,
                      reference){
  print(match.call())
  if(missing(substitute))
    stop("'substitute' is missing in the MeasError object")
  if(!is.vector(substitute) & !is.matrix(substitute))
    stop("'substitute' is not a vector or matrix")
  if(missing(reference))
    reference = NA
  if(all(is.na(reference)) & !is.matrix(substitute))
    stop("if there is no reference, replicate substitute measurements are needed")
  if(!all(is.na(reference)) & !is.vector(reference))
    stop("'reference' is not a vector")
  if(all(is.na(reference)) & is.matrix(substitute)){
    substitute <- data.frame(substitute)
    nrep <- ncol(substitute)
    nrep_range <- 1:nrep
    colnames(substitute) <- sapply(nrep_range, FUN = function(x) paste0("substitute", x))
    if(any(is.na(substitute$substitute1)) == TRUE){
      stop("the first replicate measure cannot contain missing values")}
    out <- list(substitute = substitute)
    input <- as.list(match.call()$substitute)[-1]
    input <- list(substitute = input)
    names(input$substitute) <- colnames(substitute)
    type <- "replicate"}
  else {
    out <- data.frame(substitute = substitute, reference = reference)
    attr(out, "row.names") <- NULL
    input <- c(substitute = as.list(match.call())$substitute,
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
    cat("\nThe error-prone variable", deparse((attr(x, 'input')$substitute)),
        "is correctly measured by",   deparse((attr(x, 'input')$reference)))
  }
  if(attr(x, "type") == "replicate"){
    cat("\n", deparse((attr(x, 'input')$substitute$substitute1)), " and ",
              deparse((attr(x, 'input')$substitute$substitute2)),
              " are replicate measures with classical measurement error", sep = "")
  }
  invisible(x)
}
