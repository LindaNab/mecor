#' @title Create a Measurement Error Object
#'
#' @description
#' This function creates a measurement error object, usually used as a covariate
#' or the outcome in the \code{formula} argument of \link[mecor]{mecor} if one
#' wants to correct for the measurement error in that variable using a reference
#' variable or a replicate measure.
#'
#' @param substitute a vector containing the error-prone measure
#' @param reference a vector containing the reference measure assumed without
#' measurement error
#' @param replicate a vector or matrix with replicates of the error-prone
#' measure with classical measurement error. This can either be
#' replicates obtained by using the same measurement method as the substitute
#' measure (replicates study) or replicates using a different measurement method
#' than the substitute measure (calibration study).
#' @param differential a vector containing the variable to which the measurement
#' error is differential.
#'
#' @return \code{MeasError} returns an object of \link[base]{class} "MeasError".
#'
#' An object of class \code{MeasError} is a list containing the substitute and
#' reference (and replicate or differential if applicable) variables and has
#' attributes input (the name of the substitute and reference or replicate
#' and differential (if applicable) variables) and call (the matched call).
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ## measurement error in a covariate:
#' # internal covariate-validation study
#' data(icvs)
#' with (icvs, MeasError(substitute = X_star,
#'                       reference = X))
#' # replicates study
#' data(rs)
#' with (rs, MeasError(substitute = X_star_1,
#'                     replicate = cbind(X_star_2, X_star_3)))
#' # covariate-calibration study
#' data(ccs)
#' with(ccs, MeasError(substitute = X_star,
#'                     replicate = cbind(X_star_1, X_star_2)))
#' ## measurement error in the outcome:
#' # internal outcome-validation study
#' data(iovs)
#' with(iovs, MeasError(substitute = Y_star,
#'                      reference = Y))
#' # internal outcome- validation study with differential measurement error in
#' # the dependent variable
#' data(iovs_diff)
#' with(iovs_diff, MeasError(substitute = Y_star,
#'                           reference = Y,
#'                           differential = X))
#' @export
MeasError <- function(substitute,
                      reference,
                      replicate,
                      differential) {
  check_input_MeasError(substitute,
                        reference,
                        replicate,
                        differential)
  if (!missing(replicate)) {
    # get reference from replicate (rowmeans of all non NA rows)
    reference <- get_ref_from_rep(replicate)
  }
  out <- list(substitute = substitute,
              reference = reference)
  if (!missing(replicate))
    out$replicate <- replicate
  if (!missing(differential))
    out$differential = differential
  input <- c(
    substitute = as.list(match.call())$substitute,
    reference = as.list(match.call())$reference,
    replicate = as.list(match.call())$replicate,
    differential = as.list(match.call())$differential
  )
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  class(out) <- c("MeasError", "list")
  out
}
check_input_MeasError <- function(substitute,
                                  reference,
                                  replicate,
                                  differential) {
  # checks for substitute
  if (missing(substitute))
    stop("'substitute' in the MeasError object is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the MeasError object is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the MeasError object cannot contain missing values")
  # check for reference and replicate (one of both should be non-null)
  if (!missing(reference) & !missing(replicate))
    stop("'reference' or 'replicate' in the MeasError object should be null")
  if (missing(reference) & missing(replicate))
    stop("provide a 'reference' or 'replicate' variable in the MeasError object")
  # checks for reference (internal validation study)
  if (!missing(reference)) {
    if (!is.vector(reference))
      stop("'reference' is not a vector in the MeasError object")
  }
  # checks for replicate
  if (!missing(replicate)) {
    if (!missing(differential)) {
      # differential should be missing
      stop('differential measurement error cannot be corrected in a replicates study')
    }
    if (!is.vector(replicate) & !is.matrix(replicate)) {
      stop("'replicate' is not a vector or matrix in the MeasError object")
    }
  }
}
get_ref_from_rep <- function(replicate) {
  if (is.null(ncol(replicate))) {
    # only one replicate measure provided
    reference <- replicate
  } else {
    rownums_cc <- which (stats::complete.cases(replicate))
    reference <- rowMeans(replicate)
    reference[-rownums_cc] <- NA
  }
  reference
}
