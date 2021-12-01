#' @title Create a Random Measurement Error Object
#'
#' @description
#' This function creates a random measurement error object, usually used as
#' a covariate in the \code{formula} argument of \link[mecor]{mecor} if one
#' wants to correct for random measurement error in that variable
#'
#' @param substitute a vector containing the error-prone measure
#' @param variance a numeric quantifying the assumed variance of the random measurement error
#'
#' @return \code{MeasErrorRandom} returns an object of \link[base]{class}
#' "MeasErrorRandom".
#'
#' An object of class \code{MeasErrorRandom} is a list containing the substitute
#' variable, the assumed variance of the random measurement error in that variable and, the
#' attributes input (the name of the substitute variable) and call (the matched
#' call).
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ## random measurement error in a covariate:
#' # internal covariate-validation study
#' data(bloodpressure)
#' with(bloodpressure, MeasErrorRandom(sbp30, variance = 0.25))
#' @export
MeasErrorRandom <- function(substitute,
                            variance) {
  if (missing(substitute))
    stop("'substitute' in the MeasErrorRandom object is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the MeasErrorRandom object is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the MeasErrorRandom object cannot contain missing values")
  if (!is.numeric(variance)) {
    stop("'variance' in the MeasErrorRandom object is not numeric")
  }
  out <- list(substitute = substitute,
              variance = variance)
  input <- c(substitute = as.list(match.call())$substitute)
  attr(out, "input") <- input
  attr(out, "call") <- match.call()
  class(out) <- c("MeasErrorRandom", "list")
  out
}
