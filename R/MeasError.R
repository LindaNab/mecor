#' @title Create a Measurement Error Object
#'
#' @description
#' This function creates a measurement error object, usually used as the
#' independent or dependent variable in a regression if one wants to correct for
#' the measurement error in that variable using a reference variable.
#'
#' @param substitute a vector containing the error-prone measure
#' @param reference a vector containing the reference measure assumed without
#' measurement error
#' @param replicate a vector or matrix with replicates of the error-prone
#' measure with classical measurement error. This can either be
#' replicates obtained by using the same measurement method as the substitute
#' measure (replicates study) or replicates using a different measurement method
#' than the substitute measure (calibration study).
#' @param type a character string with the study type: 'ivs' for internal
#' validation study (default if reference not null), 'rs' for replicates study
#' (default if replicate not null), 'cs' for calibration study and 'evs' for
#' external validation study
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
#' me <- with (ivs, MeasError(X_star, reference = X))
#' # replicates study
#' data(rs)
#' me <- with (rs, MeasError(X1_star, replicate = cbind(X2_star, X3_star)))
#' # calibration study
#' data(cs)
#' me <- with(cs, MeasError(X_star, replicate = cbind(X1_star, X2_star), type = "cs"))
#' @export
MeasError <- function(substitute,
                      reference,
                      replicate,
                      type){
  # checks for substitute
  if (missing(substitute))
    stop("'substitute' in the MeasError object is missing")
  if (!is.vector(substitute))
    stop("'substitute' in the MeasError object is not a vector")
  if (any(is.na(substitute)) == TRUE)
    stop("'substitute' in the MeasError object cannot contain missing values")
  # check for reference and replicate (one of both should be non-null)
  if (!missing(reference) & !missing(replicate))
    stop("'reference' and 'replicate' in the MeasError object cannot be both non-null")
  if (missing(reference) & missing(replicate))
    stop("provide a 'reference' or 'replicate' in the MeasError object")
  # checks for type
  if (!missing(type)){
    if (length(type) != 1) # when type is a vector
      stop("variable 'type' has not length 1")
    if (!is.character(type))
      stop("variable 'type' is not a character")
    if (!type %in% c("ivs", "rs", "cs", "evs"))
      stop(paste0("study type", type, " is unknown"))
  }
  # checks for reference (internal validation study)
  if (!missing(reference)){
    if (!is.vector(reference))
      stop("'reference' is not a vector in the MeasError object")
    if (missing(type)){
      type <- "ivs"
    } else if (!missing(type) & type != "ivs"){
      warning("'reference' is non-null so type is set to default 'ivs'")
      type <- "ivs"
    }
  }
  # checks for replicate (replicates study/ calibration study)
  if (!missing(replicate)){
    if (!is.vector(replicate) & !is.matrix(replicate))
      stop("'replicate' is not a vector or matrix in the MeasError object")
    if(missing(type)){
      type <- "rs"
    } else if (!missing(type) & (type != "rs" & type != "cs")){
      warning("'replicate' is non-null so type is set to default 'rs'")
      type <- "rs"
    }
    # get reference from replicate (its rowmeans)
    reference <- mecor:::get_ref_from_rep(replicate)
  }
  out <- list(substitute = substitute,
              reference = reference)
  if(!missing(replicate)) out$replicate <- replicate
  input <- c(substitute = as.list(match.call())$substitute,
             reference = as.list(match.call())$reference,
             replicate = as.list(match.call())$replicate)
  attr(out, "input") <- input
  attr(out, "type") <- type
  attr(out, "call") <- match.call()
  class(out) <- c("MeasError", "list")
  out
}
get_ref_from_rep <- function(replicate){
  reference <- rowMeans(replicate)
  reference
}
#' @export
print.MeasError <- function(x){
  cat("\nCall:\n", deparse(attributes(x)$call), "\n", sep = "")
  if(attr(x, "type") == "ivs"){
    cat("\nThe error-prone variable", deparse((attr(x, 'input')$substitute)),
        "is correctly measured by",   deparse((attr(x, 'input')$reference)))
  }
  if(attr(x, "type") == "rs" | attr(x, "type") == "cs"){
    cat("\nThe error-prone variable ", deparse((attr(x, 'input')$substitute)),
        " has replicate measures ", deparse((attr(x, 'input')$replicate)),
              " with classical measurement error", sep = "")
  }
  invisible(x)
}
