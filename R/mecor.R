#' mecor: a measurement error correction package
#'
#' mecor provides correction methods for measurement
#' error in a continuous covariate.
#'
#' @param formula an object of class \link[stats]{formula} (or one that is
#' coerced to that class): a symbolic description of the model containing
#' a \link[mecor]{MeasError} object.
#' @param data a data.frame, list or environment (or
#' object coercible by as.data.frame to a data frame) containing
#' the variables in the model
#' @param method a character string indicating the method used to correct for
#' measurement error, either "rc" (regression calibration), "rc_pooled1" (efficient
#' regression calibration using delta variance for pooling) or "rc_pooled2" (efficient
#' regression calibration using bootstrap variance for pooling).
#' @param alpha alpha level used to construct confidence intervals
#' @param B number of bootstrap samples
#'
#' @return \code{mecor} returns an object of \link[base]{class} "mecor"
#'
#' An object of class \code{mecor} is a list containing the following components:
#'
#' \item{uncorfit}{a lm.fit object of the uncorrected fit}
#' \item{corfit}{a lm.fit object of the corrected fit (if method = "rc") and a
#' matrix containing the corrected coefficients else}
#' \item{corvar}{the corrected variance using the delta method}
#' \item{ci.fieller}{fieller confidence interval (if method = "rc") else NA}
#' \item{ci.b}{bootstrap confidence interval (if B != 0)}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @references
#' L.Nab, R.H.H. Groenwold, P.M.J. Welsing, M. van Smeden.
#' Measurement error in continuous endpoints in randomised trials: an exploration of problems and solutions
#'
#' @examples
#' # measurement error in exposure
#' data(ivs)
#' fit <- mecor(Y ~ MeasError(X_star, X) + Z, data = ivs, method = "rc", B = 666)
#' data(rs)
#' mecor(Y ~ MeasError(X1_star, replicate = cbind(X2_star, X3_star)) + Z,
#' data = rs,
#' B = 999)
#' data(cs)
#' mecor(Y ~ MeasError(X_star, replicate = cbind(X1_star, X2_star)) + Z, data = cs)
#' mecor(Y ~ MeasError(W, X) + Z, ivs, method = "rc_pooled1")
#' mecor(Y ~ MeasError(W, X) + Z, ivs, method = "rc_pooled2")
#' @export
mecor <- function(formula,
                  data,
                  method = "rc",
                  alpha = 0.05,
                  B = 0,
                  use_vcov = "default"){
  if (missing(data))
    stop("data is missing without a default")
  if ((exists(deparse(substitute(data))) && is.function(data)) |
       !exists(deparse(substitute(data))))
    stop(paste0("data argument ", deparse(substitute(data)), " not found"))
  if (!is.data.frame(data))
    tryCatch({
      data <- as.data.frame(data)
      },
      error = function(e){
        message(e)
        stop(paste0("data argument ", data, " cannot be coerced to a data.frame"))
      }
    )
  if (missing(formula))
    stop("formula not found")
  if (class(formula) != "formula")
    stop("formula is not of class 'formula'")
  if (! method %in% c("rc", "erc"))
    stop("this method is not implemented")

  # Create response, covars and me (= MeasError object)
  vars_formula <- as.list(attr(terms(formula), "variables"))[-1]
  ind_me <- grep("MeasError", vars_formula) # index of MeasError in list of variables
  ind_response <- attributes(terms(formula))$response
  if (length(ind_me) == 0){
    stop("formula should contain a MeasError object")
    } else if (length(ind_me) != 1){
    stop("formula can only contain one MeasError object")
    }
  if (ind_me == ind_response){
    type <- "dep"
  } else type <- "indep"
  vars_formula_eval <- sapply(vars_formula, eval, envir = data)
  me <- vars_formula_eval[[ind_me]]
  if (type == "indep"){
    response <- as.matrix(vars_formula_eval[[ind_response]])
    colnames(response) <- vars_formula[ind_response]
    if (!length(vars_formula_eval[-c(ind_me, ind_response)]) == 0){
      covars <- sapply(vars_formula_eval[-c(ind_me, ind_response)], cbind)
      colnames(covars) <- vars_formula[-c(ind_me, ind_response)]
    } else covars <- NULL
    uncorfit <- mecor:::uncorrected(response, covars, me, type)
    if(method == "rc"){
      corfit <- mecor:::regcal(response, covars, me, B, alpha)
    }
    if(method == "erc"){
      corfit <- mecor:::efficient_regcal(response, covars, me, B, alpha, use_vcov)
    }
  } else if (type == "dep"){
    covars <- sapply(vars_formula_eval[-ind_me], cbind)
    colnames(covars) <- vars_formula[-ind_me]
    uncorfit <- mecor:::uncorrected(response = NULL, covars, me, type)
    if(method == "rc"){
      corfit <- mecor:::regcal_o(covars, me, B, alpha)
    }
  }


  # output
  out <- list(uncorfit = uncorfit,
              corfit = corfit
              )
  class(out) <- 'mecor'
  attr(out, "call") <- match.call()
  attr(out, "B") <- B
  attr(out, "alpha") <- alpha
  return(out)
}
