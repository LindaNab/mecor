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
#' measurement error, either "rc" (regression calibration), "erc" (efficient
#' regression calibration using delta variance for pooling), "irc" (inadmissible
#' regression calibration)
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
#' fit <-
#'   mecor(Y ~ MeasError(X_star, reference = X) + Z,
#'          data = ivs,
#'          method = "erc",
#'          B = 999)
#' data(rs)
#' fit <-
#'   mecor(Y ~ MeasError(X1_star, replicate = cbind(X2_star, X3_star)) + Z1 + Z2,
#'         data = rs,
#'         method = "rc",
#'         B = 999)
#' data(cs)
#' cs_rep <- subset(cs, !is.na(X1_star) & !is.na(X2_star))
#' fit <-
#' mecor(Y ~ MeasError(X1_star, replicate = X2_star) + Z,
#'       data = cs_rep,
#'       method = "mle",
#'       B = 999)
#' fit <-
#' mecor(Y ~ MeasError(X_star, replicate = cbind(X1_star, X2_star)) + Z,
#'       data = cs,
#'       method = "erc",
#'       B = 999)
#' # measurement error in the outcome
#' data(ivs_o)
#' fit <-
#'   mecor(MeasError(Y_star, reference = Y) ~ X + Z,
#'         data = ivs_o,
#'         method = "rc",
#'         B = 999)
#' # differential measurement error in the outcome
#' fit <-
#'   mecor(MeasError(Y_star, reference = Y, differential = X) ~ X,
#'         data = ivs_diff_o,
#'         method = "rc",
#'         B = 999)
#' # external information
#' calmod_fit <- lm(X ~ X_star + Z, data = ivs)
#' fit <-
#'   mecor(Y ~ MeasErrorExt(X_star, model = calmod_fit) + Z,
#'         data = ivs,
#'         B = 999)
#' me_fit <- lm(Y_star ~ Y, data = ivs_o)
#' fit <-
#'   mecor(MeasErrorExt(Y_star, model = me_fit) ~ X + Z,
#'         data = ivs_o,
#'         B = 999)
#' fit <-
#'   mecor(MeasErrorExt(Y_star, model = list(coef = c(0, 0.5))) ~ X + Z,
#'         data = ivs_o)
#' # mle
#' data(rs)
#' fit <-
#'   mecor(Y ~ MeasError(X1_star, replicate = cbind(X2_star, X3_star)) + Z1 + Z2,
#'         data = rs,
#'         method = "mle")
#' @export
mecor <- function(formula,
                  data,
                  method = "rc",
                  alpha = 0.05,
                  B = 0){
  mecor:::check_input_mecor(formula, data, method)
  # Create response, covars and me (= MeasError(Ext) object)
  vars_formula <- as.list(attr(terms(formula), "variables"))[-1]
  ind_me <- grep("MeasError", vars_formula) # index of MeasError(Ext) in list of variables
  if (length(ind_me) == 0){
    stop("formula should contain a MeasError or MeasErrorExt object")
  } else if (length(ind_me) != 1){
    stop("formula can only contain one MeasError or MeasErrorExt object")
  }
  ind_response <- attributes(terms(formula))$response
  if (ind_me == ind_response){
      type <- "dep"
  } else type <- "indep"
  vars_formula_eval <- sapply(vars_formula, eval, envir = data)
  me <- vars_formula_eval[[ind_me]]
  mecor:::check_me(me, B, type)
  if (type == "dep" & (!is.null(me$differential) | length(me$coef) == 4 | length(me$model$coef) == 4))
    type <- "dep_diff"
  if (type == "indep"){
    response <- as.matrix(vars_formula_eval[[ind_response]])
    colnames(response) <- vars_formula[ind_response]
    if (!length(vars_formula_eval[-c(ind_me, ind_response)]) == 0){
      covars <- sapply(vars_formula_eval[-c(ind_me, ind_response)], cbind)
      colnames(covars) <- vars_formula[-c(ind_me, ind_response)]
    } else covars <- NULL
  } else if (startsWith(type, "dep")){
    covars <- sapply(vars_formula_eval[-ind_me], cbind)
    colnames(covars) <- vars_formula[-ind_me]
  }
  uncorfit <- mecor:::uncorrected(response, covars, me, type)
  corfit <- switch(method,
                   "rc" = mecor:::regcal(response, covars, me, B, alpha, type),
                   "erc" = mecor:::efficient_regcal(response, covars, me, B, alpha, type),
                   "irc" = mecor:::inadmissible_regcal(response, covars, me, B, alpha, type),
                   "mle" = mecor:::mle(response, covars, me, B, alpha, type))
  # output
  out <- list(uncorfit = uncorfit,
              corfit = corfit)
  class(out) <- 'mecor'
  attr(out, "call") <- match.call()
  attr(out, "B") <- B
  attr(out, "alpha") <- alpha
  return(out)
}
check_input_mecor <- function(formula,
                              data,
                              method){
  if (missing(data))
    stop("data is missing without a default")
  if ((exists(deparse(substitute(data))) && is.function(data)) |
      !exists(deparse(substitute(data))))
    stop(paste0("data argument ", deparse(substitute(data)), " not found"))
  if (!is.data.frame(data))
    tryCatch({
      data <- as.data.frame(data)
    },
    error = function(e) {
      message(e)
      stop(paste0("data argument ", data, " cannot be coerced to a data.frame"))
    })
  if (missing(formula))
    stop("formula not found")
  if (class(formula) != "formula")
    stop("formula is not of class 'formula'")
  if (!method %in% c("rc", "erc", "irc", "mle"))
    stop("this method is not implemented"
    )
}
check_me <- function(me, B, type){
  if(class(me)[1] == "MeasErrorExt" && length(grep("MeasErrorExt.list", attributes(me)$call)) != 0){
    if (B != 0){
      B <- 0
      warning("B set to 0 since bootstrap cannot be used if the class of 'model' in the MeasErrorExt object is list")
    }
  }
  if (type == "indep" & !is.null(me$differential))
    stop("Differential measurement error is only supported in the dependent variable")
}
