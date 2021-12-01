#' mecor: a Measurement Error Correction Package
#'
#' mecor provides correction methods for measurement error in a continuous
#' covariate or outcome in linear regression models with a continuous outcome
#'
#' @param formula an object of class \link[stats]{formula} (or one that is
#' coerced to that class): a symbolic description of the regression model
#' containing a \link[mecor]{MeasError}, \link[mecor]{MeasErrorExt} or
#' \link[mecor]{MeasErrorRandom} object in one of the covariates or the outcome.
#' @param data a data.frame, list or environment (or object coercible by
#' as.data.frame to a data frame) containing the variables in the model
#' specified in \code{formula}.
#' @param method a character string indicating the method used to correct for
#' the measurement error, either "standard" (regression calibration for
#' covariate measurement error and method of moments for outcome measurement
#' error), "efficient" (efficient regression calibration for covariate
#' measurement error and efficient method of moments for outcome measurement
#' error), "valregcal" (validation regression calibration) or "mle" (maximum
#' likelihood estimation). Defaults to "standard".
#' @param B number of bootstrap samples, defaults to 0.
#'
#' @return \code{mecor} returns an object of \link[base]{class} "mecor".
#'
#' An object of class \code{mecor} is a list containing the following components:
#'
#' \item{corfit}{a list containing the corrected fit, including the coefficients
#' of the corrected fit (\code{coef}) and the variance--covariance matrix of the
#' coefficients of the corrected fit obtained by the delta method (\code{vcov}),
#' and more depending on the method used.}
#' \item{uncorfit}{an \link[stats]{lm.fit} object of the uncorrected fit.}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @references
#' L. Nab, R.H.H. Groenwold, P.M.J. Welsing, and  M. van Smeden.
#' Measurement error in continuous endpoints in randomised trials: problems and
#' solutions
#'
#' L. Nab, M. van Smeden, R.H. Keogh, and R.H.H. Groenwold.
#' mecor: an R package for measurement error correction in linear models with
#' continuous outcomes
#'
#' @examples
#' ## measurement error in a covariate/outcome:
#' # internal covariate-validation study
#' data(vat)
#' out <-
#' mecor(ir_ln ~ MeasError(wc, reference = vat) + sex + age + tbf,
#'       data = vat,
#'       method = "standard",
#'       B = 999)
#' # replicates study
#' data(bloodpressure)
#' mecor(creatinine ~ MeasError(sbp30, replicate = cbind(sbp60, sbp120)) + age,
#'       data = bloodpressure,
#'       method = "mle")
#' # outcome-calibration study
#' data(sodium)
#' mecor(MeasError(recall, replicate = cbind(urinary1, urinary2)) ~ diet,
#'       data = sodium,
#'       method = "efficient")
#' # external outcome-validation study
#' data(haemoglobin_ext)
#' calmod_fit <- lm(capillary ~ venous, data = haemoglobin_ext)
#' data(haemoglobin) # suppose reference venous is not available
#' mecor(MeasErrorExt(capillary, model = calmod_fit) ~ supplement,
#'       data = haemoglobin)
#' # sensitivity analyses
#' data(vat) # suppose reference vat is not available
#' # guesstimate the coefficients of the calibration model:
#' mecor(ir_ln ~ MeasErrorExt(wc, model = list(coef = c(0.2, 0.5, -1.3, 0, 0.6))) + sex + age + tbf,
#'       data = vat)
#' # assume random measurement error in wc of magnitude 0.25:
#' mecor(ir_ln ~ MeasErrorRandom(wc, variance = 0.25) + sex + age + tbf,
#'       data = vat)
#' data(bloodpressure) # suppose replicates sbp60 and sbp60 are not available
#' mecor(creatinine ~ MeasErrorRandom(sbp30, variance = 25) + age,
#'       data = bloodpressure)
#'
#' ## differential measurement error in the outcome:
#' # internal outcome-validation study
#' mecor(MeasError(capillary, reference = venous, differential = supplement) ~ supplement,
#'       data = haemoglobin,
#'       method = "standard")
#' @export
mecor <- function(formula,
                  data,
                  method = "standard",
                  B = 0) {
  check_input_mecor(formula,
                    data,
                    method)
  # Create response, covars and me (= MeasError(Ext/Random) object)
  vars_formula <- as.list(attr(stats::terms(formula), "variables"))[-1]
  vars_formula_eval <- sapply(vars_formula, eval, envir = data, enclos = parent.frame())
  if (any(unlist(sapply(vars_formula_eval, FUN = class)) ==  "MeasError")){
    ind_me <- which(unlist(sapply(vars_formula_eval, FUN = class)) ==  "MeasError")
  } else if (any(unlist(sapply(vars_formula_eval, FUN = class)) ==  "MeasErrorRandom")){
    ind_me <- which(unlist(sapply(vars_formula_eval, FUN = class)) ==  "MeasErrorRandom")
  } else if (any(unlist(sapply(vars_formula_eval, FUN = class)) ==  "MeasErrorExt")){
    ind_me <- which(unlist(sapply(vars_formula_eval, FUN = class)) ==  "MeasErrorExt")
  } else ind_me <- NA
  # list of variables
  check_ind_me(ind_me)
  ind_response <- attributes(stats::terms(formula))$response
  type <- get_me_type(ind_me, ind_response)
  me <- vars_formula_eval[[ind_me]]
  B <- check_me(me, B, type, method)
  if (type == "dep" & (!is.null(me$differential) | # MeasError
                       length(me$coef) == 4 | # MeasErrorExt.list
                       length(me$model$coef) == 4)) {
    # MeasErrorExt.lm
    type <- "dep_diff"
  }
  # init response and covars
  if (type == "indep") {
    response <- as.matrix(vars_formula_eval[[ind_response]])
    colnames(response) <- vars_formula[ind_response]
    if (!length(vars_formula_eval[-c(ind_me, ind_response)]) == 0) {
      covars <- sapply(vars_formula_eval[-c(ind_me, ind_response)], cbind)
      colnames(covars) <- vars_formula[-c(ind_me, ind_response)]
    } else
      covars <- NULL
  } else if (startsWith(type, "dep")) {
    covars <- sapply(vars_formula_eval[-ind_me], cbind)
    colnames(covars) <- vars_formula[-ind_me]
  }
  # corrected fit
  corfit <- switch(
    method,
    "standard" = standard(response,
                          covars,
                          me,
                          B,
                          type),
    "efficient" = efficient(response,
                            covars,
                            me,
                            B,
                            type),
    "valregcal" = valregcal(response,
                            covars,
                            me,
                            B,
                            type),
    "mle" = mle(response,
                covars,
                me,
                B,
                type)
  )
  # uncorrected fit
  uncorfit <- uncorrected(response,
                          covars,
                          me,
                          type)
  # output
  out <- list(corfit = corfit,
              uncorfit = uncorfit)
  class(out) <- 'mecor'
  attr(out, "call") <- match.call()
  attr(out, "B") <- B
  return(out)
}
check_input_mecor <- function(formula,
                              data,
                              method) {
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
  if (!method %in% c("standard", "efficient", "valregcal", "mle"))
    stop("this method is not implemented")
}

check_ind_me <- function(ind_me) {
  if (is.na(ind_me)) {
    stop("formula should contain a MeasError(Ext/Random) object")
  } else if (length(ind_me) != 1) {
    stop("formula can only contain one MeasError(Ext/Random) object")
  }
}
get_me_type <- function(ind_me,
                        ind_response) {
  if (ind_me == ind_response) {
    type <- "dep"
  } else
    type <- "indep"
}
check_me <- function(me,
                     B,
                     type,
                     method) {
  if (class(me)[1] == "MeasErrorExt" &&
      length(grep("MeasErrorExt.list", attributes(me)$call)) != 0) {
    if (B != 0) {
      B <- 0
      warning(
        "B set to 0 since bootstrap cannot be used if the class of 'model' in the MeasErrorExt object is of type list"
      )
    }
  }
  if (type == "dep" & !is.null(me$replicate)){
    if(is.null(dim(me$replicate)))
      stop("More than one replicate measure is needed for measurement error correction in the dependent variable")
  }
  if (type == "indep" & !is.null(me$differential))
    stop("Differential measurement error is only supported in the dependent variable")
  if (type == "dep" & class(me)[1] == "MeasErrorRandom") {
    stop(
      "Random measurement error in the dependent variable won't introduce bias in the fitted model, correction is not needed"
    )
  }
  if (class(me)[1] == "MeasErrorRandom") {
    if (method != "standard") {
      stop("methods different from 'standard' are not supported for MeasErrorRandom objects")
    }
  }
  if (method == "mle") {
    if (startsWith(type, "dep")) {
      stop(
        "The maximum likelihood estimator does not accommodate correction of measurement error in the dependent variable"
      )
    }
    if (class(me)[1] == "MeasErrorExt" || class(me)[1] == "MeasErrorRandom")
      stop(
        "The maximum likelihood estimator does not accommodate measurement error correction using a 'MeasErrorExt' or 'MeasErrorRandom' object"
      )
    if (class(me)[1] == "MeasError" & is.null(me$replicate)) {
      stop(
        "Replicates measures of the substitute measure in 'MeasError' are needed for maximum likelihood estimation"
      )
    }
  }
  if (method == "valregcal") {
    if (startsWith(type, "dep")) {
      stop(
        "Validation regression calibration is not suitable for measurement error in the dependent variable"
      )
    }
    if (class(me)[1] == "MeasErrorExt") {
      stop("Validation regression calbration is not suitable for external designs")
    } else if (class(me)[1] == "MeasError" &&
               (type == "indep" & !is.null(me$replicate))) {
      stop("Validation regression calibration is not suitable for a design with replicates")
    }
  }
  B
}
