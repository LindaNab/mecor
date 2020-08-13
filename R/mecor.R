#' mecor: a Measurement Error Correction Package
#'
#' mecor provides correction methods for measurement error in a continuous
#' covariate or outcome in linear regression models with a continuous outcome
#'
#' @param formula an object of class \link[stats]{formula} (or one that is
#' coerced to that class): a symbolic description of the regressino model
#' containing a \link[mecor]{MeasError}, \link[mecor]{MeasErrorExt} or
#' \link[mecor]{MeasErrorRandom} object in one of the covariates or the outcome.
#' @param data a data.frame, list or environment (or object coercible by
#' as.data.frame to a data frame) containing the variables in the model
#' specified in \code{formula}.
#' @param method a character string indicating the method used to correct for
#' the measurement error, either "rc" (regression calibration), "erc" (efficient
#' regression calibration), "irc" (inadmissible regression calibration) or "mle"
#' (maximum likelihood estimation). Defaults to "rc".
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
#' Measurement error in continuous endpoints in randomised trials: problems and solutions
#'
#' L. Nab, M. van Smeden, R.H. Keogh, and R.H.H. Groenwold.
#' mecor: an R package for measurement error correction
#'
#' @examples
#' ## measurement error in a covariate:
#' # internal covariate-validation study
#' data(icvs)
#' mecor(Y ~ MeasError(X_star, reference = X) + Z,
#'       data = icvs,
#'       method = "erc")
#' # replicates study
#' data(rs)
#' mecor(Y ~ MeasError(X1_star, replicate = cbind(X2_star, X3_star)) + Z1 + Z2,
#'       data = rs,
#'       method = "mle")
#' # covariate-calibration study
#' data(ccs)
#' mecor(Y ~ MeasError(X_star, replicate = cbind(X1_star, X2_star)) + Z,
#'       data = ccs,
#'       method = "erc")
#' # external covariate-validation study
#' data(ecvs)
#' calmod_fit <- lm(X ~ X_star + Z, data = ecvs)
#' data(icvs) # suppose reference X is not available
#' mecor(Y ~ MeasErrorExt(X_star, model = calmod_fit) + Z,
#'       data = icvs)
#' # sensitivity analyses
#' data(icvs) # suppose reference X is not available
#' # guesstimate the coefficients of the calibration model:
#' mecor(Y ~ MeasErrorExt(X_star, model = list(coef = c(0, 0.9, 0.2))) + Z,
#'       data = icvs)
#' # assume random measurement error in X_star of magnitude 0.25:
#' mecor(Y ~ MeasErrorRandom(X_star, error = 0.25) + Z,
#'       data = icvs)
#' data(rs) # suppose replicates X2_star and X3_star are not available
#' mecor(Y ~ MeasErrorRandom(X1_star, error = 0.25) + Z1 + Z2,
#'       data = rs)
#'
#' ## measurement error in the outcome:
#' # internal outcome-validation study
#' data(iovs)
#' mecor(MeasError(Y_star, reference = Y) ~ X + Z,
#'       data = iovs,
#'       method = "rc")
#' # external outcome-validation study
#' data(eovs)
#' memod_fit <- lm(Y_star ~ Y, data = eovs)
#' data(iovs) # suppose reference Y is not available
#' mecor(MeasErrorExt(Y_star, model = memod_fit) ~ X + Z,
#'       data = iovs,
#'       method = "rc")
#' # sensitivity analyses
#' data(iovs) # suppose reference Y is not available
#' # guesstimate the coefficients of the measurement error model:
#' mecor(MeasErrorExt(Y_star, model = list(coef = c(0, 0.5))) ~ X + Z,
#'       data = iovs,
#'       method = "rc")
#'
#' ## differential measurement error in the outcome:
#' # internal outcome-validation study
#' data(iovs_diff)
#' mecor(MeasError(Y_star, reference = Y, differential = X) ~ X,
#'       data = iovs_diff,
#'       method = "rc")
#' # sensitivity analysis
#' data(iovs_diff) # suppose reference Y is not available
#' # guesstimate the coefficients of the measurement error model:
#' mecor(MeasErrorExt(Y_star, model = list(coef = c(0, 0.5, 1, 1))) ~ X,
#'       data = iovs_diff,
#'       method = "rc")
#' @export
mecor <- function(formula,
                  data,
                  method = "rc",
                  B = 0){
  mecor:::check_input_mecor(formula, data, method)
  # Create response, covars and me (= MeasError(Ext/Random) object)
  vars_formula <- as.list(attr(terms(formula), "variables"))[-1]
  ind_me <- grep("MeasError", vars_formula) # index of MeasError(Ext/Random) in
                                            # list of variables
  mecor:::check_ind_me(ind_me)
  ind_response <- attributes(terms(formula))$response
  type <- mecor:::get_me_type(ind_me, ind_response)
  vars_formula_eval <- sapply(vars_formula, eval, envir = data)
  me <- vars_formula_eval[[ind_me]]
  B <- mecor:::check_me(me, B, type, method)
  if (type == "dep" & (!is.null(me$differential) | # MeasError
                       length(me$coef) == 4 | # MeasErrorExt.list
                       length(me$model$coef) == 4)){ # MeasErrorExt.lm
    type <- "dep_diff"
  }
  # init response and covars
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
  # corrected fit
  corfit <- switch(method,
                   "rc" = mecor:::regcal(response,
                                         covars,
                                         me,
                                         B,
                                         type),
                   "erc" = mecor:::efficient_regcal(response,
                                                    covars,
                                                    me,
                                                    B,
                                                    type),
                   "irc" = mecor:::inadmissible_regcal(response,
                                                       covars,
                                                       me,
                                                       B,
                                                       type),
                   "mle" = mecor:::mle(response,
                                       covars,
                                       me,
                                       B,
                                       type))
  # uncorrected fit
  uncorfit <- mecor:::uncorrected(response,
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
check_ind_me <- function(ind_me){
  if (length(ind_me) == 0){
    stop("formula should contain a MeasError(Ext/Random) object")
  } else if (length(ind_me) != 1){
    stop("formula can only contain one MeasError(Ext/Random) object")
  }
}
get_me_type <- function(ind_me,
                        ind_response){
  if (ind_me == ind_response){
    type <- "dep"
  } else type <- "indep"
}
check_me <- function(me,
                     B,
                     type,
                     method){
  if(class(me)[1] == "MeasErrorExt" && length(grep("MeasErrorExt.list", attributes(me)$call)) != 0){
    if (B != 0){
      B <- 0
      warning("B set to 0 since bootstrap cannot be used if the class of 'model' in the MeasErrorExt object is of type list")
    }
  }
  if (type == "indep" & !is.null(me$differential))
    stop("Differential measurement error is only supported in the dependent variable")
  if (type == "dep" & class(me)[1] == "MeasErrorRandom"){
    stop("Random measurement error in the dependent variable won't introduce bias in the fitted model, correction is not needed")
  }
  if (class(me)[1] == "MeasErrorRandom"){
    if (method != "rc"){
      stop("methods different than 'rc' are not supported for MeasErrorRandom objects")
    }
  }
  if (method == "mle"){
    if (startsWith(type, "dep")){
      stop("The maximum likelihood estimator does not accommodate correction of measurement error in the dependent variable")
    }
    if (class(me) == "MeasErrorExt" || class(me) == "MeasErrorRandom")
      stop("The maximum likelihood estimator does not accommodate measurement error correction using a 'MeasErrorExt' or 'MeasErrorRandom' object")
    if (class(me)[1] == "MeasError" & is.null(me$replicate)){
      stop("Replicates measures of the substitute measure in 'MeasError' are needed for maximum likelihood estimation")
    }
  }
  if (method == "irc"){
    if (startsWith(type, "dep")){
      stop("Inadmissible regression calibration is not suitable for measurement error in the dependent variable")
    }
    if (class(me)[1] == "MeasErrorExt"){
      stop("Inaddmissible regression calbration is not suitable for external designs")
    } else if (class(me)[1] == "MeasError" && (type == "indep" & !is.null(me$replicate))){
      stop("Inadmissible regression calibration is not suitable for a design with replicates")
    }
  }
  B
}
