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
#' @param erc_B number of bootstrap samples used by efficient regression
#' calibration to estimate the bootstrap vcov, used for pooling
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
#'          method = "irc",
#'          B = 999,
#'          erc_B = 0)
#' data(rs)
#' mecor(Y ~ MeasError(X1_star, replicate = cbind(X2_star, X3_star)) + Z1 + Z2,
#'       data = rs,
#'       method = "rc",
#'       B = 999)
#' data(cs)
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
#'         method = "mle",
#'         B = 999)
#' @export
mecor <- function(formula,
                  data,
                  method = "rc",
                  alpha = 0.05,
                  B = 0,
                  erc_B = 0){
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
  if (! method %in% c("rc", "erc", "irc", "mle"))
    stop("this method is not implemented")

  # Create response, covars and me (= MeasError object)
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
  if (type == "indep" & !is.null(me$differential)){
    stop("Differential measurement error is only supported in the dependent variable")
  } else if (type == "dep" & (!is.null(me$differential) | length(me$coef) == 4 | length(me$model$coef) == 4)){
    type <- "dep_diff"
  }
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
    if (method == "rc"){
      corfit <-
        mecor:::regcal(response, covars, me, B, alpha, type)
    } else if (method == "erc"){
      if (B != 0 & erc_B != 0){
        menu(c("yes", "no"),
             title = "You're about to bootstrap the vcov used by efficient regression calibration and bootstrap the variance of that estimator. This may take a while. Do you want to proceed?")
      }
      corfit <-
        mecor:::efficient_regcal(response, covars, me, B, alpha, type, erc_B)
    } else if (method == "irc"){
      if (startsWith(type, "dep")){
        stop("Inadmissible regression calibration is not suitable for measurement error in the dependent variable")
      }
      if (class(me)[1] == "MeasErrorExt"){
        stop("Inaddmissible regression calbration is not suitable for external designs")
      } else if (class(me)[1] == "MeasError" && (type == "indep" & !is.null(me$replicate))){
        stop("Inadmissible regression calibration is not suitable for a design with replicates")
      }
      corfit <-
        mecor:::inadmissible_regcal(response, covars, me, B, alpha, type)
    } else if (method == "mle"){
      if (startsWith(type, "dep")){
        stop("The maximum likelihood estimator is not suitabel for measurement error in the dependent variable")
      }
      corfit <-
        mecor:::mle(response, covars, me, B, alpha)
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
