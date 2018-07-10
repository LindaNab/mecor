#' @title Fitting Calibration Model
#'
#' @description
#' This function fits a linear regression to your calibration dataset
#' in order to extract the measurement error parameters of the measurement
#' error model from the external calibration set.
#'
#' @param formula an object of class \link[stats]{formula} (or one that is
#' coerced to that class): a symbolic description of the
#' model to be fitted in an \link[stats]{lm} model. If \code{memodel} is "systematic",
#' formula should be of the form 'V ~ Y'. If \code{memodel} is "differential,
#' \code{formula} should be of the form 'V ~ Y * X'.
#' @param data an optional data frame, list or environment (or
#' object coercible by as.data.frame to a data frame) containing
#' the variables in the model. If not found in \code{data}, the
#' variables are taken from \code{environment(formula)}, typically
#' the enviroment from which \code{mefit} is called.
#' @param mestructure a character string indicating the underlying structure of the measurement errors
#' in your data, i.e. "classical", "systematic" or "differential".
#' @param difvar variable indicating the grouping variable in formula if \code{memodel} is "differential".
#' @param robust indicates whether robust standard errors need to be calculated, the "HC3"
#' robust standard errors from \link[sandwich]{vcovHC} are used for the heteroskedasticity-consistent
#' estimation of the covariance matrix of the coefficient estimates.
#'
#' @return \code{mefit} returns an object of \link[base]{class} "mefit".
#'
#' An object of class \code{mefit} is a list containing the following components:
#'
#' \item{coefficients}{a named vector containing the coefficients of the fitted calibration model}
#' \item{vcov}{the variance-covariance matrix of the parameters of the fitted calibration model, in case
#' \code{robust} is TRUE, the heteroskedasticity-consistent covariance matrix is returned}
#' \item{size}{the size of the calibration sample}
#' \item{meanx}{mean of explanatory variable of the fitted calibration model}
#' \item{call}{the matched call}
#' \item{mestructure}{the assumed structure of the measurement errors}
#' \item{difvar}{the used grouping variable}
#' \item{diflevels}{the levels of the grouping variable 'difvar'}
#' \item{rdf}{the residual degrees of freedom}
#' \item{r.squared}{R^2, the 'fraction of variance explained by the model', see \link[stats]{lm}
#' for calculation}
#' \item{sigma}{the square root of the estimated variance of the random error, see \link[stats]{lm}
#' for calculation}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @references
#'
#' @examples
#' Xcal <- c(rep(0,500),rep(1,500))
#' Ycal <- Xcal + rnorm(1000,0,3)
#'
#' ##systematic measurement error (sme)
#' Vcal_sme <- 1 + 2 * Ycal + rnorm(1000,0,3)
#' fit_sme <- mefit(formula = Vcal_sme ~ Ycal, mestructure = "systematic")
#'
#' ##differential measurement error (dme)
#' Vcal_dme <- 2 + 2 * Xcal + 3 * Ycal + 2 * Xcal * Ycal + rnorm(1000,0,3)
#' fit_dme <- mefit(formula = Vcal_dme ~ Ycal * Xcal, mestructure = "differential", difvar = "Xcal")
#'
#' @importFrom sandwich vcovHC
#' @export
mefit <- function(formula,
                  data = NULL,
                  mestructure = "systematic",
                  difvar = NULL,
                  robust = FALSE){
  if(mestructure == "systematic"){
    if(attr(terms(formula), "variables")[4] != "NULL()"){
      stop("length of 'formula' too long, should be of form 'V ~ Y' for systematic mestructures")}
    if(!is.null(difvar)){
      warning("'mestructure' is systematic so variable difvar is set null")
      id <- NULL}
    }
  if(mestructure == "differential"){
    if(is.null(difvar)){
      stop("'mestructure' is differential but there is no difvar")  }
    if(attr(terms(formula), "variables")[5] != "NULL()"){
      stop("length of 'formula' too long, should be of form 'V ~ Y * X' for differential mestructures")  }
    if(!grepl(":", attr(terms(formula), "term.labels")[3])){
      stop("'mestructure' is differential, so there should be an interaction term in 'formula'")  }
    if(NROW(unique(get(difvar))) == 1){
      warning("'number of levels of 'difvar' is 1, so 'mestructure' is set 'systematic'")
      mestructure = "systematic"}
    if(NROW(unique(get(difvar))) != 2){
      stop("number of levels of difvar does not equal 2, this functionality is not supported by 'mecor'")  }
    }
  model <- lm(formula = formula, data = data, x = TRUE)
  if(robust == TRUE){
    vcov <- vcovHC(model) }
  else vcov <- vcov(model)
  if(mestructure == "differential") diflevels = unique(get(difvar)) else diflevels = NA
  out <- list(coefficients = model$coefficients,
              vcov = vcov,
              size = nrow(model$x),
              meanx = mean(model$x),
              call = match.call(),
              mestructure = mestructure,
              difvar = difvar,
              diflevels = diflevels,
              rdf = model$df.residual,
              r.squared = summary(model)$r.squared,
              sigma = summary(model)$sigma)
  class(out) <- 'mefit'
  return(out)
}

delta <- function(nm,
                  coef.cm,
                  mefit,
                  alpha){
  t1 <- coef(mefit)[2]
  b <- coef(nm)[2]
  vart1 <- summary(mefit)$coef[2,2]
  varnb <- summary(nm)$coef[2,2]
  varcb <- 1 / t1 ^ 2 * (varnb + b ^ 2 * vart1)
  tq <- qt((1 - alpha / 2), nm$df.residual)

  ##X <- model_naive$model[,as.character(model_naive$terms[[3]])]
  ##s_xx <- sum((X - mean(X)) ^ 2)
  ##t_q <- qt((1 - systematic$alpha / 2),(NROW(X) - 2))
  ##varbeta <- 1 / (systematic$coefficients["theta1_hat"]) ^ 2 *
  ##              (summary(model_naive)$sigma ^ 2 / s_xx +
  ##              (model_corrected$coefficients[2] ^ 2 * systematic$coefficients["t"] ^ 2) /
  ##              systematic$coefficients["s_yy"] )
  ci <- unname(c(coef.cm[2] - tq * varcb, coef.cm[2] + tq * varcb))
  return(ci)
}

fieller <- function(nm,
                    mefit,
                    alpha){
  t1 <- coef(mefit)[2]
  b <- coef(nm)[2]
  vart1 <- summary(mefit)$coef[2,2]
  varnb <- summary(nm)$coef[2,2]
  tq <- qt((1 - alpha / 2), nm$df.residual)
  v1 <- - 1 * (b * t1)
  v2 <- vart1 * tq ^ 2 - t1 ^2
  v3 <- varnb * tq ^ 2 - b ^ 2
  D <- v1 ^ 2 - v2 * v3
  if(v2 < 0 & D > 0){
    l1 <- unname((v1 - sqrt(D)) / v2)
    l2 <- unname((v1 + sqrt(D)) / v2)
    ci <- c(lower = min(l1, l2),
            upper = max(l1, l2)
    )
  }
  else ci <- c(lower = NA, upper = NA)
  return(ci)

  #X <- model_naive$model[,as.character(model_naive$terms[[3]])]
  #s_xx <- sum((X - mean(X))^2)
  #t_q <- qt((1 - systematic$alpha / 2), (NROW(X) - 2))
  #v1 <- - 1 *  (model_naive$coefficients[2] * systematic$coefficients["theta1_hat"])
  #v2 <- (systematic$coefficients["t"] ^ 2 / systematic$coefficients["s_yy"]) * t_q ^ 2 - systematic$coefficients["theta1_hat"] ^ 2
  #v3 <- (summary(model_naive)$sigma ^ 2 / s_xx) * t_q ^ 2 - model_naive$coefficients[2]^2
}

