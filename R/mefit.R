#' @title Fitting Calibration Model
#'
#' @description
#' This function fits a linear regression to your calibration dataset
#' in order to extract the measurement error parameters of the measurement
#' error model from the calibration set.
#'
#' @param formula an object of class \link[stats]{formula} (or one that is
#' coerced to that class): a symbolic description of the
#' model to be fitted in an \link[stats]{lm} model. If \code{memodel} is "systematic",
#' formula should be of the form 'V ~ Y'. If \code{me.structure} is "differential,
#' \code{formula} should be of the form 'V ~ Y * X'.
#' @param data an optional data frame, list or environment (or
#' object coercible by as.data.frame to a data frame) containing
#' the variables in the model. If not found in \code{data}, the
#' variables are taken from \code{environment(formula)}, typically
#' the enviroment from which \code{mefit} is called.
#' @param me.structure a character string indicating the underlying structure of the
#' measurement errors in your data, i.e. "classical", "systematic" or "differential".
#' @param dif.var a non-empty character string specifying the variable indicating the
#' grouping variable in formula if \code{memodel} is "differential".
#' @param robust indicates whether robust standard errors need to be calculated, the
#' "HC3" robust standard errors from \link[sandwich]{vcovHC} are used for the
#' heteroskedasticity-consistent estimation of the covariance matrix of the coefficient
#' estimates.
#'
#' @return \code{mefit} returns an object of \link[base]{class} "mefit".
#'
#' An object of class \code{mefit} is a list containing the following components:
#'
#' \item{coefficients}{a named vector containing the coefficients of the fitted calibration model}
#' \item{vcov}{the variance-covariance matrix of the parameters of the fitted calibration model, in case
#' \code{robust} is TRUE, the heteroskedasticity-consistent covariance matrix is returned}
#' \item{size}{the size of the calibration sample}
#' \item{call}{the matched call}
#' \item{me.structure}{the assumed structure of the measurement errors}
#' \item{dif.var}{the used grouping variable}
#' \item{diflevels}{the levels of the grouping variable 'dif.var'}
#' \item{rdf}{the residual degrees of freedom}
#' \item{r.squared}{R^2, the 'fraction of variance explained by the model', see \link[stats]{lm}
#' for calculation}
#' \item{sigma}{the square root of the estimated variance of the random error, see \link[stats]{lm}
#' for calculation}
#' \item{model}{the model frame used}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ##generation of external calibration set
#' Xcal <- c(rep(0, 500), rep(1, 500))
#' Ycal <- Xcal + rnorm(1000, 0, 1)
#' Vcal_cme <- Ycal + rnorm(1000, 0, 3) #classical measurement error (cme)
#' Vcal_sme <- 1 + 2 * Ycal + rnorm(1000, 0, 3) #systematic measurement error (sme)
#' Vcal_dme <- 2 + 2 * Xcal + 3 * Ycal + 2 * Xcal * Ycal + rnorm(1000, 0, 3 * (1 - Xcal) + 2 * Xcal) #differential measurement error (dme)
#'
#' fit_cme <- mefit(formula = Vcal_cme ~ Ycal - 1)
#' fit_sme <- mefit(formula = Vcal_sme ~ Ycal, me.structure = "systematic")
#' fit_dme <- mefit(formula = Vcal_dme ~ Ycal * Xcal, me.structure = "differential", dif.var = "Xcal", robust = T)
#'
#' @importFrom sandwich vcovHC
#' @import boot
#' @export
mefit <- function(formula,
                  data,
                  me.structure = "classical",
                  dif.var,
                  robust = FALSE){
  if(missing(data)) data = NULL
  else if(!is.data.frame(data)) data <- as.data.frame(data)
  if(missing(dif.var)) dif.var = NULL
  if(me.structure == "classical"){
    diflevels <- NA
    if(NROW(all.vars(formula)) != 2){
      stop("length of 'formula' wrong, should be of form 'V ~ Y' for classical me.structures")}
    if(attr(terms(formula), "intercept") == 1){
      warning("'formula' had an implied intercept term and is removed because 'me.structure' is classical")
      formula <- update(formula, . ~ . - 1)}
  }
  if(me.structure == "systematic"){
    diflevels <- NA
    if(NROW(all.vars(formula)) != 2){
      stop("length of 'formula' wrong, should be of form 'V ~ Y' for systematic me.structures")}
    if(!is.null(dif.var)){
      warning("'me.structure' is systematic so variable dif.var is set null")
      id = NULL}
    }
  if(me.structure == "differential"){
    if(is.null(dif.var)){
      stop("'me.structure' is differential but dif.var is null")  }
    if(NROW(all.vars(formula)) > 3){
      stop("length of 'formula' wrong, should be of form 'V ~ Y * X' for differential me.structures")  }
    if(!grepl(":", attr(terms(formula), "term.labels")[3])){
      stop("'me.structure' is differential, so there should be an interaction term in 'formula'")  }
    if(!is.null(data)) getdif.var <- data[dif.var]
      else if(grepl("$", dif.var, fixed = TRUE)) stop("use the 'data' argument to specify the dataframe where 'dif.var' can be found, use of $ is not supported here")
        else if(exists(dif.var)) getdif.var <- get(dif.var) else stop("'dif.var' does not exist in environment")
    if(!is.numeric(unique(getdif.var))) diflevels <- as.numeric(unlist(unique(getdif.var))) else diflevels <- unique(getdif.var)
    if(NROW(diflevels) == 1){
      warning("'number of levels of 'dif.var' is 1, so 'me.structure' is set 'systematic'")
      me.structure = "systematic"}
    if(NROW(diflevels) != 2){
      stop("number of levels of dif.var does not equal 2, this functionality is not supported by 'mecor'")  }
  }
  model <- lm(formula = formula, data = data)
  if(all(model$residuals == 0)){
    stop("all residuals of the measurement error fit are 0")  }
  if(robust == TRUE){
    vcov <- vcovHC(model) }
  else vcov <- vcov(model)
  out <- list(coefficients = model$coefficients,
              vcov = vcov,
              size = nrow(model$model),
              call = match.call(),
              me.structure = me.structure,
              dif.var = dif.var,
              diflevels = diflevels,
              rdf = model$df.residual,
              r.squared = summary(model)$r.squared,
              sigma = summary(model)$sigma,
              model = model$model)
  class(out) <- 'mefit'
  return(out)
}

delta.sme <- function(nm,
                  coef.cm,
                  mefit,
                  alpha){
  t1 <- coef(mefit)[2]
  b <- coef(nm)[2]
  vart1 <- summary(mefit)$coef[2,2] ^ 2
  varnb <- summary(nm)$coef[2,2] ^ 2
  varcb <- 1 / t1 ^ 2 * (varnb + b ^ 2 * vart1)
  tq <- qt((1 - alpha / 2), nm$df.residual)
  ci <- unname(c(coef.cm[2] - tq * sqrt(varcb), coef.cm[2] + tq * sqrt(varcb)))
  return(ci)
}

delta.dme <- function(nm,
                      coef.cm,
                      mefit,
                      alpha){
  #coefficients
  nb <- coef(nm)[2]
  na <- coef(nm)[1]
  t00 <- coef(mefit)[1]
  t10 <- ifelse(names(coef(mefit)[2]) == mefit$dif.var, coef(mefit)[3], coef(mefit)[2])
  t01 <- coef(mefit)[mefit$dif.var] + t00
  t11 <- coef(mefit)[4] + t10
  #variances
  varnb <- vcovHC(nm)[2,2]
  varna <- vcovHC(nm)[1,1]
  covnanb <- vcovHC(nm)[1,2]
  covm <- mefit$vcov
  if(names(coef(mefit)[2]) == mefit$dif.var){ #change order of covm if needed
    covm <- cbind(covm[,1], covm[,3], covm[,2], covm[,4])}
  vart00 <- covm[1,1]
  vart10 <- covm[2,2]
  vart01 <- covm[3,3] + covm[1,1] + 2 * covm[3,1]
  vart11 <- covm[4,4] + covm[2,2] + 2 * covm[4,2]
  covt10t00 <-  covm[2,1]
  covt11t01 <- covm[3,4] + covm[3,2] + covm[1,4] + covm[1,2]
  varca <- 1 / t10 ^ 2 * (varna + coef.cm[1] ^ 2 * vart10 + vart00 + 2 * coef.cm[1] * covt10t00)
  varcb <- 1 / t11 ^ 2 * ((coef.cm[2] + coef.cm[1]) ^ 2 * vart11 + varnb + varna + 2 * covnanb +
                            vart01 + 2 * (coef.cm[2] + coef.cm[1]) * covt11t01 ) + varca
  tq <- qt((1 - alpha / 2), nm$df.residual)
  ci <- unname(c(coef.cm[2] - tq * sqrt(varcb), coef.cm[2] + tq * sqrt(varcb)))
  return(ci)
}

fieller <- function(nm,
                    mefit,
                    alpha){
  t1 <- coef(mefit)[2]
  b <- coef(nm)[2]
  vart1 <- summary(mefit)$coef[2,2] ^ 2
  varnb <- summary(nm)$coef[2,2] ^ 2
  tq <- qt((1 - alpha / 2), nm$df.residual)
  v1 <- - 1 * (b * t1)
  v2 <- vart1 * tq ^ 2 - t1 ^ 2
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
}

bootstrap.sme <- function(nm, mefit, alpha, B){
  statme.sme <- function(data, indices){
    if(count == 1) {
      t1 <- coef(mefit)[2]}
    else {
      ids <- sample(1:nrow(mefit$model), replace = T)
      calsample <- mefit$model[ids,]
      t1 <- coef(lm(formula = formula(mefit), data = calsample))[2]}
    count <<- count + 1
    d <- data[indices,]
    return(coef(lm(formula = formula(nm), data = d))[2] / t1)
  }
  count <- 1
  meboot <- boot(data = nm$model, statistic = statme.sme,
                 R = B)
  ci <- boot.ci(boot.out = meboot, conf = (1 - alpha), type = "perc")
  return(ci)
}

bootstrap.dme <- function(nm, mefit, alpha, B){
  statme.dme <- function(data, indices){
    if(count == 1) {
      t00 <- coef(mefit)[1]
      t10 <- ifelse(names(coef(mefit)[2]) == mefit$dif.var, coef(mefit)[3], coef(mefit)[2])
      t01 <- coef(mefit)[mefit$dif.var] + t00
      t11 <- coef(mefit)[4] + t10}
    else {
      ids <- sample(1:nrow(mefit$model), replace = T)
      calsample <- mefit$model[ids,]
      calfit <- lm(formula = formula(mefit), data = calsample)
      t00 <- coef(calfit)[1]
      t10 <- ifelse(names(coef(calfit)[2]) == mefit$dif.var, coef(calfit)[3], coef(calfit)[2])
      t01 <- coef(calfit)[mefit$dif.var] + t00
      t11 <- coef(calfit)[4] + t10}
    count <<- count + 1
    d <- data[indices,]
    fit <- lm(formula = formula(nm), data = d)
    return( (coef(fit)[2] + coef(fit)[1] - t01)/ t11 - (coef(fit)[1] - t00) / t10)
  }
  count <- 1
  meboot <- boot(data = nm$model, statistic = statme.dme,
                 R = B)
  ci <- boot.ci(boot.out = meboot, conf = (1 - alpha), type = "perc")
  return(ci)
}

