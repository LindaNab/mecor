#' @title Fitting Measurement Error Model
#'
#' @description
#' This function fits a linear model to model the relationship between the test
#' measurement (error-prone measurement) and the reference measurement (true measurement)
#' in order to extract the measurement error parameters of the measurement
#' error model from a validation set.
#'
#' @param formula formula of the measurement error model
#' @param data optional data frame
#' @param robust boolean indicating whether robust standard errors need to be calculated, the
#' "HC3" robust standard errors from \link[sandwich]{vcovHC} are used for the
#' heteroskedasticity-consistent estimation of the covariance matrix of the coefficient
#' estimates.
#'
#' @return \code{mefit} returns an object of \link[base]{class} "mefit".
#'
#' An object of class \code{mefit} is lm object (the fit of the measurement error model) and
#' contains additionally, the following components:
#' \item{call}{the matched call}
#' \item{vcov}{the variance-covariance matrix of the parameters of the fitted measurement error model, in case
#' \item{me.model}{the assumed measurement error model}
#' \code{robust} is TRUE, the heteroskedasticity-consistent covariance matrix is returned}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ##generation of external validation set
#' Xval <- c(rep(0, 500), rep(1, 500))
#' Yval <- Xval + rnorm(1000, 0, 1)
#' Vval_cme <- Yval + rnorm(1000, 0, 3) #classical measurement error (cme)
#' Vval_sme <- 1 + 2 * Yval + rnorm(1000, 0, 3) #systematic measurement error (sme)
#' Vval_dme <- 2 + 2 * Xval + 3 * Yval + 2 * Xval * Yval + rnorm(1000, 0, 3 * (1 - Xval) + 2 * Xval) #differential measurement error (dme)
#'
#' #classical
#' fitCme <- mefit(MeasError(Vval_cme, Yval) ~ - 1)
#' #systematic
#' fitSme <- mefit(MeasError(Vval_sme, Yval) ~ 1)
#' #differential
#' fitDme <- mefit(MeasError(Vval_dme, Yval) ~ Xval)
#'
#' @importFrom sandwich vcovHC
#' @import boot
#' @export
mefit <- function(formula,
                  data,
                  robust = FALSE){
  if(missing(data)) data = NULL
  else if(!is.data.frame(data)) data <- as.data.frame(data)
  f <- mefit.formula(formula, data)
  fit <- lm(f)
  if(all(fit$residuals == 0)){
    stop("all residuals of the measurement error model fit are 0")  }
  if(robust == TRUE){
    vcov <- vcovHC(fit) }
  else vcov <- vcov(fit)
  out <- fit
  out$call <- match.call()
  out$vcov <- vcov
  out$me.model <- attr(f, "model")
  out$robust <- robust
  class(out) <- c("mefit", "lm")
  out
}

#'@export
mefit.formula <- function(formula, data)
{
  if(missing(formula)){
    stop("a formula argument is needed")}
  t <- terms(formula)
  if(charmatch(c("MeasError"), t) != 2L){
    stop("Response must be a MeasError object")}
  f <- eval.parent(t[[2L]])
  if(t[[3L]] == -1){
    f <- update(f, ~ . - 1)
    model <- "classical"}
  else if(t[[3L]] != 1){
    factor <- eval.parent(t[[3L]])
    diflevels <- unique(factor)
    if(NROW(diflevels) == 1){
      warning("number of levels of the explanatory variable in formula is 1, so formula is updated to non-differential systematic")
      model <- "systematic"}
    else {
      f <- update(f, ~ . * factor)
      model <- "differential"}}
  else {model <- "systematic"}
  out <- f
  attr(out, "call") <- match.call()
  attr(out, "model") <- model
  out
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
  t10 <- ifelse(names(coef(mefit)[2]) == mefit$differential, coef(mefit)[3], coef(mefit)[2])
  t01 <- coef(mefit)[mefit$differential] + t00
  t11 <- coef(mefit)[4] + t10
  #variances
  varnb <- vcovHC(nm)[2,2]
  varna <- vcovHC(nm)[1,1]
  covnanb <- vcovHC(nm)[1,2]
  covm <- mefit$vcov
  if(names(coef(mefit)[2]) == mefit$differential){ #change order of covm if needed
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
      t10 <- ifelse(names(coef(mefit)[2]) == mefit$differential, coef(mefit)[3], coef(mefit)[2])
      t01 <- coef(mefit)[mefit$differential] + t00
      t11 <- coef(mefit)[4] + t10}
    else {
      ids <- sample(1:nrow(mefit$model), replace = T)
      calsample <- mefit$model[ids,]
      calfit <- lm(formula = formula(mefit), data = calsample)
      t00 <- coef(calfit)[1]
      t10 <- ifelse(names(coef(calfit)[2]) == mefit$differential, coef(calfit)[3], coef(calfit)[2])
      t01 <- coef(calfit)[mefit$differential] + t00
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

