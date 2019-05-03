#' @title Fitting Measurement Error Model
#'
#' @description
#' This function fits several measurement error models to model the relationship between the
#' test measurement (error-prone measurement) and the reference measurement (true measurement)
#' in a validation set.
#'
#' @param MeasError object of class MeasError
#' @param differential vector containing the differential measurement error variable
#' @param model optional character string or vector of character strings indicating
#' the assumed measurement error model: "all" (default), "cme" (classical measurement
#' error), "sme1" (systematic measurement error with zero intercept), "sme2" (systematic
#' measurement errror with non-zero intercept) and "dme" (differential measurement error).
#' @param robust boolean indicating whether robust standard errors need to be calculated, the
#' "HC3" robust standard errors from \link[sandwich]{vcovHC} are used for the
#' heteroskedasticity-consistent estimation of the covariance matrix of the coefficient
#' estimates.
#' @param plot boolean indicating whether one wants plots of the residuals (one plot for each
#' tested model)
#'
#' @return \code{mefit} returns an object of \link[base]{class} "mefit".
#'
#' An object of class \code{mefit} is a list containing, depending on the models tested for, the following components:
#' \item{cme}{a list containing the coefficients, sigma, df and model of the classical
#' measurement error model}
#' \item{sme1}{a list containing the coefficients, sigma, df and model of the systematic
#' measurement error model with zero intercept}
#' \item{sme2}{a list containing the coefficients, sigma, df and model of the systematic
#' measurement error model with non-zero intercept}
#' \item{dme}{a list containing the coefficients, sigma, df and model of the differential
#' measurement error model}
#' \item{lrtest1}{a list containing the likelihood ratio test of, depending on which models are tested, cme vs sme1 vs sme2;
#' cme vs sme1; cme vs sme2; sme1 vs sme2}
#' \item{lrtest2}{a list containing the likelihood ratio test of, depending on whether one tests for differential measurement
#' error, sme2 vs dme}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @examples
#' ##measurement error in endpoint
#' X <- c(rep(0, 50), rep(1, 50))
#' Y <- X + rnorm(100, 0, 1)
#' Vcme <- Y + rnorm(100, 0, 3) #classical measurement error (cme)
#' Vsme <- 1 + 2*Y + rnorm(100, 0, 3) #systematic measurement error (sme)
#' Vdme <- 2 + 2*X + 3*Y + 2*X*Y + rnorm(10, 0, 3*(1-X) + 2*X) #systematic differential measurement error (dme)
#'
#' ##measurement error in exposure
#' X <- rnorm(100, 0, 1)
#' Y <- 1.5*X + rnorm(100, 0, 0.5)
#' W1 <- X + rnorm(100, 0, 0.5)
#' W2 <- X + rnorm(100, 0, 0.5)
#'
#' fit1 <- mefit(MeasError(Vcme, Y), plot = T)
#' fit2 <- mefit(MeasError(Vsme, Y), plot = T)
#' fit3 <- mefit(MeasError(Vdme, Y), X, robust = T, plot = T)
#' fit4 <- mefit(MeasError(W1, X), model = "all")
#' #fit5 <- mefit(MeasError(cbind(W1, W2), NA)) #not yet supported
#'
#' @importFrom sandwich vcovHC
#' @importFrom lmtest lrtest
#' @export
mefit <- function(MeasError,
                  differential,
                  model = "all",
                  robust = FALSE,
                  plot = TRUE){
  if(missing(differential)) differential <- NULL
  if(any(class(MeasError) != c("MeasError", "data.frame")))
    stop("MeasError is not a MeasError object")
  if(model == "all"){
    if(missing(differential)) model <- c("cme", "sme1", "sme2")
    else model <- c("cme", "sme1", "sme2", "dme")}
  reference <- MeasError$reference
  if(is.matrix(MeasError$test)){
    stop("replicate measurements are not yet supported in mefit")
    if(!all(is.na(reference))) stop("if test contains multiple measurements, reference should be NA")
    else if(all(model == c("cme", "sme1", "sme2", "dme"))){
      model = "cme"
      warning("as there are multiple measurements, the classical measurement error model is assumed")}}
  else {test <- MeasError$test}
  opar <- par(no.readonly = T)
  par(mfrow = c(2,2))
  out <- list()
  #classical measurement error
  if(any(model == "cme")){
    fitcme <- stats::lm(test ~ offset(reference) -1)
    #t <- t.test(test, reference, paired = T)
    #coef <- cbind("Estimate" = unname(t$estimate),
                  #"Std. Error" = sd(test-reference)/sqrt(NROW(test)),
                  #"t value" = unname(t$statistic), "Pr(>|t|)" = t$p.value)
    cme <- list(coefficients = "no coefficients", Sigma = summary(fitcme)$sigma,
                df = summary(fitcme)$df[2], model = "classical")
    if(plot == TRUE){
      plot(fitcme$model$`offset(reference)`, fitcme$residuals,
           xlab = "reference", ylab = "residuals", main = "Classical")
      abline(a = 0, b = 0)}
    out$cme <- cme}
  #systematic measurement error
  if(any(model == "sme1")){
    fitsme1 <- stats::lm(test ~ reference - 1)
    if(robust == T){
      sandwich_se1 <- diag(sandwich::vcovHC(fitsme1))^0.5
      t_stat1 <- coef(fitsme1)/sandwich_se1
      coef1 <- cbind("Estimate" = coef(fitsme1), "Std. Error" = sandwich_se1,
                     "t value" = t_stat1, "Pr(>|t|)" = stats::pchisq(t_stat1^2, 1, lower.tail=FALSE) )}
    else{
      coef1 <- coef(summary(fitsme1))}
    sme1 <- list(coefficients = coef1, Sigma = summary(fitsme1)$sigma,
                 df = summary(fitsme1)$df[2], model = "systematic zero intercept")
    if(plot == TRUE){
      plot(fitsme1$model$reference, fitsme1$residuals,
           xlab = "reference", ylab = "residuals", main = "Sys zero int")
      abline(a = 0, b = 0)}
    out$sme1 <- sme1}
  if(any(model == "sme2")){
    fitsme2 <- stats::lm(test ~ reference)
    if(robust == T){
      sandwich_se2 <- diag(sandwich::vcovHC(fitsme1))^0.5
      t_stat2 <- coef(fitsme2)/sandwich_se2
      coef2 <- cbind("Estimate" = coef(fitsme2), "Std. Error" = sandwich_se2,
                     "t value" = t_stat2, "Pr(>|t|)" = stats::pchisq(t_stat2^2, 1, lower.tail=FALSE) )}
    else{
      coef2 <- coef(summary(fitsme2))}
    sme2 <- list(coefficients = coef2, Sigma = summary(fitsme2)$sigma,
                 df = summary(fitsme2)$df[2], model = "systematic non-zero intercept")
    if(plot == TRUE){
      plot(fitsme2$model$reference, fitsme2$residuals,
           xlab = "reference", ylab = "residuals", main = "Sys non-zero int")
      abline(a = 0 , b = 0)}
    out$sme2 <- sme2}
  if(exists("fitcme")){
    if(exists("fitsme1") & exists("fitsme2")){
      lrtest1 <- lmtest:::lrtest(fitcme, fitsme1, fitsme2)}
    else if(exists("fitsme1")){
      lrtest1 <- lmtest:::lrtest(fitcme, fitsme1)}
    else if(exists("fitsme2")){
      lrtest1 <- lmtest:::lrtest(fitcme, fitsme2)}
  }
  else if(exists("fitsme1" & exists("fitsme2"))){
    lrtest1 <- lmtest:::lrtest(fitsme1, fitsme2)
  }
  else lrtest1 <- NULL
  #differential measurement error
  if(any(model == "dme")){
    if(!is.null(differential)){
      fitdme <- stats::lm(test ~ reference * differential)
      if(exists("fitdme") & exists("fitsme2")){
        lrtest2 <- lmtest::lrtest(fitsme2, fitdme)}
      if(robust == T){
        sandwich_se <- diag(vcovHC(fitdme))^0.5
        t_stat <- coef(fitdme)/sandwich_se
        coef <- cbind("Estimate" = coef(fitdme), "Std. Error" = sandwich_se,
                       "t value" = t_stat, "Pr(>|t|)" = stats::pchisq(t_stat^2, 1, lower.tail=FALSE) )}
      else coef <- coef(summary(fitdme))
      dme <- list(coefficients = coef, Sigma = summary(fitdme)$sigma,
                df = summary(fitdme)$df[2], model = paste("differential on", deparse(substitute(differential))))
      if(plot == TRUE){
        plot(fitdme$model$reference, fitdme$residuals,
             xlab = "reference", ylab = "residuals", main = "Differential")
        abline(a = 0, b = 0)}}
    else dme <- NULL
    out$dme <- dme}
  if(!exists("lrtest2")){
    lrtest2 <- NULL}
  out$lrtest1 <- lrtest1
  out$lrtest2 <- lrtest2
  if(plot == TRUE) par(opar)
  attr(out, "call") <- match.call()
  attr(out, "robust") <- robust
  attr(out, "size") <- NROW(MeasError)
  class(out) <- c("mefit")
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
      factor <- deparse(t[[3L]])
      f <- update(f, ~ . * eval.parent(factor))
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

fieller2 <- function(nm,
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

