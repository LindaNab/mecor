#' mecor: a measurement error correction package
#'
#' mecor provides correction methods for measurement
#' errors in continuous endpoints.
#'
#' @param formula an object of class \link[stats]{formula} (or one that is
#' coerced to that class): a symbolic description of the naive
#' model to be fitted in an \link[stats]{lm} model.
#' @param data an optional data frame, list or environment (or
#' object coercible by as.data.frame to a data frame) containing
#' the variables in the model. If not found in \code{data}, the
#' variables are taken from \code{environment(formula)}, typically
#' the enviroment from which \code{mecor} is called.
#' @param cv a non-empty character string specifying the variable
#' in \code{formula} to correct
#' @param range (a vector of) character string(s) in the format "a:b" specifying
#' the range for which \code{cv} needs to be corrected
#' @param correction (a list of) object(s) of class \link[mecor]{systematic}
#' used to correct \code{cv} for corresponding \code{range}
#'
#' @return \code{mecor} returns an object of \link[base]{class} "mecor"
#'
#' An object of class \code{mecor} is a list containing the following components:
#'
#' \item{model}{an \code{lm} object of the corrected model}
#' \item{ci}{a named vector with the confidence intervals of the effect
#' estimate of the corrected model calculated with the zero variance,
#' delta and fieller method. If fieller method produces \code{NA} the fieller
#' confidence interval is undefined.}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @references
#'
#' @examples
#' X <- c(rep(0,1000),rep(1,1000))
#' Y <- X + rnorm(2000,0,3)
#' real_model <- lm(Y ~ X)
#'
#' ##external calibration set
#' Xcal <- c(rep(0,500),rep(1,500))
#' Ycal <- Xcal + rnorm(1000,0,3)
#'
#' ##systematic measurement error (sme)
#' V <- 1 + 2 * Y + rnorm(2000,0,3)
#' naive_model_sme <- lm(V ~ X)
#' Vcal <- 1 + 2 * Ycal + rnorm(1000,0,3)
#' syst <- systematic(formula = Vcal ~ Ycal)
#' corrected_model_sme <- mecor(formula = V ~ X, cv = "V", correction = syst)
#'
#' ##differential measurement error (dme)
#' V_dme0 <- 1 + 2 * Y[1:1000] + rnorm(Y[1:1000], 0, 3)
#' V_dme1 <- 2 + 3 * Y[1001:2000] + rnorm(Y[1001:2000], 0, 3)
#' V_dme <- c(V_dme0, V_dme1)
#' naive_model_dme <- lm(V_dme ~ X)
#' Vcal_dme0 <- 1 + 2 * Ycal[1:500] + rnorm(Ycal[1:500], 0, 3)
#' Vcal_dme1 <- 2 + 3 * Ycal[501:1000] + rnorm(Ycal[501:1000], 0, 3)
#' syst1 <- systematic(formula = Vcal_dme0 ~ Ycal[1:500])
#' syst2 <- systematic(formula = Vcal_dme1 ~ Ycal[501:1000])
#' corrected_model_dme <- mecor(formula = V_dme ~ X, cv = "V_dme", range = c("1:1000", "1001:2000"), correction = list(syst1, syst2))
#'
#' @export
mecor <- function(formula,
                  data = NULL,
                  cv,
                  range = NULL,
                  correction){
  if(attr(terms(formula), "variables")[4]!="NULL()"){
    stop("variable 'formula' should be a formula describing one independent and one dependent variable")}
  if(!is.character(class(cv))){
    stop("variable 'cv' should be a character string")}
  if(names(attr(terms(formula),"factors")[,1])[1]!= cv){
    stop("variable 'cv' should be the dependent variable")}
  lc <- list()
  if(class(correction)=="systematic"){#check whether input variable 'correction' is a list, if not: make list
    lc[[1]] <- correction}
  else {lc <- correction}
  model_naive <- lm(formula, data)
  adjustedvars <- data.frame(dep = model_naive$model[,cv], indep = model_naive$model[,names(attr(terms(formula),"factors")[,1])[2]])
  for(i in 1:length(lc)){
    if(class(lc[[i]]) != 'systematic'){
      stop("'correction' should be of class systematic")}
    if(is.null(range)){ #no range provided for cv, set range to 1:NROW(cv)
      range <- paste("1",toString(NROW(model_naive$model[,cv[[i]][1]])),sep=":")}
    if(length(range)!=length(lc)){
      stop("each variable 'range' should be provided with a corresponding correction object. If no correction needed for that range, omit that range in variable 'range'.")}
    seq <- seq(strsplit(range[i], ":")[[1]][1], strsplit(range[i], ":")[[1]][2]) #extract sequence from character string range
    adjustedvars$dep[seq] <- (model_naive$model[seq,cv] - lc[[i]]$coefficients["theta0_hat"]) /
                      lc[[i]]$coefficients["theta1_hat"] #adjust variables for given range
    if(i > 1 && lc[[i]]$alpha != lc[[(i-1)]]$alpha){
      stop("alpha levels of all systematic objects should be equal")}
    else {alpha <- lc[[1]]$alpha}
    }
  model_corrected <- lm(dep ~ indep, data=adjustedvars)
  s_xx <- sum((adjustedvars$indep - mean(adjustedvars$indep))^2)
  t_q <- qt((1 - alpha / 2), (NROW(X) - 2))
  ci.zv <- c(lower = unname(model_corrected$coefficients[2] - t_q * sqrt(summary(model_corrected)$sigma ^ 2 / s_xx)),
              upper = unname(model_corrected$coefficients[2] + t_q * sqrt(summary(model_corrected)$sigma ^ 2 / s_xx)))
  if(class(correction)=="systematic"){ #only if error is systematic and not differential
    ci.delta <- delta(model_naive, model_corrected, correction)
    ci.fieller <- fieller(model_naive, model_corrected, correction)}
  else {ci.delta <- NA
  ci.fieller <- NA}
  cibeta <- c(zerovariance = ci.zv,
                  delta = ci.delta,
                  fieller = ci.fieller)

  out <- list(model = model_corrected,
                ci = cibeta
    )
  class(out) <- 'mecor'
  return(out)
}


