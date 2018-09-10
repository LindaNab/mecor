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
#' @param me.var a non-empty character string specifying the variable
#' in \code{formula} with measurement error
#' @param mefit object of class \link[mecor]{mefit}
#' used to correct \code{me.var}
#' @param dif.var an optional named vector specifying the grouping variable in \code{formula}
#' that corresponds to the dif.var used in 'mefit', which is only needed if the structure
#' of 'mefit' is 'differential'. If 'dif.var' is not specified and me.structure is
#' 'differential, mecor searches the environment for the dif.var variable specified in 'mefit'
#' object.
#' @param method a character string indicating the method used to correct for
#' measurement error, only the method "rc" (regression calibration) is implemented
#' @param alpha alpha level used to construct confidence intervals
#' @param B number of bootstrap samples
#'
#' @return \code{mecor} returns an object of \link[base]{class} "mecor"
#'
#' An object of class \code{mecor} is a list containing the following components:
#'
#' \item{coefficients}{a named vector containing the coefficients of the corrected model}
#' \item{stderr}{zero variance standard errors of coefficients of the corrected model}
#' \item{coefficients.nm}{a named vector containing the coefficients of the naive (uncorrected) model}
#' \item{rdf}{the residual degrees of freedom of the corrected model}
#' \item{call}{matched call}
#' \item{ci}{a named matrix containing the confidence intervals for the effect estimate of the corrected
#' model}
#' \item{dif.var}{a named vector linking the grouping variable in formula to the grouping variable in
#' the used mefit object (used for differential measurement error models)}
#' \item{Rbtstrp}{number of bootstrap replicates used to construct Bootstrap Confidence Interval}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @references
#' L.Nab, R.H.H. Groenwold, P.M.J. Welsing, M. van Smeden.
#' Measurement error in continuous endpoints in randomised trials: an exploration of problems and solutions
#'
#' @examples
#' ##data generation
#' X <- c(rep(0, 1000), rep(1, 1000))
#' Y <- X + rnorm(2000, 0, 1)
#' V_sme <- 1 + 2 * Y + rnorm(2000, 0, 3) #systematic measurement error (sme)
#' V_dme <- 2 + 2 * X + 3 * Y + 2 * X * Y + rnorm(2000, 0, 3 * (1 - X) + 2 * X) #differential measurement error (dme)
#' rm <- lm(Y ~ X) #real model (rm)
#'
#' ##generation of external calibration set
#' Xcal <- c(rep(0, 500), rep(1, 500))
#' Ycal <- Xcal + rnorm(1000, 0, 1)
#' Vcal_sme <- 1 + 2 * Ycal + rnorm(1000, 0, 3)
#' Vcal_dme <- 2 + 2 * Xcal + 3 * Ycal + 2 * Xcal * Ycal + rnorm(1000, 0, 3 * (1 - Xcal) + 2 * Xcal)
#'
#' ##solve systematic measurement error (sme)
#' nm_sme <- lm(V_sme ~ X) #compare naive model (nm) with rm
#' fit_sme <- mefit(formula = Vcal_sme ~ Ycal, me.structure = "systematic")
#' cm_sme <- mecor(formula = V_sme ~ X, me.var = "V_sme", mefit = fit_sme, method = "rc", B = 999) #compare with nm and rm
#'
#' ##solve differential measurement error (dme)
#' nm_dme <- lm(V_dme ~ X) ##compare with rm
#' fit_dme <- mefit(formula = Vcal_dme ~ Ycal * Xcal, data = caldata, me.structure = "differential", dif.var = "Xcal", robust = TRUE)
#' cm_dme <- mecor(formula = V_dme ~ X, data = data, me.var = "V_dme", mefit = fit_dme, dif.var = c("X" = "Xcal"), method = "rc", robust = T, B = 999)
#'
#' @export
mecor <- function(formula,
                  data,
                  me.var,
                  mefit,
                  dif.var,
                  method = "rc",
                  robust = FALSE,
                  alpha = 0.05,
                  B = 0){
  if(missing(data)) data = NULL
  else if(!is.data.frame(data)) data <- as.data.frame(data)
  if(missing(dif.var) && mefit$me.structure == "differential"){
    warning("dif.var is not specified, so dif.var specified in 'mefit' is used")
    if(!any(grepl(mefit$dif.var, formula))) stop("dif.var specified in 'mefit' does not occur in formula")
    dif.var <- mefit$dif.var
    names(dif.var) <- mefit$dif.var}
  else if(missing(dif.var)) dif.var = NULL
  if(NROW(all.vars(formula)) != 2){
    stop("variable 'formula' should be a formula describing one independent and one dependent variable")}
  if(!is.character(class(me.var))){
    stop("variable 'me.var' should be a character string")}
  if(names(attr(terms(formula),"factors")[,1])[1] == me.var){
    me.var <- cbind(me.var, "dep")}
    else stop("variable 'me.var' should be the dependent variable")
  if(class(mefit) != "mefit"){
    stop("variable 'mefit' should be of class 'mefit'")}
  if(mefit$me.structure == "classical" && me.var[2] == "dep"){
    stop("me.structure of 'mefit' is classical and 'me.var' is the dependent variable so there is nothing to correct")}
  if(mefit$me.structure == "differential") {
    if(!is.null(data)) {getdif.var <- data[names(dif.var)]} else getdif.var <- get(names(dif.var))
    if(!is.numeric(unique(getdif.var))) diflevels <- as.numeric(unlist(unique(getdif.var))) else diflevels <- unique(getdif.var)
    if(!all.equal(mefit$diflevels, diflevels, check.attributes = FALSE)){
    stop("levels of dif.var in mefit object differ from the levels of the grouping variable specified in formula")}
  }
  nm <- lm(formula, data)
  coef.nm <- summary(nm)$coef
  if(robust == TRUE){
    vcov <- vcovHC(nm) }
  else vcov <- vcov(nm)
  ci.cm <- matrix(data = NA, nrow = 2L, ncol = 2L,
                  dimnames = list(c('Zero Variance (ZV)', 'Delta'), c('Lower', 'Upper')))
  if(mefit$me.structure == "systematic" && me.var[2] == "dep"){
    t0 <- unname(coef(mefit)[1])
    t1 <- unname(coef(mefit)[2])
    int <- (unname(coef.nm[1,1]) - t0) / t1
    slope <- unname(coef.nm[2,1]) / t1
    coef.cm <- c('(Intercept)' = int, 'X' = slope)
    tq <- qt((1 - alpha / 2), nm$df.residual)
    stderr.cm <- coef.nm[,2] / t1 ^ 2
    zv.l <- coef.cm[2] - tq * stderr.cm[2]
    zv.u <- coef.cm[2] + tq * stderr.cm[2]
    ci.cm[1,] <- c(zv.l, zv.u)
    ci.cm[2,] <- delta.sme(nm, coef.cm, mefit, alpha)
    ci.cm <- rbind(ci.cm, 'Fieller' = fieller(nm, mefit, alpha))
    if(B != 0){
      bt <- bootstrap.sme(nm, mefit, alpha, B)
      ci.cm <- rbind(ci.cm, 'Bootstrap'= bt$percent[4:5])}
  }
  if(mefit$me.structure == "differential" && me.var[2] == "dep"){
    t00 <- unname(coef(mefit)[1])
    t10 <- ifelse(names(coef(mefit)[2]) == mefit$dif.var, unname(coef(mefit)[3]), unname(coef(mefit)[2]))
    t01 <- unname(coef(mefit)[mefit$dif.var]) + t00
    t11 <- unname(coef(mefit)[4]) + t10
    int <- (unname(coef.nm[1,1]) - t00) / t10
    slope <- (unname(coef.nm[2,1]) + unname(coef.nm[1,1]) - t01) / t11 - int
    coef.cm <- c('(Intercept)' = int, 'X' = slope)
    stderr.cm <- c(sqrt(vcov[1,1] / t10^2),
                   sqrt(vcov[2,2]/ t11^2 - vcov[1,1] / t11^2 + vcov[1,1] / t10^2))
    tq <- qt((1 - alpha / 2), nm$df.residual)
    zv.l <- coef.cm[2] - tq * stderr.cm[2]
    zv.u <- coef.cm[2] + tq * stderr.cm[2]
    ci.cm[1,] <- c(zv.l, zv.u)
    ci.cm[2,] <- delta.dme(nm, coef.cm, mefit, alpha)
    if(B != 0){
      bt <- bootstrap.dme(nm, mefit, alpha, B)
      ci.cm <- rbind(ci.cm, 'Bootstrap'= bt$percent[4:5])}
  }
  out <- list(coefficients = coef.cm,
              stderr = stderr.cm,
              coefficients.nm = coef.nm,
              rdf = nm$df.residual,
              call = match.call(),
              ci = ci.cm,
              dif.var = dif.var,
              Rbtstrp = ifelse(exists("bt"), bt$R, NA))
  class(out) <- 'mecor'
  return(out)
}


