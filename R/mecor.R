#' mecor: a measurement error correction package
#'
#' mecor provides correction methods for measurement
#' errors in a continuous covariate.
#'
#' @param formula an object of class \link[stats]{formula} (or one that is
#' coerced to that class): a symbolic description of the model containing
#' a \link[mecor]{MeasError} object.
#' @param data a data.frame, list or environment (or
#' object coercible by as.data.frame to a data frame) containing
#' the variables in the model
#' @param method a character string indicating the method used to correct for
#' measurement error, either "rc" (regression calibration), "rc_pooled1" (efficient
#' regression calibration using delta variance for pooling) or "rc_pooled2" (efficient
#' regression calibration using bootstrap variance for pooling).
#' @param alpha alpha level used to construct confidence intervals
#' @param B number of bootstrap samples
#'
#' @return \code{mecor} returns an object of \link[base]{class} "mecor"
#'
#' An object of class \code{mecor} is a list containing the following components:
#'
#' \item{naivefit}{a lm.fit object of the uncorrected fit}
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
#' ##data generation
#' #measurement error in exposure
#' nobs <- 1e3
#' Z <- rnorm(nobs, 0, 1)
#' X <- Z + rnorm(nobs, 0, 1)
#' Y <- 0.5 * X + 2 * Z + rnorm(nobs, 0, 1)
#' W <- X + rnorm(nobs, 0, 0.5)
#' X <- ifelse(rbinom(nobs, 0, 0.9) == 1, NA, X)
#' data <- data.frame(Z, X, W, Y)
#' W2 <- X + rnorm(nobs, 0, 0.5)
#' W2 <- ifelse(rbinom(nobs, 0, 0.8) == 1, NA, W2)
#' data2 <- data.frame(Z, W, W2, Y)
#'
#' mecor(Y ~ MeasError(W, X) + Z, data)
#' mecor(Y ~ MeasError(W, X) + Z, data, method = "rc_pooled1")
#' mecor(Y ~ MeasError(W, X) + Z, data, method = "rc_pooled2")
#' mecor(Y ~ MeasError(cbind(W, W2), NA) + Z, data2)
#' @import boot
#' @export
mecor <- function(formula,
                  data,
                  method = "rc",
                  alpha = 0.05,
                  B = 0){
  if(missing(data)) stop("data not found") #data = NULL
  else if(!is.data.frame(data)) data <- as.data.frame(data)
  if(missing(formula)) stop("formula not found")
  if(! method %in% c("rc", "rc_pooled1", "rc_pooled2")) stop("this method is not implemented")

  #Create MeasError object
  l <- as.list(attr(terms(formula), "variables"))[-1]
  indx <- grep("MeasError", l)
  if(length(indx) == 0){
    stop("formula should contain a MeasError object")}
  else if(length(indx) != 1){
    stop("formula can only contain one MeasError object")}
  if(indx == 1) mevar <- "dep"
  else mevar <- "indep"

  temp <- sapply(l, eval, envir = data)
  me <- temp[[indx]]
  vars <- sapply(temp[-indx], cbind)
  colnames(vars) <- l[-indx]

  if(mevar == "indep" & {vtp <- attributes(me)$type} == "internal"){
    mlist <- mecor:::rcm(vars, me)
    naivefit <- stats::lm.fit(mlist$x, mlist$y)
    if(method == "rc"){
      res <- mecor:::regcal(mlist, naivefit, B = B, alpha = alpha)}
    if(method == "rc_pooled1"){
      res <- mecor:::regcal_pooled(mlist, naivefit, pooled.var = "delta", B = B, alpha = alpha)}
    if(method == "rc_pooled2"){ #uses B = 999 for the bootvar used to pool the estimates
      res <- mecor:::regcal_pooled(mlist, naivefit, pooled.var = "bootstrap", B = B, alpha = alpha)}
  }
  else if(mevar == "indep" && vtp == "replicate"){
    mlist <- mecor:::rcm(vars, me)
    naivefit <- stats::lm.fit(mlist$x, mlist$y)
    if(method == "rc"){
      res <- mecor:::regcal(mlist, naivefit, B = B, alpha = alpha)}
    if(method == "rc_pooled1" | method == "rc_pooled2"){
      stop("mecor is currently not able to do a pooled regression calibration
           in case of replicate data")
    }
  }
  else if(mevar == "dep"){
    y <- me$test
    x <- cbind(1, vars[,2:ncol(vars)]) #design matrix
    lc <- l[-indx]
    colnames(x) <- c("(Intercept)", lc)
    stop("mecor cannot correct for measurment error in the dependent variable")}

  #MECORS output
  out <- list(naivefit = naivefit,
              corfit = res$corfit,
              corvar = res$corvar,
              ci = res$ci
              )
  class(out) <- 'mecor'
  attr(out, "call") <- match.call()
  attr(out, "B") <- B
  attr(out, "alpha") <- alpha
  return(out)
}


