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
#' data2 <- data.frame(Z, W, W2, Y)
#'
#' mecor(Y ~ MeasError(W, X) + Z, data)
#' mecor(Y ~ MeasError(W, X) + Z, data, method = "rc_pooled1")$corfit
#' mecor(Y ~ MeasError(W, X) + Z, data, method = "rc_pooled2")$corfit
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

  if(mevar == "indep"){
    mlist <- mecor:::rcm(vars, me)
    naivefit <- stats::lm.fit(mlist$x, mlist$y)
    if(method == "rc"){
      resrc <- mecor:::regcal(mlist, naivefit, B = B, alpha = alpha)
      corfit <- resrc$corfit
      corvar <- resrc$corvar
      ci.fieller <- resrc$ci.fieller
      if(B!=0) ci.b <- resrc$ci.b}
    if(method == "rc_pooled1"){
      resrc_p <- mecor:::regcal_pooled(mlist, naivefit, pooled.var = "delta", B = B, alpha = alpha)
      corfit <- resrc_p$corfit
      corvar <- resrc_p$var
      if(B!=0) ci.b <- resrc_p$ci.b
    }
    if(method == "rc_pooled2"){
      resrc_p <- mecor:::regcal_pooled(mlist, naivefit, pooled.var = "bootstrap", B = B, alpha = alpha)
      corfit <- resrc_p$corfit
      corvar <- resrc_p$var
      if(B!=0) ci.b <- resrc_p$ci.b
    }
  }
  if(mevar == "dep"){
    y <- me$test
    x <- cbind(1, vars[,2:ncol(vars)]) #design matrix
    lc <- l[-indx]
    colnames(x) <- c("(Intercept)", lc)}

  if({vtp <- attributes(me)$type} == "internal"){
    #xint <- cbind(1, me$reference)
    #if(ncol(vars) > 1) xint <- cbind(xint, vars[,2:ncol(vars)]) #design matrix only with internal valdata
    #if(mevar == "indep"){
      #y <- vars[,1] #vector containing the outcomes
      #x <- cbind(1, me$test)
      #if(ncol(vars) > 1) x <- cbind(x, vars[,2:ncol(vars)]) #design matrix with measurement error
      #lc <- l[-c(1,indx)]
      #cnx <- c("(Intercept)", attributes(me)$input$test)
      #if(length(lc)!=0) cnx <- c(cnx, lc)
      #colnames(x) <- cnx
      #cnxint <- c("(Intercept)", attributes(me)$input$reference)
      #if(length(lc)!=0) cnxint <- c(cnxint, lc)
      #colnames(xint) <- cnxint}
    if(mevar == "dep"){
      y <- me$test
      x <- cbind(1, vars[,2:ncol(vars)]) #design matrix
      lc <- l[-indx]
      colnames(x) <- c("(Intercept)", lc)}
  }

  if(vtp == "replicate"){
    #if(mevar == "indep"){
      #y <- vars[,1]
      #x <- cbind(1, me$test$test1)
      #if(ncol(vars) > 1) x <- cbind(x, vars[,2:ncol(vars)]) #design matrix with measurement error
      #lc <- l[-c(1,indx)]
      #cnx <- c("(Intercept)", attributes(me)$input$test$test1)
      #if(length(lc)!=0) cnx <- c(cnx, lc)
      #colnames(x) <- cnx}
    if(mevar == "dep") stop("mecor is currently not developed to correct for measurement error in
              the dependent var using replicate measurements")
  }

  # if(vtp == "internal"){ #internal validation data
  #   if(mevar == "indep"){ #measurement error in the independent variable
  #     naivefit <- stats::lm.fit(mlist$x, mlist$y)
  #     if(method == "rc" | method == "rc.pooled"){
  #       corfit <- rc(mlist, naivefit)
  #       calfit <- stats::lm.fit(x[!is.na(me$reference),], me$reference[!is.na(me$reference)])
  #       corx <- x
  #       e <- calfit$coef%*%t(x) #predicted values
  #       corx[,2] <- e
  #       colnames(corx)[2] <- as.character(attributes(me)$input$reference)
  #       corfit <- stats::lm.fit(corx, y)
  #       if(method == "rc"){
  #         deltavar <- mecor:::vardelta(naivefit, calfit, names(corfit$coef))
  #         corvar <- list(deltavar = deltavar)
  #         #ci.fieller <- mecor:::fieller(naivefit, calfit, alpha)
  #         if(B != 0){
  #           bd <- data.frame(y, me$reference, x)
  #           colnames(bd) <- c("Y", as.character(attributes(me)$input$reference), colnames(x))
  #           ci.b <- mecor:::boot_rc(bd, me, alpha, B)
  #         }
  #       }
  #       if(method == "rc.pooled"){
  #         intfit <- stats::lm.fit(xint[!is.na(me$reference),], y[!is.na(me$reference)])
  #         intvar <- diag(mecor:::vcovfromfit(intfit))
  #         wrc <- 1/deltavar * (1/(1/deltavar + 1/intvar))
  #         corfit.rcp <- wrc * corfit$coef + (1 - wrc) * intfit$coef
  #         if(B != 0){
  #         }
  #       }
  #     }
  #   }
  #   if(mevar == "dep"){
  #     if(method == "rc"){
  #     }
  #   }
  # }

  # #replicate validation data
  # if(vtp == "replicate"){
  #   if(mevar == "indep"){
  #     naivefit <- stats::lm.fit(x, y)
  #     if(method == "rc"){
  #       calfit <- stats::lm.fit(x[!is.na(me$test$test2),], me$test$test2[!is.na(me$test$test2)])
  #       corx <- x
  #       e <- calfit$coef%*%t(x) #predicted values
  #       corx[,2] <- e
  #       colnames(corx)[2] <- "cor"
  #       corfit <- stats::lm.fit(corx, y)
  #       deltavar <- mecor:::vardelta(naivefit, calfit, names(corfit$coef))
  #       corvar <- list(deltavar = deltavar)
  #       if(B != 0){
  #         bd <- data.frame(y, me$test$test2, x)
  #         colnames(bd) <- c("Y", as.character(attributes(me)$input$test$test2), colnames(x))
  #         ci.b <- mecor:::boot_rep_rc(bd, me, alpha, B)
  #       }
  #     }
  #   }
  # }

  #MECORS output
  out <- list(naivefit = naivefit,
              corfit = corfit,
              corvar = {if(exists("corvar")) corvar else NA},
              ci.fieller = {if(exists("ci.fieller")) ci.fieller else NA},
              ci.b = {if(B != 0) ci.b[,1:2] else NA}
              )
  class(out) <- 'mecor'
  attr(out, "call") <- match.call()
  attr(out, "B") <- B
  attr(out, "alpha") <- alpha
  return(out)
}


