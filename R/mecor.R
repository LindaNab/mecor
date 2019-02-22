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
#' @param true.var a non empty character string specifying the variable
#' that corresponds to the correctly measured me.var. Required if data.type
#' in mefit object is 'internal'.
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
#' cm <- mecor(formula = Y ~ MeasError(W1, X) + Z, rc())
#' X ~ W + Z
#' cm <- mecor(Y ~ MeasError(W, X) + Z, mime())
#'
#' ##data generation
#' X <- c(rep(0, 1000), rep(1, 1000))
#' Y <- X + rnorm(2000, 0, 1)
#' V_sme <- 1 + 2 * Y + rnorm(2000, 0, 3) #systematic measurement error (sme)
#' V_dme <- 2 + 2 * X + 3 * Y + 2 * X * Y + rnorm(2000, 0, 3 * (1 - X) + 2 * X) #differential measurement error (dme)
#' rm <- lm(Y ~ X) #real model (rm)
#'
#' ##solve systematic measurement error (sme)
#' nm_sme <- lm(V_sme ~ X) #compare naive model (nm) with rm
#' cm_sme <- mecor(MeasError(Vsme, NA) ~ X, method = rc(), B = 999) #compare with nm and rm
#'
#'
#' ##solve differential measurement error (dme)
#' #nm_dme <- lm(V_dme ~ X) ##compare with rm
#' #fit_dme <- mefit(formula = Vcal_dme ~ Ycal * Xcal, data = caldata, me.structure = "differential", dif.var = "Xcal", robust = TRUE)
#' #cm_dme <- mecor(formula = V_dme ~ X, data = data, me.var = "V_dme", mefit = fit_dme, dif.var = c("X" = "Xcal"), method = "rc", robust = T, B = 999)
#'
#' @import boot
#' @export
mecor <- function(formula,
                  data,
                  method = rc(),
                  robust = FALSE,
                  alpha = 0.05,
                  B = 0){
  if(missing(data)) stop("data not found")#data = NULL
  else if(!is.data.frame(data)) data <- as.data.frame(data)
  if(missing(formula)) stop("formula not found")

  l <- as.list(attr(terms(formula), "variables"))[-1]
  indx <- grep("MeasError", l)
  if(indx == 1) type <- "outcome"
  else type <- "expl"
  me <- l[indx]
  if(length(me) == 0){
    stop("formula should contain a MeasError object")}
  else if (length(me) != 1){
    stop("formula can only contain one MeasError object")}
  me <- me[[1]]
  evalme <- eval.parent(me)
  test <- evalme$test
  reference <- evalme$reference
  if(is.null(reference)){
    warning("reference is null, so mecor assumes external validation data")
    valdata <- "external"}
  else valdata <- "internal"

  #internal validation data
  if(valdata == "internal"){
    if(type == "expl"){
      pf <- update(formula, bquote(. ~ . - .(me))) #plain
      nf <- update(pf, . ~ . + test) #naive
      nm <- lm(nf, data)
      if(method == rc()){
        #rc.int(formula, data, )
        cf <- update(nf, reference ~ .) #cal
        cm <- lm(cf, data)
        calval <- predict(cm, data)
        cf <- update(pf, . ~ . + calval)
        lm(cf)
        calval[!is.na(reference)] <- reference[!is.na(reference)]
        lm(cf)
      }
    }
    if(type == "outcome"){
      nf <- update(formula, bquote(test ~ .))
      nm <- lm(nf, data)
      if(method == rc()){
      }
    }
  }

  #external validation data
  if(valdata == "external"){
  }

  #MECORS output
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

rc <- function(formula, data, keep.original = TRUE){
}

