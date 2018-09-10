#' @title Summarizing Measurement Error Correction
#'
#' @description
#' \code{summary} method for class "mecor"
#'
#' @param object an object of class "mecor", a result of a call to \link[mecor]{mefit}
#'
#' @return
#' The function \code{summary.mecor} returns a list of summary statistics of the fitted
#' calibration model given in \code{object} using the components \code{"call"}, \code{"size"},
#' \code{"rdf"}, \code{"r.squared"} and \code{"sigma"} from its argument plus
#'
#' \item{coefficients}{a px4 matrix with columns for the estimated coefficient, its standard error,
#' t-statistic and corresponding (two-sided) p-value.}
#'
#' @seealso
#' The calibration model fitting function \link[mecor]{mefit}, \link[base]{summary}
#'
#' Function \link[stats]{coef} will extract the matrix of coefficients with standard errors,
#' t-statistics and p-values
#'
#' @examples
#' ## Continuing the mecor() example:
#' coef(fit_sme)
#' summary(cm_sme)
#'
#' @export summary.mecor
#' @export
summary.mecor <- function(object){
  z <- object
  est <- z$coefficients
  se <- z$stderr
  t <- est/se
  rdf <- z$rdf
  out <- z[c("call")]
  out$coef.nm <- z$coefficients.nm
  out$coefficients <- cbind('Estimate' = est,
                       'Std. Error (ZV)' = se,
                       't value' = t,
                       `Pr(>|t|)` = 2 * pt(abs(t), rdf,
                                           lower.tail = FALSE))
  out$rdf <- rdf
  out$ci <- z$ci
  out$est <- est
  out$dif.var <- z$dif.var
  out$Rbtstrp <- z$Rbtstrp
  class(out) <- "summary.mecor"
  out
}

#' @export
print.summary.mecor <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  if(!is.null(x$dif.var)){
    cat("\nDifferential measurement error model:\n", paste(names(x$dif.var), "specified in formula, is linked to", x$dif.var, "specified in mefit object\n"), sep = "")}
  cat("\nCoefficients Uncorrected Model:\n")
  printCoefmat(x$coef.nm, signif.stars = F)
  cat("\nCoefficients Corrected Model:\n")
  printCoefmat(x$coefficients, signif.stars = F)
  cat("\nConfidence Intervals for", paste(names(x$est)[2], ":", sep = ""), "\n")
  print(x$ci)
  if(!is.na(x$call$B) && x$call$B != 0){
    cat("Bootstrap Confidence Interval is based on", x$Rbtstrp, "bootstrap replicates")  }
  invisible(x)
}

#' @export
print.mecor <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  if(!is.null(x$dif.var)){
    cat("\nDifferential measurement error model:\n", paste(names(x$dif.var), "specified in formula, is linked to", x$dif.var, "specified in mefit object\n"), sep = "")}
  if (length(coef(x))) {
    cat("\nCoefficients corrected model:\n")
    print(x$coefficients)
  }
  else cat("\nNo coefficients\n")
  cat("\n")
  invisible(x)
}
