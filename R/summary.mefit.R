#' @title Summarizing Calibration Model Fits
#'
#' @description
#' \code{summary} method for class "mefit"
#'
#' @param object an object of class "mefit", a result of a call to \link[mecor]{mefit}
#'
#' @return
#' The function \code{summary.mefit} returns a list of summary statistics of the fitted
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
#' ## Continuing the mefit() example:
#' coef(fit_sme)
#' summary(fit_sme)
#'
#' @export summary.mefit
#' @export
summary.mefit <- function(object){
  z <- object
  est <- z$coefficients
  se <- sqrt(diag(z$vcov))
  t <- est/se
  rdf <- z$df.residual
  out <- z[c("call")]
  out$coefficients <- cbind('Estimate' = est,
                            'Std. Error' = se,
                            't' = t,
                            `Pr(>|t|)` = 2 * pt(abs(t), rdf,
                                                lower.tail = FALSE))
  out$size <- nrow(z$model)
  out$me.structure <- z$me.structure
  out$dif.var <- z$dif.var
  out$rdf <- rdf
  out$r.squared <- summary.lm(z)$r.squared
  out$sigma <- summary.lm(z)$sigma
  class(out) <- "summary.mefit"
  out
}

#' @export
print.summary.mefit <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
  cat("\nCharacteristics Calibration Data:")
  if(x$me.structure == "differential")
    cat("\nThe assumed measurement error structure is:", x$me.structure, "grouped by", x$dif.var)
  else cat("\nThe assumed measurement error structure is:", x$me.structure)
  cat("\nSize of calibration data set:", x$size)
  cat("\nResidual standard error:", x$sigma, "on", x$rdf, "degrees of freedom")
  cat("\nMultiple R-squared:", x$r.squared)
  invisible(x)
}

#' @export
print.mefit <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print(x$coefficients)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}
