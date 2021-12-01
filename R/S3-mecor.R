#' @title Summarizing Measurement Error Correction
#'
#' @description
#' \code{summary} method for class "mecor"
#'
#' @param object an object of class "mecor", a result of a call to
#' \link[mecor]{mecor}.
#' @param alpha probability of obtaining a type II error.
#' @param zerovar a boolean indicating whether standard errors and confidence
#' intervals using the zerovariance method must be added to the summary object.
#' @param fieller a boolean indicating whether confidence intervals using the
#' fieller method must be added to the summary object.
#' @param ... additional arguments affecting the summary produced
#'
#' @return
#' The function \code{summary.mecor} returns a list of summary statistics of the
#' fitted corrected model and fitted uncorrected model.
#'
#' \item{call}{the matched call}
#' \item{c}{summary of the corrected fit}
#' \item{uc}{summary of the uncorrected fit}
#' \item{B}{number of bootstrap replicates used}
#' \item{alpha}{alpha level used}
#'
#' @seealso
#' The model fitting function \link[mecor]{mecor}, \link[base]{summary}
#'
#' @examples
#' ## measurement error in a covariate:
#' # internal covariate-validation study
#' data(vat)
#' mecor_fit <- mecor(ir_ln ~ MeasError(wc, reference = vat) + sex + age + tbf,
#'                    data = vat,
#'                    method = "standard")
#' summary(mecor_fit)
#' summary(mecor_fit, zerovar = TRUE, fieller = TRUE)
#' summary(mecor_fit, alpha = 0.10)
#'
#' @export
summary.mecor <- function(object,
                          alpha = 0.05,
                          zerovar = FALSE,
                          fieller = FALSE,
                          ...) {
  z <- object
  z1 <- z$uncorfit
  z2 <- z$corfit
  if (zerovar == TRUE && is.null(z2$zerovar_vcov)) {
    warning("there is no 'zerovar_vcov' object, zerovar set to FALSE")
    zerovar <- FALSE
  }
  if (fieller == TRUE && is.null(z2$fieller)) {
    warning("there is no 'fieller' object, fieller set to FALSE")
    fieller <- FALSE
  }
  # uncorrected
  rdf1 <- z1$df.residual
  rss1 <- sum(z1$residuals ^ 2)
  tq <- stats::qt((1 - alpha / 2), rdf1)
  uc <- list(
    residuals = z1$residuals,
    rdf = rdf1,
    sigma = sqrt(rss1 / rdf1)
  )
  uc$coefficients <- cbind(
    Estimate = (coef1 <- z1$coef),
    'Std. Error' = (se1 <-
                      sqrt(diag(
                        vcovfromfit(z1)
                      ))),
    't value' = (t1 <- coef1 / se1),
    'Pr(>|t|)' = 2 * stats::pt(abs(t1), rdf1, lower.tail = FALSE)
  )
  uc$ci <- cbind(
    Estimate = coef1,
    'LCI' = coef1 - tq * se1,
    'UCI' = coef1 + tq * se1
  )
  uc$ci <- round(uc$ci, 6)
  # corrected
  coefficients <- cbind(Estimate = (coef2 <- z2$coef))
  ci <- cbind(Estimate = coef2)
  zq <- stats::qnorm((1 - alpha / 2))
  if (!is.null(z2$vcov)) {
    SE <- {
      se2 <- sqrt(diag(z2$vcov))
    }
    LCI <- coef2 - zq * se2
    UCI <- coef2 + zq * se2
    coefficients <- cbind(coefficients, SE)
    ci <- cbind(ci, LCI, UCI)
  }
  if ({
    B <- attr(z, "B")
  } != 0) {
    boot_vcov <- stats::cov(z2$boot$coef)
    SE_btstr <- sqrt(diag(boot_vcov))
    boot_ci <- apply(z2$boot$coef,
                     2,
                     FUN = stats::quantile,
                     probs = c(alpha / 2, 1 - alpha / 2))
    LCI_btstr <- boot_ci[1,]
    UCI_btstr <- boot_ci[2,]
    coefficients <- cbind(coefficients,
                          'SE (btstr)' = SE_btstr)
    ci <- cbind(ci,
                'LCI (btstr)' = LCI_btstr,
                'UCI (btstr)' = UCI_btstr)
  }
  if (zerovar == T) {
    SE_zerovar <- {
      se2_zv <- sqrt(diag(z2$zerovar_vcov))
    }
    LCI_zerovar <- coef2 - zq * se2_zv
    UCI_zerovar <- coef2 + zq * se2_zv
    coefficients <- cbind(coefficients, 'SE (zerovar)' = SE_zerovar)
    ci <- cbind(ci,
                'LCI (zerovar)' = LCI_zerovar,
                'UCI (zerovar)' = UCI_zerovar)
  }
  if (fieller == TRUE) {
    fieller_ci <- calc_fieller_ci(
      z2$fieller$lambda1,
      z2$fieller$var_lambda1,
      z2$fieller$phi_star,
      z2$fieller$var_phi_star,
      alpha
    )
    fieller_ci_temp <- matrix(nrow = nrow(ci), ncol = 2)
    if (nrow(fieller_ci) == 1) {
      fieller_ci_temp[2,] <- fieller_ci
      fieller_ci <- fieller_ci_temp
    } else if (nrow(fieller_ci) > 1) {
      fieller_ci_temp[-1,] <- fieller_ci
      fieller_ci <- fieller_ci_temp
    }
    ci <- cbind(ci,
                'LCI (fieller)' = fieller_ci[, 1],
                'UCI (fieller)' = fieller_ci[, 2])
  }
  c <- list(coefficients = round(coefficients, 6))
  if (dim(ci)[2] > 1) {
    c$ci <- round(ci, 6)
  }
  c$method <- z2$method
  out <- list(call = attr(z, "call"))
  out$c <- c
  out$uc <- uc
  out$B <- B
  out$alpha <- alpha
  class(out) <- "summary.mecor"
  out
}
#' @export
print.summary.mecor <- function(x, ...) {
  cat("\nCall:\n",
      paste(deparse(x$call),
            sep = "\n",
            collapse = "\n"),
      "\n",
      sep = "")
  cat("\nCoefficients Corrected Model:\n")
  print(x$c$coefficients)
  if (!is.null(x$c$ci)) {
    cat("\n",
        paste((1 - x$alpha) * 100,
              "%",
              sep = ""),
        " Confidence Intervals:\n",
        sep = "")
    print(x$c$ci)
  }
  if (x$B != 0) {
    cat(
      "Bootstrap Confidence Intervals are based on",
      x$B,
      "bootstrap replicates using percentiles \n"
    )
  }
  cat("\nThe measurement error is corrected for by application of",
      x$c$method,
      "\n")
  if (!is.null(x$c$matrix)) {
    cat("\nModel Matrix:\n")
    print(x$c$matrix)
  }
  cat("\nCoefficients Uncorrected Model:\n")
  stats::printCoefmat(x$uc$coefficients,
                      signif.stars = F)
  cat("\n",
      paste((1 - x$alpha) * 100,
            "%",
            sep = ""),
      " Confidence Intervals:\n",
      sep = "")
  print(x$uc$ci)
  cat("\nResidual standard error:",
      x$uc$sigma,
      "on",
      x$uc$rdf,
      "degrees of freedom\n")
  invisible(x)
}
#' @export
print.mecor <- function(x, ...) {
  cat("\nCall:\n",
      paste(deparse(attr(x, "call")),
            sep = "\n",
            collapse = "\n"),
      "\n",
      sep = "")
  if (length(x$corfit$coef)) {
    cat("\nCoefficients Corrected Model:\n")
    print(x$corfit$coef)
  }
  if (length(x$uncorfit$coef)) {
    cat("\nCoefficients Uncorrected Model:\n")
    print(x$uncorfit$coef)
  }
  cat("\n")
  invisible(x)
}
