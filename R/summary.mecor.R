#' @title Summarizing Measurement Error Correction
#'
#' @description
#' \code{summary} method for class "mecor"
#'
#' @param object an object of class "mecor", a result of a call to \link[mecor]{mecor}
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
#' #coef(fitSme)
#' #summary(cm_sme)
#'
#' @export summary.mecor
#' @export
#'
summary.mecor <- function(object){
  z <- object
  z1 <- z$naivefit
  z2 <- z$corfit
  alpha <- attr(z, "alpha")
  #uncorrected
  rdf1 <- z1$df.residual
  rss1 <- sum(z1$residuals^2)
  tq <- qt((1 - alpha / 2), rdf1)
  uc <- list(residuals = z1$residuals, rdf = rdf1, sigma = sqrt(rss1/rdf1))
  uc$coefficients <- cbind(Estimate = (coef1 <- z1$coef),
                 'Std. Error' = (se1 <- sqrt(diag(mecor:::vcovfromfit(z1)))),
                 't value' = (t1 <- coef1/se1),
                 'Pr(>|t|)' = 2 * pt(abs(t1), rdf1, lower.tail = FALSE))
  uc$ci <- cbind(Estimate = coef1,
                 'Lower CI' = coef1 - tq * se1,
                 'Upper CI' = coef1 + tq * se1)
  uc$ci <- round(uc$ci, 6)
  #corrected
  if(length({q <- attributes(z2)$type}) == 0){
    coefficients <- cbind(Estimate = (coef2 <- z2$coefficients),
              SE = (se2 <- sqrt(z$corvar$var)),
              't value' = (t2 <- coef2/se2),
              'Pr(>|t|)' = 2 * pt(abs(t2), rdf1, lower.tail = FALSE)) #rdf unknown?
    c <- list(coefficients = round(coefficients, 6))
    c$ci <- cbind(Estimate = coef2,
                  'LCI' = coef2 - tq * se2,
                  'UCI' = coef2 + tq * se2,
                  'LCI (btstr)'= (if((B <- attr(z, "B")) != 0) z$ci$bootci[,1] else NA),
                  'UCI (btstr)'= (if(B != 0) z$ci$bootci[,2] else NA))
    c$ci <- round(c$ci, 6)
  }
  else if(q == "lm.fit"){
    rdf2 <- z2$df.residual
    rss2 <- sum(z2$residuals^2)
    c <- list(residuals = z2$residuals, rdf = rdf2, sigma = sqrt(rss2/rdf2))
    c$coefficients <- cbind(Estimate    = (coef2 <- z2$coef),
                           'SE (uncor.)'    = (se2 <- sqrt(diag(mecor:::vcovfromfit(z2)))),
                           't value'  = (t2 <- coef2/se2),
                           'Pr(>|t|)' = 2 * pt(abs(t2), rdf2, lower.tail = FALSE),
                           'SE (delta)' = dse <- sqrt(z$corvar$deltavar),
                           't value'  = dt <- coef2/dse,
                           'Pr(>|t|)' = 2 * pt(abs(dt), rdf2, lower.tail = FALSE))
    c$coefficients <- round(c$coefficients, 6)
    c$ci <- cbind(Estimate = coef2,
                 'LCI (uncor)' = coef2 - tq * se2,
                 'UCI (uncor)' = coef2 + tq * se2,
                 'LCI (delta)' = coef2 - tq * dse,
                 'UCI (delta)' = coef2 + tq * dse,
                 'LCI (btstr)'= (if((B <- attr(z, "B")) != 0) z$ci$bootci[,1] else NA),
                 'UCI (btstr)'= (if(B != 0) z$ci$bootci[,2] else NA),
                 'LCI (fllr)' = z$ci$fiellerci[,1],
                 'UCI (fllr)' = z$ci$fiellerci[,2])
    c$ci <- round(c$ci, 6)}
  out <- list(call = attr(z, "call"), B = B)
  out$uc <- uc
  out$c <- c
  out$alpha <- alpha
  class(out) <- "summary.mecor"
  out
}

#' @export
print.summary.mecor <- function(x){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nCoefficients Corrected Model:\n")
  print(x$c$coefficients)
  cat("\n", paste((1-x$alpha)*100, "%", sep =""), " Confidence Intervals:\n", sep = "")
  print(x$c$ci)
  if(x$B != 0){
    cat("Bootstrap Confidence Intervals are based on", x$B, "bootstrap replicates using percentiles \n")  }
  if(length(x$c$sigma) == 0 && length(x$c$rdf) == 0){
  cat("\nResidual standard error: unknown\n")}
  else{
    cat("\nResidual standard error:", x$c$sigma, "on", x$c$rdf, "degrees of freedom\n")}
  cat("\nCoefficients Uncorrected Model:\n")
  printCoefmat(x$uc$coefficients, signif.stars = F)
  cat("\n", paste((1-x$alpha)*100, "%", sep =""), " Confidence Intervals:\n", sep = "")
  print(x$uc$ci)
  cat("\nResidual standard error:", x$uc$sigma, "on", x$uc$rdf, "degrees of freedom")
  invisible(x)
}

#' @export
print.mecor <- function(x){
  cat("\nCall:\n", paste(deparse(attr(x, "call")), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  if(length(x$corfit$coef)) {
    cat("\nCoefficients Corrected Model:\n")
    print(x$corfit$coef)
  }
  if(length(x$naivefit$coef)) {
    cat("\nCoefficients Uncorrected Model:\n")
    print(x$naivefit$coef)
  }
  cat("\n")
  invisible(x)
}
