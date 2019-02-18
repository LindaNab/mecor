#' @title Summarizing Measurement Error Model fits
#'
#' @description
#' \code{summary} method for class "mefit"
#'
#' @param object an object of class "mefit", a result of a call to \link[mecor]{mefit}
#'
#' @return
#' The function \code{summary.mefit} returns a list of summary statistics of the fitted
#' measurement erro models given in \code{object}
#'
#' @seealso
#' The measurement error model fitting function \link[mecor]{mefit}, \link[base]{summary}
#'
#' @examples
#' ## Continuing the mefit() example:
#' #coef(fit1)
#' #summary(fit1)
#'
#' @export summary.mefit
#' @export
summary.mefit <- function(object){
  out <- object
  class(out) <- "summary.mefit"
  out
}

#' @export
print.summary.mefit <- function(x){
  cat("\nCall:\n", paste(deparse(attr(x, "call")), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  if(length(coef(x$cme))){
    cat("\nClassical Measurement Error:")
    cat("\nCoefficients:\n")
    print(coef(x$cme))
    cat("\nResidual standard error:", x$cme$Sigma, "on", x$cme$df, "degrees of freedom\n")}
  if(length(coef(x$sme1))){
    cat("\nSystematic Measurement Error (zero intercept):\n")
    print(coef(x$sme1))
    cat("\nResidual standard error:", x$sme1$Sigma, "on", x$sme1$df, "degrees of freedom\n")}
  if(length(coef(x$sme2))){
    cat("\nSystematic Measurement Error (non-zero intercept):\n")
    print(coef(x$sme2))
    cat("\nResidual standard error:", x$sme2$Sigma, "on", x$sme2$df, "degrees of freedom\n")}
  if(length(coef(x$dme))){
    cat("\nDifferential Measurement Error:\n")
    print(coef(x$dme))
    cat("\nResidual standard error:", x$dme$Sigma, "on", x$dme$df, "degrees of freedom\n")}
  cat("\nCharacteristics Calibration Data:")
  cat("\nSize of calibration data set:", attr(x, "size"))
  invisible(x)
}

#' @export
print.mefit <- function(x){
  cat("\nCall:\n", paste(deparse(attr(x, "call")), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  if (length(coef(x$cme))) {
    cat("\nClassical Measurement Error:\n")
    coef <- c(coef(x$cme)[,1])
    names(coef) <- rownames(coef(x$cme))
    print(coef)}
  if (length(coef(x$sme1))) {
    cat("\nSystematic Measurement Error (zero intercept):\n")
    coef2 <- c(coef(x$sme1)[,1])
    names(coef2) <- rownames(coef(x$sme1))
    print(coef2)}
  if (length(coef(x$sme2))) {
    cat("\nSystematic Measurement Error (non-zero intercept):\n")
    print(coef(x$sme2)[,1])}
  if (length(coef(x$dme))) {
    cat("\nDifferential Measurement Error:\n")
    print(coef(x$dme)[,1])}
  cat("\n")
  invisible(x)
}
