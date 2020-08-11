#' @export
print.MeasErrorRandom <- function(x){
  cat("\nCall:\n", deparse(attributes(x)$call), "\n", sep = "")
  cat("\nThe error-prone variable", deparse((attr(x, 'input')$substitute)),
      "is assumed measured with random measurement error of magnitude", x$error)
  invisible(x)
}
