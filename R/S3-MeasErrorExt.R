#' @export
print.MeasErrorExt <- function(x, ...) {
  cat("\nCall:\n", deparse(attributes(x)$call), "\n", sep = "")
  cat(
    "\nThe error-prone variable",
    deparse((attr(x, 'input')$substitute)),
    "will be corrected using the following model:\n"
  )
  print(x$model$coef)
  invisible(x)
}
