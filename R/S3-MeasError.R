#' @export
print.MeasError <- function(x){
  cat("\nCall:\n", deparse(attributes(x)$call), "\n", sep = "")
  if (is.null(x$replicate)){
    cat("\nThe error-prone variable", deparse((attr(x, 'input')$substitute)),
        "is correctly measured by",   deparse((attr(x, 'input')$reference)))
  } else{
    cat("\nThe error-prone variable ", deparse((attr(x, 'input')$substitute)),
        " has replicate measures ", deparse((attr(x, 'input')$replicate)),
        " with classical measurement error", sep = "")
  }
  invisible(x)
}
