#' Run linear regression with correction method
#'
#' This is a function that runs a linear regression with a correction method.
#' You can fill your function using the \code{x}.
#'
#' @param x A vector with the exposure variable
#' @param y A vector with the outcome variable
#' @param w A vector with the surrogate outcome variable
#' @param dataset A dataframe with the dataset
#' @param calset A dataframe with the calibration dataset
#' @param form_of_me Form of measurement error (default is classical)
#' @param correction_method Correction method (default is zero_variance)
#'
#' @return This function returns (..).
#'
#' @examples
#' lm_me(x, y, w)
#' lm_me(x, y, w, form_of_me = dme)
#'
#' @export
lm_me <- function(x,
                  y,
                  w,
                  dataset,
                  calset,
                  form_of_me = cme,
                  correction_method = zero_variance){
  print(x)
}
