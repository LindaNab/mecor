#' Data of a randomised controlled trial
#'
#' A dataset containing data of a randomised controlled trial
#' of 100 observations in two exposure groups and measures of the
#' true endpoint and surrogate endpoint with systematic measurement
#' error.
#'
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{X}{Exposure}
#'   \item{W}{Surrogate endpoint}
#' }
"uaetrial"

#' Calibration set of a randomised controlled trial
#'
#' A dataset containing data of a randomised controlled trial
#' of 20 observations in one exposure group and measures of the
#' true endpoint and surrogate endpoint with systematic measurement
#' error.
#'
#' @format A data frame with 20 rows and 3 variables:
#' \describe{
#'   \item{Xcal}{Exposure}
#'   \item{Ycal}{True endpoint}
#'   \item{Wcal}{Surrogate endpoint}
#' }
"uaetrial_cal"
