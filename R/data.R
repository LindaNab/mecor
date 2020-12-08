#' Internal Covariate-Validation Study
#'
#' A simulated dataset containing 1000 observations of the outcome Y, the
#' error-prone exposure X_star and the covariate Z. The gold standard exposure
#' X is observed in approximately 25\% of the study population
#'
#' @format A data frame with 1000 rows and 4 variables:
#' \describe{
#'   \item{Y}{outcome, continuous}
#'   \item{X_star}{error-prone exposure, continuous}
#'   \item{Z}{covariate, continuous}
#'   \item{X}{gold standard exposure, continuous}
#' }
#' @examples
#' data("icvs", package = "mecor")
"icvs"
#' Internal Outcome-Validation Study
#'
#' A simulated dataset containing 1000 observations of the error-prone outcome
#' Y_star, the exposure X and the covariate Z. The gold standard outcome Y is
#' observed in approximately 25\% of the study population
#'
#' @format A data frame with 1000 rows and 4 variables:
#' \describe{
#'   \item{Y_star}{error-prone outcome, continuous}
#'   \item{X}{exposure, continuous}
#'   \item{Z}{covariate, continuous}
#'   \item{Y}{gold standard outcome, continuous}
#' }
#' #' @examples
#' data("iovs", package = "mecor")
"iovs"
#' Internal Outcome Validation Study
#'
#' A simulated dataset containing 1000 observations of the error-prone outcome
#' Y_star with differential measurement error, and the binary exposure X. The
#' gold standard  outcome Y is observed in approximately 25\% of the study
#' population
#'
#' @format A data frame with 1000 rows and 3 variables:
#' \describe{
#'   \item{Y_star}{error-prone outcome, continuous}
#'   \item{X}{exposure, binary}
#'   \item{Y}{gold standard outcome, continuous}
#' }
#' @examples
#' data("iovs_diff", package = "mecor")
"iovs_diff"
#' Replicates Study
#'
#' A simulated dataset containing 1000 observations of the outcome Y, three
#' replicate measures of the error-prone exposure X1_star, X2_star and X3_star
#' and two covariates Z1 and Z2.
#'
#' @format A data frame with 1000 rows and 6 variables:
#' \describe{
#'   \item{Y}{outcome, continuous}
#'   \item{X_star_1}{first replicate of error-prone exposure, continuous}
#'   \item{X_star_2}{second replicate of error-prone exposure, continuous}
#'   \item{X_star_3}{third replicate of error-prone exposure, continuous}
#'   \item{Z1}{covariate, continuous}
#'   \item{Z2}{covariate, continuous}
#' }
#' @examples
#' data("rs", package = "mecor")
"rs"
#' Covariate-Calibration Study
#'
#' A simulated dataset containing 1000 observations of the outcome Y, the
#' error-prone exposure X_star with systematic measurement error, the covariate
#' Z and two replicates measures X1_star and X2_star of the exposure with
#' classical measurement error. The two replicates are observed in the first
#' 500 study participants.
#'
#' @format A data frame with 1000 rows and 5 variables:
#' \describe{
#'   \item{Y}{outcome, continuous}
#'   \item{X_star}{error-prone exposure with systematic measurement error,
#'   continuous}
#'   \item{Z}{covariate, continuous}
#'   \item{X_star_1}{first replicate of error-prone exposure with classical
#'   measurement error, continuous}
#'   \item{X_star_2}{second replicate of error-prone exposure with classical
#'   measurement error, continuous}
#' }
#' @examples
#' data("ccs", package = "mecor")
"ccs"
#' External Covariate-Validation Study
#'
#' A simulated dataset containing 100 observations of the gold standard
#' measurement of the exposure X, the error-prone exposure X_star and the
#' covariate Z. The outcome Y is not observed in the study population.
#' To be used in combination with the dataset \link[mecor]{icvs}.
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{X}{gold standard exposure, continuous}
#'   \item{X_star}{error-prone exposure, continuous}
#'   \item{Z}{covariate, continuous}
#' }
#' @examples
#' data("ecvs", package = "mecor")
"ecvs"
#' External Outcome-Validation Study
#'
#' A simulated dataset containing 100 observations of the gold standard
#' measurement of outcome Y and the error-prone outcome Y_star. The covariates X
#' and Z are not observed in the study population. To be used in combination
#' with the dataset \link[mecor]{iovs}.
#'
#' @format A data frame with 100 rows and 2 variables:
#' \describe{
#'   \item{Y}{gold standard outcome, continuous}
#'   \item{Y_star}{error-prone outcome, continuous}
#' }
#' @examples
#' data("eovs", package = "mecor")
"eovs"
#' Simulated dataset for the \link[mecor]{ipwm} function
#'
#' A simulated dataset containing 5000 observations of the covariates L1-L10,
#' the true exposure A and true outcome Y, and the misclassified exposure B and
#' misclassified outcome Z.
#'
#' @format A data frame with 5000 rows and 14 variables:
#' \describe{
#'   \item{L1}{covariate, binary}
#'   \item{L2}{covariate, continuous}
#'   \item{L3}{covariate, binary}
#'   \item{L4}{covariate, continuous}
#'   \item{L5}{covariate, binary}
#'   \item{L6}{covariate, binary}
#'   \item{L7}{covariate, continuous}
#'   \item{L8}{covariate, binary}
#'   \item{L9}{covariate, binary}
#'   \item{L10}{covariate, continuous}
#'   \item{A}{exposure, binary}
#'   \item{Y}{outcome, binary}
#'   \item{B}{misclassified exposure, binary}
#'   \item{Z}{misclassified outcome, binary}
#' }
#' @examples
#' data("sim", package = "mecor")
"sim"
