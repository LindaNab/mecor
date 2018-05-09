#' Run calibration
#'
#' This is a function that fits a linear regression to extract the
#' measurement error parameters of the measurement error model
#' from the external validation set in case you assume that the underlying
#' measurement error model is of the systematic form.
#'
#' @param formula an object of class \link[stats]{formula} (or one that is
#' coerced to that class): a symbolic description of the naive
#' model to be fitted in an \link[stats]{lm} model.
#' @param data an optional data frame, list or environment (or
#' object coercible by as.data.frame to a data frame) containing
#' the variables in the model. If not found in \code{data}, the
#' variables are taken from \code{environment(formula)}, typically
#' the enviroment from which \code{systematic} is called.
#' @param alpha the alpha level of which the confidence intervals
#' are based.
#'
#' @return \code{systematic} returns an object of \link[base]{class} "systematic".
#'
#' An object of class \code{systematic} is a list containing the following components:
#'
#' \item{coefficients}{a named vector of the coefficients of the calibration model}
#' \item{alpha}{the alpha level}
#' \item{model}{an \code{lm} object}
#'
#' @author Linda Nab, \email{l.nab@lumc.nl}
#'
#' @references
#'
#' @examples
#' Xcal <- c(rep(0,500),rep(1,500))
#' Ycal <- Xcal + rnorm(1000,0,3)
#' Vcal <- 1 + 2 * Ycal + rnorm(1000,0,3)
#' syst <- systematic(formula = Vcal ~ Ycal)
#'
#' Vcal_dme0 <- 1 + 2 * Ycal[1:500] + rnorm(Ycal[1:500], 0, 3)
#' Vcal_dme1 <- 2 + 3 * Ycal[501:1000] + rnorm(Ycal[501:1000], 0, 3)
#' formula0 <- Vcal_dme0 ~ Ycal[1:500]
#' formula1 <- Vcal_dme1 ~ Ycal[501:1000]
#' syst1 <- systematic(formula = formula0)
#' syst2 <- systematic(formula = formula1)
#'
#' @export
systematic <- function(formula,
                  data = NULL,
                  alpha = 0.05)
  {if( attr(terms(formula), "variables")[4]!="NULL()"){
        stop("'formula' should be a formula of a simple linear model")
  }
    model <- lm(formula,data)
    summary_cal <- summary(model)
    theta0_hat <- summary_cal$coef[1,1]
    theta1_hat <- summary_cal$coef[2,1]
    t <-summary_cal$sigma
    K <- NROW(model$model)
    s_yy <- sum((model$model$Ycal-mean(model$model$Ycal))^2)
    mean_Ycal <- mean(model$model[,names(attr(terms(formula),"factors")[,1])[2]])
    coefs <- c(theta0_hat=theta0_hat, theta1_hat=theta1_hat,
               t=t, K=K, mean_Ycal=mean_Ycal, s_yy=s_yy)

    out <- list(coefficients = coefs, # list of coefficients of measurement error model
      alpha = alpha,
      model = model # the calibration model
    )
    class(out) <- 'systematic'
    return(out)
  }

delta <- function(model_naive,
                  model_corrected,
                  systematic){
  X <- model_naive$model[,as.character(model_naive$terms[[3]])]
  s_xx <- sum((X - mean(X)) ^ 2)
  t_q <- qt((1 - systematic$alpha / 2),(NROW(X) - 2))
  varbeta <- 1 / (systematic$coefficients["theta1_hat"]) ^ 2 *
                (summary(model_naive)$sigma ^ 2 / s_xx +
                (model_corrected$coefficients[2] ^ 2 * systematic$coefficients["t"] ^ 2) /
                systematic$coefficients["s_yy"] )
  ci <- c(lower = unname(model_corrected$coefficients[2] - t_q * sqrt(varbeta)),
          upper = unname(model_corrected$coefficients[2] + t_q * sqrt(varbeta))
          )
  return(ci)
}

fieller <- function(model_naive,
                    model_corrected,
                    systematic){
  X <- model_naive$model[,as.character(model_naive$terms[[3]])]
  s_xx <- sum((X - mean(X))^2)
  t_q <- qt((1 - systematic$alpha / 2), (NROW(X) - 2))
  v1 <- - 1 *  (model_naive$coefficients[2] * systematic$coefficients["theta1_hat"])
  v2 <- (systematic$coefficients["t"] ^ 2 / systematic$coefficients["s_yy"]) * t_q ^ 2 - systematic$coefficients["theta1_hat"] ^ 2
  v3 <- (summary(model_naive)$sigma ^ 2 / s_xx) * t_q ^ 2 - model_naive$coefficients[2]^2
  D <- v1 ^ 2 - v2 * v3
  if(v2 < 0 & D > 0){
    l1 <- unname((v1 - sqrt(D)) / v2)
    l2 <- unname((v1 + sqrt(D)) / v2)
    ci <- c(lower = min(l1, l2),
          upper = max(l1, l2)
          )
    }
  else ci <- c(lower = NA, upper = NA)
  return(ci)
}

