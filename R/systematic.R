#' Run calibration
#'
#' This is a function that runs a calibration to substract
#' the measurement error parameters of the measurement error model
#' from the external validation set.
#'
#' @param formula formula the calibration model
#' @param data data
#' @param type form of measurement error,
#' \code{systematic} and \code{differential} are implemented
#' @param method method to estimate standard errors, \code{zero-variance},
#' \code{delta}, \code{fieller} and \code{bootstrap} are implemented
#'
#' @return This function returns the calibration parameters.
#'
#' @examples
#' Xcal <- c(rep(0,50),rep(1,50))
#' Ycal <- Xcal + rnorm(100,0,3)
#' Vcal <- 1 + 2*Ycal + rnorm(100,0,3)
#' syst <- systematic(formula = Vcal ~ Ycal)
#'
#' Vcal_dme <- 1 + (2 - 1) * Xcal + (2 - 2 * Xcal + 3 * Xcal) * Ycal + rnorm(100,0,3)
#' model_calibration1 <- lm(Vcal_dme[1:50] ~ Ycal[0:50])
#' model_calibration2 <- lm(Vcal_dme[51:100] ~ Y+cal[51:100])
#' calibration(model = list(model_calibration1,model_calibration2), form = "differential")
#'
#' @export
systematic <- function(formula,
                  data=NULL,
                  alpha = 0.05)
  { if( length(formula[[2]])!=1 || length(formula[[3]])!=1 ){
        stop("'formula' should be a simple linear model")
    }
    model <- lm(formula,data)
    summary_cal <- summary(model)
    theta0_hat <- summary_cal$coef[1,1]
    theta1_hat <- summary_cal$coef[2,1]
    t <-summary_cal$sigma
    K <- NROW(model$model)
    s_yy <- sum((model$model$Ycal-mean(model$model$Ycal))^2)
    mean_Ycal <- mean(model$model[,as.character(formula[[3]])])
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
  v1 <- model_naive$coefficients[2] * systematic$coefficients["theta1_hat"]
  v2 <- (systematic$coefficients["t"] ^ 2 / systematic$coefficients["s_yy"]) * t_q ^ 2 - systematic$coefficients["theta1_hat"] ^ 2
  v3 <- (summary(model_naive)$sigma ^ 2 / s_xx) * t_q ^ 2 - model_naive$coefficients[2]^2
  D <- v1 ^ 2 - v2 * v3
  if(D > 0){
  ci <- c(lower = unname((v1 - sqrt(D)) / v2),
          upper = unname((v1 + sqrt(D)) / v2)
          )
  }
  else ci <- c(lower = NA, upper = NA)
  return(ci)
}

