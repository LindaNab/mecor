#' mecor: a measurement error correction package
#'
#' mecor provides correction methods for measurement
#' errors in continuous endpoints.
#'
#' @param formula formula of the naive regression model
#' @param data data
#' @param cv a non-empty character string specifying the variable
#' in \code{formula} to correct
#' @param correction the method to correct \code{cv},
#' a object of type \code{systematic} is implemented
#'
#' @return This function returns the corrected effect estimator.
#'
#' @examples
#' Xcal <- c(rep(0,50),rep(1,50))
#' Ycal <- Xcal + rnorm(100,0,3)
#' Vcal <- 1 + 2*Ycal + rnorm(100,0,3)
#' syst <- systematic(formula = Vcal ~ Ycal)
#'
#' X <- c(rep(0,100),rep(1,100))
#' Y <- X + rnorm(200,0,3)
#' V <- 1 + 2*Y + rnorm(100,0,3)
#' model_mecor <- mecor(formula = V ~ X, cv = "V", correction = syst)
#'
#'
#' @export
mecor <- function(formula,
                  data = NULL,
                  cv,
                  correction){
  model_naive <- lm(formula, data)
  if(length(formula[[2]])!=1 || length(formula[[3]])!=1 ){
    stop("'formula' should be a simple linear model")
  }
  if(formula[[2]] != cv){
    stop("variable 'cv' should be the dependent variable")
  }
  if(class(correction) != 'systematic'){
    stop("'correction' should be of class systematic")
  }
  if(class(correction) == "systematic"){
    adjustedvar <- (model_naive$model[,cv] - correction$coefficients["theta0_hat"]) /
                    correction$coefficients["theta1_hat"]
    data$adjustedvar <- adjustedvar
    model_corrected <- lm(formula = update.formula(formula,adjustedvar ~ .), data = data)
    X <- model_naive$model[,as.character(formula[[3]])]
    s_xx <- sum((X - mean(X))^2)
    t_q <- qt((1 - correction$alpha / 2), (NROW(X) - 2))
    ci.zv <- c(lower = unname(model_corrected$coefficients[2] - t_q * sqrt(summary(model_corrected)$sigma ^ 2 / s_xx)),
               upper = unname(model_corrected$coefficients[2] + t_q * sqrt(summary(model_corrected)$sigma ^ 2 / s_xx)))
    ci.delta <- delta(model_naive, model_corrected, correction)
    ci.fieller <- fieller(model_naive, model_corrected, correction)
    cibeta <- c(zerovariance = ci.zv,
                   delta = ci.delta,
                   fieller = ci.fieller)
  }

  out <- list(model = model_corrected,
              ci = cibeta
  )
  class(out) <- 'mecor'
  return(out)
}


