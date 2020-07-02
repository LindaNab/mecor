inadmissible_regcal <- function(response, covars, me, B, alpha, type){
  mecor:::check_input_inadmissible_regcal(me, type)
  dm_inadm <- mecor:::get_dm_inadm(covars, me, type)
  inadm_fit <- stats:::lm.fit(dm_inadm, response)
  out <- list(coef = coef(inadm_fit),
              vcov = mecor:::vcovfromfit(inadm_fit))
  if (B != 0){
    boot <-
      mecor:::analysis_boot(response, covars, me,
                            B = B, alpha = alpha,
                            type = type, method = "irc")
    out$boot <- list(ci = boot$ci,
                     vcov = boot$vcov)
  }
  out
}
check_input_inadmissible_regcal <- function(me, type){
  if (startsWith(type, "dep")){
    stop("Inadmissible regression calibration is not suitable for measurement error in the dependent variable")
  }
  if (class(me)[1] == "MeasErrorExt"){
    stop("Inaddmissible regression calbration is not suitable for external designs")
  } else if (class(me)[1] == "MeasError" && (type == "indep" & !is.null(me$replicate))){
    stop("Inadmissible regression calibration is not suitable for a design with replicates")
  }
}
