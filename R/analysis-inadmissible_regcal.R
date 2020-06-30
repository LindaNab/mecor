inadmissible_regcal <- function(response, covars, me, B = 0, alpha = 0.05){
  dm_inadm <- mecor:::get_dm_inadm(covars, me)
  inadm_fit <- stats:::lm.fit(dm_inadm, response)
  out <- list(coef = coef(inadm_fit),
              vcov = mecor:::vcovfromfit(inadm_fit))
  if (B !=0){
    boot <-
      mecor:::analysis_boot(response, covars, me,
                            B = B, alpha = alpha, type = "indep", method = "irc")
    out$boot <- list(ci = boot$ci,
                     vcov = boot$vcov)
  }
  out
}
