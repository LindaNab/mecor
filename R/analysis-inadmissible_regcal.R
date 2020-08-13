inadmissible_regcal <- function(response,
                                covars,
                                me,
                                B,
                                type){
  dm_inadm <- mecor:::get_dm_inadm(covars,
                                   me,
                                   type) # design matrix
  inadm_fit <- stats:::lm.fit(dm_inadm,
                              response)
  out <- list(coef = coef(inadm_fit),
              vcov = mecor:::vcovfromfit(inadm_fit))
  if (B != 0){
    boot <- mecor:::analysis_boot(response,
                                 covars,
                                 me,
                                 B = B,
                                 type = type,
                                 method = "irc")
    colnames(boot$coef) <- names(out$coef)
    out$boot <- boot
  }
  out
}
