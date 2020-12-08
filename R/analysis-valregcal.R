valregcal <- function(response,
                      covars,
                      me,
                      B,
                      type) {
  dm_vrc <- get_dm_vrc(covars,
                       me,
                       type) # design matrix
  vrc_fit <- stats::lm.fit(dm_vrc,
                           response)
  out <- list(coef = stats::coef(vrc_fit),
              vcov = vcovfromfit(vrc_fit))
  out$method <- "validation regression calibration"
  if (B != 0) {
    boot <- analysis_boot(
      response,
      covars,
      me,
      B = B,
      type = type,
      method = "valregcal"
    )
    colnames(boot$coef) <- names(out$coef)
    out$boot <- boot
  }
  out
}
