#efficient regression calibration
efficient_regcal <- function(response, covars, me, B, alpha, type, erc_B){
  # complete case fit
  cc_fit <- mecor:::complete_case(response, covars, me, type)
  # reg_cal fit
  rc_fit <- mecor:::regcal(response, covars, me, B = erc_B, type = type)
  cc_fit_coef <- mecor:::get_coefs(cc_fit, rev = F)
  vcov_cc_fit <- mecor:::get_vcov(cc_fit, rev = F)
  inv_vcov_cc_fit <- solve(vcov_cc_fit)
  if (erc_B == 0){
    inv_vcov_rc_fit <- solve(rc_fit$vcov)
  } else if (erc_B != 0){
    inv_vcov_rc_fit <- solve(rc_fit$boot$vcov)
  }
  beta <- solve(inv_vcov_cc_fit + inv_vcov_rc_fit) %*%
    (inv_vcov_cc_fit %*% cc_fit_coef +
       inv_vcov_rc_fit %*% rc_fit$coef)
  vcov_beta <- solve(inv_vcov_rc_fit + inv_vcov_cc_fit)
  out <- list(coef = t(beta)[1, ],
              vcov = vcov_beta)
  if (B != 0){
    boot <-
      mecor:::analysis_boot(response, covars, me,
                            B = B, alpha = alpha, type = type,
                            method = "erc", erc_B = erc_B)
    out$boot <- list(ci = boot$ci,
                     vcov = boot$vcov)
  }
  out
}
