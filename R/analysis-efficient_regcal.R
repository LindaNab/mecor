#efficient regression calibration
efficient_regcal <- function(response, covars, me, B, alpha, type, calc_vcov = T){
  if (!is.null(me$replicate)){
    rownums_cc <- which (!is.na(me$reference))
    n_rep <- ncol(me$replicate)
    new_me <- me
    new_me$substitute <- me$replicate[rownums_cc, 1]
    new_me$replicate <- me$replicate[rownums_cc, 2:n_rep, drop = F]
    cc_response <- response[rownums_cc, , drop = F]
    cc_covars <- covars[rownums_cc, , drop = F]
    cc_fit <- mecor:::mle(cc_response, cc_covars, new_me, B = 0, alpha, type)
    vcov_cc_fit <- cc_fit$vcov
  } else {
    cc_fit <- mecor:::complete_case(response, covars, me, type) # complete case fit
    vcov_cc_fit <- mecor:::get_vcov(cc_fit, rev = F)
  }
  # reg_cal fit
  rc_fit <- mecor:::regcal(response, covars, me, B = 0, type = type)
  cc_fit_coef <- mecor:::get_coefs(cc_fit, rev = F)
  inv_vcov_cc_fit <- solve(vcov_cc_fit)
  inv_vcov_rc_fit <- solve(rc_fit$vcov)
  beta <- solve(inv_vcov_cc_fit + inv_vcov_rc_fit) %*%
    (inv_vcov_cc_fit %*% cc_fit_coef +
       inv_vcov_rc_fit %*% rc_fit$coef)
  out <- list(coef = t(beta)[1, ])
  if (calc_vcov){
    vcov_beta <- solve(inv_vcov_rc_fit + inv_vcov_cc_fit)
    out$vcov <- vcov_beta
  }
  if (B != 0){
    boot <- mecor:::analysis_boot(response,
                                  covars,
                                  me,
                                  B = B,
                                  type = type,
                                  method = "erc")
    colnames(boot$coef) <- names(out$coef)
    out$boot <- boot
  }
  out
}
