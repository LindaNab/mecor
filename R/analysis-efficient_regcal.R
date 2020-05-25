#efficient regression calibration
efficient_regcal <- function(response, covars, me, B, alpha, use_vcov){
  # complete case fit
  cc_fit <- mecor:::complete_case(response, covars, me)
  # reg_cal fit
  if (use_vcov == "default"){
    rc_fit <- mecor:::regcal(response, covars, me)
  } else if (use_vcov == "bootstrap"){
    rc_fit <- mecor:::regcal(response, covars, me, B = 999) # default B
  }
  cc_fit_coef <- mecor:::change_order_coefs(cc_fit$coef)
  vcov_cc_fit <- mecor:::get_vcov(cc_fit)
  inv_vcov_cc_fit <- solve(vcov_cc_fit)
  if (use_vcov == "default"){
    inv_vcov_rc_fit <- solve(rc_fit$vcov)
  } else if (use_vcov == "bootstrap"){
    inv_vcov_rc_fit <- solve(rc_fit$boot$vcov)
  }
  beta <- solve(inv_vcov_cc_fit + inv_vcov_rc_fit) %*%
    (inv_vcov_cc_fit %*% cc_fit_coef +
       inv_vcov_rc_fit %*% rc_fit$coef)
  vcov_beta <- solve(inv_vcov_rc_fit + inv_vcov_cc_fit)
  out <- list(coef = t(beta)[1, ],
              vcov = vcov_beta)
  out
}
