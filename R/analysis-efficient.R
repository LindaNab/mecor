# efficient regcal (cov me) and efficient mm (outcome me)
efficient <- function(response,
                      covars,
                      me,
                      B,
                      type,
                      calc_vcov = T) {

  # if the MeasError object contains replicates, mle will be used to obtain a
  # internal fit using only the replicates
  # when the MeasError object does not contain replicates, a complete case will
  # be used (using the reference)
  if (!is.null(me$replicate)){
    rownums_cc <- which (!is.na(me$reference))
    n_rep <- ncol(me$replicate)
    new_me <- me
    new_me$substitute <- me$replicate[rownums_cc, 1]
    new_me$replicate <- me$replicate[rownums_cc, 2:n_rep, drop = F]
    cc_response <- response[rownums_cc, , drop = F]
    cc_covars <- covars[rownums_cc, , drop = F]
    cc_fit <- mle(cc_response,
                  cc_covars,
                  new_me,
                  B = 0,
                  type) # mle fit
    vcov_cc_fit <- cc_fit$vcov
  } else {
    cc_fit <- complete_case(response,
                            covars,
                            me,
                            type) # complete case fit
    vcov_cc_fit <- get_vcov(cc_fit,
                            rev = F)
  }
  # reg_cal fit
  rc_fit <- standard(response,
                     covars,
                     me,
                     B = 0,
                     type = type)
  cc_fit_coef <- get_coefs(cc_fit,
                           rev = F)
  inv_vcov_cc_fit <- erc_solve(vcov_cc_fit, "complete case")
  inv_vcov_rc_fit <- erc_solve(rc_fit$vcov, "regcal")
  beta <- solve(inv_vcov_cc_fit + inv_vcov_rc_fit) %*%
    (inv_vcov_cc_fit %*% cc_fit_coef +
       inv_vcov_rc_fit %*% rc_fit$coef)
  out <- list(coef = t(beta)[1, ])
  if (calc_vcov){
    vcov_beta <- solve(inv_vcov_rc_fit + inv_vcov_cc_fit)
    out$vcov <- vcov_beta
  }
  out$method <- efficient_get_method(type)
  if (B != 0){
    boot <- analysis_boot(response,
                          covars,
                          me,
                          B = B,
                          type = type,
                          method = "efficient")
    colnames(boot$coef) <- names(out$coef)
    out$boot <- boot
  }
  out
}
# Catches errors when matrices are singular (so cannot be solved)
erc_solve <- function(vcov,
                      type) {
  out <- tryCatch(
    {
      solve(vcov)
    },
    error = function(cond){
      warning(paste0("The vcov matrix of ", type, " is computationally singular. The ", type, " estimate got weight 0."),
              call. = FALSE)
      return(matrix(0,
                    nrow = nrow(vcov),
                    ncol = ncol(vcov)))
    }
  )
  return(out)
}

efficient_get_method <- function(type) {
  if (startsWith(type, "dep")) {
    method <- "efficient method of moments"
  } else if (type == "indep") {
    method <- "efficient regression calibration"
  }
  method
}
