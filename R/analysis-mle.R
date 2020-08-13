#' @importFrom lme4 'fixef'
#' @importFrom lme4 'VarCorr'
#' @importFrom lme4 'lmer'
#' @importFrom lmerTest 'as_lmerModLmerTest'
mle <- function(response,
                covars,
                me,
                B,
                type,
                calc_vcov = T){
  if (is.vector(me$replicate)){
    n_rep <- 1
  } else n_rep <- ncol(me$replicate)
  data_wide <- cbind.data.frame(me$substitute,
                                me$replicate,
                                covars,
                                response)
  lm_formula <- stats:::as.formula(paste0(colnames(response),
                                          " ~ ",
                                          paste(colnames(covars),
                                                collapse = " + ")))
  lm_fit <- stats:::lm(formula = lm_formula,
                       data = data_wide)
  data_long <- stats:::reshape(data_wide,
                               varying = 1:{n_rep + 1},
                               v.names = c("X_star"),
                               idvar = "id",
                               direction = "long")
  terms_formula_dep <- paste(colnames(response),
                             paste(colnames(covars), collapse = " + "),
                             sep = " + ")
  lmm_formula <- as.formula(paste0("X_star ~ ",
                                   terms_formula_dep,
                                   " + (1|id)"))
  lmm_fit <- lme4::lmer(lmm_formula,
                        data = data_long,
                        control = lme4::lmerControl(optimizer = "bobyqa",
                                                    calc.derivs = FALSE,
                                                    optCtrl = list(maxfun =
                                                                     2e5)))
  beta <- mecor:::mle_get_coef(lm_fit,
                               lmm_fit)
  beta <- mecor:::change_names(beta,
                               me)
  # change order of coef
  out <- list(coef = mecor:::change_order_coefs(beta))
  if (calc_vcov == T){ # if mle is used in a bootstrap loop
    vcov_beta <- mecor:::mle_get_vcov(lm_fit,
                                      lmm_fit)
    vcov_beta <- mecor:::change_names(vcov_beta,
                                      me)
    out$vcov <- mecor:::change_order_vcov(vcov_beta)
  }
  if (B != 0){
    boot <- mecor:::analysis_boot(response,
                                  covars,
                                  me,
                                  B = B,
                                  type = "indep",
                                  method = "mle")
    colnames(boot$coef) <- names(out$coef)
    out$boot <- boot
  }
  out
}
mle_get_coef <- function(lm_fit, lmm_fit){
  coef_lm_fit <- stats::coef(lm_fit)
  fixed_ef_lmm_fit <- lme4::fixef(lmm_fit)
  sigma_sq <- summary(lm_fit)$sigma^2
  random_int <- attributes(lme4::VarCorr(lmm_fit)$id)$stddev^2
  beta <- mecor:::mle_coef(coef_lm_fit,
                           sigma_sq,
                           fixed_ef_lmm_fit,
                           random_int)
  beta
}
mle_coef <- function(coef_lm_fit,
                     sigma_sq,
                     fixed_ef_lmm_fit,
                     random_int){
  phi <- fixed_ef_lmm_fit[2] *
    sigma_sq / (random_int + fixed_ef_lmm_fit[2]^2 * sigma_sq)
  alpha <- coef_lm_fit[1] -
    phi * (fixed_ef_lmm_fit[1] +
             fixed_ef_lmm_fit[2] * coef_lm_fit[1])
  n_covars <- length(coef_lm_fit) - 1
  if (n_covars != 0){
    gamma <- coef_lm_fit[2:{1 + n_covars}] -
      phi * (fixed_ef_lmm_fit[3:{2 + n_covars}] +
               fixed_ef_lmm_fit[2] * coef_lm_fit[2:{1 + n_covars}])
  }
  beta <- c(phi,
            alpha,
            gamma)
  beta
}
mle_get_vcov <- function(lm_fit, lmm_fit){
  coef_lm_fit <- stats::coef(lm_fit)
  fixed_ef_lmm_fit <- lme4::fixef(lmm_fit)
  sigma_sq <- summary(lm_fit)$sigma^2
  random_int <- attributes(lme4::VarCorr(lmm_fit)$id)$stddev^2
  vec <- c(coef_lm_fit, sigma_sq, fixed_ef_lmm_fit, random_int)
  vcov_vec <- mecor:::mle_vcov_vec(lm_fit, lmm_fit)
  vcov_beta <- mecor:::deltamethod(mecor:::mle_get_coef_using_vec,
                                   vec,
                                   vcov_vec)
  names_vcov_beta <- names(fixed_ef_lmm_fit)
  names_vcov_beta[1:2] <- rev(names(fixed_ef_lmm_fit)[1:2])
  dimnames(vcov_beta) <- list(names_vcov_beta,
                              names_vcov_beta)
  vcov_beta
}
mle_vcov_vec <- function(lm_fit, lmm_fit){
  coef_lm_fit <- stats::coef(lm_fit)
  fixed_ef_lmm_fit <- lme4::fixef(lmm_fit)
  sigma_sq <- summary(lm_fit)$sigma^2
  random_int <- attributes(lme4::VarCorr(lmm_fit)$id)$stddev^2
  n_delta <- length(coef_lm_fit)
  n_kappa <- length(fixed_ef_lmm_fit)
  vcov <- matrix(nrow = {n_delta + n_kappa + 2},
                 ncol = {n_delta + n_kappa + 2}, 0)
  vcov_lm_fit <- vcov(lm_fit)
  vcov_fixed_ef <- vcov(lmm_fit)
  var_sigma_sq <- (2 * sigma_sq ^ 4) / (length(lm_fit$residuals) - 1)
  var_random_int <- lmerTest::as_lmerModLmerTest(lmm_fit)@vcov_varpar[1, 1]
  vcov[1:n_delta, 1:n_delta] <- vcov_lm_fit
  vcov[{n_delta + 1}, {n_delta + 1}] <- var_sigma_sq
  vcov[{n_delta + 2}:{n_delta + 1 + n_kappa},
       {n_delta + 2}:{n_delta + 1 + n_kappa}] <- as.matrix(vcov_fixed_ef)
  vcov[{n_delta + n_kappa + 2}:{n_delta + n_kappa + 2},
       {n_delta + n_kappa + 2}:{n_delta + n_kappa + 2}] <- var_random_int
  vcov
}
mle_get_coef_using_vec <- function(vec){
  n_delta <- floor((length(vec) - 2) / 2)
  n_kappa <- length(vec) - 2 - n_delta
  coef_lm_fit <- vec[1:n_delta]
  sigma_sq <- vec[n_delta + 1]
  fixed_ef_lmm_fit <- vec[{n_delta + 2}:{n_delta + 1 + n_kappa}]
  random_int <- vec[length(vec)]
  beta <- mecor:::mle_coef(coef_lm_fit,
                           sigma_sq,
                           fixed_ef_lmm_fit,
                           random_int)
}
