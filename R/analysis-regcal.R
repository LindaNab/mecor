regcal <- function(response,
                   covars,
                   me,
                   B = 0,
                   type,
                   calc_vcov = T){
  # estimate beta_star (uncor) and its vcov
  uncor_fit <- mecor:::uncorrected(response, covars, me, type)
  beta_star <- mecor:::get_coefs(uncor_fit)
  if (calc_vcov && type == "dep_diff"){
    dm_uncor <- mecor:::get_dm_uncor(covars, me, type)
    vcov_beta_star <- mecor:::get_vcovHC3(uncor_fit, dm_uncor)
  } else if (calc_vcov) vcov_beta_star <- mecor:::get_vcov(uncor_fit)
  # estimate calibration model
  calmod <- mecor:::calmod(response, covars, me, type, calc_vcov)
  coef_calmod <- calmod$coef
  if (!is.null(calmod$vcov)){
    vcov_calmod <- calmod$vcov
  }
  # estimate beta (cor) and its vcov
  beta <- mecor:::regcal_get_coef(beta_star, coef_calmod, type)
  if (type == "indep") beta <- mecor:::change_names(beta, me)
  out <- list(coef = mecor:::change_order_coefs(beta))
  n <- NROW(beta_star)
  calmod_matrix <- mecor:::regcal_get_calmod_matrix(coef_calmod, n, type)
  out$matrix <- calmod_matrix
  if (calc_vcov){
    if (exists("vcov_calmod")){
      # delta method
      vcov_beta <- mecor:::regcal_get_vcov(beta_star,
                                           coef_calmod,
                                           vcov_beta_star,
                                           vcov_calmod,
                                           type)
      if (type == "indep") vcov_beta <- mecor:::change_names(vcov_beta, me)
      out$vcov <- mecor:::change_order_vcov(vcov_beta)
    }
    # vcov matrix obtained by ignoring uncertainty in coef_calmod
    zerovar_vcov <- mecor:::regcal_zerovar(vcov_beta_star,
                                           calmod_matrix,
                                           type)
    if (type == "indep") zerovar_vcov <- mecor:::change_names(zerovar_vcov, me)
    out$zerovar_vcov <- mecor:::change_order_vcov(zerovar_vcov)
    # confidence intervals obtained by fieller method
    if (exists("vcov_calmod") && type != "dep_diff"){
      fieller <- mecor:::regcal_fieller(beta_star,
                                        coef_calmod,
                                        vcov_beta_star,
                                        vcov_calmod,
                                        type)
      out$fieller <- fieller
    }
  }
  # bootstrap
  if (B != 0){
    boot <- mecor:::analysis_boot(response,
                                  covars,
                                  me,
                                  B = B,
                                  type = type,
                                  method = "rc")
    colnames(boot$coef) <- names(out$coef)
    out$boot <- boot
  }
  out
}
calmod <- function(response,
                   covars,
                   me,
                   type,
                   calc_vcov = T){
  UseMethod("calmod", me)
}
calmod.MeasError <- function(response,
                             covars,
                             me,
                             type,
                             calc_vcov = T){
  calmod_fit <- mecor:::regcal_get_calmod(covars, me, type)
  coef_calmod <- mecor:::get_coefs(calmod_fit, type == "indep")
  out <- list(coef = coef_calmod)
  if (calc_vcov){
    vcov_calmod <- mecor:::get_vcov(calmod_fit, type == "indep")
    out$vcov <- vcov_calmod
  }
  out
}
calmod.MeasErrorExt <- function(response,
                                covars,
                                me,
                                type,
                                calc_vcov = T){
  coef_calmod <- mecor:::get_coefs(me$model, type == "indep")
  out <- list(coef = coef_calmod)
  if (calc_vcov){
    if(length(grep("MeasErrorExt.lm", attributes(me)$call)) != 0){
      vcov_calmod <- mecor:::get_vcov(me$model, type == "indep")
      out$vcov <- vcov_calmod
    } else if(length(grep("MeasErrorExt.list", attributes(me)$call)) != 0){
      if (!is.null(me$model$vcov)){
        if (type == "indep"){
          vcov_calmod <- mecor:::change_order_vcov(me$vcov)
        } else vcov_calmod <- me$vcov
        out$vcov <- vcov_calmod
      }
    }
  }
  out
}
calmod.MeasErrorRandom <- function(response,
                                   covars,
                                   me,
                                   type,
                                   calc_vcov = T){
  if(length(me$substitute) != nrow(covars)){
    stop("substitute and other covariates need to be of same length")
  }
  Q <- scale(cbind(me$substitute, covars), scale = F)
  matrix <- t(Q) %*% Q / (length(me$substitute) - 1)
  matrix1 <- matrix
  matrix1[1, 1] <- matrix1[1, 1] - me$error
  calmod_matrix <- solve(matrix1) %*% matrix
  # (1/lambda1        0
  #  -lambda2/lambda1 1)
  n_calmod_matrix <- nrow(calmod_matrix)
  lambda1 <- 1 / calmod_matrix[1, 1]
  lambda2 <- calmod_matrix[2:n_calmod_matrix, 1] * - lambda1
  lambda0 <- mean(me$substitute) - lambda1 * mean(me$substitute) -
    t(lambda2) %*% colMeans(covars)
  out <- list(coef = c(lambda1, lambda0, lambda2))
  out
}
# estimate calibration model
regcal_get_calmod <- function(covars, me, type){
  # get design matrix of calibration model
  dm_cm <- mecor:::get_dm_cm(covars, me, type) # design matrix calibration model
  if (startsWith(type, "dep")){
    y <- me$substitute
  } else if (type == "indep"){
    y <- me$reference
  }
  y <- y[!is.na(me$reference)]
  dm_cm <- dm_cm[!is.na(me$reference), ] # select complete cases
  calmod_fit <- stats::lm.fit(dm_cm, y)
  # output
  calmod_fit
}
# corrected coefs using regression calibration
regcal_get_coef <- function(beta_star, coef_calmod, type){
  n <- NROW(beta_star)
  calmod_matrix <- mecor:::regcal_get_calmod_matrix(coef_calmod, n, type)
  # estimate beta (cor)
  if (startsWith(type, "dep")){
    beta_star <- c(beta_star, 1)
  }
  beta <- as.numeric(beta_star %*% solve(calmod_matrix))
  # add names
  names(beta) <- names(beta_star)
  # output
  beta[1:n]
}
# create measurement error matrix
regcal_get_calmod_matrix <- function(coef_calmod, n, type){
  if (startsWith(type, "dep")){
    calmod_matrix <- matrix(0, nrow = n + 1, ncol = n + 1)
    if (type == "dep"){
      diag(calmod_matrix) <- coef_calmod[2]
      calmod_matrix[n + 1, 2] <- coef_calmod[1]
    } else if (type == "dep_diff"){
      calmod_matrix[2, 1] <- coef_calmod[4]
      calmod_matrix[3, 1] <- coef_calmod[3]
      calmod_matrix[2, 2] <- coef_calmod[2]
      calmod_matrix[3, 2] <- coef_calmod[1]
      calmod_matrix[1, 1] <- coef_calmod[4] + coef_calmod[2]
    }
    calmod_matrix[n + 1, n + 1] <- 1
    dimnames_calmod_matrix <- paste0("Theta", 1:{n + 1})
    dimnames(calmod_matrix) <- list(dimnames_calmod_matrix,
                                    dimnames_calmod_matrix)
  } else if (type == "indep"){
    calmod_matrix <- diag(n)
    calmod_matrix[1,] <- coef_calmod
    dimnames_calmod_matrix <- paste0("Lambda", 1:n)
    dimnames(calmod_matrix) <- list(dimnames_calmod_matrix,
                                    dimnames_calmod_matrix)
  }
  calmod_matrix
}
change_names <- function(x, me){
  if(class(me)[1] == "MeasError" && is.null(me$replicate)){
    name_reference <- as.character(attributes(me)$input$reference)
  } else {
    name_reference <- paste0("cor_", attributes(me)$input$substitute)
  }
  if (is.matrix(x)){
    colnames(x)[1] <- name_reference
    rownames(x)[1] <- name_reference
  } else{
    names(x)[1] <- name_reference
  }
  # output
  x
}
