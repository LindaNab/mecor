regcal <- function(response, covars, me, B = 0 , alpha = 0.05, type, formula_type){
  # estimate beta_star (uncor) and its vcov
  uncor_fit <- mecor:::uncorrected(response, covars, me, type)
  beta_star <- mecor:::get_coefs(uncor_fit)
  if (type == "dep_diff"){
    dm_uncor <- mecor:::get_dm_uncor(covars, me, type)
    vcov_beta_star <- mecor:::get_vcovHC3(uncor_fit, dm_uncor)
  } else vcov_beta_star <- mecor:::get_vcov(uncor_fit)
  # estimate calibration model
  if (class(me)[1] == "MeasError"){
    calmod_fit <- mecor:::regcal_get_calmod(covars, me, type)
    coef_calmod <- mecor:::get_coefs(calmod_fit, type == "indep")
    vcov_calmod <- mecor:::get_vcov(calmod_fit, type == "indep")
  } else if (class(me)[1] == "MeasErrorExt"){
      coef_calmod <- mecor:::get_coefs(me$model, type == "indep")
      if(length(grep("MeasErrorExt.lm", attributes(me)$call)) != 0){
        vcov_calmod <- mecor:::get_vcov(me$model, type == "indep")
      } else if(length(grep("MeasErrorExt.list", attributes(me)$call)) != 0){
        if (!is.null(me$model$vcov)){
          if (type == "indep") {
            vcov_calmod <- mecor:::change_order_vcov(me$vcov)
          } else vcov_calmod <- me$vcov
        }
        if (B != 0){
          B <- 0
          warning("B set to 0 since bootstrap cannot be used if the class of 'model' in the MeasErrorExt object is list")
        }
      }
  }
  # estimate beta (cor) and its vcov
  beta <- mecor:::regcal_get_coef(beta_star, coef_calmod, type)
  if (exists("vcov_calmod")){
    vcov_beta <- mecor:::regcal_get_vcov(beta_star, coef_calmod,
                                         vcov_beta_star, vcov_calmod, type)
  }
  # change names of beta and its vcov
  if (type == "indep"){
    beta <- mecor:::change_names(beta, me)
    if (exists("vcov_beta")){
      vcov_beta <- mecor:::change_names(vcov_beta, me)
    }
  }
  # change order of coef and vcov matrix
  out <- list(coef = mecor:::change_order_coefs(beta))
  if (exists("vcov_beta")){
    out$vcov <- mecor:::change_order_vcov(vcov_beta)
  }
  n <- NROW(beta_star)
  calmod_matrix <- mecor:::regcal_get_calmod_matrix(coef_calmod, n, type)
  out$matrix <- calmod_matrix
  # bootstrap functionality
  if (B != 0){
    boot <-
      mecor:::analysis_boot(response, covars, me,
                            B = B, alpha = alpha, type = type, method = "rc")
    out$boot <- list(ci = boot$ci,
                     vcov = boot$vcov)
  }
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
      calmod_matrix[n + 1, 2] <- coef_calmod[1]}
    else if (type == "dep_diff"){
      calmod_matrix[2, 1] <- coef_calmod[4]
      calmod_matrix[3, 1] <- coef_calmod[3]
      calmod_matrix[2, 2] <- coef_calmod[2]
      calmod_matrix[3, 2] <- coef_calmod[1]
      calmod_matrix[1, 1] <- coef_calmod[4] + coef_calmod[2]
    }
    calmod_matrix[n + 1, n + 1] <- 1
  } else if (type == "indep"){
    calmod_matrix <- diag(n)
    calmod_matrix[1,] <- coef_calmod
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
