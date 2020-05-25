regcal <- function(response, covars, me, B = 0 , alpha = 0.05, type){
  # estimate beta_star (uncor) and its vcov
  uncor_fit <- mecor:::uncorrected(response, covars, me, type)
  beta_star <- mecor:::get_coefs(uncor_fit)
  vcov_beta_star <- mecor:::get_vcov(uncor_fit)
  # estimate calibration model
  calmod_fit <- mecor:::regcal_get_calmod(response, covars, me, type)
  calmod_coefs <- mecor:::get_coefs(calmod_fit, type == "indep")
  vcov_calmod <- mecor:::get_vcov(calmod_fit, type == "indep")
  # estimate beta (cor) and its vcov
  beta <- mecor:::regcal_get_coef(beta_star, calmod_coefs, type)
  vcov_beta <- mecor:::regcal_get_vcov(beta_star, calmod_coefs,
                                       vcov_beta_star, vcov_calmod, type)
  # change names of beta and its vcov
  beta <- mecor:::change_names(beta, me)
  vcov_beta <- mecor:::change_names(vcov_beta, me)
  # change order of coef and vcov matrix
  out <- list(coef = mecor:::change_order_coefs(beta),
              vcov = mecor:::change_order_vcov(vcov_beta))
  # bootstrap functionality
  if (B != 0){
    boot <- mecor:::regcal_boot(response, covars, me, B = B, alpha = alpha)
    out$boot <- list(ci = boot$ci,
                     vcov = boot$vcov)
  }
  out
}
# corrected coefs using regression calibration
regcal_get_coef <- function(beta_star, calmod_coefs, type){
  n <- NROW(beta_star)
  calmod_matrix <- mecor:::regcal_get_calmod_matrix(calmod_coefs, n, type)
  # estimate beta (cor)
  if (type == "dep"){
    beta_star <- c(beta_star, 1)
  }
  beta <- as.numeric(beta_star %*% solve(calmod_matrix))
  # add names
  names(beta) <- names(beta_star)
  # output
  beta[1:n]
}
# estimate calibration model
regcal_get_calmod <- function(response, covars, me, type){
  # get design matrix of calibration model
  dm_cm <- mecor:::get_dm_cm(covars, me, type) # design matrix calibration model
  if (type == "dep"){
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
# create measurement error matrix
regcal_get_calmod_matrix <- function(calmod_coefs, n, type){
  if (type == "dep"){
    calmod_matrix <- matrix(0, nrow = n + 1, ncol = n + 1)
    diag(calmod_matrix) <- calmod_coefs[2]
    calmod_matrix[n + 1, 2] <- calmod_coefs[1]
    calmod_matrix[n + 1, n + 1] <- 1
  } else if (type == "indep"){
    calmod_matrix <- diag(n)
    calmod_matrix[1,] <- calmod_coefs
  }
  calmod_matrix
}
change_names <- function(x, me){
  if (is.null(me$replicate)){
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
