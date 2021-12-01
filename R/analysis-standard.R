standard <- function(response,
                     covars,
                     me,
                     B = 0,
                     type,
                     calc_vcov = TRUE) {
  # turned off when used in bootstrap
  # estimate beta_star (uncor) and its vcov
  uncor_fit <- uncorrected(response,
                           covars,
                           me,
                           type)
  beta_star <- get_coefs(uncor_fit)
  if (calc_vcov && type == "dep_diff") {
    dm_uncor <- get_dm_uncor(covars,
                             me,
                             type)
    vcov_beta_star <- get_vcovHC3(uncor_fit,
                                  dm_uncor)
  } else if (calc_vcov) vcov_beta_star <- get_vcov(uncor_fit)
  # estimate calibration model (cov-me) or measurement error model (out-me)
  model <- model(response,
                 covars,
                 me,
                 type,
                 calc_vcov)
  coef_model <- model$coef
  if (!is.null(model$vcov)) {
    vcov_model <- model$vcov
  }
  # estimate beta (cor) and its vcov
  beta <- standard_get_coef(beta_star,
                            coef_model,
                            type)
  if (type == "indep") {
    beta <- change_names(beta,
                         me)
  }
  out <- list(coef = change_order_coefs(beta))
  n <- NROW(beta_star)
  model_matrix <- standard_get_model_matrix(coef_model,
                                            n,
                                            type)
  out$matrix <- model_matrix
  if (calc_vcov) {
    if (exists("vcov_model")) {
      # delta method
      vcov_beta <- standard_get_vcov(beta_star,
                                     coef_model,
                                     vcov_beta_star,
                                     vcov_model,
                                     type)
      if (type == "indep")
        vcov_beta <- change_names(vcov_beta,
                                  me)
      out$vcov <- change_order_vcov(vcov_beta)
    }
    # vcov matrix obtained by ignoring uncertainty in coefs of model_matrix
    zerovar_vcov <- standard_zerovar(vcov_beta_star,
                                     model_matrix,
                                     type)
    if (type == "indep")
      zerovar_vcov <- change_names(zerovar_vcov,
                                   me)
    out$zerovar_vcov <- change_order_vcov(zerovar_vcov)
    # elements to build fieller ci's
    if (exists("vcov_model") && type != "dep_diff") {
      fieller <- standard_fieller(beta_star,
                                  coef_model,
                                  vcov_beta_star,
                                  vcov_model,
                                  type)
      out$fieller <- fieller
    }
  }
  out$method <- standard_get_method(type)
  # bootstrap
  if (B != 0) {
    boot <- analysis_boot(response,
                          covars,
                          me,
                          B = B,
                          type = type,
                          method = "standard")
    colnames(boot$coef) <- names(out$coef)
    out$boot <- boot
  }
  out
}
model <- function(response, covars, me, ...) {
  UseMethod("model", me)
}
model.MeasError <- function(response,
                            covars,
                            me,
                            type,
                            calc_vcov = TRUE) {
  if (type == "dep" & !is.null(me$replicate)){ # outcome calibration study
    cc <- which(!is.na(me$reference)) # select the individuals of whom all k replicate measures are available
    me2 <- MeasError(substitute = me$replicate[cc, 1], replicate = me$replicate[cc, -1]) # has to be >1, which is checked earlier in check_me
    corfit <- standard(
      response = me$substitute[cc],
      covars = NULL,
      me2,
      B = 0,
      type = "indep",
      calc_vcov
    ) # correct the regression of Y*,s ~ Y*,c with standard regression calibration to get correct model matrix
    coef_model <- corfit$coef
    out <- list(coef = coef_model)
    if (calc_vcov){
      vcov_model <- corfit$vcov
      out$vcov <- vcov_model
    }
  } else {
    model_fit <- standard_get_model(covars,
                                  me,
                                  type)
    coef_model <- get_coefs(model_fit,
                            type == "indep")
    out <- list(coef = coef_model)
    if (calc_vcov){
      vcov_model <- get_vcov(model_fit,
                             type == "indep")
      out$vcov <- vcov_model
    }
  }
  out
}
model.MeasErrorExt <- function(response,
                               covars,
                               me,
                               type,
                               calc_vcov = TRUE) {
  coef_model <- get_coefs(me$model,
                          type == "indep")
  out <- list(coef = coef_model)
  if (calc_vcov) {
    if (length(grep("MeasErrorExt.lm", attributes(me)$call)) != 0) {
      vcov_model <- get_vcov(me$model,
                             type == "indep")
      out$vcov <- vcov_model
    } else if (length(grep("MeasErrorExt.list", attributes(me)$call)) != 0) {
      if (!is.null(me$model$vcov)) {
        if (type == "indep") {
          vcov_model <- change_order_vcov(me$vcov)
        } else
          vcov_model <- me$vcov
        out$vcov <- vcov_model
      }
    }
  }
  out
}
model.MeasErrorRandom <- function(response,
                                  covars,
                                  me,
                                  type,
                                  calc_vcov = TRUE) {
  if (!is.null(covars)) {
    if (length(me$substitute) != nrow(covars)) {
      stop("substitute and other covariates need to be of same length")
    } else q <- cbind(me$substitute, covars)
  } else q <- me$substitute
  Q <- scale(q, scale = F)
  matrix <- t(Q) %*% Q / (length(me$substitute) - 1)
  matrix1 <- matrix
  matrix1[1, 1] <- matrix1[1, 1] - me$variance
  model_matrix <- solve(matrix1) %*% matrix
  # (1/lambda1        0
  #  -lambda2/lambda1 1)
  n_model_matrix <- nrow(model_matrix)
  lambda1 <- 1 / model_matrix[1, 1]
  if (!is.null(covars)) {
    lambda2 <- model_matrix[2:n_model_matrix, 1] * - lambda1
    lambda0 <- mean(me$substitute) - lambda1 * mean(me$substitute) -
      t(lambda2) %*% colMeans(covars)
    out <- list(coef = c(lambda1, lambda0, lambda2))
  } else {
    lambda0 <- mean(me$substitute) - lambda1 * mean(me$substitute)
    out <- list(coef = c(lambda1, lambda0))
  }
  out
}
# estimate calibration model/measurement error model
standard_get_model <- function(covars,
                               me,
                               type) {
  dm_cm <- get_dm_cm(covars,
                     me,
                     type) # design matrix calibration model (cov-me)
  # or measurement error model (outcome-me)
  if (startsWith(type, "dep")) {
    y <- me$substitute
  } else if (type == "indep") {
    y <- me$reference
  }
  y <- y[!is.na(me$reference)]
  dm_cm <- dm_cm[!is.na(me$reference),] # select complete cases
  model_fit <- stats::lm.fit(dm_cm, y)
  # output
  model_fit
}
# corrected coefs using regression calibration
standard_get_coef <- function(beta_star,
                              coef_model,
                              type) {
  n <- NROW(beta_star)
  model_matrix <- standard_get_model_matrix(coef_model,
                                            n,
                                            type)
  # estimate beta (cor)
  if (startsWith(type, "dep")) {
    beta_star <- c(beta_star, 1)
  }
  beta <- as.numeric(beta_star %*% solve(model_matrix))
  # add names
  names(beta) <- names(beta_star)
  # output
  beta[1:n]
}
# create calibration model matrix Lambda (cov-me)
# or measurement error model matrix Theta (out-me)
standard_get_model_matrix <- function(coef_model,
                                      n,
                                      type) {
  if (startsWith(type, "dep")) {
    model_matrix <- matrix(0, nrow = n + 1, ncol = n + 1)
    if (type == "dep") {
      diag(model_matrix) <- coef_model[2]
      model_matrix[n + 1, 2] <- coef_model[1]
    } else if (type == "dep_diff") {
      model_matrix[2, 1] <- coef_model[4]
      model_matrix[3, 1] <- coef_model[3]
      model_matrix[2, 2] <- coef_model[2]
      model_matrix[3, 2] <- coef_model[1]
      model_matrix[1, 1] <- coef_model[4] + coef_model[2]
    }
    model_matrix[n + 1, n + 1] <- 1
    dimnames_model_matrix <- paste0("Theta", 1:{
      n + 1
    })
    dimnames(model_matrix) <- list(dimnames_model_matrix,
                                   dimnames_model_matrix)
  } else if (type == "indep") {
    model_matrix <- diag(n)
    model_matrix[1, ] <- coef_model
    dimnames_model_matrix <- paste0("Lambda", 1:n)
    dimnames_model_matrix[2] <- "Lambda0"
    dimnames(model_matrix) <- list(dimnames_model_matrix,
                                   dimnames_model_matrix)
  }
  model_matrix
}
change_names <- function(x,
                         me) {
  if (class(me)[1] == "MeasError" && is.null(me$replicate)) {
    name_reference <- as.character(attributes(me)$input$reference)
  } else {
    name_reference <-
      paste0("cor_", attributes(me)$input$substitute)[1] # add cor_
  }
  if (is.matrix(x)) {
    colnames(x)[1] <- name_reference
    rownames(x)[1] <- name_reference
  } else{
    names(x)[1] <- name_reference
  }
  x
}
standard_get_method <- function(type) {
  if (startsWith(type, "dep")) {
    method <- "method of moments"
  } else if (type == "indep") {
    method <- "regression calibration"
  }
  method
}
