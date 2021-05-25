# designmatrices
# design matrix for uncorrected analysis
get_dm_uncor <- function(covars,
                         me,
                         type) {
  if (startsWith(type, "dep")) {
    dm <- cbind(1, covars)
    colnames(dm)[1] <- "(Intercept)"
  } else if (type == "indep") {
    dm <- cbind(1, me$substitute)
    colnames(dm) <-
      c("(Intercept)", attributes(me)$input$substitute)
    dm <- bind_covars(dm, covars)
  }
  dm
}
# design matrix for complete case analysis
get_dm_cc <- function(covars,
                      me,
                      type) {
  if (startsWith(type, "dep")) {
    dm <- get_dm_uncor(covars, me, type)
  } else if (type == "indep") {
    dm <- cbind(1, me$reference)
    name <- ifelse(
      !is.null(me$replicate),
      paste0("cor_", attributes(me)$input$substitute),
      attributes(me)$input$reference
    )
    colnames(dm) <- c("(Intercept)", name)
    dm <- bind_covars(dm, covars)
  }
  dm
}
# design matrix for calibration model (covariate-me) or measurement error model
# (outcome-me)
get_dm_cm <- function(covars,
                      me,
                      type) {
  if (startsWith(type, "dep")) {
    dm <- cbind(1, me$reference)
    name <- attributes(me)$input$reference
    colnames(dm) <- c("(Intercept)", name)
    if (type == "dep_diff") {
      dm <- cbind(dm, me$differential, me$reference * me$differential)
      colnames(dm)[3] <-
        as.character(attributes(me)$input$differential)
      colnames(dm)[4] <-
        paste0(colnames(dm)[2], ":", colnames(dm)[3])
    }
  } else if (type == "indep") {
    dm <- get_dm_uncor(covars,
                       me,
                       type)
  }
  dm
}
# design matrix for validation regression calibration, using the references
# measures if available and the expected values of the calibration model
# if not
get_dm_vrc <- function(covars,
                       me,
                       type) {
  dm <- get_dm_uncor(covars,
                             me,
                             type)
  model_fit <- standard_get_model(covars,
                                          me,
                                          type)
  fitted_values <- dm %*% stats::coef(model_fit)
  x <- fitted_values
  x[!is.na(me$reference)] <- me$reference[!is.na(me$reference)]
  dm[, as.character(attributes(me)$input$substitute)] <- x
  ind <-
    colnames(dm) == as.character(attributes(me)$input$substitute)
  name <- as.character(attributes(me)$input$reference)
  colnames(dm)[ind] <- name
  dm
}
bind_covars <- function(dm,
                        covars) {
  if (!is.null(covars)) {
    # if there are covars, add these to design matrix
    dm <- cbind(dm, covars)
  }
  else
    dm
}
