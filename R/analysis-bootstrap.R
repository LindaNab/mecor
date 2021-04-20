# bootstrap
analysis_boot <- function(response,
                          covars,
                          me,
                          B = 999,
                          type,
                          method) {
  # allocate memory for samples
  samples <- vector(mode = "list",
                    length = B)
  for (i in seq_along(samples)) {
    samples[[i]] <-
      get_new_sample(me,
                     response,
                     covars,
                     type,
                     method)
  }
  n_coef <- ifelse(is.null(covars), 0, ncol(covars)) +
    ifelse(type == "indep", 2, 1)
  pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
  coef <- matrix(ncol = B, nrow = n_coef)
  switch(
    method,
    "standard" = for (i in seq_along(samples)) {
      utils::setTxtProgressBar(pb, i)
      coef[, i] <-
        do.call(standard, samples[[i]])$coef
    },
    "efficient" = for (i in seq_along(samples)) {
      utils::setTxtProgressBar(pb, i)
      coef[, i] <-
        do.call(efficient, samples[[i]])$coef
    },
    "valregcal" = for (i in seq_along(samples)) {
      utils::setTxtProgressBar(pb, i)
      coef[, i] <-
        do.call(valregcal, samples[[i]])$coef
    },
    "mle" = for (i in seq_along(samples)) {
      utils::setTxtProgressBar(pb, i)
      coef[, i] <-
        do.call(mle, samples[[i]])$coef
    }
  )
  close(pb)
  out <- list(coef = t(coef))
}
get_new_sample <- function(x, ...) {
  UseMethod("get_new_sample")
}
get_new_sample.MeasError <- function(me,
                                     response,
                                     covars,
                                     type,
                                     method) {
  # sample new rownums
  rownum_filled <- which (!is.na(me$reference))
  rownum_empty <- which (is.na(me$reference))
  new_rownum_filled <- sample(rownum_filled,
                              size = NROW(rownum_filled),
                              replace = T)
  new_rownum_empty <- sample(rownum_empty,
                             size = NROW(rownum_empty),
                             replace = T)
  new_rownums <- c(new_rownum_filled, new_rownum_empty)
  # sample new variables
  new_covars <- covars[new_rownums, , drop = F]
  new_me <- me
  new_me$reference <- me$reference[new_rownums]
  new_me$substitute <- me$substitute[new_rownums]
  if (type == "dep_diff") {
    new_me$differential <- me$differential[new_rownums]
  }
  if (!is.null(me$replicate)) {
    if (is.null(dim(me$replicate))) {
      new_me$replicate <- me$replicate[new_rownums]
    } else {
      new_me$replicate <- me$replicate[new_rownums,]
    }
  }
  new_sample <- list(
    covars = new_covars,
    me = new_me,
    B = 0,
    # no bootstrapping within bootstrap
    type = type
  )
  if (type == "indep") {
    new_response <- response[new_rownums, , drop = F]
    new_sample$response = new_response
    new_sample <-
      new_sample[c("response", "covars", "me", "B", "type")]
  }
  if (method %in% c("standard", "efficient", "mle"))
    new_sample$calc_vcov = F
  new_sample
}
get_new_sample.MeasErrorRandom <- function(me,
                                           response,
                                           covars,
                                           type,
                                           method){
  new_rownums <- sample(1:NROW(me$substitute),
                        size = NROW(me$substitute),
                        replace = T)
  new_me <- me
  new_me$substitute <- new_me$substitute[new_rownums]
  new_sample <- list(
    response = response[new_rownums, , drop = F],
    covars = covars[new_rownums, , drop = F],
    me = new_me,
    B = 0, # no bootstrapping within bootstrap
    type = type
  )
  new_sample
}
get_new_sample.MeasErrorExt <- function(me,
                                        response,
                                        covars,
                                        type,
                                        method) {
  new_rownums <- sample(1:NROW(me$substitute),
                        size = NROW(me$substitute),
                        replace = T)
  new_covars <- covars[new_rownums, , drop = F]
  new_me <- me
  new_me$substitute <- new_me$substitute[new_rownums]
  new_rownums_ext_data <- sample(1:nrow(me$model$model),
                                 size = nrow(me$model$model),
                                 replace = T)
  new_ext_data <- me$model$model[new_rownums_ext_data,]
  new_me$model <- stats::lm(stats::terms(me$model), data = new_ext_data)
  new_sample <- list(
    covars = new_covars,
    me = new_me,
    B = 0, # no bootstrapping within bootstrap
    type = type
  )
  if (type == "indep") {
    new_response <- response[new_rownums, , drop = F]
    new_sample$response = new_response
    new_sample <-
      new_sample[c("response", "covars", "me", "B", "type")]
  }
  if (method %in% c("standard", "efficient", "mle"))
    new_sample$calc_vcov = F
  new_sample
}
