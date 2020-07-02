# bootstrap
analysis_boot <- function(response,
                          covars,
                          me,
                          B = 999,
                          alpha = 0.05,
                          type,
                          method){
  samples <- replicate(
    B,
    mecor:::get_new_sample(response, covars, me, type, method),
    simplify = F
  )
  n_coef <- ncol(covars) + 2
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  coef <- switch(method,
                 "rc" = vapply(seq_along(samples),
                               FUN = function(i) {
                                 setTxtProgressBar(pb, i)
                                 do.call(mecor:::regcal, samples[[i]])$coef},
                               FUN.VALUE = double(n_coef)),
                 "erc" = vapply(seq_along(samples),
                                FUN = function(i) {
                                  setTxtProgressBar(pb, i)
                                  do.call(mecor:::efficient_regcal, samples[[i]])$coef},
                                FUN.VALUE = double(n_coef)),
                 "irc" = sapply(samples,
                                FUN = function(x) do.call(mecor:::inadmissible_regcal, x)$coef
                               ),
                 "mle" = sapply(samples,
                                FUN = function(x) do.call(mecor:::mle, x)$coef
                               ))
  ci_perc <- apply(coef,
                   1,
                   FUN = quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))
  out <- list(ci = ci_perc,
              vcov = cov(t(coef)))
}
get_new_sample <- function(response, covars, me, type, method){
  UseMethod("get_new_sample", me)
}
get_new_sample.MeasError <- function(response, covars, me, type, method){
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
  if (type == "dep_diff"){
    new_me$differential <- me$differential[new_rownums]
  }
  if (!is.null(me$replicate)){
    new_me$replicate <- me$replicate[new_rownums, ]
  }
  new_sample <- list(
    covars = new_covars,
    me = new_me,
    B = 0, # no bootstrapping within bootstrap
    type = type
  )
  if (type == "indep"){
    new_response <- response[new_rownums, , drop = F]
    new_sample$response = new_response
    new_sample <- new_sample[c("response", "covars", "me", "B", "type")]
  }
  if (method %in% c("rc", "erc"))
    new_sample$calc_vcov = F
  new_sample
}
get_new_sample.MeasErrorExt <- function(response, covars, me, type){
  new_rownums <- sample(1:NROW(me$substitute),
                        size = NROW(me$substitute),
                        replace = T)
  new_covars <- covars[new_rownums, , drop = F]
  new_me <- me
  new_me$substitute <- new_me$substitute[new_rownums]
  new_rownums_ext_data <- sample(1:nrow(me$model$model),
                                 size = nrow(me$model$model),
                                 replace = T)
  new_ext_data <- me$model$model[new_rownums_ext_data, ]
  new_me$model <- lm(terms(me$model), data = new_ext_data)
  new_sample <- list(
    covars = new_covars,
    me = new_me,
    B = 0, # no bootstrapping within bootstrap
    type = type
  )
  if (type == "indep"){
    new_response <- response[new_rownums, , drop = F]
    new_sample$response = new_response
    new_sample <- new_sample[c("response", "covars", "me", "B", "type")]
  }
  if (method %in% c("rc", "erc"))
    new_sample$calc_vcov = F
  new_sample
}
