# bootstrap
analysis_boot <- function(response,
                          covars,
                          me,
                          B = 999,
                          alpha = 0.05,
                          type,
                          method,
                          erc_B = 0){
  strat_samples <- replicate(
    B,
    mecor:::get_strat_sample(response, covars, me, type, method, erc_B),
    simplify = F
  )
  if (method == "rc"){
    coef <- sapply(
      strat_samples,
      FUN = function(x) do.call(mecor:::regcal, x)$coef
    )
  } else if (method == "erc"){
    coef <- sapply(
      strat_samples,
      FUN = function(x) do.call(mecor:::efficient_regcal, x)$coef
    )
  } else if (method == "irc"){
    coef <- sapply(
      strat_samples,
      FUN = function(x) do.call(mecor:::inadmissible_regcal, x)$coef
    )
  }
  ci_perc <- apply(coef,
                   1,
                   FUN = quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))
  out <- list(ci = ci_perc,
              vcov = cov(t(coef)))
}
get_strat_sample <- function(response, covars, me, type, method, erc_B){
  rownum_filled <- which (!is.na(me$reference))
  rownum_empty <- which (is.na(me$reference))
  new_rownum_filled <- sample(rownum_filled,
                              size = NROW(rownum_filled),
                              replace = T)
  new_rownum_empty <- sample(rownum_empty,
                             size = NROW(rownum_empty),
                             replace = T)
  new_rownums <- c(new_rownum_filled, new_rownum_empty)
  if (type == "indep"){
    new_response <- response[new_rownums, , drop = F]
  }
  new_covars <- covars[new_rownums, , drop = F]
  new_me <- me
  new_me$reference <- new_me$reference[new_rownums]
  new_me$substitute <- new_me$substitute[new_rownums]
  new_sample <- list(
    covars = new_covars,
    me = new_me,
    B = 0
  ) # no bootstrapping within bootstrap
  if (exists("new_response")){
    new_sample$response = new_response
    new_sample <- new_sample[c("response", "covars", "me", "B")]
  }
  if (method == "rc" | method == "erc"){
    new_sample$type <- type
  }
  if (method == "erc"){
    new_sample$erc_B <- erc_B
  }
  new_sample
}
