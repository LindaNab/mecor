naive <- function(response, covars, me){
  dm_naive <- mecor:::get_dm_naive(covars, me)
  naive_fit <- stats::lm.fit(dm_naive, response)
  naive_fit
}
