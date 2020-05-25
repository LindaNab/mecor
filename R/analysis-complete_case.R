complete_case <- function(response, covars, me){
  dm_cc <- mecor:::get_dm_cc(covars, me)
  response <- response[!is.na(me$reference), ]
  dm_cc <- dm_cc[!is.na(me$reference), ]
  cc_fit <- stats::lm.fit(dm_cc, response)
  cc_fit
}
