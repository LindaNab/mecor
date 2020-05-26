complete_case <- function(response, covars, me, type){
  dm_cc <- mecor:::get_dm_cc(covars, me, type)
  if (type == "dep"){
    y <- me$reference
  } else if (type == "indep"){
    y <- as.numeric(response)
  }
  y <- y[!is.na(me$reference)]
  dm_cc <- dm_cc[!is.na(me$reference), ]
  cc_fit <- stats::lm.fit(dm_cc, y)
  cc_fit
}
