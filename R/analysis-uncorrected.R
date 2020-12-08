uncorrected <- function(response,
                        covars,
                        me,
                        type) {
  dm_uncor <- get_dm_uncor(covars, me, type)
  if (startsWith(type, "dep")) {
    uncor_fit <- stats::lm.fit(dm_uncor, me$substitute)
  } else if (type == "indep") {
    uncor_fit <- stats::lm.fit(dm_uncor, response)
  }
  uncor_fit
}
