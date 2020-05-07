# designmatrices

# design matrix for naive analysis
get_dm_naive <- function(covars, me){
  dm <- cbind(1, me$substitute)
  colnames_dm <- c("(Intercept)", attributes(me)$input$substitute)
  if (!is.null(covars)){ # if there are covars, add these to design matrix
    dm <- cbind(dm, covars)
    colnames_dm <- c(colnames_dm, colnames(covars))
  }
  colnames(dm) <- colnames_dm
  dm
}

# design matrix for complete case analysis
get_dm_cc <- function(covars, me){
  if (attributes(me)$type == "ivs"){
    dm <- cbind(1, me$reference)
    colnames_dm <- c("(Intercept)", attributes(me)$input$reference)
    if (!is.null(covars)){ # if there are covars, add these to design matrix
      dm <- cbind(dm, covars)
      colnames_dm <- c(colnames_dm, colnames(covars))
    }
  }
  colnames(dm) <- colnames_dm
  dm
}

# design matrix for calibration model
get_dm_cm <- function(covars, me){
  dm <- get_dm_naive(covars, me)
  dm
}
