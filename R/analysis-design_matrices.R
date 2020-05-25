# designmatrices

# design matrix for uncor analysis
get_dm_uncor <- function(covars, me, type){
  if (type == "dep"){
    dm <- cbind(1, covars)
    colnames(dm)[1] <- "(Intercept)"
  } else if (type == "indep"){
    dm <- cbind(1, me$substitute)
    colnames(dm) <- c("(Intercept)", attributes(me)$input$substitute)
    dm <- mecor:::bind_covars(dm, covars)
  }
  dm
}

# design matrix for complete case analysis
get_dm_cc <- function(covars, me, type){
  if (type == "dep"){
    dm <- get_dm_uncor(covars, me, type)
  } else if (type == "indep"){
    dm <- cbind(1, me$reference)
    colnames_dm <- c("(Intercept)", attributes(me)$input$reference)
    dm <- bind_covars(dm, covars)
  }
  dm
}

# design matrix for calibration model
get_dm_cm <- function(covars, me, type){
  if (type == "dep"){
    dm <- cbind(1, me$reference)
    colnames(dm) <- c("(Intercept)", attributes(me)$input$reference)
  } else if (type == "indep"){
    dm <- get_dm_uncor(covars, me, type)
  }
  dm
}

bind_covars <- function(dm, covars){
  if (!is.null(covars)){ # if there are covars, add these to design matrix
    dm <- cbind(dm, covars)
  }
}
