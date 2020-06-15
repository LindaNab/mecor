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
    name <- ifelse(!is.null(me$replicate),
                   paste0("cor_", attributes(me)$input$substitute),
                   attributes(me)$input$reference)
    colnames(dm) <- c("(Intercept)", name)
    dm <- mecor:::bind_covars(dm, covars)
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

get_dm_inadm <- function(covars, me){
  dm <- mecor:::get_dm_uncor(covars, me, type = "indep")
  calmod_fit <- mecor:::regcal_get_calmod(covars, me, type)
  fitted_values <- dm %*% coef(calmod_fit)
  x <- fitted_values
  x[!is.na(me$reference)] <- me$reference[!is.na(me$reference)]
  dm[, as.character(attributes(me)$input$substitute)] <- x
  ind <- colnames(dm) == as.character(attributes(me)$input$substitute)
  name <- as.character(attributes(me)$input$reference)
  colnames(dm)[ind] <- name
  dm
}

bind_covars <- function(dm, covars){
  if (!is.null(covars)){ # if there are covars, add these to design matrix
    dm <- cbind(dm, covars)
  }
}
