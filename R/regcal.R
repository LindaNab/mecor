regcal <- function(mlist, naivefit, B, alpha){ #regression calibration
  m <- mlist
  calfit <- stats::lm.fit(m$x[!is.na(m$ref),], m$ref[!is.na(m$ref)]) #cal model
  corx <- m$x
  e <- calfit$coef%*%t(m$x) #predicted values
  corx[,2] <- e
  colnames(corx)[2] <- colnames(m$ref)
  corfit <- stats::lm.fit(corx, m$y)
  attributes(corfit)$type <- "lm.fit"
  deltavar <- mecor:::vardelta(naivefit, calfit, names(corfit$coef))
  ci.fieller <- mecor:::fiellerci(naivefit, calfit, alpha)
  if(B != 0){
    bd <- data.frame(m$y, m$ref, m$x)
    colnames(bd) <- c("Y", colnames(m$ref), colnames(m$x))
    ci.b <- mecor:::boot_rc(data = bd, refname = colnames(m$ref), alpha, B)}
  out <- list(corfit = corfit,
              corvar = list(deltavar = deltavar, bootvar = {if(B!=0) ci.b[,3] else NA}),
              ci = list(fiellerci = ci.fieller, bootci = {if(B!=0) ci.b[,1:2] else NA}))
  out
}

regcal_pooled <- function(mlist, naivefit, pooled.var = "delta", B, alpha){
  m <- mlist
  xint <- m$x
  xint[,2] <- m$ref
  colnames(xint)[2] <- colnames(m$ref)
  intfit <- stats::lm.fit(xint[!is.na(m$ref),], m$y[!is.na(m$ref)])
  intvar <- diag(mecor:::vcovfromfit(intfit))
  resrc <- mecor:::regcal(m, naivefit, {if(pooled.var == "bootstrap") B = 999 else 0}, alpha)
  if(pooled.var == "bootstrap"){
      wrc <- 1/resrc$corvar$bootvar * (1/(1/resrc$corvar$bootvar + 1/intvar))
      var <- 1 / ( (1 / resrc$corvar$bootvar) +
                     (1 / intvar) )}
  if(pooled.var == "delta"){
      wrc <- 1/resrc$corvar$deltavar * (1/(1/resrc$corvar$deltavar + 1/intvar))
      var <- 1 / ( (1 / resrc$corvar$deltavar) +
                     (1 / intvar) )}
  pcorfit <- wrc * resrc$corfit$coef + (1 - wrc) * intfit$coef
  if(B != 0){
    bd <- data.frame(m$y, m$ref, m$x)
    colnames(bd) <- c("Y", colnames(m$ref), colnames(m$x))
    ci.b <- mecor:::boot_rc_pooled(data = bd, refname = colnames(m$ref),
                                   naivefit = naivefit, pooled.var = "delta",
                                   alpha, B)
  }
  out <- list(corfit = list(coefficients = pcorfit),
              corvar = list(var = var, bootvar = {if(B!=0) ci.b[,3] else NA}),
              ci = list(bootci = {if(B!=0) ci.b[,1:2] else NA}))
}

# this function creates a list containing y (a vector with the outcomes),
# x (a designmatrix containing the covariates and an intercept) and
# ref (a matrix with the reference used in the calibration model)
rcm <- function(vars, me){
  y <- vars[,1] #vector containing the outcomes
  if({vtp <- attributes(me)$type} == "internal"){
    x <- cbind(1, me$test) #test var is the second entry of the x matrix
    cnx <- c("(Intercept)", attributes(me)$input$test) #colnames x
    ref <- as.matrix(me$reference) #reference measure as ref
    colnames(ref) <- as.character(attributes(me)$input$reference)}
  else if(vtp == "replicate"){
    x <- cbind(1, me$test$test1) #the first test measure is the second entry of the x matrix
    cnx <- c("(Intercept)", attributes(me)$input$test$test1) #colnames x
    ref <- as.matrix(me$test$test2) #replicate measure as ref
    colnames(ref) <- "Corrected"}
  if(ncol(vars) > 1){ #if there are more covariates, add them to the design matrix
    x <- cbind(x, vars[,2:ncol(vars)]) #design matrix with measurement error
    cnx <- c(cnx, colnames(vars)[-1])}
  colnames(x) <- cnx
  out <- list(y = y, x = x, ref = ref)
}
