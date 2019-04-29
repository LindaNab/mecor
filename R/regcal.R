regcal <- function(mlist, naivefit, B, alpha){ #regression calibration
  m <- mlist
  calfit <- stats::lm.fit(m$x[!is.na(m$ref),], m$ref[!is.na(m$ref)])
  corx <- m$x
  e <- calfit$coef%*%t(m$x) #predicted values
  corx[,2] <- e
  colnames(corx)[2] <- colnames(m$ref)
  corfit <- stats::lm.fit(corx, m$y)
  deltavar <- mecor:::vardelta(naivefit, calfit, names(corfit$coef))
  ci.fieller <- mecor:::fiellerci(naivefit, calfit, alpha)
  #if((pooled == F & B != 0) | (pooled == T & pooled.var == "bootstrap")){
  if(B != 0){
    bd <- data.frame(m$y, m$ref, m$x)
    colnames(bd) <- c("Y", colnames(m$ref), colnames(m$x))
    ci.b <- mecor:::boot_rc(data = bd, refname = colnames(m$ref), alpha, B)}
  #if(!pooled){
    out <- list(corfit = corfit, corvar = deltavar, ci.fieller = ci.fieller,
                ci.b = {if(B!=0) ci.b else NA})
  #}
  # if(pooled){
  #   intfit <- stats::lm.fit(m$xint[!is.na(m$ref),], m$y[!is.na(m$ref)])
  #   intvar <- diag(mecor:::vcovfromfit(intfit))
  #   if(pooled.var == "bootstrap"){
  #     wrc <- 1/ci.b[,"Var"] * (1/(1/ci.b[,"Var"] + 1/intvar))}
  #   if(pooled.var == "delta"){
  #     wrc <- 1/deltavar * (1/(1/deltavar + 1/intvar))}
  #   pcorfit <- wrc * corfit$coef + (1 - wrc) * intfit$coef
  #   if(B != 0){
  #     bd <- data.frame(m$y, m$ref, m$x)
  #   }
  #   out <- list(corfit = pcorfit)
  # }
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
      wrc <- 1/resrc$ci.b[,"Var"] * (1/(1/resrc$ci.b[,"Var"] + 1/intvar))
      var <- 1 / ( (1 / resrc$ci.b[,"Var"]) +
                     (1 / intvar) )}
  if(pooled.var == "delta"){
      wrc <- 1/resrc$corvar * (1/(1/resrc$corvar + 1/intvar))
      var <- 1 / ( (1 / resrc$corvar) +
                     (1 / intvar) )}
  pcorfit <- wrc * resrc$corfit$coef + (1 - wrc) * intfit$coef
  if(B != 0){
    bd <- data.frame(m$y, m$ref, m$x)
    colnames(bd) <- c("Y", colnames(m$ref), colnames(m$x))
    ci.b <- mecor:::boot_rc_pooled(data = bd, refname = colnames(m$ref),
                                   naivefit = naivefit, pooled.var = "delta",
                                   alpha, B)
  }
  var <- 1 / ( (1 / resrc$corvar) +
               (1 / intvar) )
  out <- list(corfit = pcorfit,
              var = var,
              ci.b = {if(B!=0) ci.b else NA})
}

rcm <- function(vars, me){ #this function designs the matrices needed for regcal
  y <- vars[,1] #vector containing the outcomes
  if({vtp <- attributes(me)$type} == "internal"){
    x <- cbind(1, me$test) #test var is the second entry of the x matrix
    cnx <- c("(Intercept)", attributes(me)$input$test)
    ref <- as.matrix(me$reference)
    colnames(ref) <- as.character(attributes(me)$input$reference)
    #xint <- cbind(1, me$reference)
    #if(ncol(vars) > 1) xint <- cbind(xint, vars[,2:ncol(vars)])
    }
  if(vtp == "replicate"){
    x <- cbind(1, me$test$test1)
    cnx <- c("(Intercept)", attributes(me)$input$test$test1)
    ref <- as.matrix(me$test$test2)
    colnames(ref) <- "Corrected"}
  if(ncol(vars) > 1){
    x <- cbind(x, vars[,2:ncol(vars)]) #design matrix with measurement error
    cnx <- c(cnx, colnames(vars)[-1])}
  colnames(x) <- cnx
  #if({b <- exists("xint")} == TRUE){
    #cnxint <- c("(Intercept)", attributes(me)$input$reference)
    #if(length(lv)!=0) cnxint <- c(cnxint, colnames(vars)[-1])
    #colnames(xint) <- cnxint}
  out <- list(y = y, x = x, ref = ref)
}
