regcal <- function(mlist, B, alpha, continuous = T){ #regression calibration
  m <- mlist
  # make vector containing correct names -----------------------------------
  names <- colnames(m$x)
  names[2] <- colnames(m$ref)
  # estimate lambdas -------------------------------------------------------
  calfit <- stats::lm.fit(m$x[!is.na(m$ref),], m$ref[!is.na(m$ref)])
  lambda <- unname(calfit$coefficients)
  lambda[1:2] <- rev(lambda[1:2]) #reverse order
  vcov_lambda <- unname(mecor:::vcovfromfit(calfit))
  vcov_lambda <- mecor:::change_order_vcov(vcov_lambda) #reverse order
  # create covariate-error matrix Lambda -----------------------------------
  Lambda <- diag({ncoef <- NROW(lambda)})
  Lambda[1,] <- lambda
  # create covariate-error correction matrix -------------------------------
  A <- solve(Lambda)
  # estimate beta_star -----------------------------------------------------
  misfit <- stats::lm.fit(m$x, m$y)
  beta_star <- unname(misfit$coefficients)
  beta_star[1:2] <- rev(beta_star[1:2]) #reverse order
  vcov_beta_star <- unname(mecor:::vcovfromfit(misfit))
  vcov_beta_star <- mecor:::change_order_vcov(vcov_beta_star) #reverse order
  # estimate beta ----------------------------------------------------------
  beta <- beta_star%*%A
  # estimate cov(Lambda_rs, Lambda_tu) for all r, s, t, u = 1...ncoef ------
  # and put these element together in a matrix called vcov_Lambda ----------
  vcov_Lambda <- diag(ncoef*ncoef)
  diag(vcov_Lambda) <- 0
  vcov_Lambda[1:ncoef, 1:ncoef] <- vcov_lambda
  # S_kl is a ncoefxncoef matrix whose (i,j)th element is ------------------
  # cov(Lambda_ki, Lambda_lj) ----------------------------------------------
  S <- function(k, l) {
    vcov_Lambda[(k*ncoef-(ncoef-1)):(k*ncoef),
                (l*ncoef-(ncoef-1)):(l*ncoef)]}
  # estimate cov_A_i1j1i2j2 (formula A7) -----------------------------------
  cov_A <- function(i_1, j_1, i_2, j_2){
    out <- numeric(1)
    for(r in 1:ncoef){
      for(s in 1:ncoef){
        for(t in 1:ncoef){
          for(u in 1:ncoef){
            out <- out +
              A[i_1, r]*A[s, j_1]*A[i_2, t]*A[u, j_2]*S(r, t)[s,u]}}}}
    out
  }
  # vcov_A_mn is a ncoefxncoef matrix whose (i,j)th element is -------------
  # cov(A_im, A_jn) --------------------------------------------------------
  vcov_A <- function(m, n){
    out <- matrix(nrow = ncoef, ncol = ncoef)
    for(i_1 in 1:ncoef){
      for(i_2 in 1:ncoef){
        out[i_1, i_2] <- cov_A(i_1, m, i_2, n) } }
    out
  }
  # estimate zerovar vcov matrix of beta -----------------------------------
  zerovcov_beta <- t(A)%*%vcov_beta_star%*%A
  zerovcov_beta <- mecor:::change_order_vcov(zerovcov_beta)
  # estimate vcov_beta (formula A4) ----------------------------------------
  deltavcov_beta <- function(){
    out <- matrix(nrow = ncoef, ncol = ncoef)
    vcov_beta <- function(j_1, j_2) {
      zerovcov_beta[j_1, j_2] + t(beta_star)%*%vcov_A(j_1, j_2)%*%beta_star}
    for(j_1 in 1:ncoef){
      for(j_2 in 1:ncoef){
        out[j_1,j_2] <- vcov_beta(j_1, j_2) } }
    mecor:::change_order_vcov(out)
  }
  vcov_beta <- deltavcov_beta()
  # change names -----------------------------------------------------------
  rownames(zerovcov_beta) <- names
  colnames(zerovcov_beta) <- names
  rownames(vcov_beta) <- names
  colnames(vcov_beta) <- names
  beta[1:2] <- rev(beta[1:2]) #reverse order
  colnames(beta) <- names
  # estimate fieller ci for corrected covariate ---------------------------
  ci.fieller <- mecor:::fiellerci(misfit, calfit, alpha)
  # estimate bootstrap ci for all coefficients ----------------------------
  if(B != 0){
    bd <- data.frame(m$y, m$ref, m$x)
    colnames(bd) <- c("Y", colnames(m$ref), colnames(m$x))
    ci.b <- mecor:::boot_rc(data = bd, refname = colnames(m$ref), alpha, B)}
  out <- list(beta = beta,
              corvar = list(deltavar = {if(c) deltavar else NA},
                            bootvar = {if(B!=0) ci.b[,3] else NA}),
              ci = list(fiellerci = {if(c) ci.fieller else NA},
                        bootci = {if(B!=0) ci.b[,1:2] else NA}))
  out
}

regcal2 <- function(mlist, naivefit, B, alpha){ #regression calibration
  if(B==0){
    B = 999
    warning("No closed form variance, so bootstrap variance is calculated based on 999 bootstrap replicates")}
  m <- mlist
  coefs <- naivefit$coef[-1] #coefficients naivefit (excluding intercept)
  m2 <- cov(m$x[,-1, drop = F]) #m2 matrix
  v <- (m$test$test1 - m$x[,"repmean"])^2 + (m$test$test2 - m$x[,"repmean"])^2 #within individual variance
  varme <- mean(v)/2 #variance of measurement error term
  m1 <- m2
  m1[1,1] <- m1[1,1] - varme
  corfit <- t(solve(m1)%*%m2%*%coefs)[1,] #corrected coefficients
  ncorfit <- c("Corrected", names(corfit)[-1])
  names(corfit) <- ncorfit
  m3 <- apply(m$x[,-1, drop = F], 2, mean) #mean of vars
  int <- unname(mean(m$y) - m3%*%corfit) #corrected intercept
  corfit <- c("(Intercept)" = int, corfit)
  #deltavar <- mecor:::vardelta(naivefit, calfit, names(corfit$coef))
  #ci.fieller <- mecor:::fiellerci(naivefit, calfit, alpha)
  if(B != 0){
    bd <- data.frame(m$y, m$test, m$x)
    colnames(bd) <- c("Y", colnames(m$test), colnames(m$x))
    ci.b <- mecor:::boot_rc2(bdata = bd, testnames = colnames(m$test), alpha, B)}
  out <- list(corfit = list(coefficients = corfit),
              corvar = list(var = {if(B!=0) ci.b[,3] else NA}),
              ci = list(bootci = {if(B!=0) ci.b[,1:2] else NA}))
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
# x (a designmatrix containing an intercept, the test variable and other covariates) and
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

# this function creates a list containing y (a vector with the outcomes)
# x (a design matrix containing the intercept, the mean of the rep test variables and the other covariates)
# test (a matrix containing the rep measures)
rcm2 <- function(vars, me){
  y <- vars[,1] #vector containing the outcomes
  x <- cbind(1, {m <- (me$test$test1 + me$test$test2)/2}) #mean of the two replicate measures
  cnx <- c("(Intercept)", "repmean") #colnames x
  if(ncol(vars) > 1){ #if there are more covariates, add them to the design matrix
    x <- cbind(x, vars[,2:ncol(vars)]) #design matrix with measurement error
    cnx <- c(cnx, colnames(vars)[-1])}
  colnames(x) <- cnx
  out <- list(y = y, x = x, test = me$test)
}
