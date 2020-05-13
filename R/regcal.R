regcal <- function(response, covars, me, B = 0 , alpha = 0.05){
  # measurement error contaminated fit
  naive_fit <- mecor:::naive(response, covars, me)
  # estimate beta_star
  beta_star <- mecor:::get_coefs(naive_fit)
  vcov_beta_star <- mecor:::get_vcov(naive_fit)
  # estimate calibration model
  # get design matrix of calibration model
  dm_cm <- mecor:::get_dm_cm(covars, me) # design matrix calibration model
  dm_cm <- dm_cm[!is.na(me$reference), ] # select complete cases
  reference <- me$reference[!is.na(me$reference)]
  if (attributes(me)$type == "ivs"){
    name_reference <- as.character(attributes(me)$input$reference)
  } else {
    name_reference <- paste0("cor_", attributes(me)$input$substitute)
  }
  # fit calibration model
  calmod_fit <- stats::lm.fit(dm_cm, reference)
  lambda <- mecor:::get_coefs(calmod_fit)
  vcov_lambda <- mecor:::get_vcov(calmod_fit)
  # estimate beta and its vcov
  beta <- mecor:::regcal_coefs(beta_star, lambda)
  names(beta)[1] <- name_reference
  vcov_beta <- mecor:::regcal_vcov(beta_star, lambda,
                                    vcov_beta_star, vcov_lambda)
  colnames(vcov_beta) <- names(beta)
  rownames(vcov_beta) <- names(beta)
  out <- list(coef = beta,
              vcov = vcov_beta)
  if (B != 0){
    boot <- mecor:::boot_regcal(response, covars, me, B = B, alpha = alpha)
    out$boot <- list(ci = boot$ci,
                     vcov = boot$vcov)
  }
  out
}
# corrected coefs using regression calibration
regcal_coefs <- function(beta_star, lambda){
  Lambda <- get_Lambda(lambda)
  A <- solve(Lambda)
  beta <- as.numeric(beta_star %*% A)
  names(beta) <- names(beta_star)
  beta
}
# covariance matrix of the corrected beta
regcal_vcov <- function(beta_star, lambda, vcov_beta_star, vcov_lambda){
  # take vec = (beta_star, c(A))
  # = (phi_star, alpha_star, .., gamma_k_star, 1/lambda_1, 0, 0, ..)
  # vec is of size (2+k)+(2+k)^2
  # k = number of covariates
  Lambda <- mecor:::get_Lambda(lambda)
  vec_Lambda <- c(Lambda)
  vec_A <- mecor:::get_vec_A(vec_Lambda)
  vec <- c(beta_star, vec_A)
  # The covariance matrix of vec is of size (2+k)+(2+k)^2 x (2+k)+(2+k)^2
  # We assume that there is no covariance between beta_star and vec_Lambda
  # Thus, vcov(vec) = (vcov(beta_star), 0
  #                   0,               vcov(vec_A)) # see step 2
  vcov_vec_A <- mecor:::get_vcov_vec_A(lambda, vcov_lambda)
  vcov_vec <- mecor:::get_vcov_vec(vcov_beta_star, vcov_vec_A)
  # We assume that vec is multivariate normal with mean vec and cov vcov_vec
  # Now, there is a function f: R^{(2+k)+(2+k)^2} -> R^(2+k) so that
  # f(vec) = beta
  # Then, using the delta method vcov(beta) = Jf %*% vcov_vec %*% t(Jf)
  f <- function(vec){
    # nb = number of elements in beta_star
    nb <- floor(sqrt(NROW(vec)))
    beta <- numeric(nb)
    for(i in 1:nb){
      beta[i] <- sum(vec[1:nb]*vec[(i*nb+1):((i+1)*nb)])
    }
    beta
  }
  jf <- numDeriv::jacobian(f, vec)
  vcov_beta <- jf %*% vcov_vec %*% t(jf)
  colnames(vcov_beta) <- names(beta_star)
  rownames(vcov_beta) <- names(beta_star)
  vcov_beta
}
# create measurement error matrix
get_Lambda <- function(lambda){
  Lambda <- diag({ncoef <- NROW(lambda)})
  Lambda[1,] <- lambda
  Lambda
}
# create vec(A) is the vectorised calibration matrix
get_vec_A <- function(vec_Lambda){
  Lambda <- matrix(vec_Lambda, nrow = sqrt(NROW(vec_Lambda)))
  A <- solve(Lambda)
  vec_A <- c(A)
  vec_A
}
# vcov matrix of vec = (beta, vec_A) is a block diagonal matrix
get_vcov_vec <- function(vcov_beta_star, vcov_vec_A){
  # block diagonal matrix
  vcov_vec <- rbind(
    cbind(vcov_beta_star,
          matrix(0, nrow = nrow(vcov_beta_star), ncol = nrow(vcov_vec_A))),
    cbind(matrix(0, nrow = nrow(vcov_vec_A), ncol = nrow(vcov_beta_star)),
          vcov_vec_A)
  )
  vcov_vec
}
# vcov matrix of vec_A = vec(A) using the Delta method
get_vcov_vec_A <- function(lambda, vcov_lambda){
  # vcov(vec_A): cov matrix of the vectorised inverse of Lambda
  # Lambda = (lambda_1  lambda_0  lambda_2 ..
  #           0         1         0        ..
  #           ..        ..        ..       ..)
  # lambda = (lambda_1, lambda_0, lambda_2, ..., lambda_k)
  # vcov_lambda = (var(lambda_1), cov(lambda_1, lambda_0), ..
  #                cov(lambda_0, lambda_1), var(lambda_0), ..
  #                ..                                        )
  # vec_Lambda is the vectorised measurement error matrix
  # vec_Lambda = (lambda_1, 0, 0, 0, .., lambda_0, 1, 0, 0, 0, ....)
  Lambda <- get_Lambda(lambda)
  vec_Lambda <- c(Lambda)
  # the vcov matrix of vec_Lambda is known:
  # it is a (2+k)^2 x (2+k)^2 matrix
  # vcov(vec_Lambda) = (var(lambda_1), cov(lambda_1, 0), ...,
  #                     cov(lambda_1, lambda_0), .. (row#1)
  #                     cov(0, lambda_1), var(0), ...,
  #                     cov(0, lambda_0), .. (row #2)
  #                     ..                                    )
  vcov_vec_Lambda <- matrix(0, nrow = NROW(lambda)^2,
                            ncol = NROW(lambda)^2)
  for(i in 1:{n <- nrow(vcov_lambda)}){
    for(j in 1:n)
      vcov_vec_Lambda[(1 + n * (i - 1)),
                      (1 + n * (j - 1))] <- vcov_lambda[i, j]
  }
  # Upon applying the delta method, we use that vec_Lambda is multivariate
  # normal with mean (lambda_1, 0, .., 0, lambda_0, 1, 0, .., 0, ....) and
  # variance vcov(vec_Lambda). There is a function g: R^{(2+k)x(2+k)} ->
  # R^{(2+k)x(2+k)} so that g(vec_Lambda) = vec_A.
  # Then, vcov(vec_A) = Jg %*% vcov(vec_Lambda) %*% t(Jg)
  # A is the calibration matrix, wich is the inverse of Lambda
  # get_vec_A vectorizes A, using the vectorised vec_Lambda
  # Jacobian of the vectorised calibration matrix
  J_vec_A <- numDeriv::jacobian(get_vec_A, vec_Lambda)
  # Thus, using the Delta method, vcov(vec_A) is:
  vcov_vec_A <- J_vec_A %*% vcov_vec_Lambda %*% t(J_vec_A)
  vcov_vec_A
}
boot_regcal <- function(response, covars, me, B, alpha){
  strat_samples <- replicate(
    B,
    mecor:::get_strat_sample(response, covars, me),
    simplify = F
  )
  coef <- sapply(
    strat_samples,
    FUN = function(x) do.call(mecor:::regcal, x)$coef
  )
  ci_perc <- apply(coef,
                   1,
                   FUN = quantile,
                   probs = c(alpha / 2, 1 - alpha / 2))
  out <- list(ci = ci_perc,
              vcov = cov(t(coef)))
}
get_strat_sample <- function(response, covars, me){
  rownum_filled <- which (!is.na(me$reference))
  rownum_empty <- which (is.na(me$reference))
  new_rownum_filled <- sample(rownum_filled,
                              size = NROW(rownum_filled),
                              replace = T)
  new_rownum_empty <- sample(rownum_empty,
                             size = NROW(rownum_empty),
                             replace = T)
  new_rownums <- c(new_rownum_filled, new_rownum_empty)
  new_response <- response[new_rownums, , drop = F]
  new_covars <- covars[new_rownums, , drop = F]
  new_me <- me
  new_me$reference <- new_me$reference[new_rownums]
  new_me$substitute <- new_me$substitute[new_rownums]
  out <- list(response = new_response,
              covars = new_covars,
              me = new_me,
              B = 0) # no bootstrapping within bootstrap
  out
}




#
# regcal <- function(mlist, B, alpha, continuous = T){ #regression calibration
#   m <- mlist
#   # make vector containing correct names -----------------------------------
#   names <- colnames(m$x)
#   names[2] <- colnames(m$ref)
#   # estimate lambdas -------------------------------------------------------
#   calfit <- stats::lm.fit(m$x[!is.na(m$ref),], m$ref[!is.na(m$ref)])
#   lambda <- unname(calfit$coefficients)
#   lambda[1:2] <- rev(lambda[1:2]) #reverse order
#   vcov_lambda <- unname(mecor:::vcovfromfit(calfit))
#   vcov_lambda <- mecor:::change_order_vcov(vcov_lambda) #reverse order
#   # create covariate-error matrix Lambda -----------------------------------
#   Lambda <- diag({ncoef <- NROW(lambda)})
#   Lambda[1,] <- lambda
#   # create covariate-error correction matrix -------------------------------
#   A <- solve(Lambda)
#   # estimate beta_star -----------------------------------------------------
#   misfit <- stats::lm.fit(m$x, m$y)
#   beta_star <- unname(misfit$coefficients)
#   beta_star[1:2] <- rev(beta_star[1:2]) #reverse order
#   vcov_beta_star <- unname(mecor:::vcovfromfit(misfit))
#   vcov_beta_star <- mecor:::change_order_vcov(vcov_beta_star) #reverse order
#   # estimate beta ----------------------------------------------------------
#   beta <- beta_star%*%A
#   # estimate cov(Lambda_rs, Lambda_tu) for all r, s, t, u = 1...ncoef ------
#   # and put these element together in a matrix called vcov_Lambda ----------
#   vcov_Lambda <- diag(ncoef*ncoef)
#   diag(vcov_Lambda) <- 0
#   vcov_Lambda[1:ncoef, 1:ncoef] <- vcov_lambda
#   # S_kl is a ncoefxncoef matrix whose (i,j)th element is ------------------
#   # cov(Lambda_ki, Lambda_lj) ----------------------------------------------
#   S <- function(k, l) {
#     vcov_Lambda[(k*ncoef-(ncoef-1)):(k*ncoef),
#                 (l*ncoef-(ncoef-1)):(l*ncoef)]}
#   # estimate cov_A_i1j1i2j2 (formula A7) -----------------------------------
#   cov_A <- function(i_1, j_1, i_2, j_2){
#     out <- numeric(1)
#     for(r in 1:ncoef){
#       for(s in 1:ncoef){
#         for(t in 1:ncoef){
#           for(u in 1:ncoef){
#             out <- out +
#               A[i_1, r]*A[s, j_1]*A[i_2, t]*A[u, j_2]*S(r, t)[s,u]}}}}
#     out
#   }
#   # vcov_A_mn is a ncoefxncoef matrix whose (i,j)th element is -------------
#   # cov(A_im, A_jn) --------------------------------------------------------
#   vcov_A <- function(m, n){
#     out <- matrix(nrow = ncoef, ncol = ncoef)
#     for(i_1 in 1:ncoef){
#       for(i_2 in 1:ncoef){
#         out[i_1, i_2] <- cov_A(i_1, m, i_2, n) } }
#     out
#   }
#   # estimate zerovar vcov matrix of beta -----------------------------------
#   zerovcov_beta <- t(A)%*%vcov_beta_star%*%A
#   zerovcov_beta <- mecor:::change_order_vcov(zerovcov_beta)
#   # estimate vcov_beta (formula A4) ----------------------------------------
#   deltavcov_beta <- function(){
#     out <- matrix(nrow = ncoef, ncol = ncoef)
#     vcov_beta <- function(j_1, j_2) {
#       zerovcov_beta[j_1, j_2] + t(beta_star)%*%vcov_A(j_1, j_2)%*%beta_star}
#     for(j_1 in 1:ncoef){
#       for(j_2 in 1:ncoef){
#         out[j_1,j_2] <- vcov_beta(j_1, j_2) } }
#     mecor:::change_order_vcov(out)
#   }
#   vcov_beta <- deltavcov_beta()
#   # change names -----------------------------------------------------------
#   rownames(zerovcov_beta) <- names
#   colnames(zerovcov_beta) <- names
#   rownames(vcov_beta) <- names
#   colnames(vcov_beta) <- names
#   beta[1:2] <- rev(beta[1:2]) #reverse order
#   colnames(beta) <- names
#   # estimate fieller ci for corrected covariate ---------------------------
#   ci.fieller <- mecor:::fiellerci(misfit, calfit, alpha)
#   # estimate bootstrap ci for all coefficients ----------------------------
#   if(B != 0){
#     bd <- data.frame(m$y, m$ref, m$x)
#     colnames(bd) <- c("Y", colnames(m$ref), colnames(m$x))
#     ci.b <- mecor:::boot_rc(data = bd, refname = colnames(m$ref), alpha, B)}
#   out <- list(beta = beta,
#               corvar = list(deltavar = {if(c) deltavar else NA},
#                             bootvar = {if(B!=0) ci.b[,3] else NA}),
#               ci = list(fiellerci = {if(c) ci.fieller else NA},
#                         bootci = {if(B!=0) ci.b[,1:2] else NA}))
#   out
# }

regcal2 <- function(mlist, naivefit, B, alpha){ #regression calibration
  if(B==0){
    B = 999
    warning("No closed form variance, so bootstrap variance is calculated based on 999 bootstrap replicates")}
  m <- mlist
  coefs <- naivefit$coef[-1] #coefficients naivefit (excluding intercept)
  m2 <- cov(m$x[,-1, drop = F]) #m2 matrix
  v <- (m$substitute$substitute1 - m$x[,"repmean"])^2 + (m$substitute$substitute2 - m$x[,"repmean"])^2 #within individual variance
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
    bd <- data.frame(m$y, m$substitute, m$x)
    colnames(bd) <- c("Y", colnames(m$substitute), colnames(m$x))
    ci.b <- mecor:::boot_rc2(bdata = bd, substitutenames = colnames(m$substitute), alpha, B)}
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
# x (a designmatrix containing an intercept, the substitute variable and
# the other covariates) and
# ref (a matrix with the reference used in the calibration model)
rcm <- function(vars, me){
  y <- vars[, 1] #vector containing the outcomes
  if ({vtp <- attributes(me)$type} == "internal"){
    x <- cbind(1, me$substitute)
    colnames_x <- c("(Intercept)", attributes(me)$input$substitute) #colnames x
    ref <- as.matrix(me$reference) #reference measure as ref
    colnames(ref) <- as.character(attributes(me)$input$reference)
    }
  else if (vtp == "replicate"){
    x <- cbind(1, me$substitute$substitute1)
    colnames_x <- c("(Intercept)", attributes(me)$input$substitute$substitute1)
    ref <- as.matrix(me$substitute$substitute2) #replicate measure as ref
    colnames(ref) <- "Corrected"}
  if (ncol(vars) > 1){ #if there are more covariates, add them to the design matrix
    x <- cbind(x, vars[, 2:ncol(vars)]) #design matrix with measurement error
    colnames_x <- c(colnames_x, colnames(vars)[-1])
  }
  colnames(x) <- colnames_x
  out <- list(y = y, x = x, ref = ref)
}

# this function creates a list containing y (a vector with the outcomes)
# x (a design matrix containing the intercept, the mean of the rep substitute variables and the other covariates)
# substitute (a matrix containing the rep measures)
rcm2 <- function(vars, me){
  y <- vars[,1] #vector containing the outcomes
  x <- cbind(1, {m <- (me$substitute$substitute1 + me$substitute$substitute2)/2}) #mean of the two replicate measures
  cnx <- c("(Intercept)", "repmean") #colnames x
  if(ncol(vars) > 1){ #if there are more covariates, add them to the design matrix
    x <- cbind(x, vars[,2:ncol(vars)]) #design matrix with measurement error
    cnx <- c(cnx, colnames(vars)[-1])}
  colnames(x) <- cnx
  out <- list(y = y, x = x, substitute = me$substitute)
}
