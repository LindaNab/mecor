deltamethod <- function(func, vec, vcov_vec){
  j_func <- numDeriv::jacobian(func, vec)
  vcov <- j_func %*% vcov_vec %*% t(j_func)
}

# delta method
regcal_get_vcov <- function(beta_star, coef_calmod, vcov_beta_star, vcov_calmod){
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
      beta[i] <- sum(vec[1:nb] * vec[(i * nb + 1):((i + 1) * nb)])
    }
    beta
  }
  jf <- numDeriv::jacobian(f, vec)
  vcov_beta <- jf %*% vcov_vec %*% t(jf)
  dimnames(vcov_beta) <- list(names(beta_star), names(beta_star))
  # output
  vcov_beta
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
  # Thus, using the Delta method, vcov(vec_A) is:
  vcov_vec_A <- deltamethod(get_vec_A, vec_Lambda, vcov_Lambda)
  vcov_vec_A
}
# bootstrap
regcal_boot <- function(response, covars, me, B, alpha){
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

# delta method outcome
regcal_get_vcov_o <- function(beta_star, theta, vcov_beta_star, vcov_theta){
  # take vec = (c(beta_star, 1), c(Theta))
  # = (phi_star, alpha_star, .., gamma_k_star, 1, 1/theta_1, 0, 0, ..)
  # vec is of size (3+k)+(3+k)^2
  # k = number of covariates
  Theta <- mecor:::get_Theta(theta, NROW(beta_star) + 1)
  vec_Theta <- c(Theta)
  vec_B <- mecor:::get_vec_B(vec_Theta)
  vec <- c(c(beta_star, 1), vec_B)
  # The covariance matrix of vec is of size (3+k)+(3+k)^2 x (3+k)+(3+k)^2
  # We assume that there is no covariance between beta_star and vec_Theta
  # Thus, vcov(vec) = (vcov(beta_star), 0
  #                   0,               vcov(vec_A)) # see step 2
  vcov_vec_B <- mecor:::get_vcov_vec_B(theta, vcov_theta)
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
      beta[i] <- sum(vec[1:nb] * vec[(i * nb + 1):((i + 1) * nb)])
    }
    beta
  }
  jf <- numDeriv::jacobian(f, vec)
  vcov_beta <- jf %*% vcov_vec %*% t(jf)
  dimnames(vcov_beta) <- list(names(beta_star), names(beta_star))
  # output
  vcov_beta
}
# create vec(A) is the vectorised calibration matrix
get_vec_B <- function(vec_Theta){
  Theta <- matrix(vec_Theta, nrow = sqrt(NROW(vec_Theta)))
  B <- solve(Theta)
  vec_B <- c(B)
  vec_B
}
# vcov matrix of vec = (beta, vec_A) is a block diagonal matrix
get_vcov_vec_o <- function(vcov_beta_star, vcov_vec_A){
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
get_vcov_vec_B <- function(lambda, vcov_lambda){
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
# bootstrap
regcal_boot_o <- function(response, covars, me, B, alpha){
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
get_strat_sample_o <- function(response, covars, me){
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
