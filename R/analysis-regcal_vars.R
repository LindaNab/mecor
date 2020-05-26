deltamethod <- function(func, vec, vcov_vec){
  j_func <- numDeriv::jacobian(func, vec)
  vcov <- j_func %*% vcov_vec %*% t(j_func)
}

# delta method
regcal_get_vcov <- function(beta_star, coef_calmod, vcov_beta_star, vcov_calmod, type){
  n <- NROW(beta_star)
  if (type == "dep"){
    beta_star <- c(beta_star, 1)
    vcov_temp <- matrix(0, nrow = (n + 1), ncol = (n + 1))
    vcov_temp[1:n, 1:n] <- vcov_beta_star
    vcov_beta_star <- vcov_temp
  }
  # take vec = (beta_star, c(A))
  # = (phi_star, alpha_star, .., gamma_k_star, 1/lambda_1, 0, 0, ..)
  # vec is of size (2+k)+(2+k)^2
  # k = number of covariates
  calmod_matrix <- mecor:::regcal_get_calmod_matrix(coef_calmod, n, type)
  vec_calmod_matrix <- c(calmod_matrix)
  vec_measerr_matrix <- mecor:::get_vec_measerr_matrix(vec_calmod_matrix)
  vec <- c(beta_star, vec_measerr_matrix)
  # The covariance matrix of vec is of size (2+k)+(2+k)^2 x (2+k)+(2+k)^2
  # We assume that there is no covariance between beta_star and vec_Lambda
  # Thus, vcov(vec) = (vcov(beta_star), 0
  #                   0,               vcov(vec_A)) # see step 2
  vcov_vec_measerr_matrix <-
    mecor:::get_vcov_vec_measerr_matrix(vec_calmod_matrix, vcov_calmod, n, type)
  vcov_vec <-
    mecor:::get_vcov_vec(vcov_beta_star, vcov_vec_measerr_matrix, type)
  # We assume that vec is multivariate normal with mean vec and cov vcov_vec
  # Now, there is a function f: R^{(2+k)+(2+k)^2} -> R^(2+k) so that
  # f(vec) = beta
  # Then, using the delta method vcov(beta) = Jf %*% vcov_vec %*% t(Jf)
  vcov_beta <- mecor:::deltamethod(mecor:::regcal_using_vec,
                                   vec,
                                   vcov_vec)
  dimnames(vcov_beta) <- list(names(beta_star), names(beta_star))
  # output
  vcov_beta[1:n,1:n]
}
regcal_using_vec <- function(vec){
  # nb = number of elements in beta_star
  nb <- floor(sqrt(NROW(vec)))
  beta <- numeric(nb)
  for(i in 1:nb){
    beta[i] <- sum(vec[1:nb] * vec[(i * nb + 1):((i + 1) * nb)])
  }
  beta
}
# create vec(A) is the vectorised calibration matrix
get_vec_measerr_matrix <- function(vec_calmod_matrix){
  n <- sqrt(NROW(vec_calmod_matrix))
  calmod_matrix <- matrix(vec_calmod_matrix, nrow = n)
  measerr_matrix <- solve(calmod_matrix)
  c(measerr_matrix)
}
# vcov matrix of vec = (beta, vec_A) is a block diagonal matrix
get_vcov_vec <- function(vcov_beta_star, vcov_vec_measerr_matrix, type){
  # block diagonal matrix
  vcov_vec <- rbind(
    cbind(vcov_beta_star,
          matrix(0, nrow = nrow(vcov_beta_star),
                 ncol = nrow(vcov_vec_measerr_matrix))),
    cbind(matrix(0, nrow = nrow(vcov_vec_measerr_matrix),
                 ncol = nrow(vcov_beta_star)),
          vcov_vec_measerr_matrix)
  )
}
# vcov matrix of vec_A = vec(A) using the Delta method
get_vcov_vec_measerr_matrix <- function(vec_calmod_matrix, vcov_calmod, n, type){
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
  # the vcov matrix of vec_Lambda is known:
  # it is a (2+k)^2 x (2+k)^2 matrix
  # vcov(vec_Lambda) = (var(lambda_1), cov(lambda_1, 0), ...,
  #                     cov(lambda_1, lambda_0), .. (row#1)
  #                     cov(0, lambda_1), var(0), ...,
  #                     cov(0, lambda_0), .. (row #2)
  #                     ..                                    )
  if (type == "dep"){
    vcov_vec_calmod_matrix <- matrix(0, nrow = (n + 1)^2,
                                     ncol = (n + 1)^2)
    for(i in 1:n){
      for(j in 1:n){
        vcov_vec_calmod_matrix[(n + 1) * (i - 1) + i, (n + 1) * (j - 1) + j] <-
          vcov_calmod[2, 2]
        vcov_vec_calmod_matrix[(n + 1) * (i - 1) + i, 2 * (n + 1)] <-
          vcov_calmod[2, 1]
        vcov_vec_calmod_matrix[2 * (n + 1), (n + 1) * (j - 1) + j] <-
          vcov_calmod[1, 2]
        vcov_vec_calmod_matrix[2 * (n + 1), 2 * (n + 1)] <-
          vcov_calmod[1, 1]
      }
    }
  } else if (type == "indep"){
    vcov_vec_calmod_matrix <- matrix(0, nrow = n^2,
                            ncol = n^2)
    for(i in 1:n){
      for(j in 1:n){
        vcov_vec_calmod_matrix[(1 + n * (i - 1)),
                        (1 + n * (j - 1))] <- vcov_calmod[i, j]
      }
    }
  }
  # Upon applying the delta method, we use that vec_Lambda is multivariate
  # normal with mean (lambda_1, 0, .., 0, lambda_0, 1, 0, .., 0, ....) and
  # variance vcov(vec_Lambda). There is a function g: R^{(2+k)x(2+k)} ->
  # R^{(2+k)x(2+k)} so that g(vec_Lambda) = vec_A.
  # Then, vcov(vec_A) = Jg %*% vcov(vec_Lambda) %*% t(Jg)
  # A is the calibration matrix, wich is the inverse of Lambda
  # get_vec_A vectorizes A, using the vectorised vec_Lambda
  # Thus, using the Delta method, vcov(vec_A) is:
  vcov_vec_measerr_matrix <-
    mecor:::deltamethod(mecor:::get_vec_measerr_matrix,
                        vec_calmod_matrix,
                        vcov_vec_calmod_matrix)
}
# bootstrap
regcal_boot <- function(response, covars, me, B, alpha, type){
  strat_samples <- replicate(
    B,
    mecor:::get_strat_sample(response, covars, me, type),
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
get_strat_sample <- function(response, covars, me, type){
  rownum_filled <- which (!is.na(me$reference))
  rownum_empty <- which (is.na(me$reference))
  new_rownum_filled <- sample(rownum_filled,
                              size = NROW(rownum_filled),
                              replace = T)
  new_rownum_empty <- sample(rownum_empty,
                             size = NROW(rownum_empty),
                             replace = T)
  new_rownums <- c(new_rownum_filled, new_rownum_empty)
  if (type == "indep"){
    new_response <- response[new_rownums, , drop = F]
  }
  new_covars <- covars[new_rownums, , drop = F]
  new_me <- me
  new_me$reference <- new_me$reference[new_rownums]
  new_me$substitute <- new_me$substitute[new_rownums]
  new_sample <- list(
    covars = new_covars,
    me = new_me,
    B = 0,
    type = type
  ) # no bootstrapping within bootstrap
  if (exists("new_response")){
    new_sample$response = new_response
    new_sample <- new_sample[c("response", "covars", "me", "B", "type")]
  }
  new_sample
}
