deltamethod <- function(func, vec, vcov_vec){
  j_func <- numDeriv::jacobian(func, vec)
  vcov <- j_func %*% vcov_vec %*% t(j_func)
}

# delta method
regcal_get_vcov <- function(beta_star, coef_calmod, vcov_beta_star, vcov_calmod, type){
  n <- NROW(beta_star)
  if (startsWith(type, "dep")){
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
  } else if (type == "dep_diff"){
    vcov_vec_calmod_matrix <- matrix(0, nrow = (n + 1)^2,
                                     ncol = (n + 1)^2)
    reflectn_vcov_calmod <- mecor:::reflect_vcov_2nd_diagonal(vcov_calmod)
    vcov_vec_calmod_matrix[2:3, 2:3] <- reflectn_vcov_calmod[1:2, 1:2]
    vcov_vec_calmod_matrix[2:3, 5:6] <- reflectn_vcov_calmod[1:2, 3:4]
    vcov_vec_calmod_matrix[5:6, 2:3] <- reflectn_vcov_calmod[3:4, 1:2]
    vcov_vec_calmod_matrix[5:6, 5:6] <- reflectn_vcov_calmod[3:4, 3:4]
    vcov_vec_calmod_matrix[1, ] <-
      mecor:::get_1st_row_vcov_vec_calmod_matrix_diff_outme(vcov_calmod)
    vcov_vec_calmod_matrix[, 1] <- vcov_vec_calmod_matrix[1, ]
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
reflect_vcov_2nd_diagonal <- function(vcov){
  n <- nrow(vcov)
  for (i in 1:n){
    j <- 1
    while(j < (n + 1 - i)){
      temp <- vcov[i, j]
      vcov[i, j] <- vcov[n + 1 - j, n + 1 - i]
      vcov[n + 1 - j, n + 1 - i] <- temp
      j <- j + 1
    }
  }
  dimnames(vcov) <- list(rev(colnames(vcov)), rev(rownames(vcov)))
  vcov
}
get_1st_row_vcov_vec_calmod_matrix_diff_outme <- function(vcov_calmod){
  row <- numeric(9)
  row[1] <- vcov_calmod[4, 4] + vcov_calmod[2, 2] + 2 * vcov_calmod[4, 2]
  row[2] <- vcov_calmod[4, 4] + vcov_calmod[2, 4]
  row[3] <- vcov_calmod[4, 3] + vcov_calmod[2, 3]
  row[5] <- vcov_calmod[4, 2] + vcov_calmod[2, 2]
  row[6] <- vcov_calmod[4, 1] + vcov_calmod[2, 1]
  row
}
