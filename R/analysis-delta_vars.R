#' @importFrom numDeriv 'jacobian'
deltamethod <- function(func,
                        vec,
                        vcov_vec) {
  j_func <- numDeriv::jacobian(func, vec, method = "simple")
  vcov <- j_func %*% vcov_vec %*% t(j_func)
}
# delta method
standard_get_vcov <- function(beta_star,
                              coef_model,
                              vcov_beta_star,
                              vcov_model,
                              type) {
  n <- NROW(beta_star) # n = k + 2
  if (startsWith(type, "dep")) {
    beta_star <- c(beta_star, 1)
    vcov_temp <-
      matrix(0, nrow = (n + 1), ncol = (n + 1)) # expand with 0's
    vcov_temp[1:n, 1:n] <- vcov_beta_star
    vcov_beta_star <- vcov_temp
  }
  # k = number of covariates (Z) but the first one (covariate of interest: X)
  #
  # for covariate measurement error:
  # take vec = (beta_star, c(A))
  # (A = inverse of Lamba, the calibration model matrix)
  # = (phi_star, alpha_star, .., gamma_k_star, 1/lambda_1, 0, 0, ..)
  # vec is of size n + n^2 (cov-me)
  #
  # for outcome measurement error:
  # take vec = (beta_star, 1, c(B), 0, .., 1)
  # (B = inverse of Theta, the measurement error model matrix)
  # = (phi_star, alpha_star, .., gamma_k_star, 1/theta_1, 0, 0, ..)
  # vec is of size (n+1) + (n+1)^2 (out-me)
  #
  # for differential outcome measurement error in univariable analysis:
  # take vec = (beta_star, 1, c(B), 0, .., 1)
  # = (phi_star, alpha_star, 1, theta_11, (theta11 - theta00),
  #   (theta01 - theta00), 0, theta10, theta00, 0, 0, 1)
  # vec is of size 2+1 + (2+1)^2 [n+1 + (n+1)^2 where n = 2]
  # Lambda or Theta are referred to as the model_matrix in the code
  model_matrix <- standard_get_model_matrix(coef_model,
                                            n,
                                            type)
  vec_model_matrix <- c(model_matrix)
  # A or B is referred to as the solution matrix in the code
  vec_solution_matrix <-
    get_vec_solution_matrix(vec_model_matrix)
  vec <- c(beta_star, vec_solution_matrix)
  # The covariance matrix of vec is of size size(vec)^2 x size(vec)^2
  # We assume that there is no covariance between beta_star and the solution
  # matrix,
  # Thus, vcov(vec) = (vcov(beta_star), 0
  #                   0,               vcov(vec_solution_matrix)) # see step 2
  vcov_vec_solution_matrix <-
    get_vcov_vec_solution_matrix(vec_model_matrix,
                                 vcov_model,
                                 n,
                                 type)
  vcov_vec <-
    get_vcov_vec(vcov_beta_star,
                 vcov_vec_solution_matrix,
                 type)
  # We assume that vec is multivariate normal with mean vec and cov vcov_vec
  # Now, there is a function f: R^{n + n^2} -> R^n so that
  # f(vec) = beta
  # Then, using the delta method vcov(beta) = Jf %*% vcov_vec %*% t(Jf)
  vcov_beta <- deltamethod(standard_using_vec,
                           vec,
                           vcov_vec)
  dimnames(vcov_beta) <- list(names(beta_star), names(beta_star))
  # output
  vcov_beta[1:n, 1:n]
}
standard_using_vec <- function(vec) {
  # nb = number of elements in beta_star
  nb <- floor(sqrt(NROW(vec)))
  beta <- numeric(nb)
  for (i in 1:nb) {
    beta[i] <- sum(vec[1:nb] * vec[(i * nb + 1):((i + 1) * nb)])
  }
  beta
}
get_vec_solution_matrix <- function(vec_model_matrix) {
  n <- sqrt(NROW(vec_model_matrix))
  model_matrix <- matrix(vec_model_matrix, nrow = n)
  measerr_matrix <- solve(model_matrix)
  c(measerr_matrix)
}
# vcov matrix of vec = (beta, vec_A) is a block diagonal matrix
get_vcov_vec <- function(vcov_beta_star,
                         vcov_vec_solution_matrix,
                         type) {
  # block diagonal matrix
  vcov_vec <- rbind(cbind(vcov_beta_star,
                          matrix(
                            0,
                            nrow = nrow(vcov_beta_star),
                            ncol = nrow(vcov_vec_solution_matrix)
                          )),
                    cbind(matrix(
                            0,
                            nrow = nrow(vcov_vec_solution_matrix),
                            ncol = nrow(vcov_beta_star)
                         ),
                         vcov_vec_solution_matrix
                         ))
}
# vcov matrix of vec_A = vec(A) or vec_B = vec(B) using the Delta method
get_vcov_vec_solution_matrix <- function(vec_model_matrix,
                                         vcov_model,
                                         n,
                                         type) {
  # vcov(vec_A): cov matrix of the vectorised inverse of Lambda
  # Lambda = (lambda_1  lambda_0  lambda_2 ..
  #           0         1         0        ..
  #           ..        ..        ..       ..)
  # lambda = (lambda_1, lambda_0, lambda_2, ..., lambda_k)
  # vcov_lambda = (var(lambda_1), cov(lambda_1, lambda_0), ..
  #                cov(lambda_0, lambda_1), var(lambda_0), ..
  #                ..                                        )
  # vec_Lambda is the vectorised calibration model matrix
  # vec_Lambda = (lambda_1, 0, 0, 0, .., lambda_0, 1, 0, 0, 0, ....)
  # the vcov matrix of vec_Lambda is known:
  # it is a n^2 x n^2 matrix
  # vcov(vec_Lambda) = (var(lambda_1), cov(lambda_1, 0), ...,
  #                     cov(lambda_1, lambda_0), .. (row#1)
  #                     cov(0, lambda_1), var(0), ...,
  #                     cov(0, lambda_0), .. (row #2)
  #                     ..                                    )
  if (type == "indep") {
    vcov_vec_model_matrix <- matrix(0, nrow = n^2,
                                    ncol = n^2)
    for(i in 1:n){
      for(j in 1:n){
        vcov_vec_model_matrix[(1 + n * (i - 1)),
                              (1 + n * (j - 1))] <- vcov_model[i, j]
      }
    }
  # vcov(vec_B): cov matrix of the vectorised inverse of Theta
  # Theta = (theta_1   0         0        ..
  #          0         theta_1   0        ..
  #          ..        ..        ..       ..
  #          0         theta_0   0        1)
  # theta = (theta_0, theta_1)
  # vcov_theta = (var(theta_0)            cov(theta0, theta_1)
  #               cov(theta_1 , theta_0)  var(theta_1)         )
  # vec_Theta is the vectorised calibration model matrix
  # vec_Theta = (theta_1, 0, .., 0, 0, theta_1, .., theta_0, 0, .., 1)
  # the vcov matrix of vec_Theta is known:
  # it is a (n+1)^2 x (n+1)^2 matrix
  # vcov(vec_Theta) = (var(theta_1), cov(theta_1, 0), ...,
  #                     cov(theta_1, theta_1), .. (row #1, element n+3)
  #                     cov(0, theta_1), var(0), ...,
  #                     cov(0, theta_1), .. (row #2, element n+3)
  #                     ..                                    )
  } else if (type == "dep") {
    vcov_vec_model_matrix <- matrix(0, nrow = (n + 1)^2,
                                     ncol = (n + 1)^2)
    for(i in 1:n){
      for(j in 1:n){
        vcov_vec_model_matrix[(n + 1) * (i - 1) + i, (n + 1) * (j - 1) + j] <-
          vcov_model[2, 2]
        vcov_vec_model_matrix[(n + 1) * (i - 1) + i, 2 * (n + 1)] <-
          vcov_model[2, 1]
        vcov_vec_model_matrix[2 * (n + 1), (n + 1) * (j - 1) + j] <-
          vcov_model[1, 2]
      }
    }
  vcov_vec_model_matrix[2 * (n + 1), 2 * (n + 1)] <- vcov_model[1, 1]
  # vcov(vec_B): cov matrix of the vectorised inverse of Theta
  # Theta = (theta_11             0            0
  #          theta_11 - theta_10  theta_10     0
  #          theta_01 - theta_00  theta_00     1)
  # theta = (theta_00, theta_10, (theta_01-theta_00), (theta_11-theta10))
  # E[Y*|Y,X] = theta_00+theta_10Y+(theta_01-theta_00)X+(theta_11-theta10)XY
  #           = theta1+theta2Y+theta3X+theta4XY
  # vcov(theta) = (var(theta1) cov(theta1, theta2) cov(theta1, theta3) cov(theta1, theta4)
  #                cov(theta2, theta1) var(theta2) cov(theta2, theta3) cov(theta2, theta4)
  #                .... )
  # vec_Theta is the vectorised calibration model matrix
  # vec_Theta = (theta4 + theta2, theta4, theta3, 0, theta2, theta1, 0, 0, 1)
  # the vcov matrix of vec_Theta is known:
  # it is a 3^2 x 3^2 matrix
  # vcov(vec_Theta) = (var(theta4 + theta2) cov(theta4 + theta2, theta4) ....
  #                    cov(theta4, theta4 + theta2) var(theta4, theta4) ....
  #                    cov(theta3, theta4 + theta2) cov(theta3, theta4) ....
  #                    ....
  } else if (type == "dep_diff") {
    vcov_vec_model_matrix <- matrix(0, nrow = (n + 1)^2, # n = 2
                                     ncol = (n + 1)^2)
    # first, reflect the vcov_model in its 2nd diagonal:
    reflectn_vcov_model <- reflect_vcov_2nd_diagonal(vcov_model)
    # var(theta4)          cov(theta4, theta3)
    # vcov(theta3, theta4) var(theta3):
    vcov_vec_model_matrix[2:3, 2:3] <- reflectn_vcov_model[1:2, 1:2]
    # cov(theta4, theta2) cov(theta4, theta1)
    # cov(theta3, theta2) cov(theta3, theta1):
    vcov_vec_model_matrix[2:3, 5:6] <- reflectn_vcov_model[1:2, 3:4]
    # cov(theta2, theta4) cov(theta2, theta3)
    # cov(theta1, theta4) cov(theta1, theta3):
    vcov_vec_model_matrix[5:6, 2:3] <- reflectn_vcov_model[3:4, 1:2]
    # var(theta2)         cov(theta2, theta1)
    # cov(theta1, theta2) var(theta1):
    vcov_vec_model_matrix[5:6, 5:6] <- reflectn_vcov_model[3:4, 3:4]
    vcov_vec_model_matrix[1, ] <-
      get_1st_row_vcov_vec_model_matrix_diff_outme(vcov_model)
    vcov_vec_model_matrix[, 1] <- vcov_vec_model_matrix[1, ] # first column =
                                                             # first row
  }
  # Upon applying the delta method, we use that vec_Lambda is multivariate
  # normal with mean (lambda_1, 0, .., 0, lambda_0, 1, 0, .., 0, ....) and
  # variance vcov(vec_Lambda). There is a function g: R^{nxn} ->
  # R^{nxn} so that g(vec_Lambda) = vec_A.
  # Then, vcov(vec_A) = Jg %*% vcov(vec_Lambda) %*% t(Jg)
  # A is the calibration matrix, wich is the inverse of Lambda
  # get_vec_solution_matrix vectorizes A, using the vectorised vec_Lambda
  # Thus, using the Delta method, vcov(vec_A) is:
  # [idem dito vec_Theta --> vec_B]
  vcov_vec_solution_matrix <-
    deltamethod(get_vec_solution_matrix,
                vec_model_matrix,
                vcov_vec_model_matrix)
}
# the 2nd diagonal is the one from the upper right corner to the lower left
# corner. The elements on that diagonal will stay the same and the other
# elements of the matrix will flip.
reflect_vcov_2nd_diagonal <- function(vcov) {
  n <- nrow(vcov)
  for (i in 1:n) {
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
# vec_Theta = (theta4 + theta2, theta4, theta3, 0, theta2, theta1, 0, 0, 1)
# the vcov matrix of vec_Theta is known:
# it is a 3^2 x 3^2 matrix
# the first row is:
# (var(theta4 + theta2), cov(theta4 + theta2, theta4),
# cov(theta4 + theta2, theta3), 0,
# cov(theta4 + theta2, theta2), cov(theta4 + theta2, theta1),
# 0, 0, 0)
get_1st_row_vcov_vec_model_matrix_diff_outme <- function(vcov_model) {
  row <- numeric(9)
  # var(theta4 + theta2) = var(theta4) + var(theta2) + 2*cov(theta4, theta2)
  row[1] <- vcov_model[4, 4] + vcov_model[2, 2] + 2 * vcov_model[4, 2]
  # cov(theta4 + theta2, theta4) = var(theta4) + cov(theta2, theta4)
  row[2] <- vcov_model[4, 4] + vcov_model[2, 4]
  # cov(theta4 + theta2, theta3) = cov(theta4, theta3) + cov(theta2, theta3)
  row[3] <- vcov_model[4, 3] + vcov_model[2, 3]
  # cov(theta4 + theta2, theta2) = cov(theta4, theta2) + var(theta2)
  row[5] <- vcov_model[4, 2] + vcov_model[2, 2]
  # cov(theta4 + theta2, theta1) = cov(theta4, theta1) + cov(theta2, theta1)
  row[6] <- vcov_model[4, 1] + vcov_model[2, 1]
  row
}
# coefficients needed to calculate fieller based ci
standard_fieller <- function(beta_star,
                             coef_model,
                             vcov_beta_star,
                             vcov_model,
                             type) {
  if (type == "indep") {
    lambda1 <- coef_model[1]
    var_lambda1 <- vcov_model[1, 1]
    phi_star <- beta_star[1]
    var_phi_star <- vcov_beta_star[1, 1]
  } else if (type == "dep") {
    lambda1 <- coef_model[2] # theta1
    var_lambda1 <- vcov_model[2,2] # var(theta1)
    phi_star <- beta_star[-2] # (phi*, gamma*)
    var_phi_star <- diag(vcov_beta_star)[-2] # var(phi*), var(gamma*)
  }
  out <- list(lambda1 = unname(lambda1),
              var_lambda1 = unname(var_lambda1),
              phi_star = unname(phi_star),
              var_phi_star = unname(var_phi_star))
  out
}
# fieller method
# output is a matrix ci of length(beta_star) x 2, containing LCI and UCI
# in this order: phi*, gamma*
# gives lci and uci for phi* when there is measurement error in a covariate
# gives lci and uci for phi* and gamma* when there is measrurement error in the
# outcome
# used in mecor::summary.mecor
calc_fieller_ci <- function(lambda1,
                            var_lambda1,
                            phi_star,
                            var_phi_star,
                            alpha) {
  ci <- matrix(nrow = length(phi_star), ncol = 2)
  colnames(ci) <- c('LCI', 'UCI')
  z <- stats::qnorm(1 - alpha / 2)
  f0 <- z^2 * var_phi_star - phi_star^2
  f1 <- - phi_star * lambda1
  f2 <- z^2 * var_lambda1 - lambda1^2
  D <- f1 ^ 2 - f0 * f2
  if (length(phi_star) == 1 && (f2 < 0 & D > 0)) {
    l1 <- unname((f1 - sqrt(D)) / f2)
    l2 <- unname((f1 + sqrt(D)) / f2)
    ci[1, ] <- c(min(l1, l2), max(l1, l2))
  } else if (length(phi_star) > 1) {
    if(f2 < 0 & all(D > 0)){
      l1 <- unname((f1 - sqrt(D)) / f2)
      l2 <- unname((f1 + sqrt(D)) / f2)
      ci[1:length(phi_star), ] <- c(pmin(l1, l2), pmax(l1, l2))
    }
  }
  ci
}
# zerovariance ignores uncertainty in model_matrix
# output is vcov matrix of coefficients of beta_star (in that order)
# beta* = (phi*, alpha*, gamma*)
standard_zerovar <- function(vcov_beta_star,
                           model_matrix,
                           type) {
  n <- nrow(vcov_beta_star)
  dimnames_coef <- dimnames(vcov_beta_star)
  if (startsWith(type, "dep")) {
    vcov_temp <- matrix(0, nrow = (n + 1), ncol = (n + 1))
    vcov_temp[1:n, 1:n] <- vcov_beta_star
    vcov_beta_star <- vcov_temp
  }
  vcov <- t(solve(model_matrix)) %*%
            vcov_beta_star %*%
            solve(model_matrix) # sandwich method t(Lambda)*vcov_beta_star*Lambda
  vcov <- vcov[1:n, 1:n]
  dimnames(vcov) <- dimnames_coef
  vcov
}
