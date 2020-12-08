## tools
get_coefs <- function(fit,
                      rev = TRUE) {
  coefs <- fit$coef
  if (rev == TRUE) {
    # reverse coefficients
    coefs <- change_order_coefs(fit$coef)
  }
  coefs
}
get_vcov <- function(fit,
                     rev = TRUE) {
  vcov <- vcovfromfit(fit)
  if (rev == TRUE) {
    # reverse coefficients
    vcov <- change_order_vcov(vcov)
  }
  vcov
}
get_vcovHC3 <- function(fit,
                        dm,
                        rev = TRUE) {
  vcov <- vcovHC3fromfit(fit, dm)
  if (rev == TRUE) {
    # reverse coefficients
    vcov <- change_order_vcov(vcov)
  }
  vcov
}
vcovfromfit <- function(fit) {
  p <- fit$rank
  p1 <- 1L:p
  R <- chol2inv(fit$qr$qr[p1, p1, drop = FALSE])
  r <- fit$residuals
  rss <- sum(r ^ 2)
  rdf <- fit$df.residual
  resvar <- rss / rdf
  se <- sqrt(diag(R) * resvar)
  sigma <- sqrt(resvar)
  vcov <- sigma ^ 2 * R
  dimnames(vcov) <-
    list(names(fit$coefficients), names(fit$coefficients))
  vcov
}
vcovHC3fromfit <- function(fit,
                           dm) {
  n <- NROW(dm)
  hat <- dm %*% solve(t(dm) %*% dm) %*% t(dm)
  diaghat <- diag(hat)
  df <- fit$df.residual
  res <- stats::residuals(fit)
  omega <- res ^ 2 / (1 - diaghat) ^ 2
  rval <- sqrt(omega) * dm
  rval <- crossprod(rval) / n
  meat <- rval
  p <- fit$rank
  p1 <- 1L:p
  R <- chol2inv(fit$qr$qr[p1, p1, drop = FALSE])
  bread <- R * n
  vcovHC3 <- (1 / n) * (bread %*% meat %*% bread)
  vcovHC3
}
# order coefficients
change_order_coefs <- function(coefs) {
  coefs[1:2] <- rev(coefs[1:2])
  names(coefs)[1:2] <- rev(names(coefs)[1:2])
  coefs
}
# order vcov
change_order_vcov <- function(vcov) {
  #save upper 2x2 matrix
  temp <- vcov[1:2, 1:2]
  #rearrange diagonal
  diag(temp) <- c(temp[2, 2], temp[1, 1])
  #swap first row and second row
  vcov[1:2, ] <- rbind(vcov[2, ], vcov[1, ])
  #swap first column and second column
  vcov[, 1:2] <- cbind(vcov[, 2], vcov[, 1])
  #replace upper 2x2 matrix
  vcov[1:2, 1:2] <- temp
  colnames(vcov)[1:2] <- rev(colnames(vcov)[1:2])
  rownames(vcov)[1:2] <- rev(rownames(vcov)[1:2])
  vcov
}
