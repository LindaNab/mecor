vardelta <- function(naivefit, calfit, names){
  vn <- mecor:::vcovfromfit(naivefit)
  vcal <- mecor:::vcovfromfit(calfit)
  cn <- naivefit$coef
  ccal <- calfit$coef
  varb <- 1/ccal[2]^2 * vn[2,2] + (cn[2]/ccal[2]^2)^2 * vcal[2,2]
  v <- function(i) {vn[i,i] - 2*ccal[i]/ccal[2]*vn[i,2] + (ccal[i]/ccal[2])^2*vn[2,2] +
    (cn[2]/ccal[2])^2*vcal[i,i] + (cn[2]*ccal[i]/ccal[2]^2)^2*vcal[2] -
    2*(cn[2]^2*ccal[i]^2/ccal[2]^3)*vcal[i,2]}
  if((n <- NROW(cn)) >= 3) out <- c(v(1), varb, {if(n > 3) diag(v(3:n)) else v(n)})
  else out <- c(v(1), varb)
  names(out) <- names
  out
}

vcovfromfit <- function(fit){
  p <- fit$rank
  p1 <- 1L:p
  R <- chol2inv(fit$qr$qr[p1, p1, drop = FALSE])
  r <- fit$residuals
  rss <- sum(r^2)
  rdf <- fit$df.residual
  resvar<- rss/rdf
  se <- sqrt(diag(R) * resvar)
  sigma <- sqrt(resvar)
  vcov <- sigma^2 * R
  dimnames(vcov) <- list(names(fit$coefficients), names(fit$coefficients))
  vcov
}
