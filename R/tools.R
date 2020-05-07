## tools
get_coefs <- function(fit){
  coefs <- mecor:::change_order_coefs(fit$coefficients)
  coefs
}
get_vcov <- function(fit){
  vcov <- mecor:::vcovfromfit(fit)
  vcov <- mecor:::change_order_vcov(vcov)
  vcov
}
# order coefficients
change_order_coefs <- function(coefs){
  coefs[1:2] <- rev(coefs[1:2])
  names(coefs)[1:2] <- rev(names(coefs)[1:2])
  coefs
}
