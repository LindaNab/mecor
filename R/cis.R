boot_rc <- function(data, refname, alpha, B){
  statrc <- function(data, refname, indices){
    d <- data[indices,]
    ref <- d[, refname]
    x <- as.matrix(d[,!names(data) %in% c(refname, "Y")])
    calfit <- stats::lm.fit(x[!is.na(ref),], ref[!is.na(ref)])
    corx <- x
    e <- calfit$coef%*%t(x)
    corx[,2] <- e
    colnames(corx)[2] <- refname
    corfit <- stats::lm.fit(corx, d[,"Y"])
    return(corfit$coef)
  }
  strata <- !is.na(data[, refname])
  meboot <- boot(data = data, statistic = statrc,
                 strata = strata, refname = refname, R = B)
  out <- matrix(nrow = {ids <- NROW(meboot$t0)}, ncol = 3,
                dimnames = list(names(meboot$t0), c("Lower", "Upper", "Var")))
  var <- diag(var(meboot$t))
  for(i in 1:ids){
    ci <- boot.ci(boot.out = meboot, conf = (1 - alpha), type = "perc", index = i)
    out[i,] <- c(ci$percent[,4], ci$percent[,5], var[i])
  }
  return(out)
}

boot_rc_pooled <- function(data, refname, naivefit, pooled.var, alpha, B){
  statrc_p <- function(data, refname, naivefit, pooled.var = "delta", indices){
    d <- data[indices,]
    ref <- as.matrix(d[, refname])
    colnames(ref) <- refname
    x <- as.matrix(d[,!names(d) %in% c(refname, "Y")])
    xint <- x
    xint[,2] <- ref
    colnames(xint)[2] <- refname
    m <- list(y = d$Y, x = x, ref = ref)
    resrc <- mecor:::regcal(m, naivefit, {if(pooled.var == "bootstrap") 999 else 0}, alpha)
    xint <- x
    xint[,2] <- ref
    colnames(xint)[2] <- colnames(m$ref)
    intfit <- stats::lm.fit(xint[!is.na(m$ref),], m$y[!is.na(m$ref)])
    intvar <- diag(mecor:::vcovfromfit(intfit))
    if(pooled.var == "bootstrap"){
      wrc <- 1/resrc$corvar$bootvar * (1/(1/resrc$corvar$bootvar + 1/intvar))}
    if(pooled.var == "delta"){
      wrc <- 1/resrc$corvar$deltavar * (1/(1/resrc$corvar$deltavar + 1/intvar))}
    pcorfit <- wrc * resrc$corfit$coef + (1 - wrc) * intfit$coef
    return(pcorfit)
  }
  strata <- !is.na(data[, refname])
  meboot <- boot(data = data, statistic = statrc_p,
                 strata = strata, refname = refname,
                 naivefit = naivefit, pooled.var = "delta", R = B)
  out <- matrix(nrow = {ids <- NROW(meboot$t0)}, ncol = 3)
  var <- diag(var(meboot$t))
  for(i in 1:ids){
    ci <- boot.ci(boot.out = meboot, conf = (1 - alpha), type = "perc", index = i)
    out[i,] <- c(ci$percent[,4], ci$percent[,5], var[i])
  }
  rownames(out) <- names(meboot$t0)
  colnames(out) <- c("Lower", "Upper", "Var")
  return(out)
}

fiellerci <- function(naivefit, calfit, alpha){
  out <- matrix(ncol = 2, nrow = {n<- NROW(naivefit$coef)},
                dimnames = list(names(naivefit$coef), c("Lower", "Upper")))
  ncov <- mecor:::vcovfromfit(naivefit)
  ccov <- mecor:::vcovfromfit(calfit)
  normq <- qnorm(alpha/2, lower.tail = FALSE)
  nb <- naivefit$coefficients[2]
  l1 <- calfit$coefficients[2]
  v1 <- normq^2 * ccov[2,2] - l1^2
  v2 <- l1 * nb
  v3 <- ncov[2,2] * normq^2 - nb^2
  if({D <- v2^2 - (v1 * v3)} > 0){
    l1 <- (- 1 * v2 - sqrt(D)) / v1
    l2 <- (- 1 * v2 + sqrt(D)) / v1
    ci <- c(lower = min(l1, l2),
            upper = max(l1, l2))}
  out[2,] <- ci
  out
}

# boot_rc <- function(data, rc_corf, rcf, alpha, B){
#   statrc <- function(data, rcf, rc_corf, indices){
#     d <- data[indices,]
#     rcbootm <- lm(rcf, d)
#     d$rc_reference <- predict(rcbootm, d)
#     return(coef(lm(formula = rc_corf, data = d)))
#   }
#   meboot <- boot(data = data, statistic = statrc, strata = !is.na(data$reference),
#                  rcf = rcf, rc_corf = rc_corf, R = B)
#   ids <- NROW(meboot$t0)
#   out <- matrix(nrow = ids, ncol = 2)
#   for(i in 1:ids){
#     ci <- boot.ci(boot.out = meboot, conf = (1 - alpha), type = "perc", index = i)
#     out[i,] <- c(ci$percent[,4], ci$percent[,5])
#   }
#   return(out)
# }

# var_rc <- function(data, nm, rcm, alpha){
#   nb <- coef(nm)['test']
#   l <- coef(rcm)['test']
#   varnb <- summary(nm)$coef['test',2]^2
#   varl<- summary(rcm)$coef['test',2]^2
#   cb <- nb / l
#   varcb <- varnb/(l^2) + (nb/l^2)^2*varl
#   tq <- qt((1 - alpha / 2), nm$df.residual)
#   out <- unname(c(cb - tq * sqrt(varcb), cb + tq * sqrt(varcb)))
# }
