bootstrap.rc <- function(data, rc_corf, rcf, alpha, B){
  statrc <- function(data, rcf, rc_corf, indices){
    main <- data
    valdata <- data[!is.na(data$reference),]
    if(count > 1){
      ids <- sample(1:nrow(valdata), replace = T)
      valsample <- valdata[ids,] #sample new validation dataset
      rcbootm <- lm(rcf, valsample)
      main$rc_reference <- predict(rcbootm, main)}
    count <<- count + 1
    d <- main[indices,]
    return(coef(lm(formula = rc_corf, data = d)))
  }
  count <- 1
  meboot <- boot(data = data, statistic = statrc,
                 rcf = rcf, rc_corf = rc_corf, R = B)
  ids <- NROW(meboot$t0)
  out <- matrix(nrow = ids, ncol = 2)
  for(i in 1:ids){
    ci <- boot.ci(boot.out = meboot, conf = (1 - alpha), type = "perc", index = i)
    out[i,] <- c(ci$percent[,4], ci$percent[,5])
  }
  return(out)
}

bootstrap2.rc <- function(data, rc_corf, rcf, alpha, B){
  statrc <- function(data, rcf, rc_corf, indices){
    d <- data[indices,]
    rcbootm <- lm(rcf, d)
    d$rc_reference <- predict(rcbootm, d)
    return(coef(lm(formula = rc_corf, data = d)))
  }
  meboot <- boot(data = data, statistic = statrc, strata = !is.na(data$reference),
                 rcf = rcf, rc_corf = rc_corf, R = B)
  ids <- NROW(meboot$t0)
  out <- matrix(nrow = ids, ncol = 2)
  for(i in 1:ids){
    ci <- boot.ci(boot.out = meboot, conf = (1 - alpha), type = "perc", index = i)
    out[i,] <- c(ci$percent[,4], ci$percent[,5])
  }
  return(out)
}


#ids <- sample(1:NROW(data[is.na(data$reference),]), replace = T)
#id <- sample(1:NROW(data[!is.na(data$reference),]), replace = T)
#d <- data[c(ids, id),]
