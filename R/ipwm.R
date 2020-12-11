#' Weighting for Confounding and Joint Misclassification of Exposure and Outcome
#'
#' \code{ipwm} implements a method for estimating the marginal causal odds ratio by constructing weights (modified inverse probability weights) that address both confounding and joint misclassification of exposure and outcome.
#'
#' @param formulas a list of objects of \link[base]{class} \code{\link[stats]{formula}} specifying the probability models for the stats::terms of some factorisation of the joint conditional probability function of \code{exposure_true}, \code{exposure_mis}, \code{outcome_true} and \code{outcome_mis}, given covariates
#' @param data \code{\link[base]{data.frame}} containing \code{exposure.true}, \code{exposure.mis}, \code{outcome.true}, \code{outcome.mis} and covariates. Missings (\code{NA}s) are allowed on variables \code{exposure_true} and \code{outcome_true}.
#' @param exposure_true a character string specifying the name of the true exposure variable that is free of misclassification but possibly unknown (\code{NA}) for some (but not all) subjects
#' @param exposure_mis a character string specifying the name of the counterpart of \code{exposure_true} that is available on all subjects but potentially misclassifies subjects as exposed or as non-exposed. The default (\code{exposure_mis = NULL}) indicates absence of exposure misclassification
#' @param outcome_true a character string specifying the name of the true outcome variable that is free of misclassification but possibly unknown (\code{NA}) for some (but not all) subjects
#' @param outcome_mis a character string specifying the name of the counterpart of \code{outcome_true} that is available on all subjects but potentially misclassifies subjects' outcomes. The default (\code{outcome_mis = NULL}) indicates absence of outcome misclassification
#' @param nboot number of bootstrap samples. Setting \code{nboot == 0} results in point estimation only.
#' @param conf_level the desired confidence level of the confidence interval
#' @param fix_nNAs logical indicator specifying whether or not to fix the joint distribution of \code{is.na(exposure_true)} and \code{is.na(outcome_true)}. If \code{TRUE}, stratified bootstrap sampling is done according to the missing data pattern.
#' @param semiparametric logical indicator specifying whether or not to parametrically sample \code{exposure_true}, \code{exposure_mis}, \code{outcome_true} and \code{outcome_mis}. If \code{semiparametric == TRUE}, it is assumed that the missing data pattern is conditionally independent of these variables given covariates. Provided \code{nboot > 0}, the missing data pattern and covariates are sampled nonparametrically. \code{semiparametric} is ignored if \code{nboot == 0}.
#' @param optim_args arguments passed onto \code{\link[stats]{optim}} if called. See Details below for more information.
#' @param force_optim logical indicator specifying whether or not to force the \code{\link[stats]{optim}} function to be called
#' @param sp scalar shrinkage parameter in the interval \code{(0, Inf)}. Values closer to zero result in greater shrinkage of the estimated odds ratio to unity; \code{sp == Inf} results in no shrinkage.
#' @param print logical indicator specifying whether or not to print the output.
#'
#' @details
#' This function is an implementation of the weighting method described by Penning de Vries et al. (2018).
#' The method defaults to the estimator proposed by Gravel and Platt (2018) in the absence of exposure misclassification.
#'
#' The function assumes that the exposure or the outcome has a misclassified version. An error is issued when both \code{outcome_mis} and \code{exposure_mis} are set to \code{NULL}.
#'
#' Provided \code{force_optim = FALSE}, \code{ipwm} is considerably more efficient when the \code{\link[stats]{optim}} function is not invoked; i.e., when (1) \code{exposure_mis = NULL} and the formula for \code{outcome_true} does not contain stats::terms involving \code{outcome_mis} or \code{exposure_true}, (2) \code{outcome_mis = NULL} and the formula for \code{exposure_true} does not contain stats::terms involving \code{exposure_mis} or \code{outcome_true}, or (3) \code{all(is.na(data[, exposure_true]) == is.na(data[, outcome_true]))} and the formulas for \code{exposure_true} and \code{outcome_true} do not contain stats::terms involving \code{exposure_mis} or \code{outcome_mis}. In these cases, \code{ipwm} uses iteratively reweighted least squares via the \code{\link[stats]{glm}} function for maximum likelihood estimation. In all other cases, \code{optim_args} is passed on to \code{\link[stats]{optim}} for optimisation of the joint likelihood of \code{outcome_true}, \code{outcome_mis}, \code{exposure_true} and \code{exposure_mis}.
#'
#' @return \code{ipwm} returns an object of \link[base]{class} \code{ipwm}.
#' The returned object is a list containing the following elements:
#'
#' \item{logOR}{the estimated log odds ratio;}
#' \item{call}{the matched function call.}
#'
#' If \code{nboot != 0}, the list also contains
#'
#' \item{SE}{a bootstrap estimate of the standard error for the estimator of the log odds ratio;}
#' \item{CI}{a bootstrap percentile confidence interval for the log odds ratio.}
#'
#' @author Bas B. L. Penning de Vries, \email{b.b.l.penning_de_vries@lumc.nl}
#'
#' @references
#' Gravel, C. A., & Platt, R. W. (2018). Weighted estimation for confounded binary outcomes subject to misclassification. \emph{Statistics in medicine}, 37(3), 425-436. https://doi.org/10.1002/sim.7522
#'
#' Penning de Vries, B. B. L., van Smeden, M., & Groenwold, R. H. H. (2020). A weighting method for simultaneous adjustment for confounding and joint exposure-outcome misclassifications. \emph{Statistical Methods in Medical Research}, 0(0), 1-15. https://doi.org/10.1177/0962280220960172
#'
#' @examples
#' data(sim) # simulated data on 10 covariates, exposure A and outcome Y.
#' formulas <- list(
#'   Y ~ A + L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + L10 + B + Z,
#'   A ~ L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + L10 + B + Z,
#'   Z ~ L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + L10 + B,
#'   B ~ L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + L10
#' )
#' \dontrun{
#' ipwm_out <- ipwm(
#'   formulas = formulas,
#'   data = sim,
#'   outcome_true = "Y",
#'   outcome_mis = "Z",
#'   exposure_true = "A",
#'   exposure_mis = "B",
#'   nboot = 200,
#'   sp = 1e6
#' )
#' ipwm_out
#' }
#'
#' @export
ipwm <- function(formulas,data,outcome_true,outcome_mis=NULL,exposure_true,exposure_mis=NULL,nboot=1000,conf_level=0.95,
                 fix_nNAs=FALSE,semiparametric=FALSE,optim_args=list(method="BFGS"),force_optim=FALSE,sp=Inf,print=TRUE
){
  n <- nrow(data)
  if(!is.numeric(sp) || length(sp)!=1L || sp<=0) stop("'sp' should be Inf or a positive number.")
  if(missing(outcome_true) || missing(exposure_true)) stop("'outcome_true' and 'exposure_true' should be specified.")
  nullZ <- is.null(outcome_mis)
  nullB <- is.null(exposure_mis)
  dep_var <- lapply(formulas,"[[",2L)
  if(any(lengths(dep_var)>1L)) stop("formulas should have untransformed binary variables as dependent variables.")
  dep_var <- as.character(dep_var)
  mt <- match(c(outcome_true,exposure_true,outcome_mis,exposure_mis),dep_var)
  if(any(is.na(mt))) stop(paste0("formula missing for variable",ifelse(sum(is.na(mt))>1,"s: ",": "),
                                 paste(c(outcome_true,exposure_true,outcome_mis,exposure_mis)[is.na(mt)],collapse=", "),"."))
  if(nullZ && nullB) stop("outcome_mis or exposure_mis should be non-NULL.")
  if(nullZ){
    if(any(is.na(data[,outcome_true]))) stop("outcome_true should have no missing values if outcome_mis==NULL.")
    Z <- data[,outcome_true]
    outcome_mis <- "Z"
    while(outcome_mis%in%colnames(data)) outcome_mis <- paste0(outcome_mis,".")
    data <- cbind(Z,data)
    colnames(data)[1L] <- outcome_mis
    formulas <- c(stats::as.formula(paste0(outcome_mis,"~-1+offset(ifelse(",outcome_true,",Inf,-Inf))")),formulas)
    dep_var <- lapply(formulas,"[[",2L)
    dep_var <- as.character(dep_var)
    mt <- match(c(outcome_true,exposure_true,outcome_mis,exposure_mis),dep_var)
  } else if(nullB){
    if(any(is.na(data[,exposure_true]))) stop("exposure_true should have no missing values if exposure_mis==NULL.")
    B <- data[,exposure_true]
    exposure_mis <- "B"
    while(exposure_mis%in%colnames(data)) exposure_mis <- paste0(exposure_mis,".")
    data <- cbind(B,data)
    colnames(data)[1L] <- exposure_mis
    formulas <- c(stats::as.formula(paste0(exposure_mis,"~-1+offset(ifelse(",exposure_true,",Inf,-Inf))")),formulas)
    dep_var <- lapply(formulas,"[[",2L)
    dep_var <- as.character(dep_var)
    mt <- match(c(outcome_true,exposure_true,outcome_mis,exposure_mis),dep_var)
  }
  formulas <- formulas[mt]
  dep_var <- dep_var[mt]
  if(any(is.na(data[,dep_var[3:4]]))) stop("outcome_mis and exposure_mis should not contain missing values.")
  ind_var <- lapply(formulas,function(x)all.vars(x[-1L]))
  wh_fo <- which(unlist(lapply(fact_order,
                               function(i)
                                 all(sapply(seq_len(4L),function(j)!dep_var[i[j]]%in%
                                              do.call(c,split(ind_var[i],factor(seq_len(4L)>j,c(TRUE,FALSE)))[[1L]])))
  )))
  if(!length(wh_fo)) stop("invalid factorisation of conditional joint model for (outcome_true, exposure_true, outcome_mis, exposure_mis) given covariates.")
  S <- 1L+is.na(data[,dep_var[1L]])*2L+is.na(data[,dep_var[2L]])
  tS <- table(factor(S,levels=1:4))
  opt <- force_optim
  if(!opt){
    if(nullB) w <- which(wh_fo%in%1:6)
    else if(nullZ) w <- which(wh_fo%in%7:12)
    else w <- which(wh_fo%in%c(1L,2L,7L,8L)[tS[c(2L,2L,3L,3L)]==0L])
    if(length(w)) fo <- fact_order[[wh_fo[w[1L]]]]
    else{
      fo <- fact_order[[wh_fo[1L]]]
      opt <- TRUE
    }
  } else fo <- fact_order[[wh_fo[1L]]]
  stat <- function(data){
    S <- 1L+is.na(data[,dep_var[1L]])*2L+is.na(data[,dep_var[2L]])
    par <- get_par(data,S)
    pB <- pbern(data[,dep_var[4L]],mean(data[,dep_var[4L]]))
    fn <- function(y,z,b) pr(data,par,y,data[,dep_var[4L]],z,b)
    v <- fn(1,0,0)+fn(1,0,1)+fn(1,1,0)+fn(1,1,1)
    w <- v+fn(0,0,0)+fn(0,0,1)+fn(0,1,0)+fn(0,1,1)
    num <- pB*v/w
    fn2 <- function(y,a) pr(data,par,y,a,1,data[,dep_var[4L]])
    denom <- fn2(0,0)+fn2(0,1)+fn2(1,0)+fn2(1,1)
    W <- num/denom
    p1 <- mean((data[,dep_var[3L]]*W)[data[,dep_var[4L]]==1])
    p0 <- mean((data[,dep_var[3L]]*W)[data[,dep_var[4L]]==0])
    if(p1>1) p1 <- 1
    if(p1<0) p1 <- 0
    if(p0>1) p0 <- 1
    if(p0<0) p0 <- 0
    p1 <- p1*(1+2/sp)^-1+(sp+2)^-1
    p0 <- p0*(1+2/sp)^-1+(sp+2)^-1
    logOR <- log(odds(p1))-log(odds(p0))
    return(logOR)
  }
  if(nboot){
    if(semiparametric) grd <- expand.grid(0:1,0:1,0:1,0:1)
    gen_boot <- function(data,par){
      if(fix_nNAs){
        wh <- split(seq_len(n),S)
        m <- lapply(wh,length)
        wh <- unlist(mapply(sample,wh,m,replace=TRUE))
      } else wh <- sample(seq_len(n),size=n,replace=TRUE)
      out <- data[wh,]
      isna <- is.na(out)
      if(all(rowSums(isna[,dep_var[1:2]])>0) && !fix_nNAs) return(gen_boot(data,par))
      if(semiparametric){
        fn3 <- function(x,y) pr(out[y,,drop=FALSE],par,grd[x,1L],grd[x,2L],grd[x,3L],grd[x,4L])
        p <- outer(seq_len(16L),seq_len(n),fn3)
        cp <- p
        for(i in 2:16) cp[i,] <- cp[i,]+cp[i-1L,]
        wh <- (stats::runif(n)<p[1L,])*1L
        for(i in 2:16){
          nwh <- !wh
          wh[nwh] <- (stats::runif(sum(nwh))<p[i,nwh]/(1-cp[i-1,nwh]))*i
        }
        out[,dep_var] <- grd[wh,]
        out[isna] <- NA
      }
      return(out)
    }
  }
  if(!opt){
    if(print) cat("Using stats::glm (IWLS) for maximum likelihood estimation ...\n",sep="")
    pr <- function(data,par,Y=data[,dep_var[1L]],A=data[,dep_var[2L]],Z=data[,dep_var[3L]],B=data[,dep_var[4L]]){
      if(!missing(Y)) data[,dep_var[1L]] <- Y
      if(!missing(A)) data[,dep_var[2L]] <- A
      if(!missing(Z)) data[,dep_var[3L]] <- Z
      if(!missing(B)) data[,dep_var[4L]] <- B
      p <- suppressWarnings(lapply(1:4,function(u)
        pbern(data[,dep_var[u]],stats::predict(par[[u]],newdata=data,type='response'))
      ))
      out <- p[[1L]]*p[[2L]]*p[[3L]]*p[[4L]]
      return(out)
    }
    get_par <- function(data,S) suppressWarnings(lapply(formulas,stats::glm,data=data,family=stats::binomial))
  } else{
    method <- optim_args$method
    if(print) cat("Using stats::optim (method = ",method,") for maximum likelihood estimation ...\n",sep="")
    pr <- function(data,par,Y=data[,dep_var[1L]],A=data[,dep_var[2L]],Z=data[,dep_var[3L]],B=data[,dep_var[4L]]){
      if(!missing(Y)) data[,dep_var[1L]] <- Y
      if(!missing(A)) data[,dep_var[2L]] <- A
      if(!missing(Z)) data[,dep_var[3L]] <- Z
      if(!missing(B)) data[,dep_var[4L]] <- B
      mm <- lapply(formulas,stats::model.matrix,data=data)
      tm <- lapply(formulas,stats::terms)
      offset <- lapply(tm,function(x) do.call(cbind,lapply(seq_len(length(attr(x,"offset"))),
                                                           function(i)eval(attr(x,"variables")[[attr(x,"offset")[[i]]+1L]],envir=data))))
      for(i in seq_len(length(offset))) if(is.null(offset[[i]])) offset[[i]] <- matrix(0,nrow=nrow(data))
      fn4 <- function(w,x,y,z) pbern(z,as.vector(expit(rowSums(w)+x%*%y)))
      out <- mapply(fn4,offset,mm,par,as.list(data[,dep_var]),SIMPLIFY=FALSE)
      out <- out[[1L]]*out[[2L]]*out[[3L]]*out[[4L]]
      return(out)
    }
    get_par <- function(data,S){
      npar <- integer(4L)
      for(i in 1:4){
        ft <- stats::terms.formula(formulas[[i]])
        npar[i] <- length(attr(ft,"term.labels"))+attr(ft,"intercept")
      }
      z <- by(data,S,list,simplify=FALSE)
      for(i in c('2','4')){
        fn5 <- function(x){
          for(j in 1:length(z[[i]])) z[[i]][[j]][,dep_var[2L]] <- x
          return(z[[i]])
        }
        if(i%in%names(z)) z[[i]] <- do.call(c,lapply(0:1,fn5))
      }
      for(i in c('3','4')){
        fn6 <- function(x){
          for(j in seq_len(length(z[[i]]))) z[[i]][[j]][,dep_var[1L]] <- x
          return(z[[i]])
        }
        if(i%in%names(z)) z[[i]] <- do.call(c,lapply(0:1,fn6))
      }
      z <- do.call(rbind,lapply(z,do.call,what=rbind))
      id <- unlist(mapply(rep,split(seq_len(n),factor(S,levels=seq_len(4L)),drop=FALSE),c(1L,2L,2L,4L)))
      v <- factor(unlist(mapply(rep,1:4,npar)),seq_len(4L))
      fn7 <- function(par=numeric(sum(npar))){
        par <- split(par,v)
        L <- pr(z,par)
        return(-sum(log(tapply(L,id,sum))))
      }
      optim_args$par <- numeric(sum(npar))
      optim_args$fn <- fn7
      optim_out <- do.call(stats::optim,optim_args)
      if(optim_out$converge) warning("Optimisation algorithm did not converge.")
      par <- split(optim_out$par,v)
      return(par)
    }
  }
  if(nboot){
    par <- get_par(data,S)
    est <- stat(data)
    boot_out <- numeric(nboot)
    for(it in seq_len(nboot)){
      bsample <- gen_boot(data,par)
      boot_out[it] <- stat(bsample)
      cat("\rBootstrapping: ",round(100*it/nboot),"% complete",sep="")
      utils::flush.console()
    }
    se <- stats::sd(boot_out)
    ci <- stats::quantile(boot_out,probs=(1-conf_level)/2+c(0,conf_level))
    out <- list(logOR=est,SE=se,CI=ci)
    if(print) cat("\n")
  } else out <- list(logOR=stat(data))
  out$call <- match.call()
  class(out) <- "ipwm"
  return(out)
}

#' @export
print.ipwm <- function(x, ...){
  object <- x
  cat("Call:\n",paste(deparse(object$call),sep="\n",collapse="\n"),"\n\n",sep="")
  OR <- formatC(exp(object$logOR),format="f",digits=3L)
  include_CI <- "CI"%in%names(object)
  if(include_CI){
    ci <- with(object,paste0(formatC(exp(CI),format="f",digits=3L)))
    CI <- paste0("(",paste0(ci,collapse=", "),")")
  } else CI <- ""
  if(include_CI){
    ci <- paste0("(",100-2*as.numeric(gsub("%$","",names(object$CI)[1L])),"% CI): ")
    cat("OR ",ci,paste(OR,CI),"\n",sep="")
  } else cat("OR: ",OR,"\n",sep="")
  invisible(object)
}

expit <- function(x) 1/(1+exp(-x))

pbern <- function(x,p) p^x*(1-p)^(1-x)

odds <- function(p) p/(1-p)

permutations <- function(n,k=n,print=TRUE){
  if(!{n >= k && k >= 0}){stop("condition (n >= k >= 0) not satisfied.")}
  if(k==1){return(as.list(1:n))}
  else{
    out <- list()
    length(out) <- m <- prod(n:(n-k+1))
    x <- 1:n
    f <- function(i){
      y <- integer(k)
      a <- m
      for(j in 1:k){
        b <- a/(n-j+1)
        d <- i%/%a-!i%%a
        p <- ceiling((i-d*a)/b)
        y[j] <- x[!x%in%y][p]
        a <- b
      }
      return(y)
    }
    if(print){
      for(i in 1:m){
        cat("\r",i,"/",m,sep=""); utils::flush.console()
        out[[i]] <- f(i)
      }; cat("\n");utils::flush.console()
    }
    else for(i in 1:m) out[[i]] <- f(i)
    return(out)
  }
}

fact_order <- permutations(4L,print=FALSE)
