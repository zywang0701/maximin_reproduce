library(MASS)
report.true <- function(Bs, cov.target, loading=NULL, delta.truth=-1, verbose=TRUE){
  Gamma = t(Bs) %*% cov.target %*% Bs
  if(delta.truth==-1){
    delta = decide_delta(Gamma, step_delta=0.1, MAX_iter=100, verbose=verbose)
  }else{
    delta = delta.truth
  }
  sol = opt.weight(Gamma, delta)
  weight.true = sol$weight
  beta.true = Bs %*% weight.true
  if(is.null(loading)){
    returnList <- list("delta"=delta,
                       "weight.true"=weight.true,
                       "beta.true"=beta.true)
  }else{
    true.val = sum(beta.true * loading)
    returnList <- list("delta"=delta,
                       "weight.true"=weight.true,
                       "beta.true"=beta.true,
                       "true.val"=true.val)
  }
  return(returnList)
}

report.true.reward <- function(Bs, cov.target, loading=NULL, delta.truth=-1, verbose=TRUE){
  Gamma = t(Bs) %*% cov.target %*% Bs
  if(delta.truth==-1){
    delta = decide_delta(Gamma, step_delta=0.1, MAX_iter=100, verbose=verbose)
  }else{
    delta = delta.truth
  }
  sol = opt.weight(Gamma, delta)
  weight.true = sol$weight
  beta.true = Bs %*% weight.true
  if(is.null(loading)){
    returnList <- list("delta"=delta,
                       "weight.true"=weight.true,
                       "beta.true"=beta.true,
                       "reward"=sol$reward)
  }else{
    true.val = sum(beta.true * loading)
    returnList <- list("delta"=delta,
                       "weight.true"=weight.true,
                       "beta.true"=beta.true,
                       "true.val"=true.val,
                       "reward"=sol$reward)
  }
  return(returnList)
}

########################################################################
#' For later simulations, we firstly create two functions for generating data
#'
#' @param mean.source 
#' @param cov.source 
#' @param cov.target 
#' @param Bs true parameters of groups, matrix dim of (p, L)
#' @param Ns length of each group, vector, length of L
#' @param L scalar
#' @param N.target scalar
gen.data.with.target <- function(mean.source, cov.source, mean.target, cov.target, N.target, Bs, Ns, L){
  N.source = sum(Ns)
  index.source = seq(N.source)
  X.source = mvrnorm(N.source, mu=mean.source, Sigma=cov.source)
  X.target = mvrnorm(N.target, mu=mean.target, Sigma=cov.target)
  # idx.source
  idx.source = rep(1:L, times=Ns)
  # Y.source
  Y.source = rep(0, N.source)
  for(l in 1:L){
    index.l = which(idx.source==l)
    Y.source[index.l] = X.source[index.l, ] %*% Bs[, l] + rnorm(Ns[l])
  }
  returnList <- list("X.source"=X.source,
                     "Y.source"=Y.source,
                     "X.target"=X.target,
                     "idx.source"=idx.source)
  return(returnList)
}

gen.data.only.source <- function(mean.source, cov.source, Bs, Ns, L, sd=1){
  N.source = sum(Ns)
  index.source = seq(N.source)
  X.source = mvrnorm(N.source, mu=mean.source, Sigma=cov.source)
  # idx.source
  idx.source = rep(1:L, times=Ns)
  # Y.source
  Y.source = rep(0, N.source)
  for(l in 1:L){
    index.l = which(idx.source==l)
    Y.source[index.l] = X.source[index.l, ] %*% Bs[, l] + rnorm(Ns[l],sd=sd)
  }
  returnList <- list("X.source"=X.source,
                  "Y.source"=Y.source,
                  "idx.source"=idx.source)
}

cov.gen <- function(p, rho){
  cov = matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      cov[i, j] = 0.6 ^ abs(i-j)
    }
  }
  return(cov)
}

cov.coef.loading.settings <- function(setting){
  cov.source = cov.gen(p, 0.6)
  cov.target = cov.source
  Bs = matrix(0, nrow=p, ncol=L)
  loading = rep(0, p)
  if(setting==1){
    # cov.target
    for(i in 1:p) cov.target[i, i] = 1.5
    for(i in 1:5){
      for(j in 1:5){
        if(i!=j) cov.target[i, j] = 0.9
      }
    }
    for(i in (p-1):p){
      for(j in (p-1):p){
        if(i!=j) cov.target[i, j] = 0.9
      }
    }
    # Bs
    for(j in 1:10) Bs[j, ] = j/40
    Bs[22:23, ] = 1
    Bs[p-2, 1] = 0.5
    Bs[(p-1):p, 1] = -0.5
    for(l in 2:L) Bs[p, l] = 1.4-0.2*l
    # loading
    loading[(p-2):p] = 1
  }
  if(setting==2){
    # cov.target
    for(i in 1:p) cov.target[i, i] = 1.1
    for(i in 1:6){
      for(j in 1:6){
        if(i!=j) cov.target[i, j] = 0.75
      }
    }
    # Bs
    for(j in 1:10) Bs[j, 1] = j/10
    for(j in 11:20) Bs[j, 1] = (10-j)/10
    Bs[21, 1] = 1/5
    Bs[22:23, 1] = 1
    for(l in 2:L){
      for(j in 1:10) Bs[j, l] = Bs[j, 1] + 0.1*(l-1)/sqrt(300)
      for(j in 11:20) Bs[j, l] = -0.3*(l-1)/sqrt(300)
      Bs[21, l] = 0.5*(l-1)
      for(j in 22:23) Bs[j, l] = 0.2*(j-1)
    }
    # loading
    loading[21:23] = 1
  }
  if(setting==3){
    set.seed(2021)
    # cov.target
    for(i in 1:p) cov.target[i, i] = 1.1
    for(i in 1:6){
      for(j in 1:6){
        if(i!=j) cov.target[i, j] = 0.75
      }
    }
    # Bs
    for(j in 1:10) Bs[j, 1] = j/10
    for(j in 11:20) Bs[j, 1] = (10-j)/10
    Bs[21, 1] = 1/5
    Bs[22:23, 1] = 1
    for(l in 2:2){
      for(j in 1:10) Bs[j, l] = Bs[j, 1] + 0.1*(l-1)/sqrt(300)
      for(j in 11:20) Bs[j, l] = -0.3*(l-1)/sqrt(300)
      Bs[21, l] = 0.5*(l-1)
      for(j in 22:23) Bs[j, l] = 0.2*(j-1)
    }
    for(l in 3:L) Bs[1:6, l] = rnorm(6)
    # loading
    Sigma.new = matrix(0, nrow=p, ncol=p)
    for(i in 1:p){
      for(j in 1:p){
        Sigma.new[i, j] = 0.5 ^ (1+abs(i-j)) / 25
      }
    }
    loading = mvrnorm(1, mu=rep(0, p), Sigma=Sigma.new)
  }
  if(setting==0){
    ## irregularity setting
    if(L!=2) stop("Setting 4 works only when L=2")
    # cov.target
    for(i in 1:p) cov.target[i, i] = 1.5
    for(i in 1:5){
      for(j in 1:5){
        if(i!=j) cov.target[i, j] = 0.9
      }
    }
    for(i in (p-1):p){
      for(j in (p-1):p){
        if(i!=j) cov.target[i, j] = 0.9
      }
    }
    # Bs
    for(j in 1:10) Bs[j, 1] = j/40
    for(j in 11:20) Bs[j, 1] = (10-j)/40
    Bs[21, 1] = 0.2
    Bs[22:23, 1] = 1
    for(j in 1:10) Bs[j, 2] = Bs[j, 1] + perb/sqrt(300)
    Bs[21, 2] = 0.5
    Bs[22:23, 2] = 0.2
    # loading
    for(j in 1:5) loading[j] = j/5
  }
  returnList <- list("cov.source"=cov.source,
                     "cov.target"=cov.target,
                     "Bs"=Bs,
                     "loading"=loading)
  return(returnList)
}
