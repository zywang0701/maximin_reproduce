### gen.data.with.target
### it generates data for simulations. The generated data includes target data
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


### gen.data.only source
### it generates data for simulations. The generated data only includes source data.
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

### util function to generate covariance
cov.gen <- function(p, rho){
  cov = matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      cov[i, j] = 0.6 ^ abs(i-j)
    }
  }
  return(cov)
}
