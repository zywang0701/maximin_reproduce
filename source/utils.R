# Compute Ridge-type weight vector
opt.weight<-function(Gamma,delta,report.reward=TRUE){
  L<-dim(Gamma)[2]
  opt.weight<-rep(NA, L)
  opt.reward<-NA
  # Problem definition
  v<-Variable(L)
  Diag.matrix<-diag(eigen(Gamma)$values)
  for(ind in 1:L){
    Diag.matrix[ind,ind]<-max(Diag.matrix[ind,ind],0.001)
  }
  Gamma.positive<-eigen(Gamma)$vectors%*%Diag.matrix%*%t(eigen(Gamma)$vectors)
  objective <- Minimize(quad_form(v,Gamma.positive+diag(delta,L)))
  constraints <- list(v >= 0, sum(v)== 1)
  prob.weight<- Problem(objective, constraints)
  if(is_dcp(prob.weight)){
    result<- solve(prob.weight)
    opt.status<-result$status
    opt.sol<-result$getValue(v)
    for(l in 1:L){
      opt.weight[l]<-opt.sol[l]*(abs(opt.sol[l])>10^{-8})
    }
  }
  if(report.reward){
    v<-Variable(L)
    objective<-Minimize(2*t(v)%*%Gamma.positive%*%opt.weight-t(opt.weight)%*%Gamma.positive%*%opt.weight)
    constraints<-list(v >= 0, sum(v)== 1)
    delta.optim<-Problem(objective, constraints)
    result<- solve(delta.optim)
    opt.reward<-result$value
    returnList <- list("weight" = opt.weight,
                       "reward" = opt.reward)
  }else{
    returnList <- list("weight" = opt.weight)
  }
  return(returnList)
}

# Compute Ridge-type weight vector without reward
opt.weight.2<-function(Gamma,delta){
  opt.weight<-rep(NA, 2)
  temp<-(Gamma[2,2]+delta-Gamma[1,2])/(Gamma[1,1]+Gamma[2,2]+2*delta-2*Gamma[1,2])
  opt.weight[1]<-min(max(temp,0),1)
  opt.weight[2]<-1-opt.weight[1]
  # opt.reward<-t(opt.weight)%*%Gamma%*%opt.weight
  # returnList <- list("weight" = opt.weight,
  #                    "reward" = opt.reward)
  returnList <- list("weight"=opt.weight)
  return(returnList)
}

#' Bias correction for initial estimator of Gamma target
Gamma.shift<-function(plug.in,X.l,X.k,omega.l,omega.k,Y.l,Y.k,Pred.l,Pred.k){
  u.lk<-proj.direction(X.l,omega.k)
  u.kl<-proj.direction(X.k,omega.l)
  n.k<-nrow(X.k)
  n.l<-nrow(X.l)
  prop.est<-plug.in+t(u.kl)%*%t(X.k)%*%(Y.k-Pred.k)/n.k+t(u.lk)%*%t(X.l)%*%(Y.l-Pred.l)/n.l
  returnList <- list("est" = prop.est,
                     "proj.lk" = u.lk,
                     "proj.kl" = u.kl
                     )
  return(returnList)
}

# Estimate the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
# when the covariance for the targeted distribution is unknown
cov.inner.shift<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Pred.mat.target,Proj.array){
  N<-dim(Pred.mat.target)[1]
  Sigma.est.l1<-(1/dim(X.l1)[1])*(t(X.l1)%*%X.l1)
  Sigma.est.k1<-(1/dim(X.k1)[1])*(t(X.k1)%*%X.k1)
  var1<-0
  if(l2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[l2,k2,]/dim(X.l1)[1]
  }
  if(k2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[k2,l2,]/dim(X.l1)[1]
  }
  if(l2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[l2,k2,]/dim(X.k1)[1]
  }
  if(k2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[k2,l2,]/dim(X.k1)[1]
  }
  var2<-mean((diag(Pred.mat.target[,k1]%*%t(Pred.mat.target[,l1]))-mean(Pred.mat.target[,k1]*Pred.mat.target[,l1]))*(diag(Pred.mat.target[,k2]%*%t(Pred.mat.target[,l2]))-mean(Pred.mat.target[,k2]*Pred.mat.target[,l2])))
  var<-var1+var2/N
  return((var))
}


# Estimate the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
# when the covariance for the targeted distribution is known
cov.inner.shift.known<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Proj.array){
  Sigma.est.l1<-(1/dim(X.l1)[1])*(t(X.l1)%*%X.l1)
  Sigma.est.k1<-(1/dim(X.k1)[1])*(t(X.k1)%*%X.k1)
  var1<-0
  if(l2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[l2,k2,]/dim(X.l1)[1]
  }
  if(k2==l1){
    var1<-var1+Var.vec[l1]*Proj.array[l1,k1,]%*%Sigma.est.l1%*%Proj.array[k2,l2,]/dim(X.l1)[1]
  }
  if(l2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[l2,k2,]/dim(X.k1)[1]
  }
  if(k2==k1){
    var1<-var1+Var.vec[k1]*Proj.array[k1,l1,]%*%Sigma.est.k1%*%Proj.array[k2,l2,]/dim(X.k1)[1]
  }
  var<-var1
  return((var))
}

decide_delta <- function(Gamma, step_delta=0.1, MAX_iter=100, verbose=TRUE){
  L = dim(Gamma)[1]
  min_eigen = min(eigen(Gamma)$values)
  max_eigen = max(eigen(Gamma)$values)
  solution0 = opt.weight(Gamma, 0)
  reward0 = solution0$reward
  delta_min = 0.1 * reward0 * (2 * L) / (L - 1)
  if(delta_min > 2){
    if(verbose) print(paste("delta path starts from", round(delta_min, 4), "which exceeds our maximum limit 2"))
    delta = 2
  }else{
    delta = delta_min
    solution = opt.weight(Gamma, delta)
    reward = solution$reward
    i = 1
    while(reward >= 0.9 * reward0){
      delta_new = delta + step_delta
      solution = opt.weight(Gamma, delta_new)
      reward = solution$reward
      if(reward <= 0.9 * reward0) break
      if(delta_new > 2){
        warning(paste("The picked delta is over our maximum limit 2. Early Stopping at iteration", i))
        break
      }
      delta = delta_new
      i = i+1
      if((i %% 10 == 0)&verbose){
        print(paste("Iteration ", i, "delta =", round(delta, 4), "reward ratio = ", round(reward/reward0, 4)))
      }
      if(i >= MAX_iter){
        warning("Delta searching stops, because it reach the Max iterations.")
        break
      }
    }
  }
  
  if((min_eigen + delta) < 0.5) warning("Fail to find a suitable delta, the estimator may be not stable enough.")
  
  solution = opt.weight(Gamma, delta)
  reward = solution$reward
  if(verbose){
    print(paste0("The picked delta is ", round(delta,4)))
    print(paste0("Reward Ratio is ", round(reward / reward0, 4)))
    print(paste0("Minimum Eigenvalue plus delta = ", round(min_eigen + delta, 4)))
  }
  return(delta)
}

# 
# bootstrap <- function(idx.source, m,B=1000, replace=FALSE){
#   mat.resample = matrix(NA, nrow=m, ncol=B)
#   uni_groups = sort(unique(idx.source))
#   for(index.g in 1:length(uni_groups)){
#     g = uni_groups[index.g]
#     m.g = floor(m*(table(idx.source)[index.g] / length(idx.source)))
#     # print(paste("group=",g,"m.g=",m.g))
#     idx.g = which(idx.source==uni_groups[g])
#     idx.rows = seq((index.g-1)*m.g+1, index.g*m.g)  
#     mat.resample[idx.rows,] = replicate(B, sample(idx.g, size=m.g, replace=replace))
#   }
#   return(na.omit(mat.resample))
# }

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
