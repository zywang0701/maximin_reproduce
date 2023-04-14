## Linear Function: do bias-correction
LF<-function(X,y,loading,intercept.loading=FALSE, intercept=TRUE,init.Lasso=NULL,lambda=NULL,mu=NULL,step=NULL,resol = 1.5,maxiter=10){
  ### Option 1: search tuning parameter with steps determined by the ill conditioned case (n=p/2)
  ### Option 2: search tuning parameter with maximum 10 steps.
  ### Option 3: fixed tuning parameter and this is not recommended without exploring the tuning parameter selection
  xnew<-loading
  p <- ncol(X);
  n <- nrow(X);
  n_y <- length(y)

  if(n_y!=n)
  {
    stop("Error: Check dimensions of X and y")
  }
  else
  {
    data = na.omit(data.frame(y,X))
    X <- as.matrix(data[,-1])
    y <- as.vector(data[,1])
    p <- ncol(X)
    n <- nrow(X)

    # implement the correction of the initial estimator
    # set up the randomization step
    if (intercept==TRUE){
      Xc <- cbind(rep(1,n),X);
      pp <- (p+1);
    } else {
      Xc <- X;
      pp <- p
    }
    if(is.null(init.Lasso)){
      # implement a lasso algorithm to get beta and sigma
      init.Lasso<-Initialization.step(X,y,lambda,intercept)
      htheta<-init.Lasso$lasso.est
    }else{
      htheta<-init.Lasso
    }
    sparsity <- sum(abs(htheta) > 0.001)
    sd.est <- sqrt(sum((y - Xc %*% htheta)^2) / n)

    # compute the initial estimator
    if(intercept == TRUE){
      loading <- rep(0,pp)
      if(intercept.loading == TRUE){
        loading[1] <- 1
      }
      if(intercept.loading == FALSE){
        loading[1] <- 0
      }
      loading[-1] <- xnew
    } else {
      if(intercept.loading == TRUE){
        print(paste("Setting intercept = FALSE and intercept.loading = FALSE"))
      }
      loading <- xnew
    }
    loading.norm<-sqrt(sum(loading^2))
    lasso.plugin<-sum(loading*htheta)

    #############################################
    ############# Correction step ###############

    if ((n>=6*p)){
      sigma.hat <- (1/n)*(t(Xc)%*%Xc);
      tmp <- eigen(sigma.hat)
      tmp <- min(tmp$values)/max(tmp$values)
    }else{
      tmp <- 0
    }
    sigma.hat <- (1/n)*(t(Xc)%*%Xc);
    if ((n>=6*p)&&(tmp>=1e-4)){
      direction <- solve(sigma.hat)%*%loading/loading.norm
    }else{
      if(is.null(step)){
        step.vec<-rep(NA,3)
        for(t in 1:3){
          index.sel<-sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
          Direction.Est.temp<-Direction_searchtuning_lin(Xc[index.sel,],loading,mu=NULL, resol, maxiter)
          step.vec[t]<-Direction.Est.temp$step
        }
        step<-getmode(step.vec)
      }
      Direction.Est<-Direction_fixedtuning_lin(Xc,loading,mu=sqrt(2.01*log(pp)/n)*resol^{-(step-1)})
      while(is.na(Direction.Est) || length(Direction.Est$proj)==0){
        step<-step-1
        Direction.Est <- Direction_fixedtuning_lin(Xc, loading, mu = sqrt(2.01 * log(pp) / n) * resol^{-(step - 1)})
      }
      direction<-Direction.Est$proj
    }
    correction = t(Xc%*%direction)%*%(y - Xc%*%htheta)/n;
    debias.est=lasso.plugin+correction*loading.norm
    se<-sd.est*sqrt(sum((Xc%*%direction)^2)/(n)^2)*loading.norm

    returnList <- list("prop.est" = debias.est,
                       "se" = se,
                       "proj"=direction,
                       "plug.in"=lasso.plugin
    )
    return(returnList)
  }
}


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
