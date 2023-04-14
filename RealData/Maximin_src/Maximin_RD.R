MaximinRD <- function(Xlist, Ylist, loading.idx, LFest, X.target=NULL, cov.target=NULL,
                      covariate.shift=TRUE, Coef.est){

  ############################################
  ########### Transfer Source Data ###########
  ############################################
  if((is.list(Xlist)==FALSE)||(is.list(Ylist)==FALSE)) stop("Error: check the type of Xlist and Ylist, they must be list")
  if(length(Xlist)!=length(Ylist)) stop("Error: check the length of Xlist and Ylist")
  L = length(Xlist)
  n.x.vec = rep(0, L)
  n.y.vec = rep(0, L)
  p.vec = rep(0, L)
  for(l in 1:L){
    n.x.vec[l] = dim(Xlist[[l]])[1]
    p.vec[l] = dim(Xlist[[l]])[2]
    n.y.vec[l] = length(Ylist[[l]])
    if(n.x.vec[l]!=n.y.vec[l]) stop(paste("Error: X and Y to the group",l,"has different number of samples"))
  }
  if(!all(p.vec==p.vec[1])) stop("Error: check the dimension p of each X, they must be the same")
  ## centralize Xs
  X.source = do.call(rbind, Xlist)
  Y.source = do.call(c, Ylist)
  idx.source = rep(1:L, times=n.y.vec)

  ##################################
  ########### Check Input ##########
  ##################################
  if((!is.null(X.target))&&(dim(X.source)[2]!=dim(X.target)[2])){
    stop("Error: check dimensions of Target Data")
  }
  if((!is.null(cov.target))&&((dim(X.source)[2]!=dim(cov.target)[2])||(dim(cov.target)[1]!=dim(cov.target)[2]))){
    stop("Error: check dimensions of cov.target")
  }

  if(is.null(X.target)) X.target = X.source  # the code adapts to No target setting
  n.source = nrow(X.source)
  n.target = nrow(X.target)
  p = ncol(X.source)
  L = length(unique(idx.source))
  uni_groups = sort(unique(idx.source))

  Pred.vec = rep(0, n.source)  # predicted outcome of source data
  Pred.mat.target = matrix(0, n.target, L)  # predicted outcomes of target data
  Point.vec = rep(0, L)  # debiased point estimators
  Var.vec = rep(0, L)  # variance of residual
  SE.vec = rep(0, L)  # SE of loading
  for(l in 1:L){
    ## obtain estimators of group l using Lasso
    index.set = which(idx.source==uni_groups[l])
    X = X.source[index.set, ]
    Y = Y.source[index.set]

    Pred.vec[index.set] = X%*%Coef.est[, l]
    Pred.mat.target[, l] = X.target%*%Coef.est[, l]
    ## obtain variance of residual for group l
    supp.l = which(abs(Coef.est[, l])>0.01)
    n.eff = max(0.9*nrow(X), nrow(X)-length(supp.l))
    Var.vec[l] = sum((Y - Pred.vec[index.set])^2) / n.eff
    est <- LFest[[l]]
    SE.vec[l] = est$se.vec[loading.idx]
    Point.vec[l] = est$est.debias.vec[loading.idx]
  }
  #########################################
  ######### compte Gamma.plugin ###########
  #########################################
  if(!is.null(cov.target)){
    Sigma.target.est = matrix(0, p, p)
    Sigma.target.est = cov.target
  }else{
    if(covariate.shift){
      Sigma.target.est = (1/n.target)*(t(X.target)%*%X.target)
    }else{
      Sigma.target.est = (t(X.source)%*%X.source + t(X.target)%*%X.target)/(n.source + n.target)
    }
  }
  Gamma.plugin = t(Coef.est)%*%Sigma.target.est%*%Coef.est
  Omega.est = Sigma.target.est%*%Coef.est

  ####################################################
  ##### conduct bias correction for Gamma.plugin #####
  ####################################################
  Gamma.prop = Gamma.plugin
  Proj.array = array(NA, dim=c(L, L, p))
  for(l in 1:L){
    for(k in l:L){
      index.set.l = which(idx.source==uni_groups[l])
      index.set.k = which(idx.source==uni_groups[k])
      X.l = X.source[index.set.l, ]
      X.k = X.source[index.set.k, ]
      Y.l = Y.source[index.set.l]
      Y.k = Y.source[index.set.k]
      Pred.l = Pred.vec[index.set.l]
      Pred.k = Pred.vec[index.set.k]

      if(covariate.shift){
        output <- Gamma.shift(Gamma.plugin[l, k], X.l, X.k, Omega.est[, l], Omega.est[, k],
                              Y.l, Y.k, Pred.l, Pred.k)
        Gamma.prop[l, k] = output$est
        Proj.array[l, k, ] = output$proj.lk
        Proj.array[k, l, ] = output$proj.kl
      }else{
        Gamma.prop[l, k] = Gamma.plugin[l, k]+t(Coef.est[,l])%*%t(X.k)%*%(Y.k-Pred.k)/nrow(X.k)+
          t(Coef.est[,k])%*%t(X.l)%*%(Y.l-Pred.l)/nrow(X.l)
        Proj.array[l, k, ] = Coef.est[,k]
        Proj.array[k, l, ] = Coef.est[,l]
      }
    }
  }
  # for(l in 2:L){
  #   for(k in 1:(l-1)){
  #     Gamma.prop[l, k] = Gamma.prop[k, l]
  #   }
  # }
  for(l in 1:L){
    for(k in 1:l){
      Gamma.prop[l, k] = Gamma.prop[k, l]
    }
  }
  ######################################################################
  ################## to obtain sampling materials ######################
  ## compute mean and covariance matrix for the sampling distribution ##
  ######################################################################
  gen.mu = Gamma.prop[lower.tri(Gamma.prop, diag=TRUE)]
  gen.dim = L*(L+1)/2
  gen.Cov = matrix(NA, nrow=gen.dim, ncol=gen.dim)
  for(k1 in 1:L){
    for(l1 in k1:L){
      index1 = index.map(L, l1, k1)
      for(k2 in 1:L){
        for(l2 in k2:L){
          index2 = index.map(L, l2, k2)
          index.set.l1 = which(idx.source==uni_groups[l1])
          index.set.k1 = which(idx.source==uni_groups[k1])
          index.set.l2 = which(idx.source==uni_groups[l2])
          index.set.k2 = which(idx.source==uni_groups[k2])
          X.l1 = X.source[index.set.l1, ]
          X.k1 = X.source[index.set.k1, ]
          X.l2 = X.source[index.set.l2, ]
          X.k2 = X.source[index.set.k2, ]

          if(!is.null(cov.target)){
            gen.Cov[index1, index2] <- cov.inner.shift.known(Var.vec, l1, k1, l2, k2,
                                                             X.l1, X.k1, X.l2, X.k2,
                                                             Proj.array)
          }else{
            gen.Cov[index1, index2] <- cov.inner.shift(Var.vec, l1, k1, l2, k2,
                                                       X.l1, X.k1, X.l2, X.k2,
                                                       Pred.mat.target, Proj.array)
          }
        }
      }
    }
  }
  tau = 0.2
  gen.Cov = gen.Cov + diag(max(tau*diag(gen.Cov), 1/floor(n.source/L)), dim(gen.Cov)[2])

  out = list(Gamma.prop = Gamma.prop,
             Coef.est = Coef.est,
             Point.vec = Point.vec,
             SE.vec = SE.vec,
             L = L,
             gen.mu = gen.mu,
             gen.Cov = gen.Cov)
  structure(out, class = "Maximin")
}


infer <- function(object, delta=0, gen.size=500, threshold=0, alpha=0.05, alpha.thres=0.01){
  ##############################
  ####### Maximin Effects ######
  ##############################
  # aggregated weights
  solution = opt.weight(object$Gamma.prop, delta, report.reward=FALSE)
  weight = solution$weight
  # point estimation of <loading, \beta_\delta^*>
  point = sum(object$Point.vec * weight)
  # maximin effect
  mm.effect = object$Coef.est %*% weight

  #############################
  ##### Generate Samples ######
  #############################
  gen.dim = object$L*(object$L+1)/2
  if(threshold==0){
    thres = qnorm(1-alpha.thres/(gen.dim*2))
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0 # records how many samples are picked
    while(n.picked < gen.size){
      S = mvrnorm(1, mu=rep(0, gen.dim), Sigma=object$gen.Cov)
      if(max(abs(S / sqrt(diag(object$gen.Cov)))) <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = object$gen.mu + S
      }
    }
  }
  if(threshold==1){
    # 1 stands for chi square threshold
    UDV_list = svd(object$gen.Cov)
    U = UDV_list$u
    D = UDV_list$d
    V = UDV_list$v
    D.sqrt = sqrt(D)
    gen.Cov.sqrt = U %*% diag(D.sqrt) %*% t(V)
    thres = 1.0 * qchisq((1-alpha.thres), df=gen.dim)
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0
    while(n.picked < gen.size){
      Z = mvrnorm(1, mu=rep(0, gen.dim), Sigma=diag(gen.dim))
      Z.normsq = sum(Z^2)
      if(Z.normsq <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = object$gen.mu + gen.Cov.sqrt %*% Z
      }
    }
  }
  if(threshold==2){
    gen.samples = matrix(mvrnorm(gen.size, mu=object$gen.mu, Sigma=object$gen.Cov),
                         nrow=gen.size, ncol=gen.dim)
  }

  ##################################
  ######### Construct CI ###########
  ##################################
  gen.weight.mat = matrix(NA, nrow=gen.size, ncol=object$L)
  gen.est = matrix(NA, nrow=gen.size, ncol=2)
  # construct CI
  for(g in 1:gen.size){
    gen.matrix = matrix(NA, nrow=object$L, ncol=object$L)
    gen.matrix[lower.tri(gen.matrix, diag=TRUE)] = gen.samples[g, ]

    for(l in 1:object$L){
      for(k in l:object$L){
        gen.matrix[l,k] = gen.matrix[k,l]
      }
    }
    gen.solution = opt.weight(gen.matrix, delta, report.reward=FALSE)
    gen.weight.vector = gen.solution$weight
    gen.weight.mat[g, ] = gen.weight.vector
    gen.point = sum(object$Point.vec * gen.weight.vector)
    gen.se = sqrt(sum(gen.weight.vector^2 * object$SE.vec^2))
    gen.est[g, 1] = gen.point
    gen.est[g, 2] = gen.se
  }

  CI.original = cbind(gen.est[, 1]-qnorm(1-alpha/2)*gen.est[, 2], gen.est[,1]+qnorm(1-alpha/2)*gen.est[,2])
  CI = na.omit(CI.original)
  uni = Intervals(CI)
  CI.union = as.matrix(interval_union(uni))
  colnames(CI.union) <- c("lower", "upper")

  out = list(weight = weight,
             point = point,
             mm.effect = mm.effect,
             CI = CI.union)
  return(out)
}


######################### Source ##############################
opt.weight<-function(Gamma,delta,report.reward=TRUE){
  ## Purpose: Compute Ridge-type weight vector
  ## Returns: weight:  the minimizer \eqn{\gamma}
  ##          reward:  the value of penalized reward
  ## ----------------------------------------------------
  ## Arguments: Gamma: regression covariance matrix, of dimension \eqn{L} x \eqn{L}
  ##            delta the ridge penalty level, non-positive.
  ##            report.reward the reward is computed or not (Default = `TRUE`)
  ## ----------------------------------------------------

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

Gamma.shift<-function(plug.in,X.l,X.k,omega.l,omega.k,Y.l,Y.k,Pred.l,Pred.k){
  ## Purpose: Bias correction for initial estimator of Gamma when covariate shifts
  ## Returns: est:       The proposed bias-corrected estimator of Gamma target
  ##          proj.lk:   The projection direction of index (l, k)
  ##          proj.kl:   The projection direction of index (k, l)
  ## ----------------------------------------------------------------------
  ## Arguments: plug.in: Initial estimator of Gamma target, of dimension \eqn{L} x \eqn{L}
  ##            X.l:     Design matrix of label \eqn{l} in training data, of dimension \eqn{n_l} x \eqn{p}
  ##            X.k:     Design matrix of label \eqn{k} in training data, of dimension \eqn{n_k} x \eqn{p}
  ##            omega.l: The l-th column of Omega matrix
  ##            omega.k: The k-th column of Omega matrix
  ##            Y.l:     Outcome vector of label \eqn{l} in training data, of length \eqn{n_l}
  ##            Y.k:     Outcome vector of label \eqn{k} in training data, of length \eqn{n_k}
  ##            Pred.l:  Predicted outcome vector of label \eqn{l} in training data
  ##            Pred.k:  Predicted outcome vector of label \eqn{k} in training data
  ## ----------------------------------------------------------------------
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

cov.inner.shift<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Pred.mat.target,Proj.array){
  ## Purpose: Estimate the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
  ##          when the covariance for the target distribution is unknown
  ## Returns: covariance between the pi(l1, k1) entry and pi(l2, k2) entry
  ## ----------------------------------------------------------------------
  ## Arguments: Var.vec:         Variance of residuals in groups, of length \eqn{L}
  ##            l1:              Index l1
  ##            k1:              Index k1
  ##            l2:              Index l2
  ##            k2:              Index k2
  ##            X.l1:            Design matrix of group \eqn{l1} in training data
  ##            X.k1:            Design matrix of group \eqn{k1} in training data
  ##            X.l2:            Design matrix of group \eqn{l2} in training data
  ##            X.k2:            Design matrix of group \eqn{k2} in training data
  ##            Pred.mat.target: Predicted outcome matrix for target design matrix
  ##                             L fitted coefficients, of dimension \eqn{n.target} x \eqn{L}
  ##            Proj.array:      Projection directions, of dimension \eqn{L} x \eqn{L} x \eqn{p}
  ## ----------------------------------------------------------------------
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

cov.inner.shift.known<-function(Var.vec,l1,k1,l2,k2,X.l1,X.k1,X.l2,X.k2,Proj.array){
  ## Purpose: Estimate the covariance between the pi(l1,k1) entry and pi(l2,k2) entry
  ##          when the covariance for the target distribution is known
  ## Returns: covariance between the pi(l1, k1) entry and pi(l2, k2) entry
  ## ----------------------------------------------------------------------
  ## Arguments: Var.vec:         Variance of residuals in groups, of length \eqn{L}
  ##            l1:              Index l1
  ##            k1:              Index k1
  ##            l2:              Index l2
  ##            k2:              Index k2
  ##            X.l1:            Design matrix of group \eqn{l1} in training data
  ##            X.k1:            Design matrix of group \eqn{k1} in training data
  ##            X.l2:            Design matrix of group \eqn{l2} in training data
  ##            X.k2:            Design matrix of group \eqn{k2} in training data
  ##            Proj.array:      Projection directions, of dimension \eqn{L} x \eqn{L} x \eqn{p}
  ## ----------------------------------------------------------------------
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

Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)

  htheta <- if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  }
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}

getmode <- function(v) {
  tbl <- table(v)
  if (all(tbl == 1)) {
    median(v)
  } else {
    as.numeric(names(which.max(tbl)))
  }
}

# Direction_fixedtuning_lin<-function(X,loading,mu=NULL){
#   pp<-ncol(X)
#   n<-nrow(X)
#   if(is.null(mu)){
#     mu<-sqrt(2.01*log(pp)/n)
#   }
#   loading.norm<-sqrt(sum(loading^2))
#   if (loading.norm==0){
#     H <- cbind(loading, diag(1, pp))
#   }else{
#     H <- cbind(loading / loading.norm, diag(1, pp))
#   }
#   v<-Variable(pp+1)
#   obj<-1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
#   prob<-Problem(Minimize(obj))
#   result<-solve(prob)
#   opt.sol<-result$getValue(v)
#   cvxr_status<-result$status
#   direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
#   returnList <- list("proj"=direction)
#   return(returnList)
# }
Direction_fixedtuning_lin<-function(X, loading, mu = NULL){
  pp <- ncol(X)
  n <- nrow(X)
  if(is.null(mu)){
    mu <- sqrt(2.01*log(pp)/n)
  }
  loading.norm <- sqrt(sum(loading^2))
  if (loading.norm == 0){
    H <- cbind(loading, diag(1, pp))
  }else{
    H <- cbind(loading / loading.norm, diag(1, pp))
  }
  v <- Variable(pp+1)
  obj <- 1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  if(result$status=="optimal" || result$status == "unbounded"){
    opt.sol<-result$getValue(v)
    cvxr_status<-result$status
    direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  }else{
    direction <- numeric(0)
  }
  returnList <- list("proj"=direction)
  return(returnList)
}

Direction_searchtuning_lin <- function(X, loading, mu=NULL, resol = 1.5, maxiter = 6){
  pp <- ncol(X)
  n <- nrow(X)
  tryno <- 1
  opt.sol <- rep(0,pp+1)
  lamstop <- 0
  cvxr_status <- "optimal"
  mu <- sqrt(2.01*log(pp)/n)
  while (lamstop == 0 && tryno < maxiter){
    lastv <- opt.sol
    lastresp <- cvxr_status
    loading.norm <- sqrt(sum(loading^2))
    if (loading.norm == 0){
      H <- cbind(loading, diag(1, pp))
    }else{
      H <- cbind(loading / loading.norm, diag(1, pp))
    }
    v <- Variable(pp+1)
    obj <- 1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
    prob <- Problem(Minimize(obj))
    result <- solve(prob)
    cvxr_status <- result$status
    if(tryno == 1){
      if(cvxr_status == "optimal"){
        incr <- 0
        mu <- mu/resol
        opt.sol <- result$getValue(v)
        temp.vec <- (-1)/2*(opt.sol[-1] + opt.sol[1]*loading/loading.norm)
        initial.sd <- sqrt(sum((X%*% temp.vec)^2)/(n)^2)*loading.norm
        temp.sd <- initial.sd
      }else{
        incr = 1
        mu = mu*resol
      }
    }else{
      if(incr == 1){
        if(cvxr_status == "optimal"){
          opt.sol <- result$getValue(v)
          lamstop <- 1
        }else{
          mu <- mu*resol
        }
      }else{
        if(cvxr_status == "optimal" && temp.sd < 3*initial.sd){
          mu <- mu/resol
          opt.sol <- result$getValue(v)
          temp.vec <- (-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
          temp.sd <- sqrt(sum((X%*% temp.vec)^2)/(n)^2)*loading.norm
        }else{
          mu <- mu*resol
          opt.sol <- lastv
          lamstop <- 1
          tryno <- tryno-1
        }
      }
    }
    tryno <- tryno + 1
  }
  direction <- (-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  step <- tryno-1
  returnList <- list("proj"=direction,
                     "step"=step)
  return(returnList)
}

index.map<-function(L,l,k){
  return((2*L-k)*(k-1)/2+l)
}

proj.direction<-function(Xc,loading,maxiter=6,resol=1.25){
  n<-dim(Xc)[1]
  p<-dim(Xc)[2]
  loading.norm<-sqrt(sum(loading^2))
  sigma.hat <- (1/n)*(t(Xc)%*%Xc);
  if ((n>=6*p)){
    tmp <- eigen(sigma.hat)
    tmp <- min(tmp$values)/max(tmp$values)
  }else{
    tmp <- 0
  }
  if ((n>=6*p)&&(tmp>=1e-4)){
    direction <- solve(sigma.hat)%*%loading
  }else{
    step.vec<-rep(NA,3)
    for(t in 1:3){
      index.sel<-sample(1:n,size=ceiling(0.5*min(n,p)), replace=FALSE)
      Direction.Est.temp<-Direction_searchtuning_lin(Xc[index.sel,],loading,mu=NULL, resol, maxiter)
      step.vec[t]<-Direction.Est.temp$step
    }
    step<-getmode(step.vec)
    Direction.Est<-Direction_fixedtuning_lin(Xc,loading,mu=sqrt(2.01*log(p)/n)*resol^{-(step-1)})
    while(is.na(Direction.Est) || length(Direction.Est$proj)==0){
      step<-step-1
      Direction.Est <- Direction_fixedtuning_lin(Xc, loading, mu = sqrt(2.01 * log(p) / n) * resol^{-(step - 1)})
    }

    direction<-loading.norm*Direction.Est$proj
  }
  return(direction)
}

LF <- function(X, y,loading, intercept.loading = TRUE, intercept = TRUE, init.coef = NULL, lambda = NULL, mu = NULL, step = NULL, resol = 1.5, maxiter = 6, alpha = 0.05, verbose = TRUE){
  xnew <- loading
  p <- ncol(X)
  n <- nrow(X)
  n_y <- length(y)

  if(n_y!=n)
  {
    stop("Error : Check dimensions of X and y")
  } else {
    data <- na.omit(data.frame(y, X))
    X <- as.matrix(data[,-1])
    y <- as.vector(data[,1])
    p <- ncol(X)
    n <- nrow(X)
    mean = colMeans(X)
    M = matrix(rep(mean,nrow(X)),byrow = T, nrow = nrow(X), ncol = ncol(X))
    X = X - M
    if(intercept.loading == TRUE && intercept == TRUE){
      xnew = xnew - mean
    }
    col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
    Xnor <- X %*% diag(col.norm)
    if(is.null(init.coef)){
      ####### implement a lasso algorithm to get beta and sigma
      init.coef <-  Initialization.step(X, y, lambda, intercept)
      htheta <- init.coef$lasso.est
    } else {
      htheta <- init.coef
    }
    if (intercept == TRUE){
      Xb <- cbind(rep(1,n),Xnor)
      Xc <- cbind(rep(1,n),X)
      pp <- (p+1)
    } else {
      Xb <- Xnor
      Xc <- X
      pp <- p
    }

    sparsity <- sum(abs(htheta) > 0.001)
    sd.est <- sqrt(sum((y - Xb %*% htheta)^2) / max(0.9*n, n - sparsity))

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
    loading.norm <- sqrt(sum(loading^2))
    lasso.plugin <- sum(loading*htheta)

    count <- 0
    for(i in 1:ncol(X)){
      if(length(unique(X[,i])) == 1){
        count <- count + 1
      }
    }
    if(count!=0 && intercept == TRUE)
    {
      stop("Data is singular")
    } else {
      if ((n >= 6*p)){
        sigma.hat <- (1/n)*(t(Xc)%*%Xc)
        tmp <- eigen(sigma.hat)
        tmp <- min(tmp$values)/max(tmp$values)
      } else {
        tmp <- 0
      }
      sigma.hat <- (1/n)*(t(Xc)%*%Xc)
      if ((n >= 6*p) && (tmp >= 1e-4)){
        direction <- solve(sigma.hat)%*%loading/loading.norm
      } else {
        if(is.null(step)){
          step.vec <- rep(NA,3)
          for(t in 1:3){
            index.sel <- sample(1:n, size=ceiling(0.5*min(n,p)), replace=FALSE)
            Direction.Est.temp <-  Direction_searchtuning_lin(Xc[index.sel,], loading, mu = NULL, resol, maxiter)
            step.vec[t] <- Direction.Est.temp$step
          }
          step <-  getmode(step.vec)
        }
        Direction.Est<- Direction_fixedtuning_lin(Xc, loading, mu = sqrt(2.01*log(pp)/n)*resol^{-(step-1)})
        while(is.na(Direction.Est) || length(Direction.Est$proj)==0){
          step <- step-1
          Direction.Est <-  Direction_fixedtuning_lin(Xc, loading, mu = sqrt(2.01 * log(pp) / n) * resol^{-(step - 1)})
        }
        if(verbose == TRUE){
          print(paste("step is", step))
        }
        direction <- Direction.Est$proj
      }
      correction <- t(Xc%*%direction)%*%(y - Xc%*%htheta)/n
      debias.est <- lasso.plugin + correction*loading.norm
      se <- sd.est*sqrt(sum((Xc%*%direction)^2)/(n)^2)*loading.norm
      CI <- c(debias.est - qnorm(1-alpha/2)*se, debias.est + qnorm(1-alpha/2)*se)
      if(debias.est - qnorm(1-alpha)*se > 0){
        dec <- 1
      }else{
        dec <- 0
      }
      returnList <- list("prop.est" = debias.est,
                         "se" = se,
                         "CI" = CI,
                         "decision" = dec,
                         "proj" = direction,
                         "plug.in" = lasso.plugin
      )
      return(returnList)
    }
  }
}

Initialization.step <- function(X, y, lambda = NULL, intercept = FALSE) {
  n <- nrow(X)
  col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
  Xnor <- X %*% diag(col.norm)
  htheta <- Lasso(Xnor, y, lambda = lambda, intercept = intercept)
  if (intercept == TRUE) {
    Xb <- cbind(rep(1, n), Xnor)
    col.norm <- c(1, col.norm)
  } else {
    Xb <- Xnor
  }
  htheta <- htheta * col.norm
  returnList <- list("lasso.est" = htheta)
  return(returnList)
}

measure_instability <- function(object, delta=0, gen.size=500, threshold=0, alpha.thres=0.01){
  spnorm.Gamma.diff = rep(0, gen.size)
  l2norm.gamma.diff = rep(0, gen.size)

  #############################
  ##### Generate Samples ######
  #############################
  gen.dim = object$L*(object$L+1)/2
  if(threshold==0){
    thres = qnorm(1-alpha.thres/(gen.dim*2))
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0 # records how many samples are picked
    while(n.picked < gen.size){
      S = mvrnorm(1, mu=rep(0, gen.dim), Sigma=object$gen.Cov)
      if(max(abs(S / sqrt(diag(object$gen.Cov)))) <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = object$gen.mu + S
      }
    }
  }
  if(threshold==1){
    # 1 stands for chi square threshold
    UDV_list = svd(object$gen.Cov)
    U = UDV_list$u
    D = UDV_list$d
    V = UDV_list$v
    D.sqrt = sqrt(D)
    gen.Cov.sqrt = U %*% diag(D.sqrt) %*% t(V)
    thres = 1.0 * qchisq((1-alpha.thres), df=gen.dim)
    gen.samples = matrix(0, nrow=gen.size, ncol=gen.dim)
    n.picked = 0
    while(n.picked < gen.size){
      Z = mvrnorm(1, mu=rep(0, gen.dim), Sigma=diag(gen.dim))
      Z.normsq = sum(Z^2)
      if(Z.normsq <= thres){
        n.picked = n.picked + 1
        gen.samples[n.picked, ] = object$gen.mu + gen.Cov.sqrt %*% Z
      }
    }
  }
  if(threshold==2){
    gen.samples = matrix(mvrnorm(gen.size, mu=object$gen.mu, Sigma=object$gen.Cov),
                         nrow=gen.size, ncol=gen.dim)
  }
  ################################################
  ###### Compute Gamma.diff and gamma.diff #######
  ################################################
  gen.Gamma.array = array(NA, dim=c(gen.size, object$L, object$L))
  gen.weight.mat = matrix(NA, nrow=gen.size, object$L)

  solution = opt.weight(object$Gamma.prop, delta, report.reward=FALSE)
  weight.prop = solution$weight

  for(g in 1:gen.size){
    gen.matrix = matrix(NA, nrow=object$L, ncol=object$L)
    gen.matrix[lower.tri(gen.matrix, diag=TRUE)] = gen.samples[g,]
    for(l in 1:object$L){
      for(k in l:object$L){
        gen.matrix[l, k] = gen.matrix[k, l]
      }
    }
    gen.Gamma.array[g,,] = gen.matrix
    gen.solution = opt.weight(gen.matrix, delta, report.reward=FALSE)
    gen.weight.mat[g,] = gen.solution$weight
  }
  for(g in 1:gen.size){
    Gamma.diff = gen.Gamma.array[g,,] - object$Gamma.prop
    spnorm.Gamma.diff[g] = sqrt(max((eigen(Gamma.diff)$values)^2))
    gamma.diff = gen.weight.mat[g,] - weight.prop
    l2norm.gamma.diff[g] = sqrt(sum(gamma.diff^2))
  }
  measure = mean(l2norm.gamma.diff^2 / spnorm.Gamma.diff^2)
  out <- list(measure = measure)
  return(out)
}
