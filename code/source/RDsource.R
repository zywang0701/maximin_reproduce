### MaximinRD
### The main function to prepare materials, used in Real Data application.
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

### infer
### after running "MaximinRD()", it is used to do inference.
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

### gen.prepare
### generate samples and compute the optimal weights for each sample.
gen.prepare <- function(object, delta=0, gen.size=500, threshold=0, alpha.thres=0.01){
  ##############################
  ####### Maximin Effects ######
  ##############################
  # aggregated weights
  solution = opt.weight(object$Gamma.prop, delta, report.reward=FALSE)
  weight = solution$weight
  # # point estimation of <loading, \beta_\delta^*>
  # point = sum(object$Point.vec * weight)
  # # maximin effect
  # mm.effect = object$Coef.est %*% weight
  
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
      for(k in 2:object$L){
        gen.matrix[l, k] = gen.matrix[k, l]
      }
    }
    gen.solution = opt.weight(gen.matrix, delta, report.reward=FALSE)
    gen.weight.vector = gen.solution$weight
    gen.weight.mat[g, ] = gen.weight.vector
  }
  
  out = list(gen.weight.mat = gen.weight.mat)
  return(out)
}

### pvalue.search
### In real data application, it used to search p-value for maximin effect.
pvalue.search <- function(gen.est, tol=1e-4){
  alpha0 = 0
  alpha1 = 1
  alpha_try = (alpha0+alpha1)/2
  step = 1
  while((alpha1-alpha0)>tol){
    if(pvalue.test(gen.est, alpha_try)){
      alpha0 = alpha_try
    }else{
      alpha1 = alpha_try
    }
    alpha_try = (alpha0+alpha1)/2
    step = step+1
  }
  return(list(alpha_try = alpha_try,
              alpha0 = alpha0,
              alpha1 = alpha1,
              step = step))
}

### pvalue.test
### it is a util function for pvalue.search
pvalue.test <- function(gen.est, alpha){
  CI.ori = cbind(gen.est[,1] - qnorm(1-alpha/2)*gen.est[,2], 
                 gen.est[,1] + qnorm(1-alpha/2)*gen.est[,2])
  CI = na.omit(CI.ori)
  CI = c(min(CI[,1]), max(CI[,2]))
  if((CI[1]<0)&(CI[2]>0)){
    include=TRUE
  }else{
    include=FALSE
  }
  return(include)
}