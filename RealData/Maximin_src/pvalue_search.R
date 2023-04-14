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
    # gen.point = sum(object$Point.vec * gen.weight.vector)
    # gen.se = sqrt(sum(gen.weight.vector^2 * object$SE.vec^2))
    # gen.est[g, 1] = gen.point
    # gen.est[g, 2] = gen.se
  }
  
  # CI.original = cbind(gen.est[, 1]-qnorm(1-alpha/2)*gen.est[, 2], gen.est[,1]+qnorm(1-alpha/2)*gen.est[,2])
  # CI = na.omit(CI.original)
  # uni = Intervals(CI)
  # CI.union = as.matrix(interval_union(uni))
  # colnames(CI.union) <- c("lower", "upper")

  out = list(gen.weight.mat = gen.weight.mat)
  return(out)
}
# 
# gen.est.build <- function(object, weight.mat){
#   gen.size = nrow(weight.mat)
#   gen.est = matrix(NA, nrow=gen.size, ncol=2)
#   for(i in 1:gen.size){
#     gen.weight.vector = weight.mat[i,]
#     gen.point = sum(object$Point.vec * gen.weight.vector)
#     gen.se = sqrt(sum(gen.weight.vector^2 * object$SE.vec^2))
#     gen.est[g, 1] = gen.point
#     gen.est[g, 2] = gen.se
#   }
#   return(gen.est)
# }

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