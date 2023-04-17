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

getmode <- function(v) {
  tbl <- table(v)
  if (all(tbl == 1)) {
    median(v)
  } else {
    as.numeric(names(which.max(tbl)))
  }
}

Lasso <- function(X, y, lambda = NULL, intercept = TRUE) {
  p <- ncol(X)
  n <- nrow(X)

  htheta <- if (is.null(lambda)) {
    lambda <- sqrt(qnorm(1 - (0.1 / p)) / n)
    outLas <- slim(X, y, lambda = lambda, method = "lq", q = 2,
                   verbose = FALSE)
    # Objective : sqrt(RSS/n) + lambda * penalty
    c(as.vector(outLas$intercept), as.vector(outLas$beta))
  } else if (lambda == "CV") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.1se))
  } else if (lambda == "CV.min") {
    outLas <- cv.glmnet(X, y, family = "gaussian", alpha = 1,
                        intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = outLas$lambda.min))
  } else if (lambda == "scalreg") {
    Xc <- if (intercept) {
      cbind(rep(1, n), X)
    } else {
      X
    }
    outLas <- scalreg(Xc, y)
    # return object
    if (intercept) {
      outLas$coefficients
    } else {
      # add a coefficient for the (not estimated) intercept b/c of implementation
      c(0, outLas$coefficients)
    }
  } else {
    outLas <- glmnet(X, y, family = "gaussian", alpha = 1,
                     intercept = intercept)
    # Objective : 1/2 * RSS/n + lambda * penalty
    as.vector(coef(outLas, s = lambda))
  }

  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}

Initialization.step <- function(X, y, lambda = NULL, intercept = FALSE) {
  n <- nrow(X)
  col.norm <- 1 / sqrt((1 / n) * diag(t(X) %*% X))
  Xnor <- X %*% diag(col.norm)

  ### Call Lasso
  htheta <- Lasso(Xnor, y, lambda = lambda, intercept = intercept)

  ### Calculate return quantities
  if (intercept == TRUE) {
    Xb <- cbind(rep(1, n), Xnor)
    col.norm <- c(1, col.norm)
  } else {
    Xb <- Xnor
  }
  sparsity <- sum(abs(htheta) > 0.001)
  sd.est <- sqrt(sum((y - Xb %*% htheta)^2) / n)
  htheta <- htheta * col.norm
  returnList <- list("lasso.est" = htheta,
                     "sigma" = sd.est,
                     "sparsity" = sparsity)
  return(returnList)
}

Direction_fixedtuning_lin<-function(X,loading,mu=NULL){
  pp<-ncol(X)
  n<-nrow(X)
  if(is.null(mu)){
    mu<-sqrt(2.01*log(pp)/n)
  }
  loading.norm<-sqrt(sum(loading^2))
  if (loading.norm==0){
    H <- cbind(loading, diag(1, pp))
  }else{
    H <- cbind(loading / loading.norm, diag(1, pp))
  }
  v<-Variable(pp+1)
  obj<-1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
  prob<-Problem(Minimize(obj))
  result<-solve(prob)
  opt.sol<-result$getValue(v)
  cvxr_status<-result$status
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  returnList <- list("proj"=direction)
  return(returnList)
}

Direction_searchtuning_lin<-function(X,loading,mu=NULL, resol = 1.5, maxiter = 10){
  pp<-ncol(X)
  n<-nrow(X)
  tryno = 1;
  opt.sol = rep(0,pp+1);
  lamstop = 0;
  cvxr_status = "optimal";
  mu = sqrt(2.01*log(pp)/n);
  #mu.initial= mu;
  while (lamstop == 0 && tryno < maxiter){
    ###### This iteration is to find a good tuning parameter
    lastv = opt.sol;
    lastresp = cvxr_status;
    loading.norm<-sqrt(sum(loading^2))
    if (loading.norm==0){
      H <- cbind(loading, diag(1, pp))
    }else{
      H <- cbind(loading / loading.norm, diag(1, pp))
    }
    v<-Variable(pp+1)
    obj<-1/4*sum((X%*%H%*%v)^2)/n+sum((loading/loading.norm)*(H%*%v))+mu*sum(abs(v))
    prob<-Problem(Minimize(obj))
    result<-solve(prob)
    cvxr_status<-result$status

    if(tryno==1){
      if(cvxr_status=="optimal"){
        incr = 0;
        mu=mu/resol;
        opt.sol<-result$getValue(v) ### we should move this line from above to here
        temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
        initial.sd<-sqrt(sum((X%*% temp.vec)^2)/(n)^2)*loading.norm ##what's this?
        temp.sd<-initial.sd
      }else{
        incr = 1;
        mu=mu*resol;
      }
    }else{
      if(incr == 1){ ### if the tuning parameter is increased in the last step
        if(cvxr_status=="optimal"){
          lamstop = 1;
          opt.sol<-result$getValue(v)
        }else{
          mu=mu*resol;
        }
      }else{
        if(cvxr_status=="optimal"&&temp.sd<3*initial.sd){ ##Why this condition on sd?
          mu = mu/resol;
          opt.sol<-result$getValue(v)
          temp.vec<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
          temp.sd<-sqrt(sum((X%*% temp.vec)^2)/(n)^2)*loading.norm
          #print(temp.sd)
        }else{
          mu=mu*resol;
          opt.sol=lastv;
          lamstop=1;
          tryno=tryno-1
        }
      }
    }
    tryno = tryno + 1;
  }
  direction<-(-1)/2*(opt.sol[-1]+opt.sol[1]*loading/loading.norm)
  step<-tryno-1
  returnList <- list("proj"=direction,
                     "step"=step)
  return(returnList)
}

#' Mapping the index of lower triangular matrix to its vectorized version
#'
#' @param L the 1st/2nd dimension of the matrix
#' @param l the row index of the matrix
#' @param k the column index of the matrix
#'
#' @return the index of vector after mapping
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


