opt.weight<-function(Gamma,delta=0,report.reward=TRUE){
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

compute.reward <- function(B, Sigma.Q, beta){
  temp1 = as.vector(t(B)%*%Sigma.Q%*%beta)
  temp2 = as.numeric(t(beta)%*%Sigma.Q%*%beta)
  reward = min(temp1 - temp2)
  reward
}

compute.reward2 <- function(X, y, beta){
  mean(y^2 - (y - X%*%beta)^2)
}
