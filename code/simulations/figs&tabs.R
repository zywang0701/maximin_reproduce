############################ table 1 ############################ 
## get one row (setting) in the table 1
tab1 <- function(filename){
  load(filename)
  ## column I(delta)
  measure = mean(measure.sim)
  ## column coverage-normality and length-normality
  se = sd(mat.store[,1])
  CI.oracle = cbind(mat.store[,1] - 1.96*se, mat.store[,1] + 1.96*se)
  true.val = trueList$true.val
  cov.norm = mean((true.val < CI.oracle[,2])*(true.val > CI.oracle[,1]))
  len.norm = 2*1.96*se
  ## column coverage-proposed and length-proposed
  cov.prop = mean(mat.store[,8])
  len.prop = mean(mat.store[,3])
  ## length ratio
  ratio = len.prop/len.norm
  
  return(c(measure, cov.norm, cov.prop, len.norm, len.prop, ratio))
}

############################ figure 3 ############################ 
## get left, mid, or right two sub-figures
fig3 <- function(filename, pos){
  library(ggplot2)
  load(filename)
  
  ## p1 corresponds to top subfigure; p2 corresponds to right subfigure
  
  ## left two figures
  if(pos=='left'){
    df = data.frame(matrix(mat.store[,1], ncol=1))
    colnames(df) = c('est')
    p1=ggplot(df, aes(x=est))+
      geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
      geom_density(alpha=.2, fill="#FF6666")+
      geom_vline(xintercept = trueList$true.val, color="red", size=1)+
      geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
      labs(x="est", title="Non-regularity + Instability")+
      theme(plot.title = element_text(hjust = 0.5))+
      xlim(0.3, 1.0)
    
    df = data.frame(matrix(weight.sim[,4], ncol=1))
    colnames(df) = c("weight")
    p2=ggplot(df, aes(x=weight))+
      geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
      geom_density(alpha=.2, fill="#FF6666")+
      labs(x="4th group's weight")+
      theme(plot.title = element_text(hjust = 0.5))
  }
  ## middle two figures
  if(pos=='mid'){
    df = data.frame(matrix(mat.store[,1], ncol=1))
    colnames(df) = c("est")
    p1=ggplot(df, aes(x=est))+
      geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
      geom_density(alpha=.2, fill="#FF6666")+
      geom_vline(xintercept = trueList$true.val, color="red", size=1)+
      geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
      labs(x="est", title="Non-regularity")+
      theme(plot.title = element_text(hjust = 0.5))+
      xlim(c(-0.25, 0.45))
    
    df = data.frame(matrix(weight.sim[,1], ncol=1))
    colnames(df) = c("weight")
    p2=ggplot(df, aes(x=weight))+
      geom_histogram(aes(y=..density..), colour="black",fill="white", bins=60)+
      geom_density(alpha=.2, fill="#FF6666")+
      labs(x="1st group's weight")+
      theme(plot.title = element_text(hjust = 0.5))+
      xlim(c(-0.005, 1.0))
  }
  ## right two figures
  if(pos=='right'){
    df = data.frame(matrix(mat.store[,1], ncol=1))
    colnames(df) = c("est")
    p1=ggplot(df, aes(x=est))+
      geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
      geom_density(alpha=.2, fill="#FF6666")+
      geom_vline(xintercept = trueList$true.val, color="red", size=1)+
      geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
      labs(x="est", title="Normaltiy")+
      theme(plot.title = element_text(hjust = 0.5))+
      xlim(c(-0.35, 0.35))
    
    df = data.frame(matrix(weight.sim[,1], ncol=1))
    colnames(df) = c("weight")
    p2=ggplot(df, aes(x=weight))+
      geom_histogram(aes(y=..density..), colour="black",fill="white", bins=60)+
      geom_density(alpha=.2, fill="#FF6666")+
      labs(x="1st group's weight")+
      theme(plot.title = element_text(hjust = 0.5))+
      xlim(c(0, 1))
    p2
  }
  return(list(p1=p1,
              p2=p2)) 
}

############################ table 2 ############################
## get one specific setting's coverage error and length ratio
tab2 <- function(filename){
  load(filename)
  
  delta.vec = c(0, 0.1, 0.5, 1, 2)
  err.vec = ratio.vec = rep(NA, length(delta.vec))
  for(i.delta in 1:length(delta.vec)){
    ## length normality
    points = mat.store.List[[i.delta]][,1]
    se = sd(points)
    len.norm = 2*1.96*se
    ## coverage-proposed and length-proposed
    cov.prop = mean(mat.store.List[[i.delta]][,8])
    len.prop = mean(mat.store.List[[i.delta]][,3])
    ## length ratio
    ratio = len.prop/len.norm
    
    err.vec[i.delta] = abs(cov.prop-0.95)
    ratio.vec[i.delta] = ratio
  }
  
  return(list(err = mean(err.vec),
              ratio = mean(ratio.vec)))
}

############################ figure 4 ############################
## obtain coverage, ci length, length ratio for each n in {100,200,300,500}
fig4 <- function(filename){
  load(filename)
  
  delta.vec = c(0, 0.1, 0.5, 1, 2)
  cov.vec = len.vec = ratio.vec = rep(NA, length(delta.vec))
  for(i.delta in 1:length(delta.vec)){
    ## length normality
    points = mat.store.List[[i.delta]][,1]
    se = sd(points)
    len.norm = 2*1.96*se
    ## coverage-proposed and length-proposed
    cov.prop = mean(mat.store.List[[i.delta]][,8])
    len.prop = mean(mat.store.List[[i.delta]][,3])
    ## length ratio
    ratio = len.prop/len.norm
    
    cov.vec[i.delta] = cov.prop
    len.vec[i.delta] = len.prop
    ratio.vec[i.delta] = ratio
  }
  
  return(list(cov.vec=cov.vec,
              len.vec=len.vec,
              ratio.vec=ratio.vec))
}

############################ figure 5 ############################
## obtain coverage, ci length, length ratio for method in {"CS Known", "CS", "No CS"}
fig5 <- function(filename){
  load(filename)
  
  delta.vec = c(0, 0.1, 0.5, 1, 2)
  cov.vec = len.vec = ratio.vec = rep(NA, length(delta.vec))
  for(i.delta in 1:length(delta.vec)){
    ## length normality
    points = mat.store.List[[i.delta]][,1]
    se = sd(points)
    len.norm = 2*1.96*se
    ## coverage-proposed and length-proposed
    cov.prop = mean(mat.store.List[[i.delta]][,8])
    len.prop = mean(mat.store.List[[i.delta]][,3])
    ## length ratio
    ratio = len.prop/len.norm
    
    cov.vec[i.delta] = cov.prop
    len.vec[i.delta] = len.prop
    ratio.vec[i.delta] = ratio
  }
  
  return(list(cov.vec=cov.vec,
              len.vec=len.vec,
              ratio.vec=ratio.vec))
}