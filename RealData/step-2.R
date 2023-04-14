###################
##### Step-2 ######
###################

#### Hetero Test pvalue ####
cor.setting = 2 #{1,2}
pheno.set = c(11,19,36,38)
dataList = readRDS(paste0("dataList",cor.setting,".rds"))
Cov.debias.mat.list = rep(list(NA), length(pheno.set))
Coefs.est = readRDS('LassoEst-cor2-lamMin.rds')
for(i.pheno.l in 1:length(pheno.set)){
  pheno.l = pheno.set[i.pheno.l]
  filename = paste0("Rdata/debias_result0726_cor",cor.setting,"pheno",pheno.l,".RData")
  load(filename)
  X.cent = scale(X, center = T, scale = F)
  # X.cent = cbind(1, X.cent)
  Sigma.X = t(X.cent)%*%X.cent/nrow(X.cent)
  #beta.init = as.vector(coef(cv.glmnet(X, y, intercept=F, standardize=T), s="lambda.min"))[-1]
  beta.init = Coefs.est[,pheno.l]
  pred = as.vector(X.cent%*%beta.init)
  sigma.sq = mean((y-pred)^2)
  Cov.debias.mat.list[[i.pheno.l]] = sigma.sq*t(out$proj.mat)%*%Sigma.X%*%(out$proj.mat)/nrow(X)
}

pvalue.mat = matrix(NA, nrow=length(pheno.set), ncol=length(pheno.set))
for(l1 in 1:(length(pheno.set)-1)){
  for(l2 in (l1+1):length(pheno.set)){
    pheno.l1 = pheno.set[l1]
    pheno.l2 = pheno.set[l2]
    Cov.debias.l12 = Cov.debias.mat.list[[l1]] + Cov.debias.mat.list[[l2]]
    Z.gen.l12 = MASS::mvrnorm(n=1000, mu=rep(0, ncol(Cov.debias.l12)), Sigma = Cov.debias.l12)
    T.l12 = apply(Z.gen.l12, MARGIN = 1, FUN=function(Z) max(abs(Z)))
    filename1 = paste0("Rdata/debias_result0726_cor",cor.setting,"pheno",pheno.l1,".RData")
    load(filename1)
    b.l1 = out$est.debias.vec
    filename2 = paste0("Rdata/debias_result0726_cor",cor.setting,"pheno",pheno.l2,".RData")
    load(filename2)
    b.l2 = out$est.debias.vec
    T.l12_star = max(abs(b.l1 - b.l2))
    pvalue.mat[l1,l2] = (1 + sum(T.l12 >= T.l12_star))/(1+1000)
  }
}
round(pvalue.mat,4)
for(i in 1:nrow(pvalue.mat)){
  for(j in 1:i){
    if(i!=j) pvalue.mat[i,j] = pvalue.mat[j,i]
  }
}
df = data.frame(round(pvalue.mat,4))
colnames(df) = rownames(df) = pheno.names
library(kableExtra)
kbl(df, 'latex', digits=4, align='c')%>%
  kable_styling(latex_options = c("scale_down"))