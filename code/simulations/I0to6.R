source("../../source/gendata.R")
source("../../source/MaximinInference.R")
source("../../source/LFsource.R")
source("../../source/utils.R")
library(MASS)
library(CVXR)
library(glmnet)
library(intervals)

# picked setting, setting = {0,1,2,3,4,5,6}
setting = 1
# nsim
nsim = 500

#############################Do not modify below################################
# The number of groups, L
L = 4
# sample size for each group
n = 1000
# dim
p = 500
# delta
delta = 0
# covariance matrix
cov.source = diag(p) # cov.gen(p, 0.6)
cov.target = cov.source
# generate mean for source and target data
mean.source = rep(0, p)
mean.target = rep(0, p)
# Ns
Ns = rep(n, L)

# hyperparameters for simulations
gen.size = 500

# The setting of Coefs
if(setting==0){
  seed = 0
  noise = 0
}
if(setting==1){
  seed = 42
  noise = 0.05
}else if(setting==2){
  seed = 20
  noise = 0.05
}else if(setting==3){
  seed = 36
  noise = 0.1
}else if(setting==4){
  seed = 17
  noise = 0.15
}else if(setting==5){
  seed = 12
  noise = 0.2
}else if(setting==6){
  seed = 31
  noise = 0.25
}

set.seed(seed)
Bs = matrix(0, p, L)
b = rep(0, p)
b[1:10] = seq(1:10)/20
for(l in 1:L){
  Bs[,l] = b
  Bs[1:5, l] = Bs[1:5, l] + rnorm(5, sd=noise)
}
loading = rep(0, p)
loading[1:5] = 1

trueList <- report.true(Bs, cov.target, loading, delta.truth=delta, verbose=TRUE)

set.seed(NULL)
weight.sim = matrix(NA, nrow=nsim, ncol=L)
mm.est.sim = matrix(0, nrow=nsim, ncol=p+1) # intercept=TRUE, store Maximin Effects
pick.sim = matrix(0, nrow=nsim, ncol=4) # count how many estimators picked for each threshold version
# mat.store stores the info for summary
# version1,2,3 refers to different ways to generate samples
# 1st column: point estimators
# 2nd column: length of CI - 1st version
# 3rd column: length of CI - 2nd version with alpha 0.01
# 4th column: length of CI - 2nd version with alpha 0.05
# 5th column: length of CI - 3rd version with alpha 0.01
# 6th column: length of CI - 3rd version with alpha 0.05
# 7th column: covered or not - 1st version
# 8th column: covered or not - 2nd version with alpha 0.01
# 9th column: covered or not - 2nd version with alpha 0.05
# 10th column: covered or not - 3rd version with alpha 0.01
# 11th column: covered or not - 3rd version with alpha 0.05
# 12th column: det - 1st version
# 13th column: det - 2nd version with alpha 0.01
# 14th column: det - 2nd version with alpha 0.05
# 15th column: det - 3rd version with alpha 0.01
# 16th column: det - 3rd version with alpha 0.05
mat.store = matrix(NA, nrow=nsim, ncol=16)
colnames(mat.store)<-c("point","CI0","CI1-1","CI1-5","CI2-1","CI2-5",
                       "Cover0","Cover1-1","Cover1-5","Cover2-1", "Cover2-5",
                       "det0","det1-1","det1-5","det2-1","det2-5")
spnorm.Gamma.diff.sim = matrix(NA, nrow=nsim, ncol=gen.size)
l2norm.gamma.diff.sim = replicate(length(delta.vec), matrix(NA, nrow=nsim, ncol=gen.size), simplify = FALSE)
measure.sim = rep(NA, nsim)
for(i.sim in 1:nsim){
  dataList <- gen.data.only.source(mean.source, cov.source, Bs, Ns, L)
  s1List <- mm.s1(dataList$X.source, dataList$Y.source, dataList$idx.source, NULL, loading, covariate.shift = FALSE)
  
  s3List0 <- mm.s3(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, gen.size=gen.size, threshold=0)
  pick1.001 <- pick.samples(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, set.samples=s3List0$gen.samples, threshold=1, alpha=0.01)
  pick1.005 <- pick.samples(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, set.samples=s3List0$gen.samples, threshold=1, alpha=0.05)
  pick2.001 <- pick.samples(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, set.samples=s3List0$gen.samples, threshold=2, alpha=0.01)
  pick2.005 <- pick.samples(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, set.samples=s3List0$gen.samples, threshold=2, alpha=0.05)
  pick.sim[i.sim, 1] = pick1.001$pick.size
  pick.sim[i.sim, 2] = pick1.005$pick.size
  pick.sim[i.sim, 3] = pick2.001$pick.size
  pick.sim[i.sim, 4] = pick2.005$pick.size
  
  s2List <- mm.s2(s1List$Gamma.prop, s1List$Coef.est, s1List$Point.vec, delta=delta)
  gen.Gamma.array <- mm.s3.measure(s1List$L, s3List0$gen.samples, s3List0$gen.size)$gen.Gamma.array
  for(i.gen in 1:gen.size){
    Gamma.diff = gen.Gamma.array[i.gen,,] - s1List$Gamma.prop
    spnorm.Gamma.diff.sim[i.sim, i.gen] = sqrt(max((eigen(Gamma.diff)$values)^2))
  }
  s4List0 <- mm.s4.measure(s1List$L, s3List0$gen.samples, s3List0$gen.size, delta)
  for(i.gen in 1:gen.size){
    gamma.diff = s4List0$gen.weight.mat[i.gen,] - s2List$weight.prop
    l2norm.gamma.diff.sim[[index.delta]][i.sim, i.gen] = sqrt(sum(gamma.diff^2))
  }
  measure.sim[i.sim] = mean(l2norm.gamma.diff.sim[[index.delta]][i.sim,]/spnorm.Gamma.diff.sim[i.sim,])
  
  s4List0 <- mm.s4(s1List$Point.vec, s1List$SE.vec, s1List$L, s3List0$gen.samples, s3List0$gen.size, delta=delta)
  s4List1.001 <- pick.CI(s4List0$gen.est, pick1.001$index.thres)
  s4List1.005 <- pick.CI(s4List0$gen.est, pick1.005$index.thres)
  s4List2.001 <- pick.CI(s4List0$gen.est, pick2.001$index.thres)
  s4List2.005 <- pick.CI(s4List0$gen.est, pick2.005$index.thres)

  weight.sim[i.sim, ] = s2List$weight.prop
  mm.est.sim[i.sim, ] = s2List$mm.est
  mat.store[i.sim, 1] = s2List$point
  mat.store[i.sim, 2] = s4List0$CI.length
  mat.store[i.sim, 3] = s4List1.001$CI.length
  mat.store[i.sim, 4] = s4List1.005$CI.length
  mat.store[i.sim, 5] = s4List2.001$CI.length
  mat.store[i.sim, 6] = s4List2.005$CI.length
  
  cover0 = sum((trueList$true.val < s4List0$CI.union[, 2])*
                 (trueList$true.val > s4List0$CI.union[, 1]))
  cover1.001 = sum((trueList$true.val < s4List1.001$CI.union[, 2])*
                     (trueList$true.val > s4List1.001$CI.union[, 1]))
  cover1.005 = sum((trueList$true.val < s4List1.005$CI.union[, 2])*
                     (trueList$true.val > s4List1.005$CI.union[, 1]))
  cover2.001 = sum((trueList$true.val < s4List2.001$CI.union[, 2])*
                     (trueList$true.val > s4List2.001$CI.union[, 1]))
  cover2.005 = sum((trueList$true.val < s4List2.005$CI.union[, 2])*
                     (trueList$true.val > s4List2.005$CI.union[, 1]))
  
  mat.store[i.sim, 7] = cover0
  mat.store[i.sim, 8] = cover1.001
  mat.store[i.sim, 9] = cover1.005
  mat.store[i.sim, 10] = cover2.001
  mat.store[i.sim, 11] = cover2.005
  
  det0 = !sum((0 < s4List0$CI.union[, 2])*(0 > s4List0$CI.union[, 1]))
  det1.001 = !sum((0 < s4List1.001$CI.union[, 2])*(0 > s4List1.001$CI.union[, 1]))
  det1.005 = !sum((0 < s4List1.005$CI.union[, 2])*(0 > s4List1.005$CI.union[, 1]))
  det2.001 = !sum((0 < s4List2.001$CI.union[, 2])*(0 > s4List2.001$CI.union[, 1]))
  det2.005 = !sum((0 < s4List2.005$CI.union[, 2])*(0 > s4List2.005$CI.union[, 1]))
  
  mat.store[i.sim, 12] = det0
  mat.store[i.sim, 13] = det1.001
  mat.store[i.sim, 14] = det1.005
  mat.store[i.sim, 15] = det2.001
  mat.store[i.sim, 16] = det2.005
}

filename <- paste("Setting-(I",setting,")-L",L,"-n",n,"-p",p,".RData", sep="")
save.image(filename)
