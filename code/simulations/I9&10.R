source("../../source/gendata.R")
source("../../source/MaximinInference.R")
source("../../source/LFsource.R")
source("../../source/utils.R")

library(MASS)
library(CVXR)
library(glmnet)
library(intervals)

# setting = {9,10}
setting = 9
# sample size, n = {100,200,300,500}
n = 500
# nsim
nsim = 500
###################################################

delta = 0

if(setting==9){
  # The number of groups
  L = 2
  # dimension
  p = 500
  # Coefs
  Bs = matrix(0, p, L)
  Bs[1:30,] = 1 
  Bs[1,1] = 1000^(0.1)
  Bs[1,2] = -1000^(-0.5)
  # loading
  loading = rep(0, p)
  loading[1] = 1
}
if(setting==10){
  # The number of groups
  L = 2
  # dimension
  p = 500
  # Coefs
  Bs = matrix(0, p, L)
  Bs[1:10,1] = seq(1,10)/20
  Bs[1:10,2] = -seq(1,10)/20
  # loading
  loading = rep(0, p)
  loading[1:5] = seq(1,5)/5
}

# covariance matrix
cov.source = diag(p) # cov.gen(p, 0.6)
cov.target = cov.source
# generate mean for source and target data
mean.source = rep(0, p)
mean.target = rep(0, p)
# Ns
Ns = rep(n, L)
# truth
trueList.List <- list()
for(index.delta in 1:length(delta.vec)){
  delta = delta.vec[index.delta]
  trueList <- report.true.reward(Bs, cov.target, loading, delta.truth=delta, verbose=TRUE)
  trueList.List[[index.delta]] = trueList
}

# gen.size
gen.size=500
# seed
set.seed(NULL)
min.ev.sim = rep(0, nsim)
weight.sim = matrix(0, nrow=nsim, ncol=L) # store weights
mm.est.sim = matrix(0, nrow=nsim, ncol=p+1) # intercept=TRUE, store Maximin Effects
count.sim = matrix(0, nrow=nsim, ncol=4) # count how many estimators picked for each threshold version
# mat.store stores the info for summary
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
  count.sim[i.sim, 1] = pick1.001$pick.size
  count.sim[i.sim, 2] = pick1.005$pick.size
  count.sim[i.sim, 3] = pick2.001$pick.size
  count.sim[i.sim, 4] = pick2.005$pick.size
  
  # start_time = Sys.time()
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
  # end_time = Sys.time()
  # if(i.sim%%2==1) print(paste("The sample process takes", round(end_time-start_time, 4)))
  
  weight.sim[i.sim, ] = s2List$weight.prop
  mm.est.sim[i.sim, ] = s2List$mm.est
  mat.store[i.sim, 1] = s2List$point
  mat.store[i.sim, 2] = s4List0$CI.length
  mat.store[i.sim, 3] = s4List1.001$CI.length
  mat.store[i.sim, 4] = s4List1.005$CI.length
  mat.store[i.sim, 5] = s4List2.001$CI.length
  mat.store[i.sim, 6] = s4List2.005$CI.length
  
  trueList = trueList.List[[index.delta]]
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

filename <- paste("Setting-(I",settting,")-L",L,"-n",n,"-p",p,".RData", sep="")
save.image(filename)
