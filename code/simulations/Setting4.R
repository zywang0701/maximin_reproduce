source("../../source/gendata.R")
source("../../source/MaximinInference.R")
source("../../source/LFsource.R")
source("../../source/utils.R")
library(MASS)
library(CVXR)
library(glmnet)
library(intervals)

# sample size for each group, n={100,200,300,500}
n = 500
# number of groups, L={2,5,10}
L = 2
# nsim
nsim = 500

#############################Do not modify below################################
# delta
delta.vec = c(0, 0.1, 0.5, 1, 2)
# dimension, p
p = 500
# covariance matrix
cov.source = cov.gen(p, 0.6)
cov.target = cov.source
for(i in 1:p) cov.target[i, i] = 1.5
for(i in 1:5){
  for(j in 1:5){
    if(i!=j) cov.target[i, j] = 0.9
  }
}
for(i in 499:500){
  for(j in 499:500){
    if(i!=j) cov.target[i, j] = 0.9
  }
}
# generate mean for source and target data
mean.source = rep(0, p)
mean.target = rep(0, p)
# coefs
Bs = matrix(0, p, L)
Bs[1:10, 1] = seq(1,10)/40
Bs[498:500, 1] = c(0.5, -0.5, -0.5)
for(l in 2:L){
  Bs[seq(10*l+1, 10*(l+1)), l] = Bs[1:10, 1]
  Bs[498:500,l] = 0.5^(l-1)*Bs[498:500,1]
}
# loading
loading = rep(0, p)
loading[498:500] = 1

# Ns
Ns = rep(n, L)
# number of samples for target
N.target = 2000
# truth
trueList.List <- list()
for(index.delta in 1:length(delta.vec)){
  delta = delta.vec[index.delta]
  trueList <- report.true.reward(Bs, cov.target, loading, delta.truth=delta, verbose=TRUE)
  trueList.List[[index.delta]] = trueList
}
# true.reward = rep(0, 5)
# for(i in 1:5) true.reward[i] = trueList.List[[i]]$reward
# true.reward / true.reward[1]
# 
# # weights
# do.call(cbind, lapply(trueList.List, FUN=function(out) out$weight.true))
# # Gamma.true
# Gamma.true = t(Bs)%*%cov.target%*%Bs
# Gamma.true
# min.ev = min(eigen(Gamma.true)$values)
# min.ev

# simulation 500 times
gen.size = 500

### store for one element in delta.vec
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
### store for multiple deltas in delta.vec
weight.sim.List = replicate(length(delta.vec), weight.sim, simplify=FALSE)
mm.est.sim.List = replicate(length(delta.vec), mm.est.sim, simplify=FALSE)
mat.store.List = replicate(length(delta.vec), mat.store, simplify=FALSE)

for(i.sim in 1:nsim){
  if(i.sim%%5==1) print(paste("Running Simulation --", i.sim))

  dataList <- gen.data.with.target(mean.source, cov.source, mean.target, cov.target, N.target, Bs, Ns, L)
  s1List <- mm.s1(dataList$X.source, dataList$Y.source, dataList$idx.source, dataList$X.target, loading, covariate.shift = TRUE)

  s3List0 <- mm.s3(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, gen.size=gen.size, threshold=0)
  pick1.001 <- pick.samples(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, set.samples=s3List0$gen.samples, threshold=1, alpha=0.01)
  pick1.005 <- pick.samples(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, set.samples=s3List0$gen.samples, threshold=1, alpha=0.05)
  pick2.001 <- pick.samples(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, set.samples=s3List0$gen.samples, threshold=2, alpha=0.01)
  pick2.005 <- pick.samples(s1List$gen.mu, s1List$gen.Cov, s1List$gen.dim, set.samples=s3List0$gen.samples, threshold=2, alpha=0.05)
  count.sim[i.sim, 1] = pick1.001$pick.size
  count.sim[i.sim, 2] = pick1.005$pick.size
  count.sim[i.sim, 3] = pick2.001$pick.size
  count.sim[i.sim, 4] = pick2.005$pick.size
  
  # newly added at 10.8
  min.ev.sim[i.sim] = min(eigen(s1List$Gamma.prop)$values)
  
  for(index.delta in 1:length(delta.vec)){
    delta = delta.vec[index.delta]
    # start_time = Sys.time()
    s2List <- mm.s2(s1List$Gamma.prop, s1List$Coef.est, s1List$Point.vec, delta=delta)
    s4List0 <- mm.s4(s1List$Point.vec, s1List$SE.vec, s1List$L, s3List0$gen.samples, s3List0$gen.size, delta=delta)
    s4List1.001 <- pick.CI(s4List0$gen.est, pick1.001$index.thres)
    s4List1.005 <- pick.CI(s4List0$gen.est, pick1.005$index.thres)
    s4List2.001 <- pick.CI(s4List0$gen.est, pick2.001$index.thres)
    s4List2.005 <- pick.CI(s4List0$gen.est, pick2.005$index.thres)
    # end_time = Sys.time()
    # if(i.sim%%2==1) print(paste("The sample process takes", round(end_time-start_time, 4)))
    weight.sim.List[[index.delta]][i.sim, ] = s2List$weight.prop
    mm.est.sim.List[[index.delta]][i.sim, ] = s2List$mm.est
    mat.store.List[[index.delta]][i.sim, 1] = s2List$point
    mat.store.List[[index.delta]][i.sim, 2] = s4List0$CI.length
    mat.store.List[[index.delta]][i.sim, 3] = s4List1.001$CI.length
    mat.store.List[[index.delta]][i.sim, 4] = s4List1.005$CI.length
    mat.store.List[[index.delta]][i.sim, 5] = s4List2.001$CI.length
    mat.store.List[[index.delta]][i.sim, 6] = s4List2.005$CI.length
    
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
    
    mat.store.List[[index.delta]][i.sim, 7] = cover0
    mat.store.List[[index.delta]][i.sim, 8] = cover1.001
    mat.store.List[[index.delta]][i.sim, 9] = cover1.005
    mat.store.List[[index.delta]][i.sim, 10] = cover2.001
    mat.store.List[[index.delta]][i.sim, 11] = cover2.005
    
    det0 = !sum((0 < s4List0$CI.union[, 2])*(0 > s4List0$CI.union[, 1]))
    det1.001 = !sum((0 < s4List1.001$CI.union[, 2])*(0 > s4List1.001$CI.union[, 1]))
    det1.005 = !sum((0 < s4List1.005$CI.union[, 2])*(0 > s4List1.005$CI.union[, 1]))
    det2.001 = !sum((0 < s4List2.001$CI.union[, 2])*(0 > s4List2.001$CI.union[, 1]))
    det2.005 = !sum((0 < s4List2.005$CI.union[, 2])*(0 > s4List2.005$CI.union[, 1]))
    
    mat.store.List[[index.delta]][i.sim, 12] = det0
    mat.store.List[[index.delta]][i.sim, 13] = det1.001
    mat.store.List[[index.delta]][i.sim, 14] = det1.005
    mat.store.List[[index.delta]][i.sim, 15] = det2.001
    mat.store.List[[index.delta]][i.sim, 16] = det2.005
  }
}

filename <- paste("Setting4","-n",n,"-L",L,"-delta.vec","-p",p,".RData", sep="")
save.image(filename)

######## For summarization ###############
# All materials are saved in
# - count.sim
# - truthList.List
# - weight.sim.List
# - mm.est.sim.List
# - mat.store.List
