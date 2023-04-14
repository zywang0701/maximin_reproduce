########## Figure 6 ###########

########## Top subfigure #########
library(ggplot2)
library(reshape2)
cor.setting = 2

pheno.set = c(11, 19, 36, 38); L = length(pheno.set)
dataList = readRDS('dataList2.rds')
pheno.names = colnames(dataList$pheno)[pheno.set]
pheno.names[3] = '5-Fluorouracil'
indx.plot = c(c(420, 443, 437),
              c(423, 245, 424),
              c(364, 229),
              c(6, 177),
              c(63, 84))
indx.set = indx.plot
n.indx = length(indx.set)
est.mat = se.mat = matrix(NA, nrow=length(indx.set), ncol=length(pheno.set))
colnames(est.mat) = pheno.set
for(i.pheno.l in 1:length(pheno.set)){
  pheno.l = pheno.set[i.pheno.l]
  filename = paste0("Rdata/debias_result0726_cor",cor.setting,"pheno",pheno.l,".RData")
  load(filename)
  est.mat[, i.pheno.l] = out$est.debias.vec[indx.set]
  se.mat[, i.pheno.l] = out$se.vec[indx.set]
}

point.est = est.mat
CI.min.est = est.mat - 1.96*se.mat
CI.max.est = est.mat + 1.96*se.mat
df.matrix.1 = matrix(NA, nrow=L*n.indx, ncol=5)
colnames(df.matrix.1) <- c("group","index","center","lower","upper")
i.coef.set = seq(1,n.indx)
for(i.l in 1:L){
  l = pheno.set[i.l]
  for(i.i.coef in 1:n.indx){
    row.index = (i.l-1)*n.indx+i.i.coef
    df.matrix.1[row.index, 1] = pheno.names[i.l]
    df.matrix.1[row.index, 2] = indx.set[i.i.coef]
    df.matrix.1[row.index, 3] = point.est[i.i.coef, i.l]
    df.matrix.1[row.index, 4] = CI.min.est[i.i.coef, i.l]
    df.matrix.1[row.index, 5] = CI.max.est[i.i.coef, i.l]
  }
}

df1<-data.frame(as.factor(df.matrix.1[,1]),as.numeric(df.matrix.1[,2]),as.numeric(df.matrix.1[,3]),as.numeric(df.matrix.1[,4]),as.numeric(df.matrix.1[,5]))
colnames(df1)<-c("group","index","center","lower","upper")
df1$group = factor(df1$group, levels = pheno.names)
df1$index  =factor(df1$index, levels = indx.plot)
pd <- position_dodge(0.6)
title = "CIs for the regression coefficients for different growth media"
p1 <- ggplot(df1, aes(x=factor(index), y=center, colour=group, group=group))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.5, position = pd,linetype="solid")+
  geom_point(position=pd, shape=17, size=2) +
  theme(plot.title = element_text(hjust = 0.45), legend.position = "bottom")+
  labs(title=title, x="SNP Index", y="CI", color="Growth media")+
  geom_hline(yintercept = 0, linetype="dotted")
p1 

#################### Bottom subfigure ####################
source("Maximin_src/Maximin_RD.R")
source("Maximin_src/pvalue_search.R")
library(CVXR)
library(MASS)
library(ggplot2)
library(reshape2)

pheno.set = c(11, 19, 36, 38); L = length(pheno.set)
cor.setting = 2
data.list = readRDS(paste0('dataList',cor.setting,'.rds'))
Coefs.est = readRDS("LassoEst-cor2-lamMin.rds")
L = length(pheno.set)
Xlist = Ylist = LFest = rep(list(NA), L)
for(l in 1:L){
  X = data.list$snp_data
  Y = data.list$pheno[,pheno.set[l]]
  na.set = is.na(Y)
  Xlist[[l]] = X[!na.set,]
  Ylist[[l]] = Y[!na.set]
  pheno.l = pheno.set[l]
  filename = paste0("Rdata/debias_result0726_cor",cor.setting,"pheno",pheno.l,".RData")
  load(filename)
  LFest[[l]] = out
}

mm.rd.1 = MaximinRD(Xlist, Ylist, loading.idx = 420, LFest=LFest,
                    covariate.shift = F, Coef.est = Coefs.est[,pheno.set])
saveRDS(mm.rd.1, 'result_0811/mmobject-shiftF.rds')
cor.setting = 2
delta.set = c(0, 0.2, 0.5)
pheno.set = c(11, 19, 36, 38); L = length(pheno.set)
snp.idx = 1:513
summary.point = summary.CI.lo = summary.CI.up = matrix(NA, nrow=513, ncol=length(delta.set))
mm = readRDS('result_0811/mmobject-shiftF.rds') # we save shift-F already in part1-2
for(i.delta in 1:length(delta.set)){
  delta = delta.set[i.delta]
  weight = opt.weight(mm$Gamma.prop, delta = delta)$weight
  weight.mat = gen.prepare(mm, delta=delta)$gen.weight.mat
  for(coef in snp.idx){
    print(coef)
    Point.vec = rep(0, L)
    SE.vec = rep(0, L)
    for(l in 1:L){
      pheno.l = pheno.set[l]
      load(paste0("Rdata/debias_result0726_cor",cor.setting,"pheno",pheno.l,".RData"))
      Point.vec[l] = out$est.debias.vec[coef]
      SE.vec[l] = out$se.vec[coef]
    }
    gen.size = nrow(weight.mat)
    gen.output = matrix(NA, nrow=gen.size, ncol=2)
    for(g in 1:gen.size){
      gen.weight.vector = weight.mat[g,]
      gen.est = sum(Point.vec * gen.weight.vector)
      gen.se = sqrt(sum(gen.weight.vector^2 * SE.vec^2))
      gen.output[g, 1] = gen.est
      gen.output[g, 2] = gen.se
    }
    summary.point[coef, i.delta] = sum(Point.vec * weight)
    alpha = 0.05
    CI.ori = cbind(gen.output[,1] - qnorm(1-alpha/2)*gen.output[,2],
                   gen.output[,1] + qnorm(1-alpha/2)*gen.output[,2])
    CI = na.omit(CI.ori)
    summary.CI.lo[coef, i.delta] = min(CI[,1])
    summary.CI.up[coef, i.delta] = max(CI[,2])
  }
}
summary.result = list(summary.point = summary.point,
                      summary.CI.lo = summary.CI.lo,
                      summary.CI.up = summary.CI.up)
saveRDS(summary.result, 'result_0811/summary_result_513SNP_shiftF.rds')
indx.plot = c(c(420, 443, 437),
              c(423, 245, 424),
              c(364, 229),
              c(6, 177),
              c(63, 84))
point.est = summary.result$summary.point[indx.plot,]
CI.min.est = summary.result$summary.CI.lo[indx.plot,]
CI.max.est = summary.result$summary.CI.up[indx.plot,]
delta.set = c(0, 0.2, 0.5)
coef.index.set = indx.plot
df.matrix.2 = matrix(NA, nrow=length(delta.set)*length(coef.index.set), ncol=5)
colnames(df.matrix.2) <- c("delta","index","center","lower","upper")
for(i.delta in 1:length(delta.set)){
  delta = delta.set[i.delta]
  for(i.coef.index in 1:length(coef.index.set)){
    coef.index = coef.index.set[i.coef.index]
    row.index = (i.delta-1)*length(coef.index.set)+i.coef.index
    df.matrix.2[row.index, 1] = delta
    df.matrix.2[row.index, 2] = coef.index
    df.matrix.2[row.index, 3] = point.est[i.coef.index, i.delta]
    df.matrix.2[row.index, 4] = CI.min.est[i.coef.index, i.delta]
    df.matrix.2[row.index, 5] = CI.max.est[i.coef.index, i.delta]
  }
}

df2<-data.frame(as.factor(df.matrix.2[,1]),df.matrix.2[,2],df.matrix.2[,3],df.matrix.2[,4],df.matrix.2[,5])
colnames(df2)<-c("delta","index","center","lower","upper")
df2$index  =factor(df2$index, levels = indx.plot)
library(ggplot2) 
library("gridExtra")
pd <- position_dodge(0.5)
title<-"CIs for the maximin effect"
p2<-ggplot(df2, aes(x=factor(index), y=center, colour=delta, group=delta)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.5,position=pd) +
  geom_point(position=pd, size=2)+theme(plot.title = element_text(hjust = 0.45), legend.position = "bottom")+
  labs(title=title,x="SNP Index",y="CI",color="delta")+ 
  geom_hline(yintercept=0,linetype="dotted")
p2
library(ggpubr)
ggarrange(p1, p2, nrow=2, heights=c(1,1))