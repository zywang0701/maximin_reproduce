library(MASS)
library(glmnet)
library(CVXR)
library(ggplot2)
library(reshape2)
library(ggpubr)
source('../source/RDsource.R')
source("../../source/LFsource.R")
source("../../source/utils.R")
###################################################
######### Step1: pre-process original data ########
###################################################

## read data
geno.data<- read.table("../data/TransData.txt", header = TRUE)
ind.data<- read.table("../data/ReduceInd.txt", header = FALSE)
ind.vector<-as.matrix(ind.data)
geno_red<-geno.data[,ind.vector] ## geno_red is the covariates dataset we study
pheno.data<-read.table("../data/BYxRM_PhenoData.txt",header=TRUE)
pheno<-as.matrix(pheno.data)  ## pheno is the outcome dataset we study

## remove snps with absolute correlation above 0.5
cor_thres = 0.85
cor_mat = cor(geno_red)
indx_keep = c(3611, 3638, 652) # these three signals based on our previous analysis
for(j in 1:ncol(cor_mat)){
  if((!(j%in%indx_keep))& (max(abs(cor_mat[j, indx_keep])) < cor_thres)){
    indx_keep = c(indx_keep, j)
  }
}
indx_keep = sort(indx_keep)
snp_data = geno_red[,indx_keep] # 1008*513
snp_data = scale(snp_data, center=T, scale = T)
pheno = scale(pheno, center = T, scale = T)

## we are interested in four medias ##
# media 11: Ethanol; media 19: Lactose; media 36: 5-Fluorouracil; media 38: Xylose
pheno.set = c(11,19,36,38)
pheno.names = colnames(pheno)[pheno.set]
L = length(pheno.set)

## obtain initial estimators for regression coefficients
## also apply LF to obtain bias-corrected ones
Coefs.est = matrix(0, nrow=ncol(snp_data), ncol=L)
LF.outs = rep(list(NA), L)
for(l in 1:L){
  na.set = as.vector(is.na(pheno[,pheno.set[l]]))
  X = as.matrix(snp_data[!na.set, ])
  y = as.vector(pheno[!na.set, pheno.set[l]])
  # initial
  out = cv.glmnet(X, y, intercept=F, standardize=T)
  Coefs.est[, l] = as.vector(coef(out, s="lambda.min"))[-1]
  # LF
  loading.mat = diag(ncol(X))
  out = LF(X, y, loading.mat, intercept=F ,intercept.loading = F)
  LF.outs[[l]] = out
}

##############################################
######### Step2: Preliminary Analysis ########
##############################################

## covariance matrix for each media
Cov.debias.list = rep(list(NA), L)
for(l in 1:L){
  na.set = as.vector(is.na(pheno[,pheno.set[l]]))
  X = as.matrix(snp_data[!na.set, ])
  y = as.vector(pheno[!na.set, pheno.set[l]])
  Sigma.X = t(X)%*%X/nrow(X)
  pred = X%*%Coefs.est[,l]
  sigma.sq = mean((y-pred)^2)
  proj.mat = LF.outs[[l]]$proj.mat
  Cov.debias.mat.list[[l]] = sigma.sq * t(proj.mat)%*%Sigma.X%*%proj.mat / nrow(X)
}

## hetero test by boostrap
pvalue.mat = matrix(NA, nrow=L, ncol=L)
for(l1 in 1:(L-1)){
  for(l2 in (l1+1):L){
    Cov.debias = Cov.debias.list[[l1]] + Cov.debias.list[[l2]]
    Z.gen = mvrnorm(n=1000, mu=rep(0, ncol(Cov.debias)), Sigma=Cov.debias)
    T.gen = apply(Z.gen, MARGIN=1, FUN=function(Z) max(abs(Z)))
    b.l1 = LF.outs[[l1]]$est.debias.vec
    b.l2 = LF.outs[[l2]]$est.debias.vec
    T.obs = max(abs(b.l1 - b.l2))
    pvalue.mat[l1,l2] = (1+sum(T.gen >= T.obs))/(1+1000)
  }
}
for(i in 1:L){
  for(j in 1:i){
    if(i!=j) pvalue.mat[i,j] = pvalue.mat[j,i]
  }
}
## the above pvalue.mat reproduce the Table-3 in the paper

#########################################
######### Step3: Maximin Effects ########
#########################################
## In this section, we are going to reproduce the figure 6 in the paper

###### top-subfigure #######
## the snps index we're going to plot
indx.set = c(c(420, 443, 437),
             c(423, 245, 424),
             c(364, 229),
             c(6, 177),
             c(63, 84))
n.indx = length(indx.set)
est.mat = se.mat = matrix(NA, nrow=n.indx, L)
for(l in 1:L){
  est.mat[,l] = LF.outs[[l]]$est.debias.vec[indx.set]
  se.mat[,l] = LF.outs[[l]]$se.vec[indx.set]
}
point.est = est.mat
CI.min.est = est.mat - 1.96*se.mat
CI.max.est = est.mat + 1.96*se.mat
df.matrix.1 = matrix(NA, nrow=L*n.indx, ncol=5)
colnames(df.matrix.1) <- c("group","index","center","lower","upper")
for(l in 1:L){
  for(i.coef in 1:n.indx){
    row.index = (l-1)*n.indx+i.coef
    df.matrix.1[row.index, 1] = pheno.names[l]
    df.matrix.1[row.index, 2] = indx.set[i.coef]
    df.matrix.1[row.index, 3] = point.est[i.coef, l]
    df.matrix.1[row.index, 4] = CI.min.est[i.coef, l]
    df.matrix.1[row.index, 5] = CI.max.est[i.coef, l]
  }
}
df1<-data.frame(as.factor(df.matrix.1[,1]),as.numeric(df.matrix.1[,2]),as.numeric(df.matrix.1[,3]),as.numeric(df.matrix.1[,4]),as.numeric(df.matrix.1[,5]))
colnames(df1)<-c("group","index","center","lower","upper")
df1$group = factor(df1$group, levels = pheno.names)
df1$index = factor(df1$index, levels = indx.set)
pd <- position_dodge(0.6)
title = "CIs for the regression coefficients for different growth media"
p1 <- ggplot(df1, aes(x=factor(index), y=center, colour=group, group=group))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.5, position = pd,linetype="solid")+
  geom_point(position=pd, shape=17, size=2) +
  theme(plot.title = element_text(hjust = 0.45), legend.position = "bottom")+
  labs(title=title, x="SNP Index", y="CI", color="Growth media")+
  geom_hline(yintercept = 0, linetype="dotted")

##### bot-subfigure #####
Xlist = Ylist = rep(list(NA), L)
for(l in 1:L){
  na.set = as.vector(is.na(pheno[,pheno.set[l]]))
  X = as.matrix(snp_data[!na.set, ])
  y = as.vector(pheno[!na.set, pheno.set[l]])
  Xlist[[l]] = X
  Ylist[[l]] = y
}
## apply Maximin to obtain how each model is going to aggregate
mm = MaximinRD(Xlist, Ylist, loading.idx=420, LFest=LF.outs, covariate.shift=F, Coef.est=Coefs.est)
## ridge penalty
delta.set = c(0, 0.2, 0.5)

summary.point = summary.CI.lo = summary.CI.up = matrix(NA, nrow=513, ncol=length(delta.set))
for(i.delta in 1:length(delta.set)){
  delta = delta.set[i.delta]
  weight = opt.weight(mm$Gamma.prop, delta=delta)$weight
  weight.mat = gen.prepare(mm, delta=delta)$gen.weight.mat
  for(coef in 1:513){
    print(coef)
    Point.vec = rep(0, L)
    SE.vec = rep(0, L)
    for(l in 1:L){
      Point.vec[l] = LF.outs[[l]]$est.debias.vec[coef]
      SE.vec[l] = LF.outs[[l]]$se.vec[coef]
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

indx.set = c(c(420, 443, 437),
              c(423, 245, 424),
              c(364, 229),
              c(6, 177),
              c(63, 84))
point.est = summary.result$summary.point[indx.set,]
CI.min.est = summary.result$summary.CI.lo[indx.set,]
CI.max.est = summary.result$summary.CI.up[indx.set,]
coef.index.set = indx.plot
df.matrix.2 = matrix(NA, nrow=length(delta.set)*length(indx.set), ncol=5)
colnames(df.matrix.2) <- c("delta","index","center","lower","upper")
for(i.delta in 1:length(delta.set)){
  delta = delta.set[i.delta]
  for(i.coef in 1:length(indx.set)){
    row.index = (i.delta-1)*length(indx.set)+i.coef
    df.matrix.2[row.index, 1] = delta
    df.matrix.2[row.index, 2] = indx.set[i.coef]
    df.matrix.2[row.index, 3] = point.est[i.coef, i.delta]
    df.matrix.2[row.index, 4] = CI.min.est[i.coef, i.delta]
    df.matrix.2[row.index, 5] = CI.max.est[i.coef, i.delta]
  }
}

df2<-data.frame(as.factor(df.matrix.2[,1]),df.matrix.2[,2],df.matrix.2[,3],df.matrix.2[,4],df.matrix.2[,5])
colnames(df2)<-c("delta","index","center","lower","upper")
df2$index  =factor(df2$index, levels = indx.set)
pd <- position_dodge(0.5)
title<-"CIs for the maximin effect"
p2<-ggplot(df2, aes(x=factor(index), y=center, colour=delta, group=delta)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.5,position=pd) +
  geom_point(position=pd, size=2)+theme(plot.title = element_text(hjust = 0.45), legend.position = "bottom")+
  labs(title=title,x="SNP Index",y="CI",color="delta")+ 
  geom_hline(yintercept=0,linetype="dotted")

##########################################
######### Step4: Significant Snps ########
##########################################
## apart from 11, 19, 36, 38; we further consider medias:
## media 18, 29, 31, 37, 39, 42, 45 with the names
## [1] "Lactate" "SDS" "Trehalose" "x6.Azauracil" "YNB" "YPD" "YPD.4C" 
medias = c(c(11, 19, 36, 38),
           c(18, 29, 31, 37, 39, 42, 45))

## for the additional medias we need to conduct LF as well
Coefs.est.add = matrix(0, nrow=ncol(snp_data), ncol=length(medias))
LF.outs.add = rep(list(NA), length(medias))
for(l in 1:length(medias)){
  na.set = as.vector(is.na(pheno[,medias[l]]))
  X = as.matrix(snp_data[!na.set, ])
  y = as.vector(pheno[!na.set, medias[l]])
  # initial
  out = cv.glmnet(X, y, intercept=F, standardize=T)
  Coefs.est.add[, l] = as.vector(coef(out, s="lambda.min"))[-1]
  # LF
  loading.mat = diag(ncol(X))
  out = LF(X, y, loading.mat, intercept=F ,intercept.loading = F)
  LF.outs.add[[l]] = out
}

## snps with FDR below 0.1
snp.union.list = rep(list(NA), length(medias))
for(l in 1:length(medias)){
  pvalue.vec = pnorm(abs(LF.outs.add[[l]]$est.debias.vec / LF.outs.add[[l]]$ se.vec), lower.tail = F)*2
  padjust.vec = p.adjust(pvalue.vec, method='BH')
  snp.union.list[[l]] = which(padjust.vec <= 0.1)
}

## these snps reproduce Table 4 in the paper

##########################################
######### Step5: Generalizability ########
##########################################

## find snps that are maximin siginificant under each penalty delta
delta.set = c(0, 0.2, 0.5)
snp.mm.list = rep(list(NA), length(delta.set))
pheno.set = c(11, 19, 36, 38)
snp.idx = 1:513
for(i.delta in length(delta.set)){
  delta = delta.set[i.delta]
  weight.mat = gen.prepare(mm, delta=delta)$gen.weight.mat
  pval.vec = rep(0, length(snp.idx))
  for(coef in 1:513){
    Point.vec = SE.vec= rep(0, L)
    for(l in 1:L){
      Point.vec[l] = LF.outs[[l]]$est.debias.vec[coef]
      SE.vec[l] = LF.outs[[l]]$se.vec[coef]
    }
    gen.size = nrow(weight.mat)
    gen.est = matrix(NA, nrow=gen.size, ncol=2)
    for(g in 1:gen.size){
      gen.weight.vector = weight.mat[g,]
      gen.point = sum(Point.vec * gen.weight.vector)
      gen.se = sqrt(sum(gen.weight.vector^2 * SE.vec^2))
      gen.est[g, 1] = gen.point
      gen.est[g, 2] = gen.se
    }
    out.pval = pvalue.search(gen.est)
    pval.vec[i.coef] = out.pval$alpha_try
  }
  padj.vec = p.adjust(pval.vec, method='BH')
  snp.mm = which(padj.vec < 0.1)
  snp.mm.list[[i.delta]] = snp.mm
}

## all snps that at least 1 significant
snp.union = sort(Reduce(union, snp.union.list))
## count how many medias each snp is significant at
snp_in_media.union = rep(list(NA), length(snp.union))
names(snp_in_media.union) = snp.union
for(i.snp in 1:length(snp.union)){
  snp = snp.union[i.snp]
  snp_in_media = c()
  for(i.media in 1:length(media.step2)){
    if(snp %in% snp.union.list[[i.media]]) snp_in_media = c(snp_in_media, media.step2[i.media])
  }
  snp_in_media.union[[i.snp]] = snp_in_media
}
num_media.union = sapply(snp_in_media.union, length)
names.order = names(num_media.union)

ps = rep(list(NA), 3)
for(i.delta in 1:3){
  delta = delta.set[i.delta]
  snp.mm = snp.mm.list[[i.delta]]
  snp.mm = snp.mm[snp.mm%in%snp.union]
  
  ## count how many medias each mm-significant snp is significant at
  snp_in_media.mm = rep(list(NA, length(snp.mm)))
  names(snp_in_media.mm) = snp.mm
  for(i.snp in 1:length(snp.mm)){
    snp = snp.mm[i.snp]
    snp_in_media = c()
    for(i.media in 1:length(medias)){
      if(snp %in% snp.union.list[[i.media]]) snp_in_media = c(snp_in_media, medias[i.media])
    }
    snp_in_media.mm[[i.snp]] = snp_in_media
  }
  num_media.mm = sapply(snp_in_media.mm, length)
  
  ## count how may medias other snp is significant at
  snp.others = sort(setdiff(snp.union, snp.mm))
  snp_in_media.others = rep(list(NA, length(snp.others)))
  names(snp_in_media.others) = snp.others
  for(i.snp in 1:length(snp.others)){
    snp = snp.others[i.snp]
    snp_in_media = c()
    for(i.media in 1:length(medias)){
      if(snp %in% snp.union.list[[i.media]]) snp_in_media = c(snp_in_media, medias[i.media])
    }
    snp_in_media.others[[i.snp]] = snp_in_media
  }
  num_media.others = sapply(snp_in_media.others, length)
  
  ####################### Draw Figures #####################
  df = data.frame(matrix(NA, nrow=length(snp.union), ncol=4))
  df[,1] = c(names(num_media.mm), names(num_media.others))
  df[,2] = c(num_media.mm, num_media.others)
  df[,3] = c(rep("Maximin Significant", length(num_media.mm)), 
             rep("Not Maximin Significant", length(num_media.others)))
  for(i in 1:nrow(df)){
    if(df[i,2] >=3){
      df[i,4] = df[i,1]
    }else{
      df[i,4] = ""
    }
  }
  colnames(df) = c("snp", "times","method", "label")
  
  df$snp = factor(df$snp, levels=names.order)
  title = paste0('delta = ', delta)
  if(i.delta <=2) legend.pos = 'none' else legend.pos = 'bottom'
  ps[[i.delta]] = ggplot(df, aes(snp, times, label=label))+
    geom_point(aes(color=method), size=2.2)+
    geom_text(aes(color=method), size=3., nudge_y = 1.2, nudge_x = 0, show.legend = F)+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          legend.position = legend.pos, plot.title = element_text(hjust = 0.5),
          legend.title=element_blank())+
    labs(title=title, x="SNP Index", y="num of media")+
    scale_y_continuous(breaks = c(1,3,5,7,9,11))
}
ggarrange(ps[[1]], ps[[2]], ps[[3]], nrow=3, heights = c(1,1,1), common.legend = TRUE,legend = 'bottom')
