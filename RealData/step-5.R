source("Maximin_src/Maximin_RD.R")
source("Maximin_src/pvalue_search.R")
library(CVXR)
library(MASS)
library(intervals)
################# maximin pvalue search ################
delta.set = c(0, 0.2, 0.5)
shift.set = c(FALSE, TRUE)
pheno.set = c(11, 19, 36, 38); L = length(pheno.set)
snp.idx = 1:513

## save generated weights at first 
for(delta in delta.set){
  for(shift in shift.set){
    print(delta)
    print(shift)
    if(shift==FALSE){
      mm = readRDS('result_0811/mmobject-shiftF.rds')
      weight.mat = gen.prepare(mm, delta=delta)$gen.weight.mat
      saveRDS(weight.mat, paste0('result_0811/GenWeightMat-delta',delta,'-shiftF.rds'))
    }else if(shift==TRUE){
      mm = readRDS('result_0811/mmobject-shiftT.rds')
      weight.mat = gen.prepare(mm, delta=delta)$gen.weight.mat
      saveRDS(weight.mat, paste0('result_0811/GenWeightMat-delta',delta,'-shiftT.rds'))
    }
  }
}

## pvalue search 
for(delta in delta.set){
  for(shift in shift.set){
    print(delta)
    print(shift)
    if(shift==FALSE){
      weight.mat = readRDS(paste0('result_0811/GenWeightMat-delta',delta,'-shiftF.rds'))
    }else{
      weight.mat = readRDS(paste0('result_0811/GenWeightMat-delta',delta,'-shiftT.rds'))
    }
    
    pval.vec = rep(0, length(snp.idx))
    for(i.coef in 1:length(snp.idx)){
      coef = snp.idx[i.coef]
      print(coef)
      Point.vec = rep(0, L)
      SE.vec = rep(0, L)
      for(l in 1:L){
        pheno.l = pheno.set[l]
        load(paste0("Rdata/debias_result0726_cor",2,"pheno",pheno.l,".RData"))
        Point.vec[l] = out$est.debias.vec[coef]
        SE.vec[l] = out$se.vec[coef]
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
    p.vec.list = list(pval.vec = pval.vec,
                      padj.vec = padj.vec)
    if(shift) TorF = 'T' else TorF = 'F'
    saveRDS(p.vec.list, paste0('result_0811/mmPvals-delta',delta,'shift',TorF,'.rds'))
  }
}

############## Replicability in medias ##############
load('result_0811/image-part2-1.RData')
snp.union = sort(Reduce(union, snp.union.list))
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
# o.union = order(num_media.union, decreasing = T)
# num_media.union = num_media.union[o.union]
names.order = names(num_media.union)

#### Without Covariate Shift ####
ps = rep(list(NA), 3)
for(i.delta in 1:3){
  delta = c(0, 0.2, 0.5)[i.delta]
  shiftTorF = 'F'
  p.vec.list = readRDS(paste0('result_0811/mmPvals-delta',delta,'shift',shiftTorF,'.rds'))
  padj.vec = p.vec.list$padj.vec
  snp.mm = which(padj.vec < 0.1)
  snp.mm = snp.mm[snp.mm%in%snp.union]
  
  snp_in_media.mm = rep(list(NA), length(snp.mm))
  names(snp_in_media.mm) = snp.mm
  for(i.snp in 1:length(snp.mm)){
    snp = snp.mm[i.snp]
    snp_in_media = c()
    for(i.media in 1:length(media.step2)){
      if(snp %in% snp.union.list[[i.media]]) snp_in_media = c(snp_in_media, media.step2[i.media])
    }
    snp_in_media.mm[[i.snp]] = snp_in_media
  }
  
  snp.others = sort(setdiff(snp.union, snp.mm))
  snp_in_media.others = rep(list(NA), length(snp.others))
  names(snp_in_media.others) = snp.others
  for(i.snp in 1:length(snp.others)){
    snp = snp.others[i.snp]
    snp_in_media = c()
    for(i.media in 1:length(media.step2)){
      if(snp %in% snp.union.list[[i.media]]) snp_in_media = c(snp_in_media, media.step2[i.media])
    }
    snp_in_media.others[[i.snp]] = snp_in_media
  }
  
  num_media.mm = sapply(snp_in_media.mm, length)
  num_media.others = sapply(snp_in_media.others, length)
  
  # ## re-order num_media
  # o.mm = order(num_media.mm, decreasing = T)
  # o.others = order(num_media.others, decreasing = T)
  # num_media.mm = num_media.mm[o.mm]
  # num_media.others = num_media.others[o.others]
  
  ####################### Draw Figures #####################
  library(ggplot2)
  df = data.frame(matrix(NA, nrow=length(snp.union), ncol=4))
  df[,1] = c(names(num_media.mm), names(num_media.others))
  df[,2] = c(num_media.mm, num_media.others)
  df[,3] = c(rep("Maximin Significant", length(num_media.mm)), 
             rep("Not Maximin Significant", length(num_media.others)))
  # for(i in 1:nrow(df)){
  #   if(df[i,2] >= 3){
  #     df[i,4] = df[i,1]
  #   }else if(df[i,3]=="Maximin"){
  #     df[i,4] = df[i,1]
  #   }else{
  #     df[i,4] = ""
  #   }
  # }
  for(i in 1:nrow(df)){
    if(df[i,2] >=3){
      df[i,4] = df[i,1]
    }else{
      df[i,4] = ""
    }
  }
  colnames(df) = c("snp", "times","method", "label")
  
  df$snp = factor(df$snp, levels=names.order)
  #title = paste0('Times of Passing FDR of each SNP across 11 Media - delta',delta,'-Shift',shiftTorF)
  title = paste0('delta = ', delta)
  if(i.delta <=2) legend.pos = 'none' else legend.pos = 'bottom'
  ps[[i.delta]] = ggplot(df, aes(snp, times, label=label))+
    geom_point(aes(color=method), size=2.2)+
    #geom_text(aes(color=method), position = position_dodge(width = 1))+
    geom_text(aes(color=method), size=3., nudge_y = 1.2, nudge_x = 0, show.legend = F)+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
          legend.position = legend.pos, plot.title = element_text(hjust = 0.5),
          legend.title=element_blank())+
    labs(title=title, x="SNP Index", y="num of media")+
    scale_y_continuous(breaks = c(1,3,5,7,9,11))
}

library(ggpubr)
ggarrange(ps[[1]], ps[[2]], ps[[3]], nrow=3, heights = c(1,1,1), common.legend = TRUE,legend = 'bottom')