### Table: each media, which SNP is fdr significant
media.step2 = c(c(11, 19, 36, 38),
                c(18, 29, 31, 37, 39, 42, 45))
pheno.set = media.step2
snp.union.list = rep(list(NA), length(pheno.set))
for(i.pheno.l in 1:length(pheno.set)){
  pheno.l = pheno.set[i.pheno.l]
  filename = paste0('RData/debias_result0726_cor',2,'pheno',pheno.l,'.RData')
  load(filename)
  
  pvalue.vec = pnorm(abs(out$est.debias.vec / out$se.vec), lower.tail = F)*2
  padjust.vec = p.adjust(pvalue.vec, method='BH')
  snp.union.list[[i.pheno.l]] = which(padjust.vec <= 0.1)
}
snp.union.list

dataList = readRDS('dataList2.rds')
pheno.names = colnames(dataList$pheno)
library(kableExtra)
df = data.frame(matrix(NA, nrow=length(media.step2), ncol=3))
colnames(df) = c('media','count','SNP index')
for(i in 1:length(pheno.set)){
  df[i,1] = pheno.names[pheno.set[i]]
  df[i,2] = length(snp.union.list[[i]])
  df[i,3] = paste(snp.union.list[[i]], collapse = ',')
}
kbl(df, 'latex', align='l')%>%
  kable_styling(latex_options = c("scale_down"))

save.image('result_0811/image-part2-1.RData')