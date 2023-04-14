############################
##### Data Preparation #####
############################
geno.data<- read.table("original_data/TransData.txt", header = TRUE)
ind.data<- read.table("original_data/ReduceInd.txt", header = FALSE)
ind.vector<-as.matrix(ind.data)
geno_red<-geno.data[,ind.vector] ## geno_red is the covariates dataset we study
pheno.data<-read.table("original_data/BYxRM_PhenoData.txt",header=TRUE)
pheno<-as.matrix(pheno.data)  ## pheno is the outcome dataset we study

### filtering-1 ###
freq = colSums(geno_red==1) / (2 * nrow(geno_red))
maf = ifelse(freq > 0.5, 1 - freq, freq)
snp_keep = which(maf > 0.1) # 4410 all kept

### filtering-2 ###
cor.setting = 2
if(cor.setting==1){
  cor_thres = 0.8
}else if(cor.setting==2){
  cor_thres = 0.85
}else if(cor.setting==3){
  cor_thres = 0.9
}
cor_mat = cor(geno_red)
indx_keep = c(3611, 3638, 652)

for(j in 1:ncol(cor_mat)){
  if((!(j%in%indx_keep))& (max(abs(cor_mat[j, indx_keep])) < cor_thres)){
    indx_keep = c(indx_keep, j)
  }
}
indx_keep = sort(indx_keep)
snp_data = geno_red[,indx_keep] # 1008*513
snp_data = scale(snp_data, center=T, scale = T)
pheno = scale(pheno, center = T, scale = T)
data.list = list("snp_data"=snp_data,
                 "pheno"=pheno,
                 "indx"=indx_keep)
saveRDS(data.list, paste0('dataList',cor.setting,'.rds'))

### Obtain Initial Estimators ###
library(glmnet)
Coefs.est = matrix(0, nrow=ncol(snp_data), ncol=ncol(pheno))
for(l in 1:ncol(pheno)){
    na.set = as.vector(is.na(pheno[,l]))
    X = as.matrix(snp_data[!na.set, ])
    y = as.vector(pheno[!na.set, l])
    out = cv.glmnet(X, y, intercept=F, standardize=T)
    Coefs.est[, l] = as.vector(coef(out, s="lambda.min"))[-1]
}
saveRDS(Coefs.est, paste0('LassoEst-cor',cor.setting,'-lamMin.rds'))

### Apply SIHR LF function to every phenos ###
library(SIHR)
for(pheno.l in 1:ncol(pheno)){
  snp_data = data.list$snp_data
  pheno = data.list$pheno
  na.set = is.na(pheno[,pheno.l])
  X = snp_data[!na.set, ]
  y = pheno[!na.set, pheno.l]
  loading.mat = diag(ncol(X))
  out = LF(X, y, loading.mat, intercept=F ,intercept.loading = F)
  
  rm(data.list, snp_data, pheno)
  filename = paste0("debias_result0726_cor",cor.setting,"pheno",pheno.l,".RData")
  save.image(filename)
}
