library(data.table)
library(Matrix)
library(lme4)
setwd("/path/dir")

#raw_exp data
rawexp_rep1<-fread('shared_exp1XX.txt', data.table=F)
rawexp_rep2<-fread('shared_exp2XX.txt', data.table=F)

#meta files
meta_shared1<-fread("meta_shared1XX.txt", data.table=F)
meta_shared2<-fread("meta_shared2XX.txt", data.table=F)


############################################################################
#Genotypes and replicateson transcript abundances modeled as random effects#
############################################################################

blup_expr1 <- matrix(0, nrow = ncol(rawexp_rep1[,-1]), ncol = nrow(rawexp_rep1) )
blup_expr1 <- data.frame(blup_expr1)
for (i in 1:nrow(rawexp_rep1)){
  gene_exp<- t(data.frame(rawexp_rep1[i,-1], rawexp_rep2[i,-1]))
  rep<-as.factor(rep(c(1,2), times=c(ncol(rawexp_rep1[,-1]),ncol(rawexp_rep1[,-1]))))
  #rep<-rep(c(1,2), times=c(ncol(rawexp_rep1[,-1]),ncol(rawexp_rep1[,-1])))
  geno<-as.factor(c(colnames(rawexp_rep1)[-1], colnames(rawexp_rep2)[-1]))
  block<-as.factor(c(meta_shared1$Column, meta_shared2$Column))
  plot<-as.factor(c(meta_shared1$Block, meta_shared2$Block))
  date<-as.factor(c(meta_shared1$Date_Sampled, meta_shared2$Date_Sampled))
  gene_exp<-data.frame(gene_exp, geno, rep, block, plot, date)
  colnames(gene_exp)[1]<- rawexp_rep1[,1][i]
  
  mmodel <- lmer(gene_exp[,1] ~ (1|geno) + (1|rep) , data=gene_exp, REML=F)
  
  global_mean <- coef(summary(mmodel))
  blup_expr1[,i]<-ranef(mmodel)$geno  + fixef(mmodel)
}

blup_expr1 <- as.data.frame(t(blup_expr1))
colnames(blup_expr1)<-colnames(rawexp_rep1)[-1]
blup_expr1 <- data.frame(rawexp_rep1[,1], blup_expr1)
colnames(blup_expr1)[1]<-'gene_id'
head(blup_expr1)[,1:4]
write.table(blup_expr1, "blup_XXfil_date_norandom.txt", quote=F, sep='\t', row.names = F )
