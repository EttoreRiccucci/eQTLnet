# ============================================================
# 04_estimate_expression_blups.R
#
# Purpose:
# Estimate genotype BLUPs for each gene using a linear mixed
# model with genotype and replicate included as random effects.
#
# Inputs:
# - shared_exp1XX.txt
# - shared_exp2XX.txt
# - meta_rep1XX.txt
# - meta_rep2XX.txt
#
# Output:
# - blup_XXfil_date_norandom.txt
# ============================================================

library(data.table)
library(lme4)

# ------------------------------------------------------------
# Load expression matrices
# ------------------------------------------------------------

expr_rep1 <- fread("shared_exp1XX.txt", data.table = FALSE)

expr_rep2 <- fread("shared_exp2XX.txt", data.table = FALSE)

# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

meta_rep1 <- fread("meta_shared1XX.txt", data.table = FALSE)

meta_rep2 <- fread("meta_shared2XX.txt", data.table = FALSE)

# ============================================================
# Estimate genotype BLUPs for each gene
# Genotype and replicate are modeled as random effects
# ============================================================

# ============================================================
# Estimate genotype BLUPs for each gene
# Genotype and replicate are modeled as random effects
# ============================================================

blup_matrix <- matrix(0, nrow = ncol(expr_rep1[,-1]), ncol = nrow(expr_rep1) )
blup_matrix <- data.frame(blup_matrix)
for (i in 1:nrow(expr_rep1)){
  gene_exp<- t(data.frame(expr_rep1[i,-1], expr_rep2[i,-1]))
  rep<-as.factor(rep(c(1,2), times=c(ncol(expr_rep1[,-1]),ncol(expr_rep1[,-1]))))
  #rep<-rep(c(1,2), times=c(ncol(expr_rep1[,-1]),ncol(expr_rep1[,-1])))
  geno<-as.factor(c(colnames(expr_rep1)[-1], colnames(expr_rep2)[-1]))
  block<-as.factor(c(meta_rep1$Column, meta_rep2$Column))
  plot<-as.factor(c(meta_rep1$Block, meta_rep2$Block))
  date<-as.factor(c(meta_rep1$Date_Sampled, meta_rep2$Date_Sampled))
  gene_exp<-data.frame(gene_exp, geno, rep, block, plot, date)
  colnames(gene_exp)[1]<- expr_rep1[,1][i]
  
  mmodel <- lmer(gene_exp[,1] ~ (1|geno) + (1|rep) , data=gene_exp, REML=F)
  
  global_mean <- coef(summary(mmodel))
  blup_matrix[,i]<-ranef(mmodel)$geno  + fixef(mmodel)
}

# ------------------------------------------------------------
# Format output
# ------------------------------------------------------------

blup_matrix<-t(blup_matrix)
colnames(blup_matrix) <- colnames(expr_rep1)[-1]

blup_matrix <- data.frame(gene_id = expr_rep1[, 1], blup_matrix)

# ------------------------------------------------------------
# Save BLUP matrix
# ------------------------------------------------------------

write.table(blup_matrix, "blup_XXfil_date_norandom.txt", sep = "\t", quote = FALSE, row.names = FALSE)


