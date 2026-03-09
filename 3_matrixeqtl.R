library('data.table')
library("factoextra")
#install.packages('ggfortify')
library('ggfortify')
setwd("/path/dir")
library(ggplot2)

gene_exp <- fread("blup_XXfil_date_norandom.txt", data.table = F)
gene_id<-gene_exp$gene_id
genotypes_nofounders <- fread("geno_eqtl.txt" ,data.table=F)

gene_exp<-gene_exp[,-1][,which(colnames(gene_exp[,-1]) %in% colnames(genotypes_nofounders))]
gene_exp<-cbind(gene_id, gene_exp)

explev<-data.frame(gene_exp[,-1])
explev<-apply(explev, 2, as.numeric)

explevlog<-apply(explev, 2, function(x) log(x+1))
explevlog<-as.data.frame(explevlog)
explevlog<-as.data.frame(cbind(gene_exp[,1], explevlog))
colnames(explevlog)[1]<-'gene_id'
write.table(explevlog, "log_blups_XXnodaterandom.txt", quote=F, sep="\t", row.names=FALSE)


##################################
#MatrixEQTL
##################################
#install.packages('MatrixEQTL')
library('MatrixEQTL')
# Linear model 
useModel = modelLINEAR

#set paramenters
workdir= "/path/dir"
SNP_file_name = "geno_eqtl.txt"

expression_file_name = "log_blups_XXnodaterandom.txt"

top_geno <- t(geno_PCA_prcomp$x[,c(1:10)])
id<-paste0("Geno_",colnames(geno_PCA_prcomp$x)[1:10])
top_geno<-data.frame(id, top_geno)
write.table(top_geno, "PCA_top10_geno.txt",quote=F,sep="\t",row.names=FALSE)
#separate file with extra covariates
covariates_file_name = "PCA_top10_geno.txt"

# Output file name
output_file_name = paste0(workdir, "/", "eqtl_10genoPC_logblupsXXnodate.txt");

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t"; # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 2000; # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"; # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
if(length(covariates_file_name)>0){
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  #pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

