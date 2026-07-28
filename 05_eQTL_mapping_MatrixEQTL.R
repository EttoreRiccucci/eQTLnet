# ==========================================================
# 05_eQTL_mapping_MatrixEQTL
# ==========================================================

library(data.table)
library(MatrixEQTL)
library(factoextra)
library(ggfortify)
library(ggplot2)

# ----------------------------------------------------------
# Load data
# ----------------------------------------------------------

gene_exp <- fread("blup_XXfil_date_norandom.txt", data.table = FALSE)

genotypes_nofounders <- fread("geno_eqtl.txt", data.table = FALSE)

# ----------------------------------------------------------
# Keep only common genotypes
# ----------------------------------------------------------

gene_id <- gene_exp$gene_id

common_cols <- intersect(colnames(gene_exp)[-1], colnames(genotypes_nofounders))

gene_exp <- gene_exp[, c("gene_id", common_cols)]

# ----------------------------------------------------------
# Log-transform expression values
# ----------------------------------------------------------

explev <- gene_exp[, -1]
explev <- data.frame(lapply(explev, as.numeric))

explevlog <- log(explev + 1)

explevlog <- cbind(gene_id = gene_exp$gene_id, explevlog)

write.table(explevlog, "log_blups_XXnodaterandom.txt", quote = FALSE, sep = "\t",row.names = FALSE)

# ----------------------------------------------------------
# PCA on genotype matrix
# ----------------------------------------------------------

geno_numeric <- genotypes_nofounders

geno_numeric <- geno_numeric[, -1]

geno_numeric <- as.data.frame(lapply(geno_numeric, as.numeric))

geno_PCA_prcomp <- prcomp(t(geno_numeric), center = TRUE, scale. = TRUE)

top_geno <- t(geno_PCA_prcomp$x[, 1:10])

id <- paste0("Geno_PC", 1:10)

top_geno <- data.frame(id,top_geno, check.names = FALSE)

covariate_output <-  "PCA_top10_geno.txt"

write.table(top_geno, covariate_output, quote = FALSE, sep = "\t", row.names = FALSE)

# ----------------------------------------------------------
# MatrixEQTL setup
# ----------------------------------------------------------

useModel <- modelLINEAR
errorCovariance <- numeric()

# ----------------------------------------------------------
# Load SNPs
# ----------------------------------------------------------

snps <- SlicedData$new()

snps$fileDelimiter <- "\t"
snps$fileOmitCharacters <- "NA"
snps$fileSkipRows <- 1
snps$fileSkipColumns <- 1
snps$fileSliceSize <- 2000

snps$LoadFile(genotype_input)

# ----------------------------------------------------------
# Load expression
# ----------------------------------------------------------

gene <- SlicedData$new()

gene$fileDelimiter <- "\t"
gene$fileOmitCharacters <- "NA"
gene$fileSkipRows <- 1
gene$fileSkipColumns <- 1
gene$fileSliceSize <- 2000

gene$LoadFile(expression_output)

# ----------------------------------------------------------
# Load covariates
# ----------------------------------------------------------

cvrt <- SlicedData$new()

cvrt$fileDelimiter <- "\t"
cvrt$fileOmitCharacters <- "NA"
cvrt$fileSkipRows <- 1
cvrt$fileSkipColumns <- 1

cvrt$LoadFile(covariate_output)

# ----------------------------------------------------------
# Run MatrixEQTL
# ----------------------------------------------------------

eqtl_output <-  "eqtl_10genoPC_logblupsXXnodate.txt"

me <- Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = eqtl_output,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
)

cat("\nAnalysis completed.\n")
cat("Results written to:\n")
cat(eqtl_output, "\n")