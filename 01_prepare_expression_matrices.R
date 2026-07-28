# ============================================================
# 01_prepare_expression_matrices.R
#
# Create replicate-specific expression matrices and retain only
# lines present in the SPET genotyping dataset.
# ============================================================

library(data.table)

# ----------------------------
# Load data
# ----------------------------

expression <- read.csv( "Field2023_TPM.csv",  stringsAsFactors = FALSE)

metadata <- read.csv("Field_2023_TPM_Metadata.csv", stringsAsFactors = FALSE)

spet <- fread("spet_geno.txt", data.table = FALSE)

# ----------------------------
# Split metadata by replicate
# ----------------------------

meta_rep1 <- subset(metadata, Replicate == 1)
meta_rep2 <- subset(metadata, Replicate == 2)

# retain lines shared with genotype data

meta_rep1 <- meta_rep1[meta_rep1$Line %in% colnames(spet)[-1],]

meta_rep2 <- meta_rep2[meta_rep2$Line %in% colnames(spet)[-1],]

# ----------------------------
# Create expression matrices
# ----------------------------

expr_rep1 <- expression[,colnames(expression) %in% meta_rep1$V]

expr_rep2 <- expression[, colnames(expression) %in% meta_rep2$V]

rownames(expr_rep1) <- expression[,1]
rownames(expr_rep2) <- expression[,1]

# rename columns using SSA line names

colnames(expr_rep1) <- meta_rep1$Line[match(colnames(expr_rep1), meta_rep1$V)]

colnames(expr_rep2) <- meta_rep2$Line[match(colnames(expr_rep2), meta_rep2$V)]

expr_rep1 <- expr_rep1[, order(colnames(expr_rep1))]
expr_rep2 <- expr_rep2[, order(colnames(expr_rep2))]

# ----------------------------
# Save outputs
# ----------------------------

write.table(expr_rep1,"exp_tot_r1.txt", sep = "\t", quote = FALSE,row.names = FALSE)

write.table(expr_rep2,"exp_tot_r2.txt", sep = "\t", quote = FALSE, row.names = FALSE)
