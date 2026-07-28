# ============================================================
# 02_filter_low_expression_genes.R
#
# Retain genes with <=20% zero values across lines in both
# replicates.
# ============================================================

library(data.table)

expr_rep1 <- fread("exp_tot_r1.txt", data.table = FALSE)

expr_rep2 <- fread("exp_tot_r2.txt", data.table = FALSE)

# ----------------------------
# Count zeros
# ----------------------------

zero_count_r1 <- rowSums(expr_rep1 == 0)
zero_count_r2 <- rowSums(expr_rep2 == 0)

genes_r1 <- data.frame(gene_id = rownames(expr_rep1), zero_count = zero_count_r1)

genes_r2 <- data.frame(gene_id = rownames(expr_rep2), zero_count = zero_count_r2)

threshold_r1 <- ncol(expr_rep1) * 0.20
threshold_r2 <- ncol(expr_rep2) * 0.20

genes_pass_r1 <- genes_r1$gene_id[genes_r1$zero_count <= threshold_r1]

genes_pass_r2 <- genes_r2$gene_id[genes_r2$zero_count <= threshold_r2]

genes_nozero_r1 <- genes_r1$gene_id[genes_r1$zero_count == 0]

genes_nozero_r2 <- genes_r2$gene_id[genes_r2$zero_count == 0]

genes_20pct <- intersect(genes_pass_r1, genes_pass_r2)

genes_nozero <- intersect(genes_nozero_r1, genes_nozero_r2)

expr_rep1 <- data.frame(gene_id = rownames(expr_rep1), expr_rep1)

expr_rep2 <- data.frame(gene_id = rownames(expr_rep2), expr_rep2)

expr_20pct_r1 <- subset(expr_rep1, gene_id %in% genes_20pct)

expr_20pct_r2 <- subset(expr_rep2, gene_id %in% genes_20pct)

expr_nozero_r1 <- subset(expr_rep1, gene_id %in% genes_nozero)

expr_nozero_r2 <- subset(expr_rep2, gene_id %in% genes_nozero)

write.table(expr_20pct_r1, "exp_XX_r1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(expr_20pct_r2, "exp_XX_r2.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(expr_nozero_r1, "exp_noz_r1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(expr_nozero_r2, "exp_noz_r2.txt", sep = "\t", quote = FALSE, row.names = FALSE)


