# ============================================================
# 03_remove_replicate_specific_zeros.R
#
# Remove genes showing zero expression in the same line in both
# biological replicates.
# ============================================================

library(data.table)

expr_rep1 <- fread("exp_XX_r1.txt", data.table = FALSE)

expr_rep2 <- fread("exp_XX_r2.txt", data.table = FALSE)

# ----------------------------
# Identify genes to remove
# ----------------------------

expr1 <- as.matrix(expr_rep1[, -1])
expr2 <- as.matrix(expr_rep2[, -1])

genes_to_remove <- which(rowSums((expr1 == 0) & (expr2 == 0)) > 0)

cat(length(genes_to_remove), "genes removed\n")

expr_filtered_r1 <- expr_rep1[-genes_to_remove, ]
expr_filtered_r2 <- expr_rep2[-genes_to_remove, ]

# ----------------------------
# Save outputs
# ----------------------------

write.table(expr_filtered_r1, "exp_XXfilt_r1.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(expr_filtered_r2, "exp_XXfilt_r2.txt", sep = "\t", quote = FALSE, row.names = FALSE)

