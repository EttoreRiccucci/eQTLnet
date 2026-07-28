###############################################################
# 07_eQTL_Gene_Set
#
# Description:
#   1. Load TFs, gene annotations and gold-standard GRN.
#   2. Define comparative and discovery gene sets.
#   3. Annotate all genes with their corresponding set.
#   4. Assign gene-set labels to filtered eQTL results.
#   5. Export annotated metadata and eQTL tables.
#
# Inputs (data/):
#   - TF_list_update.txt
#   - Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1_ab.1_xref_gene_IDs.txt
#   - AM_Final_Maize_GS.txt
#   - eqtl_clean.txt
#
# Outputs (results/):
#   - comp_noTF.txt
#   - comp_TFregulators.txt
#   - comp_TFnoregulators.txt
#   - disc_TF.txt
#   - disc_noTF.txt
#   - v5_meta_annotation.txt
#   - eqtl_cleanwithset.txt
###############################################################

############################
# Libraries
############################

library(data.table)
library(dplyr)

############################
# Load data
############################

tfs <- scan("TF_list_update.txt", what = "character")

tfs <- unique(tfs)

meta_annotation <- fread("Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1_ab.1_xref_gene_IDs.txt", data.table = FALSE)

gs <- fread("AM_Final_Maize_GS.txt", data.table = FALSE)

colnames(gs) <- c("from", "to", "weight")

############################
# Comparative sets
############################

gs_from <- intersect(unique(gs$from), tfs)

gs_to <- unique(gs$to)

regulated <- unique(gs_to)

regulated_tf <- intersect(regulated, tfs)

regulated_notf <- setdiff(regulated, regulated_tf)

regulators <- gs_from

tf_noregulators <- setdiff( regulated_tf, regulators)

comp_noTF <- regulated_notf
comp_TFregulators <- regulators
comp_TFnoregulators <- tf_noregulators

write(comp_noTF, "comp_noTF.txt")

write(comp_TFregulators, "comp_TFregulators.txt")


write(comp_TFnoregulators, "comp_TFnoregulators.txt")

############################
# Discovery sets
############################

disc_genes <- setdiff(meta_annotation$new_gene_model_ID, regulated)

disc_TF <- intersect(disc_genes, tfs)

disc_noTF <- setdiff(disc_genes,disc_TF)

write(disc_TF,  "disc_TF.txt")


write(disc_noTF, "disc_noTF.txt")


############################
# Annotate genes
############################

meta_annotation$gene_set <- "?"

meta_annotation$gene_set[
  meta_annotation$new_gene_model_ID %in% comp_noTF
] <- "comp_noTF"

meta_annotation$gene_set[
  meta_annotation$new_gene_model_ID %in% comp_TFregulators
] <- "comp_TFregulators"

meta_annotation$gene_set[
  meta_annotation$new_gene_model_ID %in% comp_TFnoregulators
] <- "comp_TFnoregulators"

meta_annotation$gene_set[meta_annotation$new_gene_model_ID %in% disc_noTF] <- "disc_noTF"

meta_annotation$gene_set[meta_annotation$new_gene_model_ID %in% disc_TF] <- "disc_TF"

write.table(meta_annotation,"v5_meta_annotation.txt",sep = "\t", quote = FALSE, row.names = FALSE)

cat("Unassigned genes:", sum(meta_annotation$gene_set == "?"), "\n")

############################
# Annotate eQTL results
############################

eqtl <- fread("eqtl_clean.txt", data.table = FALSE)

meta_eqtl <- meta_annotation[meta_annotation$new_gene_model_ID %in% eqtl$eqtl_gene,]

eqtl <- eqtl[eqtl$eqtl_gene %in% meta_eqtl$new_gene_model_ID,]

gene_set <- data.frame(eqtl_gene = meta_eqtl$new_gene_model_ID, gene_set = meta_eqtl$gene_set)

eqtl <- merge(eqtl, gene_set, by = "eqtl_gene")

write.table(eqtl,"eqtl_cleanwithset.txt", sep = "\t", quote = FALSE, row.names = FALSE)

############################
# Summary statistics
############################

cat("\nUnique eQTL genes:",length(unique(eqtl$eqtl_gene)),"\n")

eqtl_split <- split(eqtl,f = eqtl$gene_set)

cat("\nGene-set summary\n")
cat("----------------\n")

for (set_name in names(eqtl_split)) {
  
  cat(
    set_name,
    ":",
    nrow(eqtl_split[[set_name]]),
    "eQTLs |",
    length(unique(eqtl_split[[set_name]]$eqtl_gene)),
    "genes\n"
  )
}

############################
# Save split object
############################

save(eqtl_split, results_dir,"eqtl_by_gene_set.RData")

