###############################################################
# 09_Prepare_eQTLnet_inputs
#
# Description:
#   1. Load inferred GRN adjacency matrix.
#   2. Load gene expression data.
#   3. Convert adjacency matrix into edge list format.
#   4. Obtain:
#      - regulator expression matrix (X)
#      - regulated gene expression list (Y)
#      - regulator-target weights
#   5. Save objects for downstream analyses.
#
# Inputs (data/):
#   - complete_mx.RData
#   - AM_Final_Maize_GS.txt
#   - log_blups_XXnodaterandom.txt
#   - TF_list_update.txt
#
# Outputs (results/):
#   - pvalue_completenopv.RData
#   - regulators_completenopv.RData
#   - regulated_completenopv.RData
#
###############################################################

############################
# Libraries
############################

library(data.table)
library(tidyr)

############################
# Load GRN matrix
############################

load( "complete_mx.RData")

############################
# Load gold-standard network
############################

gs <- fread("AM_Final_Maize_GS.txt", data.table = FALSE)

colnames(gs) <- c(
  "from",
  "to",
  "weight"
)

############################
# Load expression matrix
############################

gene_exp <- fread("log_blups_XXnodaterandom.txt", data.table = FALSE)

rownames(gene_exp) <- gene_exp[, 1]

gene_exp <- gene_exp[, -1]

gene_exp <- t(gene_exp)


############################
# Load TFs
############################

tfs <- scan( "TF_list_update.txt", what = "character")

tfs <- unique(tfs)

############################
# Convert matrix to edge list
############################

complete_df <- as.data.frame(
  t(complete_mx)
)

complete_df <- cbind(
  from = rownames(complete_df),
  complete_df
)

complete_df <- complete_df %>%
  pivot_longer(
    cols = 2:ncol(.),
    names_to = "to",
    values_to = "weight"
  )

complete_df <- as.data.frame(complete_df)

############################
# Convert weights
#
# Original logic:
# Any value different from 1
# becomes 1
############################

complete_df$weight[complete_df$weight != 1] <- 1

############################
# Split by regulated gene
############################

by_regulated <- split(complete_df,f = complete_df$to)

############################
# Y matrix
# Expression of regulated genes
############################

regulated <- names(by_regulated)

exp_regulated <- as.data.frame(gene_exp[,colnames(gene_exp) %in% regulated,drop = FALSE])

exp_regulated <- as.list(exp_regulated)

############################
# X matrix
# Expression of all regulators
############################

exp_regulators <- gene_exp[,colnames(gene_exp) %in% colnames(complete_mx),  drop = FALSE]


############################
# Regulatory weights
############################

pvalue <- vector("list",length(by_regulated))

for (i in seq_along(by_regulated)) {
  
  by_regulated[[i]] <-by_regulated[[i]][order(by_regulated[[i]]$from),]
  
  pvalue[[i]] <-by_regulated[[i]]$weight
  
  names(pvalue[[i]]) <- by_regulated[[i]]$from
}

names(pvalue) <- names(by_regulated)

############################
# Save outputs
############################

save(pvalue, file =  "pvalue_completenopv.RData")

save(exp_regulators, file ="regulators_completenopv.RData")

save(exp_regulated,file ="regulated_completenopv.RData")

