###############################################################
# 06_eQTL_SNP_Prioritization
#
# Description:
#   1. Load eQTL associations
#   2. Group nearby SNPs (150 bp window)
#   3. Remove redundant SNP clusters
#   4. Select the most significant SNP per cluster
#   5. Generate filtered eQTL table
#   6. Annotate SNPs with chromosome-specific LD decay
#
# Input files (expected in data/):
#   - half_decay_LD.csv
#   - TF_list_update.txt
#   - eqtl_10genoPC_logblupsXXnodate.txt
#
# Output files (written to results/):
#   - top_snps_ineqtl.RData
#   - eqtl_clean.txt
#   - snp_meta.txt
###############################################################

############################
# Libraries
############################

library(data.table)
library(vctrs)
library(dplyr)
library(stringr)

############################
# Load data
############################

hd_LD <- read.csv2( "half_decay_LD.csv")

meanLD <- mean(hd_LD$Mb)

eqtl <- fread( "eqtl_10genoPC_logblupsXXnodate.txt", data.table = FALSE)

tfs <- scan( "TF_list_update.txt", what = "character")

############################
# TF and non-TF eQTL genes
############################

tf_id <- tfs

tf_eqtl <- unique(eqtl$gene)[unique(eqtl$gene) %in% tf_id]

notf_eqtl <- unique(eqtl$gene)[!unique(eqtl$gene) %in% tf_id]

############################
# SNP metadata
############################

snp_name <- unique(eqtl$SNP)

write(snp_name, file.path(results_dir, "snp_name_logblupsXXnodate.txt"))

snp_split <- str_split_fixed(snp_name, ":", 2)

snp_meta_tosplit <- data.frame(
  snp_name = snp_name,
  snp_loc = as.numeric(snp_split[, 2]),
  snp_chr = as.numeric(snp_split[, 1])
)

write(unique(snp_meta_tosplit$snp_chr), file.path(results_dir, "snp_chr_logblupsXXnodate.txt"))

write(unique(snp_meta_tosplit$snp_loc), file.path(results_dir, "snp_loc_logblupsXXnodate.txt"))

############################
# SNPs associated per gene
############################

eqtl_genes <- unique(eqtl$gene)

snps <- vector("list", length(eqtl_genes))

for (i in seq_along(eqtl_genes)) {
  snps[[i]] <- eqtl[eqtl$gene %in% eqtl_genes[i],][, 1]
  names(snps)[i] <- eqtl_genes[i]
}

############################
# Group SNPs within 150 bp
############################

gene_snps <- list()

for (i in seq_along(snps)) {
  snp_meta <- snp_meta_tosplit[snp_meta_tosplit$snp_name %in% snps[[i]],]
  snp_meta <- split(snp_meta, f = snp_meta$snp_chr)
  snp_meta_nozeros <- list_drop_empty(snp_meta)
  chr_snps <- list()
  
  for (j in seq_along(snp_meta_nozeros)) {
    chr_snp <- snp_meta_nozeros[[j]]
    chr_snp_pos <- sort(chr_snp$snp_loc)
    chr_int <- list()
    
    for (k in seq_along(chr_snp_pos)) {
      snp_int <- chr_snp_pos[k]:(chr_snp_pos[k] + 150)
      chr_int[[k]] <- intersect(snp_int, chr_snp_pos)
    }
    
    chr_snps[[j]] <- chr_int
  }
  
  gene_snps[[i]] <- chr_snps
}

############################
# Remove redundant clusters
############################

gene_snps_nodup <- list()

for (i in seq_along(gene_snps)) {
  list_un <- list()
  
  for (j in seq_along(gene_snps[[i]])) {
    list_genesnps <- gene_snps[[i]][[j]]
    un <- unique(unlist(list_genesnps))
    dup <- unique(unlist(list_genesnps)[duplicated(unlist(list_genesnps))])
    not_rep <- setdiff(un, dup)
    to_keep <- list()
    
    if (length(dup) > 0) {
      for (l in seq_along(dup)) {
        list_obj <- which(sapply(list_genesnps,function(x) dup[l] %in% x))
        max_len <- max(lengths(list_genesnps[list_obj]))

        to_keep[[l]] <-list_genesnps[list_obj][
            which(lengths(list_genesnps[list_obj]) == max_len)]
      }
    } else {
      to_keep <- list_genesnps
    }
    to_keep <- list_drop_empty(to_keep)
    not_rep_tokeep <- as.list(not_rep[!not_rep %in% unlist(to_keep)]
    )

    list_un[[j]] <- c(to_keep, not_rep_tokeep)
  }
  gene_snps_nodup[[i]] <- list_un
  names(gene_snps_nodup)[i] <- eqtl_genes[i]
}

############################
# Select best SNP
# (lowest p-value)
############################

snps_gene_tokeep <- list()

for (i in seq_along(gene_snps_nodup)) {

  snp_meta <- snp_meta_tosplit[snp_meta_tosplit$snp_name %in% snps[[i]],]
  snp_meta <- split(snp_meta, f = snp_meta$snp_chr)

  eqtl_gene <- names(gene_snps_nodup[i])

  by_gene <- eqtl[eqtl$gene == eqtl_gene,]

  snps_gene_chr <- list()

  for (j in seq_along(gene_snps_nodup[[i]])) {
    list_genesnps <- gene_snps_nodup[[i]][[j]]
    best_snps <- c()
    for (k in seq_along(list_genesnps)) {
      pvalue <- c()
      snp_name_tmp <- c()
      for (l in seq_along(
        unlist(list_genesnps[[k]]))) {
        snp_gene_name <-snp_meta[[j]]$snp_name[
            snp_meta[[j]]$snp_loc == unlist(list_genesnps[[k]])[l]]

        pvalue[l] <- by_gene$`p-value`[
            by_gene$SNP == snp_gene_name[1]]
        snp_name_tmp[l] <- snp_gene_name[1]
      }

      snp_df <- data.frame(snp_name = snp_name_tmp, pvalue = pvalue)
      pvalue_tokeep <- min(pvalue)
      best_snps[k] <-snp_df$snp_name[which(snp_df$pvalue == pvalue_tokeep)[1]]
    }
    snps_gene_chr[[j]] <- best_snps
  }

  snps_gene_tokeep[[i]] <- snps_gene_chr

  names(snps_gene_tokeep)[i] <- eqtl_genes[i]
}

for (i in seq_along(snps_gene_tokeep)) {

  snps_gene_tokeep[[i]] <-
    unlist(snps_gene_tokeep[[i]])
}

############################
# Save representative SNPs
############################

save(
  snps_gene_tokeep,
  file = file.path(
    results_dir,
    "top_snps_ineqtl.RData"
  )
)

############################
# Create cleaned eQTL table
############################

eqtl_clean <- list()

for (i in seq_along(snps_gene_tokeep)) {

  eqtl_gene <- rep(
    names(snps_gene_tokeep[i]),
    length(snps_gene_tokeep[[i]])
  )

  SNP <- unname(
    unlist(snps_gene_tokeep[i])
  )

  eqtl_stats <- eqtl[eqtl$gene == unique(eqtl_gene),]

  eqtl_stats <- eqtl_stats[eqtl_stats$SNP %in% SNP,]

  eqtl_clean[[i]] <- data.frame(SNP = eqtl_stats$SNP,
    eqtl_gene = eqtl_stats$gene, pvalue = eqtl_stats$`p-value`)
}

eqtl_clean <- bind_rows(eqtl_clean)

write.table(eqtl_clean,file.path(results_dir, "eqtl_clean.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

############################
# Annotate SNPs with LD decay
############################

snp_meta_tosplit$LDdecay_int <- hd_LD$Mb[snp_meta_tosplit$snp_chr]

write.table(snp_meta_tosplit,
  file.path(results_dir, "snp_meta.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

############################
# Summary
############################

cat("\nPipeline completed successfully.\n")
cat("Number of eQTL genes:", length(eqtl_genes), "\n")
cat("Number of unique SNPs:", nrow(snp_meta_tosplit), "\n")
cat("Outputs written to:", results_dir, "\n")
