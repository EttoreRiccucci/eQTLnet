###############################################################
# 08_Identify_Genes_in_eQTL_Intervals
#
# Description:
#   1. Load filtered eQTL SNPs.
#   2. Construct LD-based genomic intervals around SNPs.
#   3. Identify genes overlapping each interval.
#   4. Generate SNP-gene overlap table.
#
# Inputs:
#   results/snp_meta.txt
#   results/v5_meta_annotation.txt
#   results/eqtl_cleanwithset.txt
#
# Outputs:
#   results/genes_on_snpint.txt
#
###############################################################

############################
# Libraries
############################

library(data.table)
library(GenomicRanges)

############################
# Load data
############################

snp_meta <- fread( "snp_meta.txt", data.table = FALSE)

meta_annotation <- fread("v5_meta_annotation.txt", data.table = FALSE)

eGWAS <- fread("eqtl_cleanwithset.txt", data.table = FALSE)

############################
# Standardize column names
############################

colnames(eGWAS) <- c(
  "egwas_gene",
  "SNP",
  "pvalue",
  "egwas_gene_set"
)

############################
# Retrieve SNP metadata
############################

assosnp_loc <- snp_meta[
  snp_meta$snp_name %in% eGWAS$SNP,
]

############################
# Define LD intervals
############################

ld_bp <- assosnp_loc$LDdecay_int * 1000000

assosnp_loc$int_start <- pmax(assosnp_loc$snp_loc - ld_bp, 0)

assosnp_loc$int_end <- assosnp_loc$snp_loc + ld_bp

############################
# Convert SNPs to GRanges
############################

assosnp_loc$snp_chr <- paste0(
  "chr",
  assosnp_loc$snp_chr
)

gr_snp <- GRanges(
  seqnames = assosnp_loc$snp_chr,
  ranges = IRanges(
    start = assosnp_loc$int_start,
    end = assosnp_loc$int_end
  )
)

############################
# Convert genes to GRanges
############################

gr_all <- GRanges(
  seqnames = meta_annotation$new_gene_model_chr,
  ranges = IRanges(
    start = meta_annotation$start,
    end = meta_annotation$end
  )
)

############################
# Find overlaps
############################

OL_hits <- findOverlaps(
  query = gr_snp,
  subject = gr_all
)


############################
# Build overlap dataframe
############################

genes_on_snpint <- cbind(
  assosnp_loc[queryHits(OL_hits), ],
  meta_annotation[subjectHits(OL_hits), ]
)

genes_on_snpint <- data.frame(
  SNP = genes_on_snpint$snp_name,
  chr_snp = genes_on_snpint$snp_chr,
  pos_snp = genes_on_snpint$snp_loc,
  int_snp_start = genes_on_snpint$int_start,
  int_snp_end = genes_on_snpint$int_end,
  gene_id = genes_on_snpint$new_gene_model_ID,
  chr_gene = genes_on_snpint$new_gene_model_chr,
  gene_start = genes_on_snpint$start,
  gene_end = genes_on_snpint$end,
  gene_set = genes_on_snpint$gene_set
)

############################
# Extract only SNPs present
# in the final eQTL table
############################

egwas_snp <- eGWAS[
  eGWAS$SNP %in% genes_on_snpint$SNP,
]

genes_on_egwas_snpint <- genes_on_snpint[
  genes_on_snpint$SNP %in% eGWAS$SNP,
]

############################
# Export results
############################

write.table(genes_on_snpint, "genes_on_snpint.txt",  sep = "\t",  quote = FALSE,  row.names = FALSE)

