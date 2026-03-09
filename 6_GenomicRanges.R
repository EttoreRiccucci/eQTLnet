setwd("/path/dir")
library('data.table')

################################################
#blupsXX
snp_meta<-fread('snp_meta.txt', data.table=F)
head(snp_meta)
meta_annotation<-fread('v5_meta_annotation.txt', data.table = F)
head(meta_annotation)

eGWAS<-fread("eqtl_cleanwithset.txt", data.table=F)
head(eGWAS)
colnames(eGWAS)<-c('egwas_gene', 'SNP', 'pvalue', 'egwas_gene_set')


assosnp_loc<-snp_meta[which(snp_meta$snp_name %in% eGWAS$SNP),]
head(assosnp_loc)
#define interval around associated snps
##########SNP centered in the middle of the interval#####
########################################################################
assosnp_loc$int_start<-c()
for (i in 1:nrow(assosnp_loc)){
  if  (assosnp_loc$snp_loc[i]-((assosnp_loc$LDdecay_int[i])*1000000) > 0){  #assosnp_loc$LDdecay_int[i]%/%2 if you want to take only ld interval, not double it before and after the snp
    assosnp_loc$int_start[i]<-assosnp_loc$snp_loc[i]-((assosnp_loc$LDdecay_int[i])*1000000) #assosnp_loc$LDdecay_int[i]%/%2
  }  else{
    assosnp_loc$int_start[i]<- c(0)
  }
}
head(assosnp_loc)
assosnp_loc$int_end<-c()
for (i in 1:nrow(assosnp_loc)){
  assosnp_loc$int_end[i]<-assosnp_loc$snp_loc[i]+((assosnp_loc$LDdecay_int[i])*1000000) #assosnp_loc$LDdecay_int[i]%/%2
}
head(assosnp_loc)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GenomicRanges")  
library(GenomicRanges)
#The GRanges class represents a collection of genomic ranges that each have a single start and end location on the genome.
assosnp_loc$snp_chr <- paste('chr', assosnp_loc$snp_chr, sep='')
gr_snp=GRanges(seqnames=assosnp_loc$snp_chr,
               ranges=IRanges(start=assosnp_loc$int_start,
                              end=assosnp_loc$int_end),
)
gr_snp

gr_all=GRanges(seqnames=meta_annotation$new_gene_model_chr,
               ranges=IRanges(start=meta_annotation$start,
                              end=meta_annotation$end),
)
gr_all

#findOverlaps takes a query and a subject as inputs and returns a Hits object containing the index pairings for the overlapping elements.
OL_hits<-findOverlaps(gr_snp, gr_all)
OL_hits

length(assosnp_loc[queryHits(OL_hits),][,1])
length(meta_annotation[subjectHits(OL_hits),][,1])

genes_on_snpint<-cbind(assosnp_loc[queryHits(OL_hits),], meta_annotation[subjectHits(OL_hits),])
head(genes_on_snpint)
genes_on_snpint<-data.frame(genes_on_snpint$snp_name, genes_on_snpint$snp_chr, genes_on_snpint$snp_loc, genes_on_snpint$int_start, genes_on_snpint$int_end, genes_on_snpint$new_gene_model_ID, genes_on_snpint$new_gene_model_chr, genes_on_snpint$start, genes_on_snpint$end, genes_on_snpint$set)
colnames(genes_on_snpint)<-c('SNP', 'chr_snp', 'pos_snp', 'int_snp_start', 'int_snp_end', 'gene_id', 'chr_gene', 'gene_start', 'gene_end', 'set')

# insert egwas_genes positions on the s2 df
egwas_snp<-eGWAS[which(eGWAS$SNP %in% genes_on_snpint$SNP),]
genes_on_egwas_snpint<-genes_on_snpint[which(genes_on_snpint$SNP %in% eGWAS$SNP),]
write.table(genes_on_snpint, "genes_on_snpint.txt", quote=F, sep='\t', row.names = F)


