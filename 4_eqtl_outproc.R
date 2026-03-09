setwd("C:/Users/e.riccucci/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/LD_interval")
hd_LD<-read.csv('half_decay_LD.csv', header=T, stringsAsFactors = F)
meanLD<-mean(hd_LD$Mb)
library(data.table)
#install.packages('leaflet')
library(vctrs)
library(leaflet)
library(dplyr)
library(stringr)

setwd("C:/Users/e.riccucci/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24")


eqtl<-fread('eqtl_10genoPC_logblupsXXnodate.txt', data.table = F)

#dentify all snps genes (TF if the egwas gene was TF or non-TF and tg is the egwas gene was tf) including the one after opening the interval with LD
tfs<-scan('TF_list_update.txt', what = 'character')
tf_id<-tfs

tf_eqtl<-unique(eqtl$gene)[which(unique(eqtl$gene) %in% tf_id)]
notf_eqtl<-unique(eqtl$gene)[which(!unique(eqtl$gene) %in% tf_id)]

snp_name<-unique(eqtl$SNP)
write(snp_name, 'snp_name_logblupsXXnodate.txt')

snp_chr<-c()
snp_loc<-c()
for (i in 1:length(snp_name)){
  snp_chr[i]<-str_split(snp_name[i], ":")[[1]][1]
  snp_loc[i]<-str_split(snp_name[i], ":")[[1]][2]
}
head(snp_loc)
head(snp_name)

write(snp_chr, 'snp_chr_logblupsXXnodate.txt')
write(snp_loc, 'snp_loc_logblupsXXnodate.txt')

snp_meta_tosplit <- data.frame(snp_name, as.numeric(snp_loc), as.numeric(snp_chr))
colnames(snp_meta_tosplit)[2]<-'snp_loc'
colnames(snp_meta_tosplit)[3]<-'snp_chr'
head(snp_meta_tosplit)
dim(snp_meta_tosplit)

eqtl_genes<-unique(eqtl$gene)
snps<-list() #snps for each eqtl gene
for (i in 1:length(eqtl_genes)){
  snps[[i]]<-eqtl[which(eqtl$gene %in% eqtl_genes[i]),][,1]
  names(snps)[i]<-eqtl_genes[i]
}

#devo dividere gli snps per cromosomi perche le posizioni sono riferite al singolo chr e non all'intero genoma
gene_snps<-list()
for (i in 1:length(snps)){ #number of genes #i=3
  #if (length(snps[[i]])>1){
  snp_meta<-snp_meta_tosplit[which(snp_meta_tosplit$snp_name %in% snps[[i]]),]
  snp_meta<-split(snp_meta, f=snp_meta$snp_chr)
  snp_meta_nozeros<-list_drop_empty(snp_meta)
  chr_snps<-list()
  for (j in 1:length(snp_meta_nozeros)){ #number of chromosomes with snps asociated to the ith gene
    chr_snp<-snp_meta_nozeros[[j]]
    chr_snp_pos<-chr_snp$snp_loc
    chr_snp_pos<-chr_snp_pos[order(chr_snp_pos)]
    chr_int<-list()
    for (k  in 1:(length(chr_snp_pos))){ #number of snps reported in chr j
      snp_int<-c(chr_snp_pos[k] : (chr_snp_pos[k]+150))
      chr_int[[k]]<-intersect(snp_int, chr_snp_pos)
    }
    chr_snps[[j]]<-chr_int #snps in chr j
  }
  gene_snps[[i]]<-chr_snps #list of lists: first level (i): genes; second (j) chromosomes
}

#hai i vari gruppi di snps racchiusi in 150 bps, adesso vuoi togliere quelli che sono da soli e sono compresi in un altro gruppo di snps
gene_snps_nodup<-list()
for (i in 1:length(gene_snps)){#i=3
  list_un<-list()
  for(j in 1:length(gene_snps[[i]])){#j=4
    list_genesnps<-gene_snps[[i]][[j]]
    un<-unique(unlist(list_genesnps))
    dup<-unique(unlist(list_genesnps)[which(duplicated(unlist(list_genesnps)))])
    not_rep<-setdiff(un, dup)
    to_keep<-list()
    if (length(dup)>0){
      for (l in 1:length(dup)){ 
        list_obj<-which(sapply(list_genesnps, FUN=function(X) dup[l] %in% X)) #which list element contains a particular value in R?
        max_len<-max(lengths(list_genesnps[list_obj]))
        to_keep[[l]]<-list_genesnps[list_obj][which(length(list_genesnps[list_obj]) == max_len)]
      }
    } else {
      to_keep<-list_genesnps
    }
    to_keep<-list_drop_empty(to_keep)
    not_rep_tokeep<-as.list(not_rep[which(!not_rep %in% unlist(to_keep))])
    list_un[[j]]<-c(to_keep, not_rep_tokeep)
    #list_un[[j]]<-list_drop_empty(list_un)
  }
  gene_snps_nodup[[i]]<-list_un
  names(gene_snps_nodup)[i]<-eqtl_genes[i]
}

gene_snps_nodup[i]
snps[67]

snps_gene_tokeep<-list()
for (i in 1:length(gene_snps_nodup)){
  snp_meta<-snp_meta_tosplit[which(snp_meta_tosplit$snp_name %in% snps[[i]]),]
  snp_meta<-split(snp_meta, f=snp_meta$snp_chr)
  eqtl_gene<-names(gene_snps_nodup[i])
  by_gene<-eqtl[which(eqtl$gene == eqtl_gene),]
  snps_gene_chr<-list()
  for(j in 1:length(gene_snps_nodup[[i]])){
    list_genesnps<-gene_snps_nodup[[i]][[j]]#snps in chr
    best_snps<-c()
    for (k in 1:length(list_genesnps)){ #number of groups of snp (150 bps) in each chr
      pvalue<-c()
      snp_name<-c()
      for (l in 1:length(unlist(list_genesnps[[k]]))){ #number of snps in each group
        #snp_gene_name<-snp_meta_tosplit$snp_name[which(snp_meta_tosplit$snp_loc == unlist(list_genesnps[[k]])[l])]
        #devi selezionare il chr oportuno, snp diversi potrebbero avere la stessa posizione su chr diversi
        snp_gene_name<-snp_meta[[j]]$snp_name[which(snp_meta[[j]]$snp_loc == unlist(list_genesnps[[k]])[l])]
        pvalue[l]<-by_gene$`p-value`[which(by_gene$SNP == snp_gene_name[1])]
        snp_name[l]<-snp_gene_name
      } 
      snp_df<-data.frame(snp_name, pvalue)
      pvalue_tokeep<-min(pvalue)
      best_snps[k]<-snp_df$snp_name[which(snp_df[,2]==pvalue_tokeep)][1]
    }
    snps_gene_chr[[j]]<-best_snps  
  }
  snps_gene_tokeep[[i]]<-snps_gene_chr
  names(snps_gene_tokeep)[i]<-eqtl_genes[i]
}



for (i in 1:length(snps_gene_tokeep)){
  snps_gene_tokeep[[i]]<-unlist(snps_gene_tokeep[[i]])
}


save(snps_gene_tokeep, file='top_snps_ineqtl.RData')

#il file "eqtl_clean.txt" è il dataframe con i geni eqtl e gli snp filtrati per reads SPET, con riportati i rispettivi pvalue

eqtl_gene<-c()
SNP<-c()
eqtl_clean<-list()
for (i in 1:length(snps_gene_tokeep)){
  eqtl_gene<-rep(names(snps_gene_tokeep[i]), length(snps_gene_tokeep[[i]]))
  SNP<-unname(unlist(snps_gene_tokeep[i]))
  eqtl_stats<-eqtl[which(eqtl$gene == unique(eqtl_gene)),]
  eqtl_stats<-eqtl_stats[which(eqtl_stats$SNP %in% SNP),]
  eqtl_clean[[i]]<-data.frame(eqtl_stats$SNP, eqtl_stats$gene, eqtl_stats$`p-value`)
  colnames(eqtl_clean[[i]])<-c('SNP', 'eqtl_gene', 'pvalue')
}

eqtl_clean<-bind_rows(eqtl_clean)
head(eqtl_clean)
head(eqtl)
write.table(eqtl_clean, "eqtl_clean.txt", quote=F, sep="\t", row.names=FALSE)

#create snp_meta with ld infeqtl_genes#create snp_meta with ld info per chr
for (i in 1:nrow(snp_meta_tosplit)){
  if (snp_meta_tosplit$snp_chr[i] == 1){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[1]
  } else if (snp_meta_tosplit$snp_chr[i] == 2){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[2]
  } else if (snp_meta_tosplit$snp_chr[i] == 3){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[3]
  } else if (snp_meta_tosplit$snp_chr[i] == 4){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[4]
  } else if (snp_meta_tosplit$snp_chr[i] == 5){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[5]
  } else if (snp_meta_tosplit$snp_chr[i] == 6){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[6]
  } else if (snp_meta_tosplit$snp_chr[i] == 7){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[7]
  } else if (snp_meta_tosplit$snp_chr[i] == 8){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[8]
  } else if (snp_meta_tosplit$snp_chr[i] == 9){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[9]
  } else if (snp_meta_tosplit$snp_chr[i] == 10){
    snp_meta_tosplit$LDdecay_int[i]<-hd_LD$Mb[10]
  }
} #$LDdecay_int reports LD interval in Mb

head(snp_meta_tosplit)
write.table(snp_meta_tosplit, "snp_meta.txt", quote=F, sep="\t", row.names=FALSE)

