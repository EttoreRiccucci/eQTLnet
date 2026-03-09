library(data.table)

setwd("/path/dir")
meta_annotation<-fread('Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1_ab.1_xref_gene_IDs.txt', data.table = F)
tfs<-scan('TF_list_update.txt', what='character')
tfs<-unique(tfs)
##GS
gs<-fread('AM_Final_Maize_GS.txt', data.table=F)
colnames(gs)<-c('from','to','weight')
exp_tot<-read.csv('Field2023_TPM.csv', header=T, stringsAsFactors = F)
# 
metadata<-read.csv('Field_2023_TPM_Metadata.csv', header=T, stringsAsFactors = F)
rils_name<-fread('INFO.9.founders.all.txt', data.table=F)
rils_name<-data.frame(rils_name$Taxa, rils_name$accession)

##geno
spet<-fread('spet_geno.txt', data.table = F)

#add SSA rils name to exp_tot
meta_rep1<-metadata[which(metadata$Replicate==1),]
meta_rep2<-metadata[which(metadata$Replicate==2),]
setdiff(meta_rep1$Line, meta_rep2$Line) #5 rils have only 1 replicate (1 or 2)

#keep only rils in rep 1 and 2 that are shared with spet
meta_rep1<-meta_rep1[which(meta_rep1$Line %in% colnames(spet[-1])),]
meta_rep2<-meta_rep2[which(meta_rep2$Line %in% colnames(spet[-1])),]

write.table(meta_rep1, 'meta_rep1.txt', quote = F, row.names = F, sep='\t')
write.table(meta_rep2, 'meta_rep2.txt', quote = F, row.names = F, sep='\t')

rawexp_rep1<-exp_tot[,which(colnames(exp_tot) %in% meta_rep1$V)] 
rownames(rawexp_rep1)<-exp_tot$X
rawexp_rep2<-exp_tot[,which(colnames(exp_tot) %in% meta_rep2$V)] 
rownames(rawexp_rep2)<-exp_tot$X

rep1_cnames<-meta_rep1$Line[which(meta_rep1$V %in% colnames(rawexp_rep1))]
rep2_cnames<-meta_rep2$Line[which(meta_rep2$V %in% colnames(rawexp_rep2))]

colnames(rawexp_rep1) <- rep1_cnames
colnames(rawexp_rep2) <- rep2_cnames

rawexp_rep1<-rawexp_rep1[,order(colnames(rawexp_rep1))]
rawexp_rep2<-rawexp_rep2[,order(colnames(rawexp_rep2))]
v_names_rep1<-meta_rep1[which(meta_rep1$Line %in% colnames(rawexp_rep1)[-1]),][,1]
v_names_rep2<-meta_rep2[which(meta_rep2$Line %in% colnames(rawexp_rep2)[-1]),][,1]

write.table(rawexp_rep1, 'exp_tot_r1.txt', quote = F, row.names = F, sep='\t')

write.table(rawexp_rep2, 'exp_tot_r2.txt', quote = F, row.names = F, sep='\t')
#devi guardare nelle due repliche quali geni hanno valore di espressione zero in entrambe le repliche--> li elimini
#prendi solo i geni che hanno al massimo 15% di zeros sulle linee -->vedi quant tf e tfin gs tieni cosi prima di fare eqtl mapping,
#confron

#remove genes with zero expression across replicates
#keep only genes with max 20 % zeros across lines

#look for the number of zeros across lines for each gene in rep1
 count_r1<-c()
 for (i in 1:nrow(rawexp_rep1)) {
   count_r1[i]<-length(which(rawexp_rep1[i,]==0))
    count_r1<-unlist(count_r1, use.names=F)
 }
 write(count_r1,'count_zeros_r1.txt')

  count<-scan('count_zeros_r1.txt', character())
#
#look for the number of zeros across lines for each gene in rep2
count_r2<-c()
for (i in 1:nrow(rawexp_rep2)) {
  count_r2[i]<-length(which(rawexp_rep2[i,]==0))
  # count_r1<-unlist(count_r1, use.names=F)
}
write(count_r2,'count_zeros_r2.txt')
count<-scan('count_zeros_r2.txt', character())

#colnames(exp)[1]<-c('gene_id')

zeros_r1<-data.frame(exp_tot[,1], count_r1)
colnames(zeros_r1)[1]<-c('gene_id')

zeros_r2<-data.frame(exp_tot[,1], count_r2)
colnames(zeros_r2)[1]<-c('gene_id')

XX_exp_genes_r1<-zeros_r1[which(zeros_r1$count_r1<= (ncol(rawexp_rep1))*0.2),]

#length(which(zeros_r2$count_r2 <= (ncol(rawexp_rep2))*0.2))#20656
XX_exp_genes_r2<-zeros_r2[which(zeros_r2$count_r2<= (ncol(rawexp_rep2))*0.2),]
nozeros_r2<-zeros_r2[which(zeros_r2$count_r2 == 0),]#12487
gene_tokeepXX<-intersect(XX_exp_genes_r1[,1], XX_exp_genes_r2[,1])#20234
gene_tokeepnoz<-intersect(nozeros_r1[,1], nozeros_r2[,1])#9984


rawexp_rep1<-data.frame(exp_tot[,1], rawexp_rep1)
colnames(rawexp_rep1)[1]<-'gene_id'
exp_tot_r1XX <- rawexp_rep1[which(rawexp_rep1$gene_id %in% gene_tokeepXX),]
exp_tot_r1noz <- rawexp_rep1[which(rawexp_rep1$gene_id %in% gene_tokeepnoz),]
write.table(exp_tot_r1XX, 'exp_XX_r1.txt', quote = F, row.names = F, sep='\t')
write.table(exp_tot_r1noz, 'exp_noz_r1.txt', quote = F, row.names = F, sep='\t')

rawexp_rep2<-data.frame(exp_tot[,1], rawexp_rep2)
colnames(rawexp_rep2)[1]<-'gene_id'
exp_tot_r2XX <- rawexp_rep2[which(rawexp_rep2$gene_id %in% gene_tokeepXX),]
exp_tot_r2noz <- rawexp_rep2[which(rawexp_rep2$gene_id %in% gene_tokeepnoz),]
write.table(exp_tot_r2XX, 'exp_XX_r2.txt', quote = F, row.names = F, sep='\t')
write.table(exp_tot_r2noz, 'exp_noz_r2.txt', quote = F, row.names = F, sep='\t')

# #how many TF among the new genes??
# length(intersect(exp_tot_r1XX$gene_id, tfs))#1149
# length(intersect(exp_tot_r1noz$gene_id, tfs))#449
# 
# #how many TF among the new genes are in gs$from??
# length(intersect(exp_tot_r1XX$gene_id, unique(gs$from)))#126
# length(intersect(exp_tot_r1noz$gene_id, unique(gs$from)))#55
# 
# #how many genes are in gs$to?
# length(intersect(exp_tot_r1XX$gene_id, unique(gs$to)))#1707
# length(intersect(exp_tot_r1noz$gene_id, unique(gs$to)))#795

#check which genes has zeros in both replicates
##XX
rawexp_rep1<-fread('exp_XX_r1.txt', data.table=F)
rawexp_rep2<-fread('exp_XX_r2.txt', data.table=F)
line_zeros<-c()
coord_1<-list()
for (i in 1:ncol(rawexp_rep1[,-1])){
  for(j in 1:nrow(rawexp_rep1[,-1])){
    if (rawexp_rep1[,-1][j,i] == 0){ #all genes with zero exp in line i
      line_zeros[j]<-j
    } else {
      line_zeros[j]<-NA
    }
    line_zeros <- line_zeros[!is.na(line_zeros)]
  }
  coord_1[[i]]<- line_zeros
}

line_zeros<-c()
coord_2<-list()
for (i in 1:ncol(rawexp_rep2[,-1])){
  for(j in 1:nrow(rawexp_rep2[,-1])){
    if (rawexp_rep2[,-1][j,i] == 0){ #all genes with zero exp in line i
      line_zeros[j]<-j
    } else {
      line_zeros[j]<-NA
    }
    line_zeros <- line_zeros[!is.na(line_zeros)]
  }
  coord_2[[i]]<- line_zeros
}

intersect(colnames(rawexp_rep1[-1]), colnames(rawexp_rep2[-1]))

gene_toremove<-list()
for (i in 1:length(coord_2)){ #look for genes with zero exp in both rep
  gene_toremove[[i]]<-intersect(coord_1[[i]], coord_2[[i]])
}
gene_toremove<-unique(unlist(gene_toremove, use.names=F)) #remove genes with zero exp in both replicates

write(gene_toremove, 'gene_toremove.txt')
gene_toremove<-scan('gene_toremove.txt')
raw_exp_filt_1<- rawexp_rep1[-gene_toremove,]
raw_exp_filt_2<- rawexp_rep2[-gene_toremove,]

# #how many TF among the new genes??
# length(intersect(raw_exp_filt_1$gene_id, tfs))#851
# 
# #how many TF among the new genes are in gs$from??
# length(intersect(raw_exp_filt_1$gene_id, unique(gs$from)))#101
# 
# #how many genes are in gs$to?
# length(intersect(raw_exp_filt_1$gene_id, unique(gs$to)))#1361

write.table(raw_exp_filt_1, 'exp_XXfilt_r1.txt', quote = F, row.names = F, sep='\t')
write.table(raw_exp_filt_2, 'exp_XXfilt_r2.txt', quote = F, row.names = F, sep='\t')

exp2<-fread('exp_XXfilt_r2.txt', data.table=F)
exp1<-fread('exp_XXfilt_r1.txt', data.table=F)


