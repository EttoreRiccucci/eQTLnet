library(data.table)
setwd("/path/dir")
tfs<-scan('TF_list_update.txt', what = 'character')
tfs<-unique(tfs)
meta_annotation<-fread('Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1_ab.1_xref_gene_IDs.txt', data.table = F)
##GS
gs<-fread('AM_Final_Maize_GS.txt', data.table=F)
colnames(gs)<-c('from','to','weight')
######COMPARATIVE analysis--> only genes in GS
gs_from<-intersect(unique(gs$from), tfs)
gs_to<-unique(gs$to)
length(intersect(gs_from, gs_to)) #same length(gs_from)

#use all genes in gs_to for comparative analysis
length(unique(gs_to))
regulated<-unique(gs_to) #all regulated (include all regulators) #31298

regulated_tf<-intersect(regulated, tfs) #1951
regulated_notf<-setdiff(regulated, regulated_tf)#1045
regulators<- gs_from #216
tf_noregulators<-setdiff(regulated_tf, regulators) #1740

#sets names: comp_noTF, comp_TFregulators, comp_TFnoregulators
comp_noTF = regulated_notf
comp_TFregulators = regulators
comp_TFnoregulators = tf_noregulators
write(comp_noTF, 'comp_noTF.txt')
write(comp_TFregulators, 'comp_TFregulators.txt')
write(comp_TFnoregulators, 'comp_TFnoregulators.txt')

#DISCOVERY --> in the final analysis I will use all the genes, here I will take all those that are non in gs
disc_genes<- setdiff(meta_annotation$new_gene_model_ID, regulated) #31260
disc_TF<-intersect(disc_genes, tfs)#226
disc_noTF<-setdiff(disc_genes, disc_TF)#31034
write(disc_TF, 'disc_TF.txt')
write(disc_noTF, 'disc_noTF.txt')
##########

meta_annotation$gene_set<-c()
for (i in 1:length(meta_annotation$new_gene_model_ID)){
  if (meta_annotation$new_gene_model_ID[i] %in% comp_noTF) {
    meta_annotation$set[i]<-c('comp_noTF')
  } else if (meta_annotation$new_gene_model_ID[i] %in% comp_TFregulators){
    meta_annotation$set[i]<-c('comp_TFregulators')
  } else if (meta_annotation$new_gene_model_ID[i] %in% comp_TFnoregulators) {
    meta_annotation$set[i]<-c('comp_TFnoregulators')
  } else if (meta_annotation$new_gene_model_ID[i] %in% disc_noTF) {
    meta_annotation$set[i]<-c('disc_noTF')
  } else if (meta_annotation$new_gene_model_ID[i] %in% disc_TF) {
    meta_annotation$set[i]<-c('disc_TF')
  }else {
    meta_annotation$set[i]<-c('?')
  }
}
head(meta_annotation)

write.table(meta_annotation, 'v5_meta_annotation.txt', sep='\t', quote= F, row.names= F)
length(which(meta_annotation$gene_set == '?'))

#################################################

#assign set to our eqtl results
eqtl<-fread("eqtl_clean.txt", data.table = F)
meta_eqtl<-meta_annotation[which(meta_annotation$new_gene_model_ID %in% eqtl$eqtl_gene),]
head(meta_eqtl)

eqtl<-eqtl[which(eqtl$eqtl_gene %in% meta_eqtl$new_gene_model_ID),]
head(eqtl)
gene_set<-data.frame(meta_eqtl$new_gene_model_ID, meta_eqtl$set)
colnames(gene_set)<-c('eqtl_gene', 'set')
eqtl<-merge(eqtl,gene_set, by='eqtl_gene')
head(eqtl)
write.table(eqtl, 'eqtl_cleanwithset.txt', sep='\t', quote= F, row.names= F)
length(unique(eqtl$eqtl_gene))

eqtl<-split(eqtl, f=eqtl$set)
names(eqtl)
dim(eqtl[[1]])#comp_noTF 669
length(unique(eqtl[[1]]$eqtl_gene)) #198
dim(eqtl[[2]])#comp_TFnoregulators 2539
length(unique(eqtl[[2]]$eqtl_gene))#230
dim(eqtl[[3]])#comp_TFregulators 93
length(unique(eqtl[[3]]$eqtl_gene))#36
dim(eqtl[[4]])#disc_noTF 16718
length(unique(eqtl[[4]]$eqtl_gene))#4892 
dim(eqtl[[5]])#disc_TF 86
length(unique(eqtl[[5]]$eqtl_gene))#31
class(eqtl[[2]])

