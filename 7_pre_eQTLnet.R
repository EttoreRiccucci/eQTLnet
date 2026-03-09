library(tidyr)
library(data.table)
setwd("C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24/app2")
load('complete_mx.Rdata')
dim(complete_mx)
setwd("C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24")
gs<-fread('AM_Final_Maize_GS.txt', data.table=F)
colnames(gs)<-c('from','to','weight')


gene_exp<-fread('log_blups_XXnodaterandom.txt', data.table = F)
row.names(gene_exp) <- gene_exp[,1]
gene_exp <- gene_exp[,-1]
gene_exp <- t(gene_exp)
head(gene_exp)[,1:4]
tfs<-scan('TF_list_update.txt', what='character')
tfs<-unique(tfs)


length(intersect(colnames(complete_mx), gs$from )) #101
length(intersect(rownames(complete_mx), gs$to )) #1361


complete_df<-as.data.frame(t(complete_mx))
complete_df<- cbind(rownames(t(complete_mx)), complete_df)
head(complete_df)[,1:4]
complete_df<-complete_df %>%
  pivot_longer(
    cols = 2:ncol(complete_df),
    names_to = 'to',
    values_to = 'weight'
  )

complete_df<- as.data.frame(complete_df)
colnames(complete_df)[1]<-'from'

dim(complete_df)
length(unique(complete_df$from))#851
length(unique(complete_df$to))#15929

tochange<-which(!complete_df$weight == 1)
complete_df[tochange,3]<-1
by_regulated<-split(complete_df, f = complete_df$to)

####y
regulated<-c()
for (i in 1:length(by_regulated)){
  regulated[i]<-names(by_regulated[i])
}
exp_regulated <- as.data.frame(gene_exp[,which(colnames(gene_exp) %in% regulated)])
class(exp_regulated)
exp_regulated<-as.list(exp_regulated) #list of regulated

###x
# regulators<-list() #regulators for each regulated
# exp_regulators <- list()
# for (i in 1:length(by_regulated)) {
#   regulators[[i]]<-by_regulated[[i]][,1]
#   exp_regulators[[i]] <- as.data.frame(gene_exp[,which(colnames(gene_exp) %in% regulators[[i]])])
#   colnames(exp_regulators[[i]])<-colnames(gene_exp)[which(colnames(gene_exp) %in% regulators[[i]])]
# }
# names(regulators) <- names(by_regulated) #i nomi di ciascun oggetto della lista corrispondono al regulated gene
# names(exp_regulators) <- names(by_regulated)
# dim(exp_regulators[[6]])
# length(exp_regulators)
# 
# 
# dim_exp_regulators<-c()
# for (i in 1:length(exp_regulators)){
#   dim_exp_regulators[i]<-length(exp_regulators[[i]])
# }
# dim_exp_regulators

#for regulators_exp I need a unique matrix with all exp values for TFs
exp_regulators<-gene_exp[,which(colnames(gene_exp) %in% colnames(complete_mx))]
dim(exp_regulators)
head(exp_regulators)[,1:4]
###pv
pvalue<-list()
for(i in 1:length(by_regulated)){
  by_regulated[[i]]<-by_regulated[[i]][order(by_regulated[[i]][,1]),]
  pvalue[[i]]<-by_regulated[[i]]$weight
  names(pvalue[[i]])<-by_regulated[[i]]$from
}

dim_pv<-c()
for (i in 1:length(pvalue)){
  dim_pv[i]<-length(pvalue[[i]])
}

dim_pv == dim_exp_regulators

dim_pv[i]
dim_exp_regulators[i]
pvalue[i]
length(pvalue[[i]])
head(exp_regulators[[i]])
by_regulated[i]
dim(by_regulated[[i]])


names(pvalue[[24]])
colnames(exp_regulators)

save(pvalue, file = 'pvalue_completenopv.Rdata')
save(exp_regulators, file = 'regulators_completenopv.Rdata')
save(exp_regulated, file = 'regulated_completenopv.Rdata')

