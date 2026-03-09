# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(PRROC)
library(scales) # Include scales for rescaling
library(precrec)

# Function to load RData file and return its content
loadRData <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load data
setwd("C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24/")
beta_df <- loadRData('beta_dfnorm.Rdata')
length(unique(beta_df$from))
length(unique(beta_df$to))

setwd("C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24/app2")
beta_mx <- loadRData('beta_mx.Rdata')
head(beta_mx)[,1:4]
dim(beta_mx)

length(which(colnames(beta_mx) %in% positive_interactions$from))
length(which(rownames(beta_mx) %in% positive_interactions$to))

setwd("C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24/")
rsq<-loadRData('rsq.Rdata')
adj_rsq<-loadRData('adj_rsq.Rdata')
rsq

##### read Gold standard network
##############################
#Because GS contains only positive interactions,
# Load the positive interactions dataframe
setwd("C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24")
positive_interactions <- fread("AM_Final_Maize_GS.txt", data.table=F )
colnames(positive_interactions) <- c( "from","to", 'weight')
# Get the list of all possible gene pairs that do not appear in the positive interactions dataframe
all_gene_pairs <- expand.grid(from = unique(positive_interactions$from),
                              to = unique(positive_interactions$to))
negative_interactions <- all_gene_pairs[!(paste(all_gene_pairs$from, all_gene_pairs$to) %in% paste(positive_interactions$from, positive_interactions$to)), ]
negative_interactions$weight <- 0

#####2-1- For the balance case you set 20 different seeds in R to sample the same 
#number of negative interactions as the positive ones. By doing this you will obtain 20 different negative 
#data frames with the same length (number of interactions) as the positive set.
#2-2- Combine the unique positive and one of the 20 negative sets every time to 
#get 20 different GS.
#2-3- For each GS, you compute the AUROC, AUPR, and F1 scores and save them as a table

# set.seed(20)
# neg_tk<-negative_interactions[sample(nrow(negative_interactions), nrow(positive_interactions)), ]

#neg_tk_2<-negative_interactions[sample(nrow(negative_interactions), nrow(positive_interactions)), ]
neg_tk<-list()
gs<-list()
for (i in 1:20){
  set.seed(i)
  neg_tk[[i]]<-negative_interactions[sample(nrow(negative_interactions), nrow(positive_interactions)),]
  gs[[i]] <- rbind(positive_interactions, neg_tk[[i]])
}

neg_tk[[14]]


####no rsq/adjrsq--> all shared interactions
# Find shared 'int' values between the two datasets

beta_df$int <- paste0(beta_df$from, beta_df$to)
shared<-list() 
common_gs<-list()
common_beta_df <- list()
curves <- list()
for(i in 1:length(gs)){
  gs[[i]]$int <- paste0(gs[[i]]$from, gs[[i]]$to)
  shared[[i]] <- intersect(beta_df$int, gs[[i]]$int)
  common_beta_df[[i]] <- beta_df[beta_df$int %in% shared[[i]], ]
  common_gs[[i]] <- gs[[i]][gs[[i]]$int %in% shared[[i]], ]
  
  # Remove duplicates from common_gs based on 'int' column
  common_gs_unique <- common_gs[[i]] %>% distinct(int, .keep_all = TRUE)
  
  # Order both datasets by 'int' column
  common_gs_unique_ordered <- common_gs_unique %>% arrange(int)
  common_beta_df_ordered <- common_beta_df[[i]] %>% arrange(int)
  
  # Ensure both datasets are aligned by the 'int' column
  if (!identical(common_gs_unique_ordered$int, common_beta_df_ordered$int)) {
    stop("The datasets are not aligned by the 'int' column")
  }
  
  # Select rows in common_beta_df_ordered where corresponding common_gs_unique_ordered$weight == 1
  positive_class <- common_beta_df_ordered %>% filter(common_gs_unique_ordered$weight == 1) # only 485 positive interactions
  
  # Select rows in common_beta_df_ordered where corresponding common_gs_unique_ordered$weight == 0
  negative_class <- common_beta_df_ordered %>% filter(common_gs_unique_ordered$weight == 0) # 2311 negative ones
  
  sscurves_nweighted_rr_pval <- evalmod(scores = common_beta_df_ordered$weight, labels = common_gs_unique_ordered$weight)
  plot(sscurves_nweighted_rr_pval)
  curves[[i]]<-auc(sscurves_nweighted_rr_pval)
  
}

curves[[20]]

# Filter both datasets to include only rows with shared 'int' values
common_beta_df <- beta_df[beta_df$int %in% shared, ]
common_gs <- gs[gs$int %in% shared, ]

length(unique(gs$int))#4423868 all the interactions are unique!
length(unique(common_gs$int))#2796 all the interactions are unique!

# gs_dup<-gs[which(duplicated(gs$int)),]
# head(gs_dup)
# gs[which(gs$int == 'Zm00001eb067270Zm00001eb065870'),]
# gs[556,]
# which(gs_dup$weight==0)
# length(unique(gs_dup$from))
# length(unique(gs_dup$to))
# length(unique(gs_dup$int))#3816

# Remove duplicates from common_gs based on 'int' column
common_gs_unique <- common_gs %>% distinct(int, .keep_all = TRUE)

# Order both datasets by 'int' column
common_gs_unique_ordered <- common_gs_unique %>% arrange(int)
common_beta_df_ordered <- common_beta_df %>% arrange(int)

# Ensure both datasets are aligned by the 'int' column
if (!identical(common_gs_unique_ordered$int, common_beta_df_ordered$int)) {
  stop("The datasets are not aligned by the 'int' column")
}

# Select rows in common_beta_df_ordered where corresponding common_gs_unique_ordered$weight == 1
positive_class <- common_beta_df_ordered %>% filter(common_gs_unique_ordered$weight == 1) # only 485 positive interactions

# Select rows in common_beta_df_ordered where corresponding common_gs_unique_ordered$weight == 0
negative_class <- common_beta_df_ordered %>% filter(common_gs_unique_ordered$weight == 0) # 2311 negative ones

sscurves_nweighted_rr_pval <- evalmod(scores = common_beta_df_ordered$weight, labels = common_gs_unique_ordered$weight)
plot(sscurves_nweighted_rr_pval)
auc(sscurves_nweighted_rr_pval)


##rsq
trs<-seq(from = 0, to = 1, by = 0.05)
roc<-rep(0, length(trs))
prc<-rep(0, length(trs))

curves<-data.frame(trs, roc, prc)
for (i in 1:length(trs)){
  rsq_tk<-rsq[which(rsq > trs[i])]
  myint_tk<-beta_df[which(beta_df$to %in% names(rsq_tk)),]
  
  shared <- intersect(myint_tk$int, gs$int)
  
  # Filter both datasets to include only rows with shared 'int' values
  common_beta_df <- beta_df[beta_df$int %in% shared, ]
  common_gs <- gs[gs$int %in% shared, ]
  
  length(unique(gs$int))#4423868 all the interactions are unique!
  length(unique(common_gs$int))#2796 all the interactions are unique!
  # Remove duplicates from common_gs based on 'int' column
  common_gs_unique <- common_gs %>% distinct(int, .keep_all = TRUE)
  
  # Order both datasets by 'int' column
  common_gs_unique_ordered <- common_gs_unique %>% arrange(int)
  common_beta_df_ordered <- common_beta_df %>% arrange(int)
  
  # Ensure both datasets are aligned by the 'int' column
  if (!identical(common_gs_unique_ordered$int, common_beta_df_ordered$int)) {
    stop("The datasets are not aligned by the 'int' column")
  }
  
  # Select rows in common_beta_df_ordered where corresponding common_gs_unique_ordered$weight == 1
  positive_class <- common_beta_df_ordered %>% filter(common_gs_unique_ordered$weight == 1) # only 485 positive interactions
  
  # Select rows in common_beta_df_ordered where corresponding common_gs_unique_ordered$weight == 0
  negative_class <- common_beta_df_ordered %>% filter(common_gs_unique_ordered$weight == 0) # 2311 negative ones
  
  sscurves_nweighted_rr_pval <- evalmod(scores = common_beta_df_ordered$weight, labels = common_gs_unique_ordered$weight)
  plot(sscurves_nweighted_rr_pval)
  auc(sscurves_nweighted_rr_pval)
  
  curves[i,2]<-auc(sscurves_nweighted_rr_pval)[1,4]
  curves[i,3]<-auc(sscurves_nweighted_rr_pval)[2,4]
}

curves_rsq<-curves
write.table(curves_rsq, "C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/curves_rsq_tnormalized_snpfrom.txt", quote=F, sep='\t', row.names = F)

####adj_rsq
trs<-seq(from = 0, to = 1, by = 0.05)
roc<-rep(0, length(trs))
prc<-rep(0, length(trs))

curves<-data.frame(trs, roc, prc)
for (i in 1:length(trs)){
  adj_rsq_tk<-adj_rsq[which(adj_rsq > trs[i])]
  myint_tk<-beta_df[which(beta_df$to %in% names(adj_rsq_tk)),]
  
  shared <- intersect(myint_tk$int, gs$int)
  
  # Filter both datasets to include only rows with shared 'int' values
  common_beta_df <- beta_df[beta_df$int %in% shared, ]
  common_gs <- gs[gs$int %in% shared, ]
  
  length(unique(gs$int))#4423868 all the interactions are unique!
  length(unique(common_gs$int))#2796 all the interactions are unique!
  # Remove duplicates from common_gs based on 'int' column
  common_gs_unique <- common_gs %>% distinct(int, .keep_all = TRUE)
  
  # Order both datasets by 'int' column
  common_gs_unique_ordered <- common_gs_unique %>% arrange(int)
  common_beta_df_ordered <- common_beta_df %>% arrange(int)
  
  # Ensure both datasets are aligned by the 'int' column
  if (!identical(common_gs_unique_ordered$int, common_beta_df_ordered$int)) {
    stop("The datasets are not aligned by the 'int' column")
  }
  
  # Select rows in common_beta_df_ordered where corresponding common_gs_unique_ordered$weight == 1
  positive_class <- common_beta_df_ordered %>% filter(common_gs_unique_ordered$weight == 1) # only 485 positive interactions
  
  # Select rows in common_beta_df_ordered where corresponding common_gs_unique_ordered$weight == 0
  negative_class <- common_beta_df_ordered %>% filter(common_gs_unique_ordered$weight == 0) # 2311 negative ones
  
  sscurves_nweighted_rr_pval <- evalmod(scores = common_beta_df_ordered$weight, labels = common_gs_unique_ordered$weight)
  plot(sscurves_nweighted_rr_pval)
  auc(sscurves_nweighted_rr_pval)
  
  curves[i,2]<-auc(sscurves_nweighted_rr_pval)[1,4]
  curves[i,3]<-auc(sscurves_nweighted_rr_pval)[2,4]
}
curves_adj_rsq<-curves
write.table(curves_adj_rsq, "C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/previenna/curves_adjrsq_normalized_snpfrom.txt", quote=F, sep='\t', row.names = F)
