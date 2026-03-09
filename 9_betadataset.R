##########average_nosets
setwd("C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24/")
library(data.table)
#install.packages('vctrs')
library(vctrs)
#install.packages('rlist')
library(rlist)
library(dplyr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

x_g <- loadRData('regulators_completenopv.Rdata')
dim(x_g)
y_g <- loadRData('regulated_completenopv.Rdata')
beta <- loadRData('beta_cv_nopv.Rdata')
length(beta[[14]])
names(beta[[12]])[1]
y_scaled <- list()
for (i in 1:length(y_g)){
  y_scaled[[i]] <- scale(y_g[[i]], center = TRUE, scale = FALSE)#and set it as as.matrix()
}
names(y_scaled)<-names(y_g)

#lookfor constant y
id<-c()
l<-c()
for ( i in 1:length(y_scaled)){
  l[i]<-length(unique(y_scaled[[i]]))
  
}

id<-which(l==1) 
#remove constant y
y_scaled<-y_scaled[-id]

from<-c()
to<-c()
weight<-c()
beta_df<-list()
for (i in 1:length(y_scaled)){
  from<- colnames(x_g)
  to<- rep(names(y_scaled)[i], ncol(x_g))
  names(beta[[i]])<-c()
  weight<-beta[[i]]
  beta_df[[i]]<-data.frame(from, to, weight)
} 

beta_df<-bind_rows(beta_df)
dim(beta_df) #13541963 3
head(beta_df)
save(beta_df, file = 'beta_df_nopv.Rdata')

####create beta_df normalized
beta_cv<-loadRData('beta_cv_nopv.Rdata')

#remove self interactions
# for (i in 1:length(beta_cv)){
#   for (j in 1:length(beta_cv[[i]])){
#     if (names(beta_cv[[i]])[j] == names(beta_cv[i])){
#       beta_cv[[i]]<-beta_cv[[i]][-j]
#     }
#   }
# }
# 
# beta_cv<-list_drop_empty(beta_cv)
#save(beta_cv, file = 'beta_cv_noself.Rdata')

min_max_normalize <- function(x) {
  abs((x - min(x)) / (max(x) - min(x)))
}

beta_norm<-list()
for (i in 1:length(beta_cv)){
  beta_norm[[i]]<-unname(beta_cv[[i]])
  l<-length(beta_cv[[i]])
  if (l==1){
    beta_norm[[i]]<-beta_cv[[i]]
  } else {
    beta_norm[[i]]<-min_max_normalize(beta_cv[[i]])
  }
}
names(beta_norm)<-names(beta_cv)

l_bn<-c()
for (i in 1:length(beta_norm)){
  l_bn[i]<-length(beta_norm[[i]])
}
length(which(l_bn==850))

from<-c()
to<-c()
weight<-c()
beta_dfnorm<-list()
for (i in 1:length(y_scaled)){
  from<- colnames(x_g)
  to<- rep(names(y_scaled)[i], ncol(x_g))
  names(beta_norm[[i]])<-c()
  weight<-beta_norm[[i]]
  beta_dfnorm[[i]]<-data.frame(from, to, weight)
} 

beta_dfnorm<-bind_rows(beta_dfnorm)
dim(beta_dfnorm)
head(beta_dfnorm)
save(beta_dfnorm, file = 'beta_dfnorm_nopv.Rdata')
length(unique(beta_dfnorm$to))

#let'create the beta_mx!!
library(tidyr)

beta_tb<-pivot_wider(beta_dfnorm, names_from = from, values_from = weight)
dim(beta_tb)
#beta_tb[is.na(beta_tb)] <- 1
#beta_tb<-t(beta_tb)
# beta_df_allint<-beta_tb %>%
#   pivot_longer(!to, names_to = "from", values_to = "weight")
# 
# beta_df_allint<-as.data.frame(beta_df_allint)
# class(beta_df_allint)
# beta_df_allint<-data.frame(beta_df_allint[,2], beta_df_allint[,1], beta_df_allint[,3])
# colnames(beta_df_allint)<-c('from','to', 'weight')
# head(beta_df_allint)
# dim(beta_df_allint)
# save(beta_df_allint, file = "C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24/app2/beta_df_allint.Rdata")

# beta_mx<-matrix(as.numeric(beta_tb), ncol = ncol(beta_tb), dimnames = list(c(beta_tb[,1]),
#                                                                            c(colnames(beta_tb))))
rnames<-unique(beta_tb[,1])[[1]]
beta_tb<-beta_tb[,-1]

beta_mx<-apply(as.matrix.noquote(beta_tb),2,as.numeric)
rownames(beta_mx)<-rnames
save(beta_mx, file = "aabeta_mx_nopv.Rdata")
head(beta_mx)[,1:5]
dim(beta_mx)



