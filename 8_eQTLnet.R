#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#' Loss Function in  P-values weighted Ridge Regression __CV
#' @author Alain Mbebi
#'
#' @param betas vector of parameter estimates
#'
#' @return A single numeric, value of the loss function.
#' @export
#'
#' @examples
pval_differentially_weighted_ridge_loss_func <- function(betas) {
  sum((y - X %*% betas)^2) + lambda * sum(TF_weights * betas^2)
}

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#' P-values weighted Ridge Regression__CV
#' This is a special case of Differentially-weighted Ridge Regression for which the weights are p-values from eGWAS
#' @param y is the vector of response variable that is the gene expression level for non-TF across accession
#' @param X is the matrix of predictor variables (TF significantly associated with y above from eGWAS). This is either from eQTL or a TF from expression matrix
#' @param lambda is the ridge regularization parameter
#' @param optim_method is the optimization method, passed to optim()
#' @param predictor_weights vector of length ncol(X) with weights corresponding to significant p-values for TF (predictor variables, X). 
#' @return A list with parameter estimates, fitted values and multiple R-squared.
#' @importFrom glmnet glmnet
#'
#' @examples

pval_differentially_weighted_ridge <- function(y, X, lambda, TF_weights, optim_method = "Nelder-Mead") {
  # Perform 10-fold cross-validation to select the best lambda 
  lambdas_to_try <- 10^seq(-3, 5, length.out = 100) ##if it takes long, you can set length.out = 50 or 30
  # Setting alpha = 0 to implement ridge regression
  ridge_cv <- cv.glmnet(X, y, alpha = 0, lambda = lambdas_to_try, standardize = TRUE, nfolds = 10)###you can see the standardized ==true I was talking about here
  
  # Best cross-validated lambda
  lambda_cv <- ridge_cv$lambda.min
  # Fit final model, get its sum of squared residuals and multiple R-squared
  #model_init <- glmnet(X, y, alpha = 0, lambda = lambda, standardize = FALSE) 
  model_cv   <- glmnet(X, y, alpha = 0, lambda = lambda_cv, standardize = TRUE)
  #betas_init <- as.vector(model_init$beta)
  betas_cv <- as.vector(model_cv$beta)
  #coef <- optim(betas_init, pval_differentially_weighted_ridge_loss_func, method = optim_method)$par
  coef <- optim(betas_cv, pval_differentially_weighted_ridge_loss_func, method = optim_method)$par
  fitted <- X %*% coef
  rsq <- cor(y, fitted)^2
  rss <- sum((fitted - y) ^ 2)
  tss <- sum((y - mean(y)) ^ 2)
  p=ncol(X)
  n=length(y)
  rsq <- 1-(rss/n)/(tss/n)
  adj_rsq <- 1-(rss/(n-p-1))/(tss/(n-1)) ##this is the adjested r2 test it also and here p is the number of TF in each model (i.e. p=ncol(X) )
  names(coef) <- colnames(X)
  output <- list("coef" = coef, "fitted" = fitted, "rsq" = rsq, "adj_rsq" = adj_rsq)
  return(output)
}
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#Here is an example run with the first target gene in your y_mg matrix, the 
#first TF matrix in your x_mg and the first vector of p-values in your pv_mg vector.
setwd("C:/Users/39334/OneDrive - Scuola Superiore Sant'Anna/Desktop/Research/Ph.D_SSSA/Work/Potsdam/GRN/eQTL_Mara/ER_GRN/Cambridge_data/exp/fresh_start/16.08.24/")
setwd('/home/e.riccucci/MAGIC_maize/eQTL/cambridge_data/exp/fresh_start/16.08.24/nopv')

load('pvalue_completenopv.Rdata')
load('regulators_completenopv.Rdata')
load('regulated_completenopv.Rdata')

# emp_vec<-c()
# for (i in 1 : length(exp_regulators)){
#   exp_regulators[[i]]<-as.data.frame(exp_regulators[[i]])
#   if (ncol(exp_regulators[[i]]) == 1){
#     emp_vec <- c(rep(0,length(exp_regulators[[i]])))
#     exp_regulators[[i]] <- cbind(emp_vec, exp_regulators[[i]][, 1, drop = FALSE])
#     #colnames(exp_regulators[[i]])[2] <- colnames(exp_regulators[[i]])
#     pvalue[[i]] <- c(0, 1/(-log10(pvalue[[i]])))
#   } else {
#     exp_regulators[[i]] <- exp_regulators[[i]]
#     pvalue[[i]] <- 1/(-log10(pvalue[[i]]))
#   }
#   exp_regulators[[i]]<-as.matrix(exp_regulators[[i]])
# }

class(exp_regulators)

# Before we used to scale everything and this is not correct (at least it conflicts with glmnet and I think this is why even with the previous data we didn't get correct results)
# THis is what you should do:
# 
# #Prep_data
# # Center y, X will be standardized in the modelling function
y_scaled <- list()
for (i in 1:length(exp_regulated)){
  y_scaled[[i]] <- scale(exp_regulated[[i]], center = TRUE, scale = FALSE)#and set it as as.matrix()
}
names(y_scaled)<-names(exp_regulated)
library(glmnet)
#lookfor constant y
id<-c()
l<-c()
for ( i in 1:length(y_scaled)){
  l[i]<-length(unique(y_scaled[[i]]))
  
}

id<-which(l==1) 
#remove constant y
y_scaled<-y_scaled[-id]
pvalue<-pvalue[-id]

test_pval_differentially_weighted_ridge <- list()
beta_gene_1_pval_weighted_ridge <- list()
rsq <- c()
adj_rsq <- c()
fitted<- list()
for (i in 1:length(y_scaled)){
  y=y_scaled[[i]][,1]
  X=as.matrix(exp_regulators)#x_g[[11]]
  lambda=.5
  TF_weights=unname(pvalue[[i]])#pv[[11]]#
  
  test_pval_differentially_weighted_ridge[[i]]  <- pval_differentially_weighted_ridge(y, X, lambda, TF_weights, optim_method = "Nelder-Mead")
  beta_gene_1_pval_weighted_ridge[[i]] = test_pval_differentially_weighted_ridge[[i]]$coef
  rsq[i] <- test_pval_differentially_weighted_ridge[[i]]$rsq
  adj_rsq[i] <- test_pval_differentially_weighted_ridge[[i]]$adj_rsq
  fitted[[i]] <- test_pval_differentially_weighted_ridge[[i]]$fitted
}

beta_gene_1_pval_weighted_ridge[[i]] 
length(beta_gene_1_pval_weighted_ridge[[14]]) 

# beta_len<-length(beta_gene_1_pval_weighted_ridge)
# for ( i in 1 : beta_len){
#   if (names(beta_gene_1_pval_weighted_ridge[[i]])[1] == 'emp_vec'){
#     beta_gene_1_pval_weighted_ridge[[i]] <- beta_gene_1_pval_weighted_ridge[[i]][-1]
#   }
# } 

exp_regulators[6]
names(exp_regulated)
names(fitted)<-names(y_scaled)
for (i in 1:length(fitted)){
  colnames(fitted[[i]])<-names(fitted[[i]])
}
names(beta_gene_1_pval_weighted_ridge)<-names(y_scaled)
names(rsq)<-names(y_scaled)
names(fitted)<-names(y_scaled)
names(adj_rsq)<-names(y_scaled)
save(beta_gene_1_pval_weighted_ridge, file = 'beta_cv_nopv.Rdata')
save(rsq, file = 'rsq_nopv.Rdata')
save(adj_rsq, file = 'adj_rsq_nopv.Rdata')
save(fitted, file = 'fitted_nopv.Rdata')

length(fitted)
