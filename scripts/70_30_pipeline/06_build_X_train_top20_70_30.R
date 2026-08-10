#!/usr/bin/env Rscript

library(methylKit)

#Loading the top 20 CpG pool

top20_dmps_70_30 <- readRDS("myDiff_train_70_30_min12_top20_df.rds")
top20_ids_70_30 <- paste(top20_dmps_70_30$chr, top20_dmps_70_30$start, top20_dmps_70_30$end, sep = ".")

#Loading the full imputated training KNN and mean beta matrices

beta_train_knn_70_30 <- readRDS("beta_matrix_train_70_30_min12_knn_imputed.rds")
beta_train_mean_70_30 <- readRDS("beta_matrix_train_70_30_min12_mean_imputed.rds")


#Subsetting the top 20 CpGs and transposing them (X)

X_train_knn_20_70_30 <- t(beta_train_knn_70_30[rownames(beta_train_knn_70_30) %in% top20_ids_70_30, ])
X_train_mean_20_70_30 <- t(beta_train_mean_70_30[rownames(beta_train_mean_70_30) %in% top20_ids_70_30, ])


#Saving them as rds objects:

saveRDS(X_train_knn_20_70_30, "X_train_70_30_knn_imputed_top20.rds")

saveRDS(X_train_mean_20_70_30, "X_train_70_30_mean_imputed_top20.rds")



