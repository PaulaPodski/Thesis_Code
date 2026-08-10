#!/usr/bin/env Rscript

# 70/30 train split building the top 100 feature matrices (the X)

# Loading the top 100 DMP pool

top100_dmps_70_30 <- readRDS("myDiff_train_70_30_min12_top100_df.rds")
top100_ids_70_30 <- paste(top100_dmps_70_30$chr, top100_dmps_70_30$start, top100_dmps_70_30$end, sep = ".")

#Loading the fully imputated KNN and mean imputated beta matrices

beta_train_knn <- readRDS("beta_matrix_train_70_30_min12_knn_imputed.rds")
beta_train_mean <- readRDS("beta_matrix_train_70_30_min12_mean_imputed.rds")

#Subsetting the top 100 DMPs and transposing

X_train_knn_100 <- t(beta_train_knn[rownames(beta_train_knn) %in% top100_ids_70_30, ])
X_train_mean_100 <- t(beta_train_mean[rownames(beta_train_mean) %in% top100_ids_70_30, ])


#Saving them as rds objects:


saveRDS(X_train_knn_100, "X_train_70_30_knn_imputed.rds")
saveRDS(X_train_mean_100, "X_train_70_30_mean_imputed.rds")


