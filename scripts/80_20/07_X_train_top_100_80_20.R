#!/usr/bin/env Rscript

#80/20 train split building the top 100 feature matrices (the X)

#Loading the top 100 DMP pool

top100_dmps_80_20 <- readRDS("myDiff_train_80_20_min12_top100_df.rds")
top100_ids_80_20 <- paste(top100_dmps_80_20$chr, top100_dmps_80_20$start, top100_dmps_80_20$end, sep = ".")

#Loading the full-pool imputed training beta matrices:

beta_train_knn <- readRDS("beta_matrix_train_80_20_min12_knn_imputed.rds")
beta_train_mean <- readRDS("beta_matrix_train_80_20_min12_mean_imputed.rds")

#Subsetting the top 100 DMPs and transposing them to samples rows

X_train_knn_100 <- t(beta_train_knn[rownames(beta_train_knn) %in% top100_ids_80_20, ])
X_train_mean_100 <- t(beta_train_mean[rownames(beta_train_mean) %in% top100_ids_80_20, ])


#Saving the rds objects

saveRDS(X_train_knn_100, "X_train_80_20_knn_imputed.rds")
saveRDS(X_train_mean_100, "X_train_80_20_mean_imputed.rds")
