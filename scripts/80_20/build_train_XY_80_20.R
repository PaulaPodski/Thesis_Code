#!/usr/bin/env Rscript
library(methylKit)

## ---- Labels ----
split <- readRDS("sample_split_80_20_min12.rds")
y_train <- as.numeric(split$train_treatment)
saveRDS(y_train, "y_train_80_20.rds")

## ---- CpG panels ----
top100_dmps <- readRDS("myDiff_train_80_20_min12_top100_df.rds")
top100_ids <- paste(top100_dmps$chr, top100_dmps$start, top100_dmps$end, sep = ".")

top20_dmps <- readRDS("myDiff_train_80_20_min12_top20_df.rds")
top20_ids <- paste(top20_dmps$chr, top20_dmps$start, top20_dmps$end, sep = ".")

## ---- KNN-imputed X matrices ----
beta_knn <- readRDS("beta_matrix_train_80_20_min12_knn_imputed.rds")
X_knn_top100 <- t(beta_knn[rownames(beta_knn) %in% top100_ids, ])
X_knn_top20  <- t(beta_knn[rownames(beta_knn) %in% top20_ids, ])
saveRDS(X_knn_top100, "X_train_80_20_knn_imputed.rds")
saveRDS(X_knn_top20,  "X_train_80_20_knn_imputed_top20.rds")

## ---- Mean-imputed X matrices ----
beta_mean <- readRDS("beta_matrix_train_80_20_min12_mean_imputed.rds")
X_mean_top100 <- t(beta_mean[rownames(beta_mean) %in% top100_ids, ])
X_mean_top20  <- t(beta_mean[rownames(beta_mean) %in% top20_ids, ])
saveRDS(X_mean_top100, "X_train_80_20_mean_imputed.rds")
saveRDS(X_mean_top20,  "X_train_80_20_mean_imputed_top20.rds")
