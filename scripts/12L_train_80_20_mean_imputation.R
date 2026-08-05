#!/usr/bin/env Rscript

library(methylKit)

beta_train <- readRDS("beta_matrix_train_80_20_min12.rds")

row_means <- rowMeans(beta_train, na.rm = TRUE)
beta_mean_imputed <- beta_train

na_idx <- which(is.na(beta_mean_imputed), arr.ind = TRUE)
beta_mean_imputed[na_idx] <- row_means[na_idx[, 1]]

saveRDS(beta_mean_imputed, "beta_matrix_train_80_20_min12_mean_imputed.rds")

