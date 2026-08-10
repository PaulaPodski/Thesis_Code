#!/usr/bin/env Rscript

library(methylKit)
library(impute)

beta_train <- readRDS("beta_matrix_train_70_30_min12.rds")

knn_result <- impute.knn(as.matrix(beta_train), k = 10, rowmax = 0.5, colmax = 0.9999, maxp = 20000, rng.seed = 362436069) # apart from colmax and k chosen, the other params are default 
beta_knn_imputed <- knn_result$data

saveRDS(beta_knn_imputed, "beta_matrix_train_70_30_min12_knn_imputed.rds")
