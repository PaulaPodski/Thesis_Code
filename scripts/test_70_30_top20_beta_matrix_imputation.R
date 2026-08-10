#!/usr/bin/env Rscript

library(methylKit)
library(impute)

#building the test beta matrix at min.per.group=4L

normalized_test_70_30 <- readRDS("normalized_test_70_30_min12.rds")

meth_test_70_30 <- unite(normalized_test_70_30, destrand = FALSE, min.per.group = 4L)
saveRDS(meth_test_70_30, "meth_test_70_30_top20_min4L.rds")

beta_test_70_30 <- percMethylation(meth_test_70_30, rowids = TRUE) / 100
saveRDS(beta_test_70_30, "beta_matrix_test_70_30_top20_min4L.rds")

top20_dmps_70_30 <- readRDS("myDiff_train_70_30_min12_top20_df.rds")
top20_ids_70_30 <- paste(top20_dmps_70_30$chr, top20_dmps_70_30$start, top20_dmps_70_30$end, sep = ".")
overlap <- sum(top20_ids_70_30 %in% rownames(beta_test_70_30))

print(overlap)

dim(beta_test_70_30)

#mean imputation

row_means <- rowMeans(beta_test_70_30, na.rm = TRUE)

beta_test_mean_imputed <- beta_test_70_30

na_idx <- which(is.na(beta_test_mean_imputed), arr.ind = TRUE)

beta_test_mean_imputed[na_idx] <- row_means[na_idx[, 1]]

#Saving it as rds object 

saveRDS(beta_test_mean_imputed, "beta_matrix_test_70_30_top20_min4L_mean_imputed.rds")

#KNN imputation (k = 10)

knn_result <- impute.knn(as.matrix(beta_test_70_30), k = 10, rowmax = 0.5, colmax = 0.9999, maxp = nrow(beta_test_70_30), rng.seed = 362436069)

beta_test_knn_imputed <- knn_result$data

#Saving it as rds object

saveRDS(beta_test_knn_imputed, "beta_matrix_test_70_30_top20_min4L_knn_imputed.rds")

