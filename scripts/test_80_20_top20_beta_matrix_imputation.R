#!/usr/bin/env Rscript
library(methylKit)
library(impute)

# ---- Step 1: Rebuild test beta matrix at min.per.group=2 ----
normalized_test_80_20 <- readRDS("normalized_test_80_20_min12.rds")

meth_test_80_20 <- unite(normalized_test_80_20, destrand = FALSE, min.per.group = 2L)
saveRDS(meth_test_80_20, "meth_test_80_20_top20_min2L.rds")

beta_test_80_20 <- percMethylation(meth_test_80_20, rowids = TRUE) / 100
saveRDS(beta_test_80_20, "beta_matrix_test_80_20_top20_min2L.rds")

top20_dmps_80_20 <- readRDS("myDiff_train_80_20_min12_top20_df.rds")
top20_ids_80_20 <- paste(top20_dmps_80_20$chr, top20_dmps_80_20$start, top20_dmps_80_20$end, sep = ".")
overlap <- sum(top20_ids_80_20 %in% rownames(beta_test_80_20))
cat("Overlap with top20:", overlap, "\n")
cat("beta_test dims:", dim(beta_test_80_20), "\n")

# ---- Step 2: Mean imputation ----
row_means <- rowMeans(beta_test_80_20, na.rm = TRUE)
beta_test_mean_imputed <- beta_test_80_20
na_idx <- which(is.na(beta_test_mean_imputed), arr.ind = TRUE)
beta_test_mean_imputed[na_idx] <- row_means[na_idx[, 1]]
saveRDS(beta_test_mean_imputed, "beta_matrix_test_80_20_top20_min2L_mean_imputed.rds")
cat("Mean imputation done. Remaining NAs:", sum(is.na(beta_test_mean_imputed)), "\n")

# ---- Step 3: KNN imputation ----
knn_result <- impute.knn(as.matrix(beta_test_80_20), k = 10, rowmax = 0.5, colmax = 0.9999, maxp = nrow(beta_test_80_20), rng.seed = 362436069)
beta_test_knn_imputed <- knn_result$data
saveRDS(beta_test_knn_imputed, "beta_matrix_test_80_20_top20_min2L_knn_imputed.rds")
cat("KNN imputation done. Remaining NAs:", sum(is.na(beta_test_knn_imputed)), "\n")

cat("80/20 test pipeline complete.\n")
