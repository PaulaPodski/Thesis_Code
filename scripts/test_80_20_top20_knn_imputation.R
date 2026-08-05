#!/usr/bin/env Rscript
library(impute)

beta_test <- readRDS("beta_matrix_test_80_20_top20_min2L.rds")

chunk_size <- 500000
n <- nrow(beta_test)
starts <- seq(1, n, by = chunk_size)

beta_knn_imputed <- beta_test

for (s in starts) {
  e <- min(s + chunk_size - 1, n)
  chunk <- beta_test[s:e, , drop = FALSE]
  result <- impute.knn(as.matrix(chunk), k = 10, rowmax = 0.5, colmax = 0.9999, maxp = chunk_size + 1, rng.seed = 362436069)
  beta_knn_imputed[s:e, ] <- result$data
}

saveRDS(beta_knn_imputed, "beta_matrix_test_80_20_top20_min2L_knn_imputed.rds")
