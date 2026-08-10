#!/usr/bin/env Rscript

library(methylKit)

normalized_test <- readRDS("normalized_test_80_20_min12.rds")

meth_test <- unite(normalized_test, destrand = FALSE, min.per.group = 5L)
saveRDS(meth_test, "meth_test_80_20_min12.rds")

beta_test <- percMethylation(meth_test, rowids = TRUE) / 100
saveRDS(beta_test, "beta_matrix_test_80_20_min12.rds")
