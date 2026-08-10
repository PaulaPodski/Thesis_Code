#!/usr/bin/env Rscript
library(methylKit)

normalized_test <- readRDS("normalized_test_80_20_min12.rds")
top100_dmps <- readRDS("myDiff_train_80_20_min12_top100_df.rds")
top100_ids <- paste(top100_dmps$chr, top100_dmps$start, top100_dmps$end, sep = ".")

meth_test <- unite(normalized_test, destrand = FALSE, min.per.group = 3L)
saveRDS(meth_test, "meth_test_80_20_min12.rds")

beta_test <- percMethylation(meth_test, rowids = TRUE) / 100
saveRDS(beta_test, "beta_matrix_test_80_20_min12.rds")

overlap <- sum(top100_ids %in% rownames(beta_test))
cat("Overlap with top100 at min.per.group=3L:", overlap, "\n")
cat("beta_test dims:", dim(beta_test), "\n")
