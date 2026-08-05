#!/usr/bin/env Rscript
library(methylKit)

normalized_test <- readRDS("normalized_test_80_20_min12.rds")
top100_dmps <- readRDS("myDiff_train_80_20_min12_top100_df.rds")
top100_ids <- paste(top100_dmps$chr, top100_dmps$start, top100_dmps$end, sep = ".")

test_meth_all <- unite(normalized_test, destrand = FALSE, min.per.group = 1L)
test_all_df <- getData(test_meth_all)
test_all_ids <- paste(test_all_df$chr, test_all_df$start, test_all_df$end, sep = ".")

overlap_count <- sum(top100_ids %in% test_all_ids)
cat("Number of top100 CpGs present in test at ANY coverage (min.per.group=1L):", overlap_count, "\n")

saveRDS(test_all_ids, "test_80_20_all_cpg_ids_min1L.rds")
