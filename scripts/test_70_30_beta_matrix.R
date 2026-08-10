#!/usr/bin/env Rscript
library(methylKit)

split          <- readRDS("sample_split_70_30_min12.rds")
test_file.list <- readRDS("test_file_list_70_30_min12.rds")

test_ids       <- split$test_ids
test_treatment <- split$test_treatment

stopifnot(length(test_file.list) == length(test_ids),
          length(test_ids) == 22,
          sum(test_treatment) == 11, sum(test_treatment == 0) == 11)

myobj_test <- methRead(
  test_file.list,
  sample.id = as.list(test_ids),
  assembly = "hg38",
  treatment = test_treatment,
  context = "CpG",
  pipeline = "bismarkCoverage",
  mincov = 10,
  dbtype = "tabix",
  dbdir = "methylDB_test_70_30_min12"
)

filtered_test   <- filterByCoverage(myobj_test, lo.count = 10, lo.perc = NULL,
                                     hi.count = NULL, hi.perc = 99.9)
normalized_test <- normalizeCoverage(filtered_test)
saveRDS(normalized_test, "normalized_test_70_30_min12.rds")

meth_test <- unite(normalized_test, destrand = FALSE, min.per.group = 12L)
saveRDS(meth_test, "meth_test_70_30_min12.rds")

beta_test <- percMethylation(meth_test, rowids = TRUE) / 100
saveRDS(beta_test, "beta_matrix_test_70_30_min12.rds")
