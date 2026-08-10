#!/usr/bin/env Rscript
library(methylKit)

normalized_test <- readRDS("normalized_test_70_30_min12.rds")
top20_dmps <- readRDS("myDiff_train_70_30_min12_top20_df.rds")
top20_ids <- paste(top20_dmps$chr, top20_dmps$start, top20_dmps$end, sep = ".")

results <- data.frame(min_per_group = integer(), n_overlap = integer(), n_total_cpgs = integer())

for (m in 11:1) {
  cat("Testing min.per.group =", m, "\n")
  
  result_row <- tryCatch({
    meth_test_m <- unite(normalized_test, destrand = FALSE, min.per.group = m)
    beta_test_m <- percMethylation(meth_test_m, rowids = TRUE) / 100
    overlap <- sum(top20_ids %in% rownames(beta_test_m))
    n_total <- nrow(beta_test_m)
    
    rm(meth_test_m, beta_test_m)
    gc()
    
    data.frame(min_per_group = m, n_overlap = overlap, n_total_cpgs = n_total)
  }, error = function(e) {
    cat("min.per.group =", m, "FAILED:", conditionMessage(e), "\n")
    data.frame(min_per_group = m, n_overlap = NA, n_total_cpgs = NA)
  })
  
  results <- rbind(results, result_row)
  saveRDS(results, "top20_thresholds_70_30_partial.rds")
}

print(results)
saveRDS(results, "top20_thresholds_70_30.rds")

