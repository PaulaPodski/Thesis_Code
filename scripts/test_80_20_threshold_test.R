#!/usr/bin/env Rscript
library(methylKit)
normalized_test <- readRDS("normalized_test_80_20_min12.rds")
top100_dmps <- readRDS("myDiff_train_80_20_min12_top100_df.rds")
top100_ids <- paste(top100_dmps$chr, top100_dmps$start, top100_dmps$end, sep = ".")
smallest_group <- 7
results <- data.frame(
  min_per_group = integer(),
  pct_of_smallest_group = numeric(),
  n_overlap = integer(),
  n_total_cpgs = integer()
)
for (m in 1:smallest_group) {
  cat("Testing min.per.group =", m, "\n")
  result_row <- tryCatch({
    meth_test_m <- unite(normalized_test, destrand = FALSE, min.per.group = m)
    beta_test_m <- percMethylation(meth_test_m, rowids = TRUE) / 100
    overlap <- sum(top100_ids %in% rownames(beta_test_m))
    n_total <- nrow(beta_test_m)
    rm(meth_test_m, beta_test_m)
    gc()
    data.frame(
      min_per_group = m,
      pct_of_smallest_group = round(m / smallest_group * 100, 1),
      n_overlap = overlap,
      n_total_cpgs = n_total
    )
  }, error = function(e) {
    cat("min.per.group =", m, "FAILED:", conditionMessage(e), "\n")
    data.frame(
      min_per_group = m,
      pct_of_smallest_group = round(m / smallest_group * 100, 1),
      n_overlap = NA,
      n_total_cpgs = NA
    )
  })
  results <- rbind(results, result_row)
  saveRDS(results, "threshold_test_80_20_partial.rds")
  write.csv(results, "threshold_test_80_20_partial.csv", row.names = FALSE)
}
print(results)
saveRDS(results, "threshold_test_80_20.rds")
write.csv(results, "threshold_test_80_20.csv", row.names = FALSE)
