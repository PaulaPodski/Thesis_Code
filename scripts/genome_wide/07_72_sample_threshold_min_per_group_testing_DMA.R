#!/usr/bin/env Rscript

#Here on the full 72 sample cohort this script tests the genome-wide DMA across different coverage (min.per.group) thresholds 

#The thresholds tested were  3,9,12,18 and 24 (which returned 0 signficant DMPs)

library(methylKit)

thresholds <- c(3, 9, 12, 18, 24)
myDiff_files <- c("myDiff_min3_alive_dead_fixed.rds", "myDiff_min9_alive_dead.rds", "myDiff_min12_alive_dead.rds", "myDiff_min18_alive_dead.rds", "myDiff_min24_alive_dead.rds") 


#for loop through all specified thresholds and calculates the total number of significant DMPs, hyper and hypomethylated DMPs that pass the required threshold of +/- 25% methylation differencr and FDR < 0.01

for (i in seq_along(thresholds)) {
  myDiff <- readRDS(myDiff_files[i])
  sig <- tryCatch({
    getMethylDiff(myDiff, difference = 25, qvalue = 0.01)
  }, error = function(e) NULL)

  if (is.null(sig)) {
    cat("min.per.group =", thresholds[i], "| Hyper: 0 | Hypo: 0 | Total sig: 0 (no significant DMPs)\n")
  } else {
    sig_df <- getData(sig)
    cat("min.per.group =", thresholds[i],
        "Hypermethylated DMPs:", sum(sig_df$meth.diff > 0),
        "Hypomethylated DMPs:", sum(sig_df$meth.diff < 0),
        "Total signigficant DMPs:", nrow(sig_df), "\n")
  }
}
