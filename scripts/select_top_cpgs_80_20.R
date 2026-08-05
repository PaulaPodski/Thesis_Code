#!/usr/bin/env Rscript
library(methylKit)

myDiff_train <- readRDS("myDiff_train_80_20_min12.rds")

# Filter to significant DMPs: FDR q < 0.01, |methylation difference| > 25%
sig_dmps <- getMethylDiff(myDiff_train, difference = 25, qvalue = 0.01)
sig_dmps_df <- getData(sig_dmps)
saveRDS(sig_dmps_df, "myDiff_train_80_20_min12_threshold_pass_df.rds")

# Rank by ascending qvalue, take top 100 and top 20
sig_dmps_df <- sig_dmps_df[order(sig_dmps_df$qvalue), ]

top100_df <- head(sig_dmps_df, 100)
top20_df  <- head(sig_dmps_df, 20)

saveRDS(top100_df, "myDiff_train_80_20_min12_top100_df.rds")
saveRDS(top20_df,  "myDiff_train_80_20_min12_top20_df.rds")

cat("Significant DMPs (q<0.01, |diff|>25):", nrow(sig_dmps_df), "\n")
cat("Top100 selected:", nrow(top100_df), "\n")
cat("Top20 selected:", nrow(top20_df), "\n")
