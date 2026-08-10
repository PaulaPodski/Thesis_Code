#!/usr/bin/env Rscript

# 70/30 train split samples: DMA significance filtering and top DMP ranking

library(methylKit)

# Loading DMA results rds object 
myDiff_train <- readRDS("myDiff_train_70_30_min12.rds")

# Filtering significant DMPs: FDR q < 0.01 and methylation difference +/- 25%     
sig_dmps <- getMethylDiff(myDiff_train, difference = 25, qvalue = 0.01)
sig_dmps_df <- getData(sig_dmps)

nrow(sig_dmps_df)               # significant DMPs
sum(sig_dmps_df$meth.diff > 0)  # Hypermethylated DMPs
sum(sig_dmps_df$meth.diff < 0)  # Hypomethylated DMPs

saveRDS(sig_dmps_df, "myDiff_train_70_30_min12_threshold_pass_df.rds")

# Building the CpG ID and rank by q value (most to least significant):    

sig_dmps_df$cpg_id <- paste(sig_dmps_df$chr, sig_dmps_df$start, sig_dmps_df$end, sep = ".")
sig_dmps_ranked <- sig_dmps_df[order(sig_dmps_df$qvalue), ]

# Extracting the top 100 and top 20 CpGs:

top100_df <- head(sig_dmps_ranked, 100)
top20_df  <- head(sig_dmps_ranked, 20)

nrow(top100_df)
nrow(top20_df)

saveRDS(top100_df, "myDiff_train_70_30_min12_top100_df.rds")
saveRDS(top20_df, "myDiff_train_70_30_min12_top20_df.rds")
