#!/usr/bin/env Rscript

library(methylKit)

meth_3 <- readRDS("meth_min3.rds")

myDiff_3 <- calculateDiffMeth(meth_3)

saveRDS(myDiff_3, "myDiff_min3_alive_dead_fixed.rds")
