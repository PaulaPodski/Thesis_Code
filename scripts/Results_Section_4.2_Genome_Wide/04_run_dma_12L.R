#!/usr/bin/env Rscript

library(methylKit)

meth_12 <- readRDS("meth_min12.rds")

myDiff_12 <- calculateDiffMeth(meth_12)

saveRDS(myDiff_12, "myDiff_min12_alive_dead.rds")
