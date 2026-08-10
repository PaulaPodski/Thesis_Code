#!/usr/bin/env Rscript

library(methylKit)

meth_24 <- readRDS("meth_min24.rds")

myDiff_24 <- calculateDiffMeth(meth_24)

saveRDS(myDiff_24, "myDiff_min24_alive_dead.rds")
