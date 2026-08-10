#!/usr/bin/env Rscript

library(methylKit)

meth_9 <- readRDS("meth_min9.rds")

myDiff_9 <- calculateDiffMeth(meth_9)

saveRDS(myDiff_9, "myDiff_min9_alive_dead.rds")
