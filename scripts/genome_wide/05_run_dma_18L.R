#!/usr/bin/env Rscript

library(methylKit)

meth_18 <- readRDS("meth_min18.rds")

myDiff_18 <- calculateDiffMeth(meth_18)

saveRDS(myDiff_18, "myDiff_min18_alive_dead.rds")
