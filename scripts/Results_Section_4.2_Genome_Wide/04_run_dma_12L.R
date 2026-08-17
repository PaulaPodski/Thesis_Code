#!/usr/bin/env Rscript

library(methylKit)

#Firstly we will load the normalised methylKit object:

normalized.myobj <- readRDS("normalized.myobj.rds")

#Then I create and save the 12L united rds object:

meth_12 <- unite(normalized.myobj, destrand = FALSE, min.per.group = 12L)

saveRDS(meth_12, "meth_min12.rds")

#Finally,I run the DMA on the united 12L rds object:

meth_12 <- readRDS("meth_min12.rds")

myDiff_12 <- calculateDiffMeth(meth_12)

saveRDS(myDiff_12, "myDiff_min12_alive_dead.rds")
