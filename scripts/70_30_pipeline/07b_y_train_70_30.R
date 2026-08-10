#!/usr/bin/env Rscript

#70/30 train split building the label vector (y)

#Load the sample split

split <- readRDS("sample_split_70_30_min12.rds")

#Extract the training labels directly (0 = alive, 1 = dead)

y_train_70_30 <- split$train_treatment

#Save it as an rds object 

saveRDS(y_train_70_30, "y_train_70_30.rds")
