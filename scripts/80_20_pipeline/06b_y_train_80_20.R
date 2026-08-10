#!/usr/bin/env Rscript

# Building the training label vector (y):

#Loading the sample split

split <- readRDS("sample_split_80_20_min12.rds")

#Extracting the training labels directly (0 = alive, 1 = dead)

y_train_80_20 <- split$train_treatment

#Saving it as an object 

saveRDS(y_train_80_20, "y_train_80_20.rds")

length(y_train_80_20)

print(table(y_train_80_20)) #Shows class balance
