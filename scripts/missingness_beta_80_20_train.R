#!/usr/bin/env Rscript

# 80/20 train split: Beta matrix missingness (before imputation)

#Loading the beta matrix object:

beta_train <- readRDS("beta_matrix_train_80_20_min12.rds")

dim(beta_train)

#Calculating the missingness of:

total_na <- sum(is.na(beta_train)) #Total missing cells
total_cells <- nrow(beta_train) * ncol(beta_train) # Sum of all cells in matrix
overall_missing_pct <- total_na / total_cells * 100 #Calculation of overall beta matrix missingness/sparsity

total_na
total_cells
round(overall_missing_pct, 2) #outputs a rounded number percentage

col_missing_pct <- colMeans(is.na(beta_train)) * 100 #Giving percentage of calculated decimal 
top10_missing <- sort(col_missing_pct, decreasing = TRUE)[1:10] #This shows the top 10 samples with the most missing data points in matrix


print(round(top10_missing, 2))

