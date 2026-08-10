#!/usr/bin/env Rscript

#70/30 train split CV deviance plots (Top 20 final alphas for KNN and mean)

library(glmnet)

#load the data

X_train_knn_20_70_30  <- readRDS("X_train_70_30_knn_imputed_top20.rds")
X_train_mean_20_70_30 <- readRDS("X_train_70_30_mean_imputed_top20.rds")
y_train_70_30          <- readRDS("y_train_70_30.rds")
y_train_70_30_factor   <- factor(y_train_70_30, levels = c(0, 1))


#setting seed 

default_fold <- 10
set.seed(42)
fold_assigned <- sample(rep(1:default_fold, length.out = length(y_train_70_30_factor)))


#fitting the models at alphas 0.7 and 0.2 knn and mean respectively

final_knn_20_70_30 <- cv.glmnet(x = X_train_knn_20_70_30, y = y_train_70_30_factor, alpha = 0.70, foldid = fold_assigned, family = "binomial", type.measure = "deviance")

final_mean_20_70_30 <- cv.glmnet(x = X_train_mean_20_70_30, y = y_train_70_30_factor, alpha = 0.20, foldid = fold_assigned,family = "binomial", type.measure = "deviance")

#Plotting the CV deviance curves 

png("cv_deviance_70_30_top20_KNN_and_Mean.png", width = 3200, height = 1800, res = 200)

par(mfrow = c(1, 2), mar = c(5, 5, 8, 2))

plot(final_mean_20_70_30,
     main = "70/30 Train Top 20, Mean Imputed (alpha=0.20): CV Deviance")

plot(final_knn_20_70_30,
     main = "70/30 Train Top 20, KNN Imputed k=10 (alpha=0.70): CV Deviance")

dev.off()




