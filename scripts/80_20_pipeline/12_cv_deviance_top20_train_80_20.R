#!/usr/bin/env Rscript

library(glmnet)

#Loading the data 

X_train_knn_20_80_20  <- readRDS("X_train_80_20_knn_imputed_top20.rds")
X_train_mean_20_80_20 <- readRDS("X_train_80_20_mean_imputed_top20.rds")
y_train_80_20          <- readRDS("y_train_80_20.rds")
y_train_80_20_factor   <- factor(y_train_80_20, levels = c(0, 1))


default_fold <- 10
set.seed(42)
fold_assigned <- sample(rep(1:default_fold, length.out = length(y_train_80_20_factor)))

# Fitting the  final models at chosen alphas (0.2 and 0.8)

final_knn_20_80_20 <- cv.glmnet(x = X_train_knn_20_80_20, y = y_train_80_20_factor,
                                  alpha = 0.20, foldid = fold_assigned,
                                  family = "binomial", type.measure = "deviance")

final_mean_20_80_20 <- cv.glmnet(x = X_train_mean_20_80_20, y = y_train_80_20_factor,
                                   alpha = 0.80, foldid = fold_assigned,
                                   family = "binomial", type.measure = "deviance")


png("cv_deviance_80_20_top20_KNN_and_Mean.png", width = 3200, height = 1800, res = 200)

par(mfrow = c(1, 2), mar = c(5, 5, 8, 2))

plot(final_mean_20_80_20,
     main = "80/20 Train Top 20, Mean Imputed (alpha=0.80): CV Deviance")

plot(final_knn_20_80_20,
     main = "80/20 Train Top 20, KNN Imputed k=10 (alpha=0.20): CV Deviance")

dev.off()

