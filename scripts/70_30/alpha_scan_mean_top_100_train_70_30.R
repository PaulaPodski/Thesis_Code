#!/usr/bin/env Rscript

#70/30 train split ELN alpha scan for the Mean imputated beta-matrix (Top 100)

library(glmnet)
library(pROC)

X_train_mean_100_70_30 <- readRDS("X_train_70_30_mean_imputed.rds")
y_train_70_30 <- readRDS("y_train_70_30.rds")
y_train_70_30_factor <- factor(y_train_70_30, levels = c(0, 1))

default_fold <- 10
set.seed(42)
fold_assigned <- sample(rep(1:default_fold, length.out = length(y_train_70_30_factor)))

alpha_table <- seq(0, 1, by = 0.05)

fit_one_alpha <- function(a, feature_matrix) {cv.glmnet(x = feature_matrix, y = y_train_70_30_factor, alpha = a, foldid = fold_assigned, family = "binomial", type.measure = "deviance")}

all_models_mean_100 <- lapply(alpha_table, fit_one_alpha, feature_matrix = X_train_mean_100_70_30)

auc_mean_100 <- sapply(all_models_mean_100, function(fit) {p <- as.numeric(predict(fit, newx = X_train_mean_100_70_30, s = "lambda.min", type = "response")); as.numeric(auc(roc(response = y_train_70_30_factor, predictor = p)))})

n_retained_mean_100 <- sapply(all_models_mean_100, function(fit) { sum(as.matrix(coef(fit, s = "lambda.min"))[-1, 1] != 0)})

deviance_mean_100 <- sapply(all_models_mean_100, function(fit) min(fit$cvm))

results_mean_100 <- data.frame(alpha = alpha_table, AUC = auc_mean_100, no_CpGs_retained = n_retained_mean_100, deviance = deviance_mean_100)

print(results_mean_100)

saveRDS(results_mean_100, "alpha_scan_70_30_mean_top100_fixed_folds.rds")
