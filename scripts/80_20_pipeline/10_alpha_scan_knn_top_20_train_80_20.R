#!/usr/bin/env Rscript

#80/20 train split ELN alpha scan for KNN imputated beta matrix (Top 20)

#Loading the relevant packages:

library(glmnet)
library(pROC)

#Loading the top 20 feature matrix (X) and outcome (y)

X_train_knn_20_80_20 <- readRDS("X_train_80_20_knn_imputed_top20.rds")

y_train_80_20 <- readRDS("y_train_80_20.rds")

#converting the y into a factor required for glmnet package for binomial models

y_train_80_20_factor <- factor(y_train_80_20, levels = c(0, 1))

#Cross validation glmnet folds - glmnet package uses 10 fold as default (folds generated once and then reused for every alpha value so all alpha values are evaluated on identical splits)

default_fold <- 10

set.seed(42) #for reproducibility

fold_assigned <- sample(rep(1:default_fold, length.out = length(y_train_80_20_factor)))

#Defining the alpha values to be scanned (0 through 1 in 0.05 increments)

alpha_table <- seq(0, 1, by = 0.05)

#creating function fitting the model for one alpha value

fit_one_alpha <- function(a, feature_matrix) {cv.glmnet(x = feature_matrix, y = y_train_80_20_factor, alpha = a, foldid = fold_assigned, family = "binomial", type.measure = "deviance")}

#Now fitting model for all alpha values. The lapply() function runs fit_one_alpha() once per alpha (21 times)
#Then putting it into a list

all_models_knn_20 <- lapply(alpha_table, fit_one_alpha, feature_matrix = X_train_knn_20_80_20)

#Extracting the AUC, non-zero coefficient CpG features (no. features retained by model per alpha value) and the min. cross validated deviance for each model:

auc_knn_20 <- sapply(all_models_knn_20, function(fit) {p <- as.numeric(predict(fit, newx = X_train_knn_20_80_20, s = "lambda.min", type = "response")); as.numeric(auc(roc(response = y_train_80_20_factor, predictor = p)))})

n_retained_knn_20 <- sapply(all_models_knn_20, function(fit) { sum(as.matrix(coef(fit, s = "lambda.min"))[-1, 1] != 0)})

deviance_knn_20 <- sapply(all_models_knn_20, function(fit) min(fit$cvm))

#Creating the final results table:

results_knn_20 <- data.frame(alpha = alpha_table, AUC = auc_knn_20, no_CpGs_retained = n_retained_knn_20, deviance = deviance_knn_20)

print(results_knn_20)

#Plotting

png("alpha_vs_apparent_AUC_knn_top20.png", width = 1800, height = 1400, res = 300)

plot(results_knn_20$alpha, results_knn_20$AUC, type = "b", pch = 19, col = "blue", xlab = "Alpha (0 = Ridge, 1 = Lasso)", ylab = "Apparent train AUC at lambda.min", main = "80/20 KNN Top 20: Alpha vs AUC", ylim = c(0.5, 1.02))

abline(h = 0.5, lty = 2, col = "grey")

dev.off()

png("alpha_scan_80_20_knn_top20_AUC.png", width = 1800, height = 1400, res = 300)
plot(results_knn_20$alpha, results_knn_20$AUC, type = "b", pch = 19, col = "blue",
     xlab = "alpha value", ylab = "AUC (lambda.min)",
     main = "Elastic Net Alpha Scan, 80/20 Train Top 20 (KNN Imputed)")
dev.off()

png("alpha_scan_80_20_knn_top20_AUC.png", width = 1800, height = 1400, res = 300)
plot(results_knn_20$alpha, results_knn_20$AUC, type = "b", pch = 19, col = "black",
     xlab = "alpha value", ylab = "AUC (lambda.min)",
     main = "Elastic Net Alpha Scan, 80/20 Train Top 20 (KNN Imputed)",
     cex.main = 0.8)
dev.off()
