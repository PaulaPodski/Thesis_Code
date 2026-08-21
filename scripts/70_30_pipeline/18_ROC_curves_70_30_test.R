#!/usr/bin/env Rscript

#70/30 test split ROC curves on the top 20 KNN and Mean imputations 

library(pROC)

#Loading the rds objects from test predictions

result_knn  <- readRDS("test_results_70_30_knn_top20_final.rds")
result_mean <- readRDS("test_results_70_30_mean_top20_final.rds")

#Plotting ROC curve for KNN imputated

roc_curve_knn <- roc(result_knn$y_test, as.numeric(result_knn$pred_prob))

png("roc_curve_70_30_knn_top20.png", width = 1600, height = 1400, res = 300)
plot(roc_curve_knn, col = "blue", main = "ROC Curve: 70/30 Split, KNN Imputation", print.auc = TRUE)
abline(a = 0, b = 1, lty = 2, col = "red")

dev.off()

#Plotting ROC curve for mean impuated 

roc_curve_mean <- roc(result_mean$y_test, as.numeric(result_mean$pred_prob))

png("roc_curve_70_30_mean_top20.png", width = 1600, height = 1400, res = 300)
plot(roc_curve_mean, col = "blue", main = "ROC Curve: 70/30 Split, Mean Imputation", print.auc = TRUE)
abline(a = 0, b = 1, lty = 2, col = "red")

dev.off()

