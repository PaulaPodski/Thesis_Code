#!/usr/bin/env Rscript

library(pROC)

#Loading the prediction test objects 
result_knn  <- readRDS("test_results_80_20_knn_top20_final.rds")
result_mean <- readRDS("test_results_80_20_mean_top20_final.rds")


#KNN ROC curve

roc_curve_knn <- roc(result_knn$y_test, as.numeric(result_knn$pred_prob))

png("roc_curve_80_20_knn_top20.png", width = 1600, height = 1400, res = 300)
plot(roc_curve_knn, col = "blue", main = "ROC Curve KNN (k = 10)", print.auc = TRUE)
abline(a = 0, b = 1, lty = 2, col = "red")
dev.off()


#Mean ROC curve

roc_curve_mean <- roc(result_mean$y_test, as.numeric(result_mean$pred_prob))

png("roc_curve_80_20_mean_top20.png", width = 1600, height = 1400, res = 300)
plot(roc_curve_mean, col = "blue", main = "ROC Curve Mean", print.auc = TRUE)
abline(a = 0, b = 1, lty = 2, col = "red")
dev.off()



