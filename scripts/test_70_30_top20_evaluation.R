#!/usr/bin/env Rscript

#70/30 test split evaluation of top 20 DMPs at train test derived alphas for KNN and mean imputation


library(glmnet)
library(pROC)


#Loading the data and firring the final models with alphas

X_train_knn_20_70_30  <- readRDS("X_train_70_30_knn_imputed_top20.rds")
X_train_mean_20_70_30 <- readRDS("X_train_70_30_mean_imputed_top20.rds")
y_train_70_30          <- readRDS("y_train_70_30.rds")
y_train_70_30_factor   <- factor(y_train_70_30, levels = c(0, 1))

default_fold <- 10
set.seed(42) #reprorucible
fold_assigned <- sample(rep(1:default_fold, length.out = length(y_train_70_30_factor)))

final_knn_20_70_30 <- cv.glmnet(x = X_train_knn_20_70_30, y = y_train_70_30_factor, alpha = 0.70, foldid = fold_assigned, family = "binomial", type.measure = "deviance")

final_mean_20_70_30 <- cv.glmnet(x = X_train_mean_20_70_30, y = y_train_70_30_factor, alpha = 0.20, foldid = fold_assigned, family = "binomial", type.measure = "deviance")

#saving these as rds objects 

saveRDS(final_knn_20_70_30, "final_eln_model_70_30_knn_imputed_top20_alpha070_FINAL.rds")
saveRDS(final_mean_20_70_30, "final_eln_model_70_30_mean_imputed_top20_alpha020_FINAL.rds")

#Loading the test data and subsetting to the available DMPs

split_70_30 <- readRDS("sample_split_70_30_min12.rds")
y_test_70_30 <- split_70_30$test_treatment
names(y_test_70_30) <- split_70_30$test_ids

beta_test_knn  <- readRDS("beta_matrix_test_70_30_top20_min4L_knn_imputed.rds")
beta_test_mean <- readRDS("beta_matrix_test_70_30_top20_min4L_mean_imputed.rds")

top20_dmps_70_30 <- readRDS("myDiff_train_70_30_min12_top20_df.rds")
top20_ids_70_30 <- paste(top20_dmps_70_30$chr, top20_dmps_70_30$start, top20_dmps_70_30$end, sep = ".")

overlap_ids_70_30 <- intersect(top20_ids_70_30, rownames(beta_test_knn))

X_test_knn_17_70_30  <- t(beta_test_knn[overlap_ids_70_30, , drop = FALSE])
X_test_mean_17_70_30 <- t(beta_test_mean[overlap_ids_70_30, , drop = FALSE])

y_test_70_30_aligned <- y_test_70_30[rownames(X_test_knn_17_70_30)]

#Here we are filling up the missing columns (because glmnet cannot work with 3 missing columns i.e. CpGs so we are adding the columns with traininf set means (Apparently this is NOT data leakage!)

model_cpgs_knn_70_30  <- rownames(coef(final_knn_20_70_30,  s = "lambda.min"))[-1]
model_cpgs_mean_70_30 <- rownames(coef(final_mean_20_70_30, s = "lambda.min"))[-1]

missing_knn_70_30  <- setdiff(model_cpgs_knn_70_30,  colnames(X_test_knn_17_70_30))
missing_mean_70_30 <- setdiff(model_cpgs_mean_70_30, colnames(X_test_mean_17_70_30))

X_test_knn_full_70_30 <- X_test_knn_17_70_30
for (cpg in missing_knn_70_30) {
  X_test_knn_full_70_30 <- cbind(X_test_knn_full_70_30, mean(X_train_knn_20_70_30[, cpg]))
  colnames(X_test_knn_full_70_30)[ncol(X_test_knn_full_70_30)] <- cpg
}
X_test_knn_full_70_30 <- X_test_knn_full_70_30[, model_cpgs_knn_70_30]

X_test_mean_full_70_30 <- X_test_mean_17_70_30
for (cpg in missing_mean_70_30) {
  X_test_mean_full_70_30 <- cbind(X_test_mean_full_70_30, mean(X_train_mean_20_70_30[, cpg]))
  colnames(X_test_mean_full_70_30)[ncol(X_test_mean_full_70_30)] <- cpg
}
X_test_mean_full_70_30 <- X_test_mean_full_70_30[, model_cpgs_mean_70_30]

#Predicting the test data now that there are 20 columns (no missing columns)

pred_prob_knn_20_70_30  <- predict(final_knn_20_70_30,  newx = X_test_knn_full_70_30,  s = "lambda.min", type = "response")
pred_prob_mean_20_70_30 <- predict(final_mean_20_70_30, newx = X_test_mean_full_70_30, s = "lambda.min", type = "response")

#Evaluating the model by extracting the AUC:

roc_knn_20_70_30  <- roc(response = y_test_70_30_aligned, predictor = as.numeric(pred_prob_knn_20_70_30))
roc_mean_20_70_30 <- roc(response = y_test_70_30_aligned, predictor = as.numeric(pred_prob_mean_20_70_30))

#Calculating the deviance using glmnet's assess.glmnet function 

dev_knn_70_30  <- assess.glmnet(final_knn_20_70_30,  newx = X_test_knn_full_70_30,  newy = y_test_70_30_aligned, family = "binomial", s = "lambda.min")$deviance

dev_mean_70_30 <- assess.glmnet(final_mean_20_70_30, newx = X_test_mean_full_70_30, newy = y_test_70_30_aligned, family = "binomial", s = "lambda.min")$deviance

#Saving results as objects

saveRDS(list(pred_prob = pred_prob_knn_20_70_30, roc = roc_knn_20_70_30, auc = auc(roc_knn_20_70_30), deviance = dev_knn_70_30, y_test = y_test_70_30_aligned), "test_results_70_30_knn_top20_final.rds")

saveRDS(list(pred_prob = pred_prob_mean_20_70_30, roc = roc_mean_20_70_30, auc = auc(roc_mean_20_70_30), deviance = dev_mean_70_30, y_test = y_test_70_30_aligned), "test_results_70_30_mean_top20_final.rds")

