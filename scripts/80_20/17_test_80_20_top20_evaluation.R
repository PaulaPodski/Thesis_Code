#!/usr/bin/env Rscript

#This script evaluates the trained split from 80/20 on final chosen alpha values of 0.2 from KNN and 0.8 from mean imputation
#There are 3 missing CpGs as earlier only 17/20 were retained from min.per.group
#There is a method of adding 3 columns so it is back to 20 (as glmnet would fail to run predict() on 17 CpGs when it was trained on 20.
#The missing columns are filled with training set means (NOT DATA LEAKAGE) because only the final model will use the 17 CpGs 

library(glmnet)
library(pROC)

#Loading the train data and fitting final models

X_train_knn_20_80_20  <- readRDS("X_train_80_20_knn_imputed_top20.rds")
X_train_mean_20_80_20 <- readRDS("X_train_80_20_mean_imputed_top20.rds")
y_train_80_20          <- readRDS("y_train_80_20.rds")
y_train_80_20_factor   <- factor(y_train_80_20, levels = c(0, 1))

set.seed(42)

fold_assigned <- sample(rep(1:10, length.out = length(y_train_80_20_factor)))

final_knn_20_80_20 <- cv.glmnet(x = X_train_knn_20_80_20, y = y_train_80_20_factor, alpha = 0.20, foldid = fold_assigned, family = "binomial", type.measure = "deviance")

final_mean_20_80_20 <- cv.glmnet(x = X_train_mean_20_80_20, y = y_train_80_20_factor, alpha = 0.80, foldid = fold_assigned, family = "binomial", type.measure = "deviance")

saveRDS(final_knn_20_80_20, "final_eln_model_80_20_knn_imputed_top20_alpha020_FINAL.rds")
saveRDS(final_mean_20_80_20, "final_eln_model_80_20_mean_imputed_top20_alpha080_FINAL.rds")


#Loading the test data and subsetting the 17 CpGs

#loaded the train/test split object containing outcome lables for 14 test samples

split_80_20 <- readRDS("sample_split_80_20_min12.rds")

y_test_80_20 <- split_80_20$test_treatment

names(y_test_80_20) <- split_80_20$test_ids

#Loaded the imputed test beta matrices

beta_test_knn  <- readRDS("beta_matrix_test_80_20_top20_min2L_knn_imputed.rds")
beta_test_mean <- readRDS("beta_matrix_test_80_20_top20_min2L_mean_imputed.rds")


#Loaded the top 20 pool from the train data 

top20_dmps_80_20 <- readRDS("myDiff_train_80_20_min12_top20_df.rds")
top20_ids_80_20 <- paste(top20_dmps_80_20$chr, top20_dmps_80_20$start, top20_dmps_80_20$end, sep = ".")


#Seeing which top 20 CpGs have coverage in test set (17 because 3 are not there)

overlap_ids_80_20 <- intersect(top20_ids_80_20, rownames(beta_test_knn))


#Subsetting test beta matrices to the same 17 CpGs from train and transposing the X 

X_test_knn_17  <- t(beta_test_knn[overlap_ids_80_20, , drop = FALSE])
X_test_mean_17 <- t(beta_test_mean[overlap_ids_80_20, , drop = FALSE])

#matching the outcome (y) to samples in test feature matrix 

y_test_80_20_aligned <- y_test_80_20[rownames(X_test_knn_17)]

# Adding the missing CpGs with training set mean

#exttracting list of 20 CpG names the model is expecting wihtouth the intercept

model_cpgs_knn  <- rownames(coef(final_knn_20_80_20,  s = "lambda.min"))[-1]
model_cpgs_mean <- rownames(coef(final_mean_20_80_20, s = "lambda.min"))[-1]

#Show which 20 CpGs are missing from test data

missing_knn  <- setdiff(model_cpgs_knn,  colnames(X_test_knn_17))
missing_mean <- setdiff(model_cpgs_mean, colnames(X_test_mean_17))

#Addinf one column per missing CpG and filling it with that CpGs mean value from train data

X_test_knn_full <- X_test_knn_17
for (cpg in missing_knn) {
  X_test_knn_full <- cbind(X_test_knn_full, mean(X_train_knn_20_80_20[, cpg]))
  colnames(X_test_knn_full)[ncol(X_test_knn_full)] <- cpg
}

#Then reordering the columns to match model order

X_test_knn_full <- X_test_knn_full[, model_cpgs_knn]

#Repeating for mean imputation

X_test_mean_full <- X_test_mean_17
for (cpg in missing_mean) {
  X_test_mean_full <- cbind(X_test_mean_full, mean(X_train_mean_20_80_20[, cpg]))
  colnames(X_test_mean_full)[ncol(X_test_mean_full)] <- cpg
}


X_test_mean_full <- X_test_mean_full[, model_cpgs_mean]

#Predicitng the test split data on both imputations

pred_prob_knn_20  <- predict(final_knn_20_80_20,  newx = X_test_knn_full,  s = "lambda.min", type = "response")

pred_prob_mean_20 <- predict(final_mean_20_80_20, newx = X_test_mean_full, s = "lambda.min", type = "response")

#Calculating the AUC for both imputations

roc_knn_20  <- roc(response = y_test_80_20_aligned, predictor = as.numeric(pred_prob_knn_20))
roc_mean_20 <- roc(response = y_test_80_20_aligned, predictor = as.numeric(pred_prob_mean_20))

auc(roc_knn_20) #AUC for KNN
auc(roc_mean_20) #AUC for mean 

#Claculating cv deviance using glmnet's assess.glmnet() function - is built specifially for test data only when calculating the deviance

dev_knn  <- assess.glmnet(final_knn_20_80_20,  newx = X_test_knn_full,  newy = y_test_80_20_aligned, family = "binomial", s = "lambda.min")$deviance

dev_mean <- assess.glmnet(final_mean_20_80_20, newx = X_test_mean_full, newy = y_test_80_20_aligned, family = "binomial", s = "lambda.min")$deviance

dev_knn #deviance for KNN
dev_mean #deviance for mean

saveRDS(list(pred_prob = pred_prob_knn_20, roc = roc_knn_20, auc = auc(roc_knn_20), deviance = dev_knn, y_test = y_test_80_20_aligned), "test_results_80_20_knn_top20_final.rds")

saveRDS(list(pred_prob = pred_prob_mean_20, roc = roc_mean_20, auc = auc(roc_mean_20), deviance = dev_mean, y_test = y_test_80_20_aligned), "test_results_80_20_mean_top20_final.rds")



