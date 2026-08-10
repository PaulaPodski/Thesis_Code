#!/usr/bin/env Rscript

#70/30 test split individual sample (mean) SHAP plots using source template

library(glmnet)
library(iml)
library(ggplot2)


#Loading the fitted Mean model and the training data it was built on and testing the predictions and labels already computed from previous steps

final_mean_20_70_30 <- readRDS("final_eln_model_70_30_mean_imputed_top20_alpha020_FINAL.rds")
X_train_mean_20_70_30 <- readRDS("X_train_70_30_mean_imputed_top20.rds")
result_mean <- readRDS("test_results_70_30_mean_top20_final.rds")

#the model would expect 20 CpG features but 17 are in the test set so we need to load the full 20, extract mean values per each column and fill in the 3 missing CpG columns

model_cpgs_mean <- rownames(coef(final_mean_20_70_30, s = "lambda.min"))[-1]
train_means_mean <- colMeans(X_train_mean_20_70_30)

#Now creating the adaptor function that iml package requires where it fills the missing CpGs with train mean values before calling predict

predict_fn_mean <- function(object, newdata) {
  newdata <- as.data.frame(newdata)
  full_input <- matrix(rep(train_means_mean, nrow(newdata)), nrow = nrow(newdata), byrow = TRUE)
  colnames(full_input) <- model_cpgs_mean
  full_input[, colnames(newdata)] <- as.matrix(newdata)
  predict(object, newx = full_input, s = "lambda.min", type = "response")[, 1]
}

#Loading the test set and subsetting to 17 needed CpGs for the test set

split_70_30 <- readRDS("sample_split_70_30_min12.rds")
beta_test_mean <- readRDS("beta_matrix_test_70_30_top20_min4L_mean_imputed.rds")
top20_dmps_70_30 <- readRDS("myDiff_train_70_30_min12_top20_df.rds")
top20_ids_70_30 <- paste(top20_dmps_70_30$chr, top20_dmps_70_30$start, top20_dmps_70_30$end, sep = ".")
overlap_ids_70_30 <- intersect(top20_ids_70_30, rownames(beta_test_mean))

X_test_mean_17_70_30 <- t(beta_test_mean[overlap_ids_70_30, , drop = FALSE])

#Building the predictor object iml package requires

predictor_mean <- Predictor$new(final_mean_20_70_30, data = as.data.frame(X_test_mean_17_70_30), y = result_mean$y_test, predict.fun = predict_fn_mean)

#Reproducibility setting the seed before SHAP analysis as iml uses internal random sampling

set.seed(42)

#creating a for loop for all 22 patients behaviour in the elastic net model to calculate and plot the SHAP plot for each patient sample

for (i in 1:nrow(X_test_mean_17_70_30)) {
  sample_id <- rownames(X_test_mean_17_70_30)[i]
  x_interest_i <- as.data.frame(X_test_mean_17_70_30[i, , drop = FALSE])
  shapley_i <- Shapley$new(predictor_mean, x.interest = x_interest_i)
  p_i <- shapley_i$plot()
  ggsave(paste0("shap_plot_70_30_mean_", sample_id, ".png"), plot = p_i, width = 7, height = 5, dpi = 300, bg = "white")
}

#second part of SHAP is creating the feature importance i.e feature contribution to the model across all patients

imp_70_30_mean <- FeatureImp$new(predictor_mean, loss = "mae")
p_imp_70_30_mean <- plot(imp_70_30_mean)
ggsave("feature_importance_70_30_mean_FeatureImp.png", plot = p_imp_70_30_mean, width = 7, height = 5, dpi = 300, bg = "white")
