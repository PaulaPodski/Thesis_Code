#!/usr/bin/env Rscript

#SHAP and Feature Importance Plots for the test samples from the 80/20 split


library(glmnet)
library(iml)
library(ggplot2)


#Load the model the train data and test outputs

final_knn_20_80_20    <- readRDS("final_eln_model_80_20_knn_imputed_top20_alpha020_FINAL.rds")
X_train_knn_20_80_20  <- readRDS("X_train_80_20_knn_imputed_top20.rds")
result_knn             <- readRDS("test_results_80_20_knn_top20_final.rds")

#Then the predict function is needed to be built and with the missing 3 CpGs (train means)

model_cpgs_knn  <- rownames(coef(final_knn_20_80_20, s = "lambda.min"))[-1]
train_means_knn <- colMeans(X_train_knn_20_80_20)

predict_fn_knn <- function(object, newdata) {
  newdata <- as.data.frame(newdata)
  full_input <- matrix(rep(train_means_knn, nrow(newdata)), nrow = nrow(newdata), byrow = TRUE)
  colnames(full_input) <- model_cpgs_knn
  full_input[, colnames(newdata)] <- as.matrix(newdata)
  predict(object, newx = full_input, s = "lambda.min", type = "response")[, 1]
}

#Loading the test data with subsetting the 17 needed CpGs

split_80_20 <- readRDS("sample_split_80_20_min12.rds")
beta_test_knn <- readRDS("beta_matrix_test_80_20_top20_min2L_knn_imputed.rds")
top20_dmps_80_20 <- readRDS("myDiff_train_80_20_min12_top20_df.rds")
top20_ids_80_20 <- paste(top20_dmps_80_20$chr, top20_dmps_80_20$start, top20_dmps_80_20$end, sep = ".")
overlap_ids_80_20 <- intersect(top20_ids_80_20, rownames(beta_test_knn))

X_test_knn_17 <- t(beta_test_knn[overlap_ids_80_20, , drop = FALSE])


#Required to build the iml Predictor Object for SHAP calculations

predictor_knn <- Predictor$new(final_knn_20_80_20, data = as.data.frame(X_test_knn_17), y = result_knn$y_test, predict.fun = predict_fn_knn)


#Then we create the SHAP for all 14 patient samples (14 plots per individual sample)

for (i in 1:nrow(X_test_knn_17)) {

  sample_id <- rownames(X_test_knn_17)[i]
  x_interest_i <- as.data.frame(X_test_knn_17[i, , drop = FALSE])

  shapley_i <- Shapley$new(predictor_knn, x.interest = x_interest_i)

  p_i <- shapley_i$plot()
  ggsave(paste0("shap_plot_80_20_knn_", sample_id, ".png"), plot = p_i,
         width = 7, height = 5, dpi = 300, bg = "white")

}

#Lastly, feature importance plot for the test split

imp_80_20_knn <- FeatureImp$new(predictor_knn, loss = "mae")
p_imp_80_20_knn <- plot(imp_80_20_knn)

ggsave("feature_importance_80_20_knn_FeatureImp.png", plot = p_imp_80_20_knn,
       width = 7, height = 5, dpi = 300, bg = "white")

#Combined SHAP summary plot across all test patients

shap_results_all_80_20_knn <- readRDS("shap_results_all_patients_80_20_knn_iml.rds")

p_summary_80_20_knn <- ggplot(shap_results_all_80_20_knn, aes(x = feature, y = phi)) +
  geom_jitter(width = 0.15, alpha = 0.6, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  labs(title = "SHAP Values Across All Test Patients (80/20 KNN)",
       x = NULL, y = "Phi (SHAP contribution)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("shap_summary_all_patients_80_20_knn.png", plot = p_summary_80_20_knn, width = 7, height = 6, dpi = 300, bg = "white")
