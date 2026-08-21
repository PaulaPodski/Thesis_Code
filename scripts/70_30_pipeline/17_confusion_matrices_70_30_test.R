#!/usr/bin/env Rscript

#Confusion Matrices (Top 20, cutoff = 0.5) on test 70/30 split
#All statistics for confusion matrices (i.e. sensitivity, specificity, etc) used caret package and plots were done via ggplot2

library(caret)
library(ggplot2)


#Loading the save rds objects from the test predictions

result_knn  <- readRDS("test_results_70_30_knn_top20_final.rds")
result_mean <- readRDS("test_results_70_30_mean_top20_final.rds")

cutoff <- 0.5


#Building the confusion matrices with caret where positive class is 1 for dead group

y_pred_class_knn  <- factor(ifelse(result_knn$pred_prob  >= cutoff, 1, 0), levels = c(0, 1))
y_pred_class_mean <- factor(ifelse(result_mean$pred_prob >= cutoff, 1, 0), levels = c(0, 1))
y_actual <- factor(result_knn$y_test, levels = c(0, 1))

cm_knn  <- confusionMatrix(y_pred_class_knn,  y_actual, positive = "1")
cm_mean <- confusionMatrix(y_pred_class_mean, y_actual, positive = "1")

#to view statistics:

print(cm_knn)
print(cm_mean)

saveRDS(cm_knn, "confusion_matrix_70_30_knn_top20.rds")
saveRDS(cm_mean, "confusion_matrix_70_30_mean_top20.rds")

#Plotting the KNN confusion matrix using source template

cm_knn_df <- as.data.frame(as.table(cm_knn$table))
colnames(cm_knn_df) <- c("Prediction", "Reference", "Freq")
cm_knn_df$Label <- ifelse(cm_knn_df$Reference == 0 & cm_knn_df$Prediction == 0, "True Negative", ifelse(cm_knn_df$Reference == 0 & cm_knn_df$Prediction == 1, "False Positive", ifelse(cm_knn_df$Reference == 1 & cm_knn_df$Prediction == 0, "False Negative", "True Positive")))

ggplot(data = cm_knn_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = paste0(Freq, "\n", Label)), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  theme_minimal() +
  labs(title = "70/30 Test Split, KNN Imputation, cutoff = 0.5", x = "Reference", y = "Prediction")

ggsave("confusion_matrix_70_30_knn_cutoff05.png", width = 5, height = 5, dpi = 300, bg = "white")

#Plotting mean confusion matrix using source template

cm_mean_df <- as.data.frame(as.table(cm_mean$table))
colnames(cm_mean_df) <- c("Prediction", "Reference", "Freq")
cm_mean_df$Label <- ifelse(cm_mean_df$Reference == 0 & cm_mean_df$Prediction == 0, "True Negative", ifelse(cm_mean_df$Reference == 0 & cm_mean_df$Prediction == 1, "False Positive", ifelse(cm_mean_df$Reference == 1 & cm_mean_df$Prediction == 0, "False Negative", "True Positive")))

ggplot(data = cm_mean_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = paste0(Freq, "\n", Label)), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  theme_minimal() +
  labs(title = "70/30 Test Split, Mean Imputation, cutoff = 0.5", x = "Reference", y = "Prediction")

ggsave("confusion_matrix_70_30_mean_cutoff05.png", width = 5, height = 5, dpi = 300, bg = "white")
