#This pipeline describes loading DMA results from train 80/20 split from min.per.group=12L and plotting results via Volcano Plot
#This code was executed interactively on Apocrita using the following command: salloc --job-name=data_analysis --cpus-per-task=4 --mem=64G --time=06:00:00

#This was also ran within an environment built before launching R - please make sure you execute this before any data loading/analysis:

source /gpfs/scratch/hfy327/miniforge3/etc/profile.d/conda.sh
conda activate methylation_env_clean

#load relevant packages:

library(methylKit)
library(ggplot2)

#Load the DMA object
myDiff_train_80_20 <- getData(readRDS("myDiff_train_80_20_min12.rds"))

#Calculating the significance thresholds for DMPs +/- 25% difference and significance

myDiff_train_80_20$sig <- "Not significant"
myDiff_train_80_20$sig[myDiff_train_80_20$qvalue < 0.01 & myDiff_train_80_20$meth.diff > 25]  <- "Hypermethylated"
myDiff_train_80_20$sig[myDiff_train_80_20$qvalue < 0.01 & myDiff_train_80_20$meth.diff < -25] <- "Hypomethylated"

#Verifying and classifying DMPs into their respective groups and showing a summary table for each:

sig_table_80_20 <- as.data.frame(table(myDiff_train_80_20$sig))
colnames(sig_table_80_20) <- c("Methylation Status", "Count")
print(sig_table_80_20)

#Plotting the volcano plot of the DMA:

volcano_plot_80_20 <- ggplot(myDiff_train_80_20, aes(x = meth.diff, y = -log10(qvalue), color = sig)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_colour_manual(values = c(
    "Hypermethylated" = "red",
    "Hypomethylated" = "blue",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-25, 25), linetype = "solid", colour = "black") +
  geom_hline(yintercept = -log10(0.01), linetype = "solid", colour = "black") +
  labs(
    x = "Methylation difference (%)",
    y = "-log10(q-value)",
    title = "80/20 Train DMA",
    colour = "Threshold"
  ) +
  theme_minimal()

ggsave("volcano_plot_80_20_train.png", plot = volcano_plot_80_20, width = 8, height = 6, dpi = 300, bg = "white")
