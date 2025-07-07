# ------------------------------------------------------
# DDS 8521 High-Dimensional Data Analysis
# RNA-Seq (HiSeq) PANCAN Dataset – Clean Workflow
# ------------------------------------------------------

# 1. Set working directory
setwd("E:/DataSets/gene+expression+cancer+rna+seq/TCGA-PANCAN-HiSeq-801x20531/TCGA-PANCAN-HiSeq-801x20531/")

# 2. Install & load packages
install.packages(c("moments", "caret"), dependencies = TRUE)
library(moments)
library(caret)
library(ggplot2)

# 3. Load your data
data <- read.csv("data.csv", header = TRUE, row.names = 1)
labels <- read.csv("labels.csv", header = TRUE, row.names = 1)

# 4. Merge labels into the data frame
data$TumorType <- labels[match(rownames(data), rownames(labels)), 1]

# 5. Prepare numeric matrix for PCA
data_numeric <- data[, -ncol(data)]               # drop TumorType
# remove columns with zero variance
data_numeric <- data_numeric[, apply(data_numeric, 2, var) != 0]

# 6. Scale and run PCA
data_scaled <- scale(data_numeric)
pca <- prcomp(data_scaled)
pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     TumorType = data$TumorType)

# 7. Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = TumorType)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of RNA-Seq Gene Expression", color = "Tumor Type")

# 8. Select a small set of genes for sampling analysis
genes_to_analyze <- colnames(data)[1:5]

# 9. Create three samples of size ~200
set.seed(1)
sample_srs   <- data[sample(nrow(data), 200), ]

set.seed(1)
idx_strat    <- createDataPartition(data$TumorType, p = 0.25, list = FALSE)
sample_strat <- data[idx_strat, ]

step         <- floor(nrow(data) / 200)
idx_sys      <- seq(1, nrow(data), by = step)
sample_sys   <- data[idx_sys, ]

# 10. Function to compute descriptive statistics
compute_stats <- function(df, vars) {
  data.frame(
    Variable = vars,
    Mean     = sapply(df[, vars], mean),
    Variance = sapply(df[, vars], var),
    Skewness = sapply(df[, vars], skewness),
    Kurtosis = sapply(df[, vars], kurtosis)
  )
}

# 11. Compute & combine stats
stats_full   <- compute_stats(data, genes_to_analyze)
stats_srs    <- compute_stats(sample_srs, genes_to_analyze)
stats_strat  <- compute_stats(sample_strat, genes_to_analyze)
stats_sys    <- compute_stats(sample_sys, genes_to_analyze)

stats_full$Sample  <- "Full"
stats_srs$Sample   <- "Simple Random"
stats_strat$Sample <- "Stratified"
stats_sys$Sample   <- "Systematic"

all_stats <- rbind(stats_full, stats_srs, stats_strat, stats_sys)
print(all_stats)

# 12. Histograms & QQ-plots for those genes
par(mfrow = c(5, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

for (g in genes_to_analyze) {
  hist(data[[g]], main = paste("Histogram –", g), xlab = g, breaks = 30, col = "gray")
  qqnorm(data[[g]], main = paste("QQ-Plot –", g))
  qqline(data[[g]], col = "red")
}
