# ==================== 设置环境和创建输出目录 ====================
# 创建输出目录
if(!dir.exists("cluster")) {
  dir.create("cluster")
}
if(!dir.exists("cluster/figures")) {
  dir.create("cluster/figures")
}
if(!dir.exists("cluster/data")) {
  dir.create("cluster/data")
}
if(!dir.exists("cluster/reports")) {
  dir.create("cluster/reports")
}

# 加载必要的包
library(mia)
library(vegan)
library(scater)
library(randomForest)
library(pheatmap)
library(ggplot2)
library(cluster)
library(RColorBrewer)
library(pROC)
library(corrplot)
library(factoextra)
library(ggpubr)
library(NbClust)
library(caret)
library(igraph)
library(dendextend)

cat("=== 开始物种聚类分析（10个聚类）===\n\n")

# ==================== 1. 计算物种相关性 ====================
cat("=== 1. 计算物种相关性 ===\n")

zero_samples <- which(colSums(assay(tse_2w, "counts")) == 0)
if(length(zero_samples) > 0) {
cat("移除以下零计数样本：", names(zero_samples), "\n")
tse_2w <- tse_2w[, -zero_samples]
}

# 提取计数数据
counts <- assay(tse_2w, "counts")
cat("数据维度:", nrow(counts), "个物种,", ncol(counts), "个样本\n")

# 计算物种间的Spearman相关性
cat("计算物种间相关性...\n")
species_cor <- cor(t(counts), method = "spearman")
cat("相关性矩阵计算完成\n")

# 保存相关性矩阵
write.csv(species_cor, "cluster/data/species_correlation_matrix.csv")
cat("相关性矩阵已保存\n\n")

# ==================== 2. 物种聚类 ====================
cat("=== 2. 物种聚类 ===\n")

# 转换为距离矩阵
dist_species <- as.dist(1 - abs(species_cor))
cat("距离矩阵计算完成\n")

# 层次聚类
hc_species <- hclust(dist_species, method = "ward.D2")
cat("层次聚类完成\n")

# 切割成10个聚类
k <- 10
species_clusters <- cutree(hc_species, k = k)
cat("聚类完成，共", k, "个聚类\n")
print(table(species_clusters))

# 保存物种聚类结果
species_cluster_df <- data.frame(
  Species = rownames(counts),
  Cluster = species_clusters,
  Cluster_Name = paste0("Cluster_", species_clusters)
)
write.csv(species_cluster_df, "cluster/data/species_clusters.csv", row.names = FALSE)
cat("物种聚类结果已保存\n\n")

# ==================== 3. 计算样本菌群组成 ====================
cat("=== 3. 计算样本菌群组成 ===\n")

# 计算每个样本的菌群丰度
cluster_abundance <- matrix(0, ncol(tse_2w), k)
colnames(cluster_abundance) <- paste0("Cluster", 1:k)
rownames(cluster_abundance) <- colnames(tse_2w)

for(i in 1:k) {
  cluster_species <- species_cluster_df$Species[species_cluster_df$Cluster == i]
  cluster_abundance[, i] <- colSums(counts[cluster_species, ])
}

# 转换为相对丰度
cluster_abundance_rel <- cluster_abundance / rowSums(cluster_abundance)

# 保存菌群丰度
write.csv(cluster_abundance_rel, "cluster/data/cluster_abundance.csv")
write.csv(data.frame(
  Sample = rownames(cluster_abundance_rel),
  Group = colData(tse_2w)$group,
  cluster_abundance_rel
), "cluster/data/cluster_abundance_with_metadata.csv", row.names = FALSE)

cat("菌群丰度计算完成\n\n")

# ==================== 4. 箱线图（修复版）====================
cat("=== 4. 绘制箱线图 ===\n")

# 方法1：使用ggplot2绘制所有箱线图在一张图上
library(ggplot2)
library(tidyr)

# 准备长格式数据
boxplot_data <- data.frame(
  Sample = rownames(cluster_abundance_rel),
  Group = colData(tse_2w)$group,
  cluster_abundance_rel
)

boxplot_long <- pivot_longer(boxplot_data, 
                              cols = starts_with("Cluster"),
                              names_to = "Cluster", 
                              values_to = "Abundance")

# 绘制所有簇的箱线图
pdf("cluster/figures/boxplot_all_clusters_ggplot.pdf", width = 14, height = 8)
p <- ggplot(boxplot_long, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~ Cluster, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = c("healthy" = "lightblue", "disease" = "lightcoral")) +
  labs(title = "Microbial Cluster Abundance: Healthy vs Disease",
       x = "", y = "Relative Abundance") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p)
dev.off()

cat("箱线图已保存到: cluster/figures/boxplot_all_clusters_ggplot.pdf\n")

# 方法2：单独绘制每个簇的箱线图
pdf("cluster/figures/boxplot_individual_clusters.pdf", width = 12, height = 10)
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
for(i in 1:k) {
  boxplot(cluster_abundance_rel[, i] ~ colData(tse_2w)$group,
          main = paste("Cluster", i),
          ylab = "Relative abundance",
          col = c("lightblue", "lightcoral"),
          names = c("Healthy", "Disease"))
}
dev.off()
cat("单独箱线图已保存到: cluster/figures/boxplot_individual_clusters.pdf\n\n")

# ==================== 5. 堆叠柱状图（修复版）====================
cat("=== 5. 绘制堆叠柱状图 ===\n")

# 计算平均组成
group_means <- aggregate(cluster_abundance_rel, 
                         by = list(Group = colData(tse_2w)$group), 
                         FUN = mean)
rownames(group_means) <- group_means$Group
group_means_matrix <- as.matrix(group_means[, -1])

# 使用ggplot2绘制堆叠柱状图
stacked_data <- data.frame(
  Group = rep(rownames(group_means_matrix), each = ncol(group_means_matrix)),
  Cluster = rep(colnames(group_means_matrix), times = nrow(group_means_matrix)),
  Abundance = as.vector(t(group_means_matrix))
)

pdf("cluster/figures/composition_stacked_bar.pdf", width = 8, height = 6)
p <- ggplot(stacked_data, aes(x = Group, y = Abundance, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = rainbow(k)) +
  labs(title = "Average Microbial Cluster Composition",
       x = "", y = "Relative Abundance") +
  theme_minimal() +
  theme(legend.position = "right")
print(p)
dev.off()

cat("堆叠柱状图已保存到: cluster/figures/composition_stacked_bar.pdf\n\n")

# ==================== 6. 丰度分布图（使用ggplot2）====================
cat("=== 6. 绘制丰度分布图 ===\n")

# 准备数据
abundance_long <- data.frame(
  Abundance = as.vector(cluster_abundance_rel),
  Cluster = rep(paste0("Cluster", 1:k), each = nrow(cluster_abundance_rel)),
  Group = rep(colData(tse_2w)$group, times = k)
)

# 箱线图（ggplot2版本）
pdf("cluster/figures/cluster_abundance_boxplot_ggplot.pdf", width = 12, height = 6)
p <- ggplot(abundance_long, aes(x = Cluster, y = Abundance, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("healthy" = "lightblue", "disease" = "lightcoral")) +
  labs(title = "Microbial Cluster Abundance Distribution",
       x = "Cluster", y = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

# 小提琴图
pdf("cluster/figures/cluster_abundance_violin_ggplot.pdf", width = 12, height = 6)
p <- ggplot(abundance_long, aes(x = Cluster, y = Abundance, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("healthy" = "lightblue", "disease" = "lightcoral")) +
  labs(title = "Microbial Cluster Abundance Distribution (Violin Plot)",
       x = "Cluster", y = "Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

cat("丰度分布图已保存\n\n")

# ==================== 7. 相关性热图（修复版）====================
cat("=== 7. 绘制相关性热图 ===\n")

# 计算聚类间平均相关性
cluster_cor_mean <- matrix(0, k, k)
for(i in 1:k) {
  for(j in 1:k) {
    species_i <- which(species_clusters == i)
    species_j <- which(species_clusters == j)
    cluster_cor_mean[i, j] <- mean(species_cor[species_i, species_j], na.rm = TRUE)
  }
}
rownames(cluster_cor_mean) <- paste("Cluster", 1:k)
colnames(cluster_cor_mean) <- paste("Cluster", 1:k)

# 保存
write.csv(cluster_cor_mean, "cluster/data/cluster_correlation_matrix.csv")

# 使用corrplot绘制热图
pdf("cluster/figures/cluster_correlation_heatmap.pdf", width = 8, height = 7)
corrplot(cluster_cor_mean, method = "color", type = "upper",
         order = "original", addCoef.col = "black",
         tl.col = "black", tl.srt = 45,
         title = "Average Correlation Between Clusters",
         mar = c(0,0,1,0))
dev.off()

# 使用pheatmap绘制（备用）
pdf("cluster/figures/cluster_correlation_pheatmap.pdf", width = 8, height = 7)
pheatmap(cluster_cor_mean,
         main = "Average Correlation Between Clusters",
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("navy", "white", "red"))(50))
dev.off()

cat("相关性热图已保存\n\n")

# ==================== 8. 差异分析和火山图 ====================
cat("=== 8. 差异分析和火山图 ===\n")

# Wilcoxon检验
diff_results <- data.frame(
  Cluster = 1:k,
  Healthy_mean = NA,
  Disease_mean = NA,
  Log2FC = NA,
  P_value = NA,
  FDR = NA
)

for(i in 1:k) {
  healthy_vals <- cluster_abundance_rel[colData(tse_2w)$group == "healthy", i]
  disease_vals <- cluster_abundance_rel[colData(tse_2w)$group == "disease", i]
  
  diff_results$Healthy_mean[i] <- mean(healthy_vals)
  diff_results$Disease_mean[i] <- mean(disease_vals)
  diff_results$Log2FC[i] <- log2((mean(disease_vals) + 0.0001) / (mean(healthy_vals) + 0.0001))
  
  test <- wilcox.test(healthy_vals, disease_vals)
  diff_results$P_value[i] <- test$p.value
}

diff_results$FDR <- p.adjust(diff_results$P_value, method = "fdr")
diff_results$Significant <- diff_results$FDR < 0.05
diff_results$Direction <- ifelse(diff_results$Log2FC > 0, "Up_in_disease", "Down_in_disease")

write.csv(diff_results, "cluster/data/differential_analysis.csv", row.names = FALSE)

# 火山图
pdf("cluster/figures/volcano_plot.pdf", width = 8, height = 6)
p <- ggplot(diff_results, aes(x = Log2FC, y = -log10(P_value), color = Significant)) +
  geom_point(size = 4) +
  geom_text(aes(label = ifelse(Significant, paste("Cluster", Cluster), "")),
            vjust = -0.5, hjust = 0.5, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  labs(title = "Volcano Plot: Disease vs Healthy",
       x = "Log2 Fold Change (Disease/Healthy)",
       y = "-log10(P-value)") +
  theme_minimal()
print(p)
dev.off()

cat("差异分析和火山图已保存\n\n")

# ==================== 9. 物种树状图 ====================
cat("=== 9. 绘制物种树状图 ===\n")

# 基础树状图
pdf("cluster/figures/species_dendrogram_basic.pdf", width = 14, height = 8)
plot(hc_species, labels = FALSE, main = "Species Clustering Dendrogram (10 clusters)",
     xlab = "Species", sub = "")
rect.hclust(hc_species, k = k, border = 2:(k+1))
dev.off()

# 简化版树状图（避免文件过大）
pdf("cluster/figures/species_dendrogram_simple.pdf", width = 10, height = 6)
plot(hc_species, labels = FALSE, main = "Species Clustering Dendrogram",
     xlab = "", sub = "", hang = -1)
dev.off()

cat("物种树状图已保存\n\n")

# ==================== 10. 生成报告 ====================
cat("=== 10. 生成报告 ===\n")

sink("cluster/reports/analysis_report.txt")

cat("========================================\n")
cat("微生物组菌群聚类分析报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 数据概况\n")
cat("   - 总样本数:", ncol(tse_2w), "\n")
cat("   - 健康组:", sum(colData(tse_2w)$group == "healthy"), "\n")
cat("   - 疾病组:", sum(colData(tse_2w)$group == "disease"), "\n")
cat("   - 物种数:", nrow(tse_2w), "\n")
cat("   - 菌群聚类数:", k, "\n\n")

cat("2. 各聚类大小\n")
for(i in 1:k) {
  cat("   - Cluster", i, ":", sum(species_clusters == i), "物种\n")
}
cat("\n")

cat("3. 差异分析结果\n")
sig_count <- sum(diff_results$Significant)
cat("   - 显著差异的菌群 (FDR < 0.05):", sig_count, "/", k, "\n")
if(sig_count > 0) {
  cat("   - 在疾病组升高的菌群:", sum(diff_results$Direction == "Up_in_disease" & diff_results$Significant), "\n")
  cat("   - 在疾病组降低的菌群:", sum(diff_results$Direction == "Down_in_disease" & diff_results$Significant), "\n")
}
cat("\n")

cat("4. 输出文件\n")
cat("   - 数据文件: cluster/data/\n")
cat("   - 图表文件: cluster/figures/\n")
cat("   - 报告文件: cluster/reports/\n")

sink()

cat("报告已保存到: cluster/reports/analysis_report.txt\n\n")

# ==================== 完成 ====================
cat("========================================\n")
cat("分析完成！\n")
cat("所有结果已保存到 'cluster' 目录\n")
cat("========================================\n")

# 列出生成的文件
cat("\n生成的文件:\n")
cat("\n数据文件:\n")
list.files("cluster/data", pattern = "\\.csv$")
cat("\n图表文件:\n")
list.files("cluster/figures", pattern = "\\.pdf$")