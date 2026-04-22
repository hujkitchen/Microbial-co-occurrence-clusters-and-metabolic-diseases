# ==================== 可视化每个Cluster的主要菌种（PDF格式）====================
cat("\n=== 可视化每个Cluster的主要菌种（PDF格式）===\n")

# 确保目录存在
if(!dir.exists("cluster/cluster_species")) {
  dir.create("cluster/cluster_species")
}
if(!dir.exists("cluster/cluster_species/figures")) {
  dir.create("cluster/cluster_species/figures")
}

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# ==================== 图1：每个Cluster的Top 10菌种条形图（PDF）====================
cat("生成每个Cluster的Top10菌种条形图（PDF）...\n")

for(i in 1:10) {
  # 获取该cluster的数据
  cluster_data <- cluster_species_list[[i]]
  
  # 按丰度排序，取前10
  top_species <- head(cluster_data[order(-cluster_data$Mean_Abundance), ], 10)
  
  # 截断过长的菌种名称
  top_species$Short_Name <- substr(top_species$Species, 1, 30)
  top_species$Short_Name <- factor(top_species$Short_Name, 
                                    levels = rev(top_species$Short_Name))
  
  # 创建条形图
  p <- ggplot(top_species, aes(x = Short_Name, y = Mean_Abundance, 
                                fill = Prevalence)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "lightblue", high = "darkblue", 
                        name = "Prevalence\n(% samples)") +
    labs(title = paste("Cluster", i, "- Top 10 Species"),
         subtitle = paste("Total species:", nrow(cluster_data)),
         x = "", y = "Mean Abundance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.y = element_text(size = 10))
  
  # 保存为PDF
  ggsave(paste0("cluster/cluster_species/figures/cluster", i, "_top10.pdf"),
         p, width = 10, height = 6)
}

cat("Top10条形图PDF已保存\n")

# ==================== 图2：所有Cluster的Top菌种热图（PDF）====================
cat("生成Top菌种热图（PDF）...\n")

# 提取每个cluster的Top5菌种
top5_all <- c()
for(i in 1:10) {
  top5 <- head(cluster_species_list[[i]]$Species[order(-cluster_species_list[[i]]$Mean_Abundance)], 5)
  top5_all <- c(top5_all, top5)
}

# 获取这些菌种的丰度数据
counts <- assay(tse_2w, "counts")
heatmap_data <- counts[top5_all, ]

# 计算每个菌种所属的cluster
species_cluster <- sapply(top5_all, function(sp) {
  for(i in 1:10) {
    if(sp %in% cluster_species_list[[i]]$Species) return(i)
  }
  return(NA)
})

# 按cluster和丰度排序
order_df <- data.frame(Species = top5_all, Cluster = species_cluster, 
                       MeanAb = rowMeans(heatmap_data))
order_df <- order_df[order(order_df$Cluster, -order_df$MeanAb), ]
heatmap_data <- heatmap_data[order_df$Species, ]
species_cluster <- species_cluster[order_df$Species]

# 采样样本（避免图形过大）
set.seed(123)
n_sample <- min(100, ncol(heatmap_data))
sample_idx <- sample(1:ncol(heatmap_data), n_sample)

# 创建热图数据
heatmap_long <- melt(log1p(heatmap_data[, sample_idx]))
colnames(heatmap_long) <- c("Species", "Sample", "LogAbundance")

# 添加cluster信息
heatmap_long$Cluster <- factor(species_cluster[heatmap_long$Species])
heatmap_long$Species <- factor(heatmap_long$Species, levels = rev(order_df$Species))

# 绘制热图
p <- ggplot(heatmap_long, aes(x = Sample, y = Species, fill = LogAbundance)) +
  geom_tile() +
  facet_grid(Cluster ~ ., scales = "free", space = "free") +
  scale_fill_gradient(low = "navy", high = "yellow", name = "Log(Count+1)") +
  labs(title = "Top 5 Species from Each Cluster",
       x = "Samples (subset)", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0, size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14),
        panel.spacing = unit(0.5, "lines"))

ggsave("cluster/cluster_species/figures/top5_species_heatmap.pdf",
       p, width = 12, height = 10)

cat("Top菌种热图PDF已保存\n")

# ==================== 图3：气泡图（PDF）====================
cat("生成气泡图（PDF）...\n")

# 合并所有cluster的数据
all_species_df <- do.call(rbind, cluster_species_list)

# 只显示每个cluster的Top15菌种
top15_all <- c()
for(i in 1:10) {
  top15 <- head(cluster_species_list[[i]]$Species[order(-cluster_species_list[[i]]$Mean_Abundance)], 15)
  top15_all <- c(top15_all, top15)
}
plot_df <- all_species_df[all_species_df$Species %in% top15_all, ]

p <- ggplot(plot_df, aes(x = Mean_Abundance, y = Prevalence, 
                          size = Mean_Abundance, color = factor(Cluster))) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = Species), size = 2.5, vjust = -0.5, hjust = 0.5, check_overlap = TRUE) +
  scale_size_continuous(range = c(2, 10)) +
  scale_color_manual(values = rainbow(10), name = "Cluster") +
  labs(title = "Species Abundance vs Prevalence by Cluster",
       x = "Mean Abundance", y = "Prevalence (% samples)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

ggsave("cluster/cluster_species/figures/species_bubble_plot.pdf",
       p, width = 14, height = 10)

cat("气泡图PDF已保存\n")

# ==================== 图4：每个Cluster的物种丰度分布箱线图（PDF）====================
cat("生成丰度分布箱线图（PDF）...\n")

all_species_df$LogAbundance <- log10(all_species_df$Mean_Abundance + 1)

p <- ggplot(all_species_df, aes(x = factor(Cluster), y = LogAbundance, fill = factor(Cluster))) +
  geom_boxplot() +
  scale_fill_manual(values = rainbow(10)) +
  labs(title = "Species Abundance Distribution by Cluster",
       x = "Cluster", y = "Log10(Mean Abundance + 1)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ggsave("cluster/cluster_species/figures/abundance_boxplot.pdf",
       p, width = 10, height = 6)

cat("丰度箱线图PDF已保存\n")

# ==================== 图5：各Cluster物种数量饼图（PDF）====================
cat("生成物种数量饼图（PDF）...\n")

cluster_sizes <- sapply(cluster_species_list, nrow)
cluster_size_df <- data.frame(
  Cluster = factor(1:10),
  Count = cluster_sizes,
  Percentage = cluster_sizes / sum(cluster_sizes) * 100
)

p <- ggplot(cluster_size_df, aes(x = "", y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%\n(n=", Count, ")")), 
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = rainbow(10)) +
  labs(title = "Species Distribution Across Clusters") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("cluster/cluster_species/figures/cluster_pie_chart.pdf",
       p, width = 8, height = 8)

cat("饼图PDF已保存\n")

# ==================== 图6：水平堆叠条形图（PDF）====================
cat("\n=== 生成水平堆叠条形图 ===\n")

# 确保目录存在
if(!dir.exists("cluster/cluster_species/figures")) {
  dir.create("cluster/cluster_species/figures", recursive = TRUE)
}

# 获取计数数据
counts <- assay(tse_2w, "counts")

# 为每个cluster选择代表性菌种（Top3）
representative_species <- c()
for(i in 1:10) {
  top3 <- head(cluster_species_list[[i]]$Species[order(-cluster_species_list[[i]]$Mean_Abundance)], 3)
  representative_species <- c(representative_species, top3)
}

cat("代表性菌种数量:", length(representative_species), "\n")

# 计算每个菌种所属的cluster
species_to_cluster <- c()
for(sp in representative_species) {
  for(i in 1:10) {
    if(sp %in% cluster_species_list[[i]]$Species) {
      species_to_cluster <- c(species_to_cluster, i)
      break
    }
  }
}

# ==================== 方法1：使用相对丰度（更稳定）====================

# 提取这些菌种的丰度数据
rep_counts <- counts[representative_species, ]

# 转换为相对丰度（每个样本内标准化）
rep_rel <- sweep(rep_counts, 2, colSums(rep_counts), "/")
rep_rel[is.na(rep_rel)] <- 0
rep_rel[is.infinite(rep_rel)] <- 0

# 按疾病状态分组
group_info <- colData(tse_2w)$condition_simplified
group_simple <- ifelse(group_info == "Healthy", "Healthy", "Disease")

# 计算每组平均
group_means <- matrix(0, nrow = length(representative_species), ncol = 2)
rownames(group_means) <- representative_species
colnames(group_means) <- c("Healthy", "Disease")

for(i in 1:length(representative_species)) {
  group_means[i, "Healthy"] <- mean(rep_rel[i, group_simple == "Healthy"], na.rm = TRUE)
  group_means[i, "Disease"] <- mean(rep_rel[i, group_simple == "Disease"], na.rm = TRUE)
}

# 转换为数据框
group_means_df <- data.frame(
  Species = rep(representative_species, 2),
  Group = rep(c("Healthy", "Disease"), each = length(representative_species)),
  Abundance = c(group_means[, "Healthy"], group_means[, "Disease"]),
  Cluster = rep(species_to_cluster, 2)
)

# 移除缺失值
group_means_df <- group_means_df[!is.na(group_means_df$Abundance) & 
                                  !is.infinite(group_means_df$Abundance), ]

# 按Cluster排序
group_means_df$Species <- factor(group_means_df$Species, 
                                  levels = representative_species)
group_means_df$Cluster <- factor(group_means_df$Cluster)

# 绘制水平堆叠条形图
p <- ggplot(group_means_df, aes(x = Species, y = Abundance, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Cluster ~ ., scales = "free", space = "free") +
  coord_flip() +
  scale_fill_manual(values = c("Healthy" = "#2E86AB", "Disease" = "#D64045")) +
  labs(title = "Representative Species Abundance: Healthy vs Disease",
       x = "", y = "Relative Abundance") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        strip.text.y = element_text(angle = 0, size = 10, face = "bold"),
        axis.text.y = element_text(size = 8))

# 保存
ggsave("cluster/cluster_species/figures/representative_species_comparison.pdf",
       p, width = 12, height = 14)

cat("水平堆叠条形图已保存\n")

# ==================== 方法2：点图（替代方案，更稳定）====================
cat("生成替代点图...\n")

# 计算效应量
effect_size <- data.frame(
  Species = representative_species,
  Cluster = species_to_cluster,
  Healthy_Mean = group_means[, "Healthy"],
  Disease_Mean = group_means[, "Disease"],
  Log2FC = log2((group_means[, "Disease"] + 0.0001) / (group_means[, "Healthy"] + 0.0001))
)

# 移除异常值
effect_size <- effect_size[is.finite(effect_size$Log2FC), ]

# 按Cluster分组
effect_size$Cluster <- factor(effect_size$Cluster)

# 绘制点图
p2 <- ggplot(effect_size, aes(x = Log2FC, y = reorder(Species, Log2FC), 
                               color = Cluster)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = rainbow(10)) +
  labs(title = "Effect Size of Representative Species (Disease vs Healthy)",
       x = "Log2 Fold Change", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.y = element_text(size = 8))

ggsave("cluster/cluster_species/figures/representative_species_effect.pdf",
       p2, width = 10, height = 12)

cat("点图已保存\n")

# ==================== 方法3：热图（每个Cluster的Top5菌种比较）====================
cat("生成Cluster Top菌种热图...\n")

# 为每个cluster选择Top5菌种
top5_per_cluster <- c()
cluster_labels <- c()
for(i in 1:10) {
  top5 <- head(cluster_species_list[[i]]$Species[order(-cluster_species_list[[i]]$Mean_Abundance)], 5)
  top5_per_cluster <- c(top5_per_cluster, top5)
  cluster_labels <- c(cluster_labels, rep(i, length(top5)))
}

# 提取丰度数据
top5_counts <- counts[top5_per_cluster, ]

# 计算健康组和疾病组的平均丰度
group_means_top5 <- matrix(0, nrow = length(top5_per_cluster), ncol = 2)
rownames(group_means_top5) <- top5_per_cluster
colnames(group_means_top5) <- c("Healthy", "Disease")

for(i in 1:length(top5_per_cluster)) {
  group_means_top5[i, "Healthy"] <- mean(top5_counts[i, group_simple == "Healthy"])
  group_means_top5[i, "Disease"] <- mean(top5_counts[i, group_simple == "Disease"])
}

# 计算Log2FC
log2fc_top5 <- log2((group_means_top5[, "Disease"] + 0.0001) / 
                     (group_means_top5[, "Healthy"] + 0.0001))

# 创建热图数据
heatmap_df <- data.frame(
  Species = top5_per_cluster,
  Cluster = factor(cluster_labels),
  Log2FC = log2fc_top5
)

# 按Cluster和Log2FC排序
heatmap_df$Species <- factor(heatmap_df$Species, 
                              levels = rev(heatmap_df$Species[order(heatmap_df$Cluster, heatmap_df$Log2FC)]))

# 绘制热图
p3 <- ggplot(heatmap_df, aes(x = Cluster, y = Species, fill = Log2FC)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Log2FC)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Log2 FC\n(Disease/Healthy)") +
  labs(title = "Top 5 Species from Each Cluster: Disease vs Healthy",
       x = "Cluster", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8))

ggsave("cluster/cluster_species/figures/top5_cluster_heatmap.pdf",
       p3, width = 10, height = 12)

cat("Cluster Top菌种热图已保存\n")

# ==================== 方法4：简单的条形图（每个Cluster的Top3）====================
cat("生成各Cluster Top3菌种条形图...\n")

for(i in 1:10) {
  # 获取该cluster的Top3菌种
  top3 <- head(cluster_species_list[[i]]$Species[order(-cluster_species_list[[i]]$Mean_Abundance)], 3)
  
  # 提取丰度数据
  top3_counts <- counts[top3, ]
  
  # 计算健康组和疾病组的平均丰度
  top3_means <- matrix(0, nrow = length(top3), ncol = 2)
  rownames(top3_means) <- top3
  colnames(top3_means) <- c("Healthy", "Disease")
  
  for(j in 1:length(top3)) {
    top3_means[j, "Healthy"] <- mean(top3_counts[j, group_simple == "Healthy"])
    top3_means[j, "Disease"] <- mean(top3_counts[j, group_simple == "Disease"])
  }
  
  # 转换为数据框
  top3_df <- data.frame(
    Species = rep(top3, 2),
    Group = rep(c("Healthy", "Disease"), each = length(top3)),
    Abundance = c(top3_means[, "Healthy"], top3_means[, "Disease"])
  )
  
  # 绘制条形图
  p4 <- ggplot(top3_df, aes(x = Species, y = Abundance, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("Healthy" = "#2E86AB", "Disease" = "#D64045")) +
    labs(title = paste("Cluster", i, "- Top 3 Species"),
         x = "", y = "Mean Abundance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0("cluster/cluster_species/figures/cluster", i, "_top3_comparison.pdf"),
         p4, width = 8, height = 4)
}

cat("各Cluster Top3菌种条形图已保存\n")

# ==================== 完成 ====================
cat("\n\n========================================\n")
cat("所有图形已保存到: cluster/cluster_species/figures/\n")
cat("生成的文件:\n")
cat("  - representative_species_comparison.pdf (水平堆叠条形图)\n")
cat("  - representative_species_effect.pdf (效应量点图)\n")
cat("  - top5_cluster_heatmap.pdf (Top5菌种热图)\n")
cat("  - cluster*_top3_comparison.pdf (各Cluster Top3菌种比较)\n")
cat("========================================\n")

# ==================== 图7：所有Cluster的Top菌种汇总图（PDF）====================
cat("生成Top菌种汇总图（PDF）...\n")

# 创建2x5的网格布局
plot_list <- list()

for(i in 1:10) {
  cluster_data <- cluster_species_list[[i]]
  top_species <- head(cluster_data[order(-cluster_data$Mean_Abundance), ], 8)
  top_species$Short_Name <- substr(top_species$Species, 1, 25)
  top_species$Short_Name <- factor(top_species$Short_Name, 
                                    levels = rev(top_species$Short_Name))
  
  p <- ggplot(top_species, aes(x = Short_Name, y = Mean_Abundance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("Cluster", i, "(", nrow(cluster_data), "species)"),
         x = "", y = "Mean Abundance") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 7))
  
  plot_list[[i]] <- p
}

# 组合图形
pdf("cluster/cluster_species/figures/all_clusters_top8.pdf", width = 15, height = 12)
grid.arrange(grobs = plot_list, ncol = 2, nrow = 5,
             top = "Top 8 Species in Each Cluster")
dev.off()

cat("汇总图PDF已保存\n")

# ==================== 图8：网络图（PDF）====================
cat("生成菌种共现网络图（PDF）...\n")

# 计算Top30菌种的相关性
all_species <- all_species_df$Species
mean_abundance_all <- all_species_df$Mean_Abundance
top30_species <- head(all_species[order(mean_abundance_all, decreasing = TRUE)], 30)

# 计算相关性
cor_matrix <- cor(t(counts[top30_species, ]), method = "spearman")

# 创建网络
library(igraph)
threshold <- 0.5
adj_matrix <- abs(cor_matrix) > threshold & cor_matrix > 0
diag(adj_matrix) <- FALSE

if(sum(adj_matrix) > 0) {
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # 添加节点属性
  V(g)$cluster <- sapply(V(g)$name, function(sp) {
    for(i in 1:10) {
      if(sp %in% cluster_species_list[[i]]$Species) return(i)
    }
    return(0)
  })
  V(g)$size <- log(rowMeans(counts[top30_species, ]) + 1) * 5
  
  # 保存为PDF
  pdf("cluster/cluster_species/figures/species_network.pdf", width = 12, height = 10)
  
  plot(g, 
       vertex.color = rainbow(10)[V(g)$cluster],
       vertex.size = V(g)$size,
       vertex.label = V(g)$name,
       vertex.label.cex = 0.6,
       vertex.label.dist = 1,
       main = "Species Co-occurrence Network (Top 30 species, r > 0.5)")
  
  legend("topright", legend = paste("Cluster", 1:10), 
         col = rainbow(10), pch = 19, cex = 0.8)
  
  dev.off()
  cat("网络图PDF已保存\n")
}

# ==================== 打印每个Cluster的菌种列表 ====================
cat("\n\n========================================\n")
cat("每个Cluster的菌种列表\n")
cat("========================================\n")

for(i in 1:10) {
  cat(sprintf("\n【Cluster %d】- 共 %d 个菌种\n", i, nrow(cluster_species_list[[i]])))
  cat("主要菌种（按丰度排序）:\n")
  
  top10 <- head(cluster_species_list[[i]][order(-cluster_species_list[[i]]$Mean_Abundance), ], 10)
  for(j in 1:nrow(top10)) {
    cat(sprintf("  %d. %s (丰度: %.4f, 出现率: %.1f%%)\n", 
                j, top10$Species[j], top10$Mean_Abundance[j], top10$Prevalence[j] * 100))
  }
  if(nrow(cluster_species_list[[i]]) > 10) {
    cat(sprintf("  ... 还有 %d 个菌种\n", nrow(cluster_species_list[[i]]) - 10))
  }
}

# ==================== 完成 ====================
cat("\n\n========================================\n")
cat("所有PDF文件已保存到: cluster/cluster_species/figures/\n")
cat("生成的文件列表:\n")
cat("  - cluster1_top10.pdf ... cluster10_top10.pdf (各Cluster Top10菌种)\n")
cat("  - top5_species_heatmap.pdf (Top5菌种热图)\n")
cat("  - species_bubble_plot.pdf (气泡图)\n")
cat("  - abundance_boxplot.pdf (丰度箱线图)\n")
cat("  - cluster_pie_chart.pdf (物种分布饼图)\n")
cat("  - representative_species_comparison.pdf (健康vs疾病比较)\n")
cat("  - all_clusters_top8.pdf (所有Cluster汇总图)\n")
cat("  - species_network.pdf (菌种共现网络图)\n")
cat("========================================\n")