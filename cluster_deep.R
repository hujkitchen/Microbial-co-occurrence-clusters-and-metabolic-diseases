# ==================== 加载必要的包 ====================
library(ggplot2)
library(pheatmap)
library(vegan)
library(randomForest)
library(pROC)
library(corrplot)
library(igraph)
library(glmnet)
library(caret)

# 创建深入分析目录
if(!dir.exists("cluster/deep_analysis")) {
  dir.create("cluster/deep_analysis")
}
if(!dir.exists("cluster/deep_analysis/figures")) {
  dir.create("cluster/deep_analysis/figures")
}
if(!dir.exists("cluster/deep_analysis/data")) {
  dir.create("cluster/deep_analysis/data")
}

cat("=== 开始深入分析 ===\n\n")

# ==================== 1. 菌群功能预测 ====================
cat("=== 1. 菌群功能预测 ===\n")

# 基于已知文献，为每个聚类分配功能
cluster_functions <- data.frame(
  Cluster = 1:10,
  Function = c(
    "Butyrate producers, anti-inflammatory",
    "Opportunistic pathogens, pro-inflammatory",
    "Mucin degraders, metabolic health",
    "Lactate producers, probiotics",
    "Succinate producers, plant polysaccharide degraders",
    "Secondary bile acid producers",
    "Mucin-degrading, IBD-associated",
    "Unclassified/novel species",
    "Propionate producers",
    "Skin/oral commensals, potential contaminants"
  ),
  Clinical_Association = c(
    "Protective",
    "Harmful",
    "Mixed",
    "Protective",
    "Protective",
    "Mixed",
    "Harmful",
    "Unknown",
    "Protective",
    "Neutral"
  )
)

write.csv(cluster_functions, "cluster/deep_analysis/data/cluster_functions.csv", row.names = FALSE)
print(cluster_functions)

# ==================== 2. 疾病风险评分 ====================
cat("\n=== 2. 疾病风险评分 ===\n")

# 基于显著差异的菌群计算风险评分
# 有害菌群（疾病组升高）: Cluster 1, 2, 3, 8
# 有益菌群（疾病组降低）: Cluster 4, 5, 6, 9, 10

harmful_clusters <- c(1, 2, 3, 8)
protective_clusters <- c(4, 5, 6, 9, 10)

# 计算风险评分
risk_score <- (rowSums(cluster_abundance_rel[, harmful_clusters]) / 
               (rowSums(cluster_abundance_rel[, protective_clusters]) + 0.0001))

# 标准化风险评分
risk_score_scaled <- scale(log1p(risk_score))

# 添加到数据
risk_df <- data.frame(
  Sample = rownames(cluster_abundance_rel),
  Group = colData(tse_2w)$group,
  Risk_Score = risk_score_scaled[,1]
)

# 保存
write.csv(risk_df, "cluster/deep_analysis/data/risk_scores.csv", row.names = FALSE)

# 可视化风险评分
pdf("cluster/deep_analysis/figures/risk_score_boxplot.pdf", width = 6, height = 5)
boxplot(Risk_Score ~ Group, data = risk_df,
        main = "Disease Risk Score Based on Microbial Clusters",
        ylab = "Standardized Risk Score",
        xlab = "",
        col = c("lightblue", "lightcoral"))
wilcox_result <- wilcox.test(Risk_Score ~ Group, data = risk_df)
legend("topright", legend = paste("p =", format(wilcox_result$p.value, scientific = TRUE, digits = 3)))
dev.off()

cat("风险评分计算完成，p值:", wilcox_result$p.value, "\n")

# ==================== 3. 机器学习预测模型优化 ====================
cat("\n=== 3. 机器学习模型优化 ===\n")

# 准备数据
set.seed(123)
X <- as.matrix(cluster_abundance_rel)
y <- ifelse(colData(tse_2w)$group == "disease", 1, 0)

# 处理不平衡数据
sample_weights <- ifelse(y == 0, 1, sum(y == 0) / sum(y == 1))

# LASSO回归
cv_lasso <- cv.glmnet(X, y, family = "binomial", 
                      alpha = 1, 
                      weights = sample_weights,
                      nfolds = 10)

# 保存模型
saveRDS(cv_lasso, "cluster/deep_analysis/data/lasso_model.rds")

# 提取重要特征
lasso_coef <- coef(cv_lasso, s = "lambda.min")
important_features <- which(as.matrix(lasso_coef)[-1, 1] != 0)
lasso_importance <- data.frame(
  Cluster = important_features,
  Coefficient = as.matrix(lasso_coef)[important_features + 1, 1]
)
lasso_importance <- lasso_importance[order(-abs(lasso_importance$Coefficient)), ]
write.csv(lasso_importance, "cluster/deep_analysis/data/lasso_importance.csv", row.names = FALSE)

cat("LASSO选出的重要特征:\n")
print(lasso_importance)

# 可视化LASSO系数
pdf("cluster/deep_analysis/figures/lasso_coefficients.pdf", width = 8, height = 5)
plot(cv_lasso, main = "LASSO Cross-validation")
dev.off()

# ==================== 4. 菌群相互作用网络 ====================
cat("\n=== 4. 菌群相互作用网络 ===\n")

# 计算疾病组和健康组的网络差异
# 健康组相关性
cor_healthy <- cor(cluster_abundance_rel[colData(tse_2w)$group == "healthy", ], method = "spearman")
# 疾病组相关性
cor_disease <- cor(cluster_abundance_rel[colData(tse_2w)$group == "disease", ], method = "spearman")

# 相关性差异
cor_diff <- cor_disease - cor_healthy

# 保存
write.csv(cor_healthy, "cluster/deep_analysis/data/healthy_correlation.csv")
write.csv(cor_disease, "cluster/deep_analysis/data/disease_correlation.csv")
write.csv(cor_diff, "cluster/deep_analysis/data/correlation_difference.csv")

# 可视化差异网络
pdf("cluster/deep_analysis/figures/network_difference.pdf", width = 10, height = 8)
corrplot(cor_diff, method = "color", type = "upper",
         title = "Correlation Difference (Disease - Healthy)",
         mar = c(0,0,2,0),
         tl.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# 构建网络图
threshold <- 0.3
adj_diff <- abs(cor_diff) > threshold
diag(adj_diff) <- FALSE

if(sum(adj_diff) > 0) {
  g_diff <- graph_from_adjacency_matrix(adj_diff, mode = "undirected")
  V(g_diff)$name <- paste("Cluster", 1:10)
  
  # 计算网络中心性
  centrality <- data.frame(
    Cluster = 1:10,
    Degree = degree(g_diff),
    Betweenness = betweenness(g_diff),
    Closeness = closeness(g_diff)
  )
  write.csv(centrality, "cluster/deep_analysis/data/network_centrality.csv", row.names = FALSE)
  
  # 可视化网络
  pdf("cluster/deep_analysis/figures/network_difference_graph.pdf", width = 8, height = 7)
  plot(g_diff, 
       vertex.size = 20,
       vertex.label = V(g_diff)$name,
       vertex.color = rainbow(10),
       main = "Network of Correlation Differences (|Δr| > 0.3)")
  dev.off()
  
  cat("网络分析完成\n")
}


# ==================== 5. 样本聚类与疾病状态 ====================
cat("\n=== 5. 基于菌群的样本聚类 ===\n")

# 使用菌群组成进行样本聚类
dist_samples <- vegdist(cluster_abundance_rel, method = "bray")
hc_samples <- hclust(dist_samples, method = "ward.D2")

# 尝试不同的聚类数
sil_width <- sapply(2:10, function(k) {
  clusters <- cutree(hc_samples, k = k)
  sil <- silhouette(clusters, dist_samples)
  mean(sil[, 3])
})

# 最佳聚类数
best_k_samples <- which.max(sil_width) + 1
cat("最佳样本聚类数:", best_k_samples, "\n")

# 切割树状图
sample_clusters <- cutree(hc_samples, k = best_k_samples)

# 检查聚类与疾病状态的一致性
contingency <- table(Predicted_Cluster = sample_clusters, True_Group = colData(tse_2w)$group)
print(contingency)

# 调整兰德指数
library(mclust)
ari <- adjustedRandIndex(colData(tse_2w)$group, sample_clusters)
cat("调整兰德指数 (ARI):", ari, "\n")

# 可视化样本聚类
pdf("cluster/deep_analysis/figures/sample_clustering.pdf", width = 10, height = 6)
plot(hc_samples, labels = FALSE, main = "Sample Clustering Based on Microbial Clusters")
rect.hclust(hc_samples, k = best_k_samples, border = 2:(best_k_samples+1))
dev.off()

# ==================== 6. 菌群多样性分析 ====================
cat("\n=== 6. 菌群多样性分析 ===\n")

# 方法：手动计算多样性指标（不依赖vegan的diversity函数）

# 确保cluster_abundance_rel是矩阵格式
abundance_matrix <- as.matrix(cluster_abundance_rel)

# 计算Shannon多样性
# H = -sum(p_i * log(p_i))
calculate_shannon <- function(x) {
  # 移除零值
  x <- x[x > 0]
  if(length(x) == 0) return(0)
  p <- x / sum(x)
  -sum(p * log(p))
}

# 计算Simpson多样性
# D = 1 - sum(p_i^2)
calculate_simpson <- function(x) {
  x <- x[x > 0]
  if(length(x) == 0) return(0)
  p <- x / sum(x)
  1 - sum(p^2)
}

# 计算Pielou均匀度
# J = H / ln(S)
calculate_pielou <- function(x) {
  x <- x[x > 0]
  if(length(x) <= 1) return(0)
  S <- length(x)
  H <- calculate_shannon(x)
  H / log(S)
}

# 计算丰富度（非零菌群数）
calculate_richness <- function(x) {
  sum(x > 0)
}

# 应用到所有样本
shannon_diversity <- apply(abundance_matrix, 1, calculate_shannon)
simpson_diversity <- apply(abundance_matrix, 1, calculate_simpson)
pielou_evenness <- apply(abundance_matrix, 1, calculate_pielou)
richness <- apply(abundance_matrix, 1, calculate_richness)

# 创建多样性数据框
diversity_df <- data.frame(
  Sample = rownames(abundance_matrix),
  Group = colData(tse_2w)$group,
  Shannon = shannon_diversity,
  Simpson = simpson_diversity,
  Pielou = pielou_evenness,
  Richness = richness
)

# 保存多样性数据
write.csv(diversity_df, "cluster/deep_analysis/data/diversity_metrics.csv", row.names = FALSE)

# 查看数据摘要
cat("多样性数据计算完成\n")
print(head(diversity_df))

# ==================== 多样性统计检验 ====================

# Shannon多样性检验
shannon_test <- wilcox.test(Shannon ~ Group, data = diversity_df)
cat("\nShannon多样性检验 p值:", shannon_test$p.value, "\n")

# Simpson多样性检验
simpson_test <- wilcox.test(Simpson ~ Group, data = diversity_df)
cat("Simpson多样性检验 p值:", simpson_test$p.value, "\n")

# Pielou均匀度检验
pielou_test <- wilcox.test(Pielou ~ Group, data = diversity_df)
cat("Pielou均匀度检验 p值:", pielou_test$p.value, "\n")

# Richness检验
richness_test <- wilcox.test(Richness ~ Group, data = diversity_df)
cat("Richness检验 p值:", richness_test$p.value, "\n")

# ==================== 效应量计算 ====================
library(effsize)

cohen_d_shannon <- cohen.d(diversity_df$Shannon[diversity_df$Group == "healthy"],
                           diversity_df$Shannon[diversity_df$Group == "disease"])
cat("\nCohen's d效应量 (Shannon):", cohen_d_shannon$estimate, 
    "(", cohen_d_shannon$magnitude, ")\n")

# ==================== 可视化 ====================

# 1. 箱线图
pdf("cluster/deep_analysis/figures/diversity_boxplot.pdf", width = 12, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

boxplot(Shannon ~ Group, data = diversity_df,
        main = paste("Shannon Diversity\np =", format(shannon_test$p.value, scientific = TRUE, digits = 3)),
        ylab = "Shannon Diversity",
        xlab = "",
        col = c("lightblue", "lightcoral"),
        names = c("Healthy", "Disease"))

boxplot(Simpson ~ Group, data = diversity_df,
        main = paste("Simpson Diversity\np =", format(simpson_test$p.value, scientific = TRUE, digits = 3)),
        ylab = "Simpson Diversity",
        xlab = "",
        col = c("lightblue", "lightcoral"),
        names = c("Healthy", "Disease"))

boxplot(Pielou ~ Group, data = diversity_df,
        main = paste("Pielou Evenness\np =", format(pielou_test$p.value, scientific = TRUE, digits = 3)),
        ylab = "Evenness",
        xlab = "",
        col = c("lightblue", "lightcoral"),
        names = c("Healthy", "Disease"))

boxplot(Richness ~ Group, data = diversity_df,
        main = paste("Richness\np =", format(richness_test$p.value, scientific = TRUE, digits = 3)),
        ylab = "Number of Clusters",
        xlab = "",
        col = c("lightblue", "lightcoral"),
        names = c("Healthy", "Disease"))

dev.off()
cat("箱线图已保存\n")

# 2. 小提琴图
library(ggplot2)

diversity_long <- data.frame(
  Value = c(diversity_df$Shannon, diversity_df$Simpson, 
            diversity_df$Pielou, diversity_df$Richness),
  Index = rep(c("Shannon", "Simpson", "Pielou", "Richness"), 
              each = nrow(diversity_df)),
  Group = rep(diversity_df$Group, 4)
)

pdf("cluster/deep_analysis/figures/diversity_violin.pdf", width = 12, height = 6)
ggplot(diversity_long, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  facet_wrap(~ Index, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("healthy" = "lightblue", "disease" = "lightcoral")) +
  labs(title = "Alpha Diversity Metrics: Healthy vs Disease",
       x = "", y = "Diversity Value") +
  theme_minimal() +
  theme(legend.position = "bottom")
dev.off()
cat("小提琴图已保存\n")

# 3. 密度图
pdf("cluster/deep_analysis/figures/diversity_density.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))

# Shannon密度
plot(density(diversity_df$Shannon[diversity_df$Group == "healthy"]), 
     main = "Shannon Diversity Distribution",
     xlab = "Shannon Diversity", col = "blue", lwd = 2)
lines(density(diversity_df$Shannon[diversity_df$Group == "disease"]), 
      col = "red", lwd = 2)
legend("topright", legend = c("Healthy", "Disease"), 
       col = c("blue", "red"), lwd = 2)

# Simpson密度
plot(density(diversity_df$Simpson[diversity_df$Group == "healthy"]), 
     main = "Simpson Diversity Distribution",
     xlab = "Simpson Diversity", col = "blue", lwd = 2)
lines(density(diversity_df$Simpson[diversity_df$Group == "disease"]), 
      col = "red", lwd = 2)

# Pielou密度
plot(density(diversity_df$Pielou[diversity_df$Group == "healthy"]), 
     main = "Pielou Evenness Distribution",
     xlab = "Evenness", col = "blue", lwd = 2)
lines(density(diversity_df$Pielou[diversity_df$Group == "disease"]), 
      col = "red", lwd = 2)

# Richness密度
plot(density(diversity_df$Richness[diversity_df$Group == "healthy"]), 
     main = "Richness Distribution",
     xlab = "Number of Clusters", col = "blue", lwd = 2)
lines(density(diversity_df$Richness[diversity_df$Group == "disease"]), 
      col = "red", lwd = 2)

dev.off()
cat("密度图已保存\n")

# ==================== 7. ROC曲线比较 ====================
cat("\n=== 7. ROC曲线比较 ===\n")

# 确保有平衡的数据集
set.seed(123)

# 准备数据
if(!exists("cluster_abundance_rel")) {
  cat("错误：cluster_abundance_rel不存在，请先运行前面的分析\n")
} else {
  
  # 创建平衡数据集
  group <- colData(tse_2w)$group
  healthy_idx <- which(group == "healthy")
  disease_idx <- which(group == "disease")
  sampled_healthy <- sample(healthy_idx, length(disease_idx))
  balanced_idx <- c(sampled_healthy, disease_idx)
  balanced_abundance <- cluster_abundance_rel[balanced_idx, ]
  balanced_group <- group[balanced_idx]
  
  # 确保group是因子并重命名级别（避免randomForest问题）
  balanced_group <- factor(balanced_group, levels = c("healthy", "disease"))
  
  cat("平衡数据集大小:", nrow(balanced_abundance), "\n")
  cat("健康组:", sum(balanced_group == "healthy"), "\n")
  cat("疾病组:", sum(balanced_group == "disease"), "\n\n")
  
  # ==================== 模型1：使用所有10个菌群 ====================
  cat("训练模型1：使用所有10个菌群...\n")
  
  # 方法1：使用randomForest（确保y是因子）
  rf_model_all <- randomForest(x = balanced_abundance,
                               y = balanced_group,
                               ntree = 500,
                               importance = TRUE)
  
  # 预测概率
  pred_all <- predict(rf_model_all, type = "prob")
  
  # 检查pred_all的结构
  if(is.null(dim(pred_all))) {
    # 如果是二分类，pred_all应该是矩阵
    pred_all_disease <- as.numeric(pred_all == "disease")
  } else {
    pred_all_disease <- pred_all[, "disease"]
  }
  
  # 计算ROC
  roc_all <- roc(balanced_group, pred_all_disease)
  auc_all <- auc(roc_all)
  
  cat("模型1 AUC:", round(auc_all, 3), "\n")
  
  # ==================== 模型2：使用LASSO选出的重要菌群 ====================
  cat("\n运行LASSO特征选择...\n")
  
  # 准备LASSO数据（需要数值型响应变量）
  X <- as.matrix(balanced_abundance)
  y_numeric <- ifelse(balanced_group == "disease", 1, 0)
  
  # LASSO回归
  library(glmnet)
  set.seed(123)
  cv_lasso <- cv.glmnet(X, y_numeric, family = "binomial", alpha = 1, nfolds = 10)
  
  # 提取重要特征
  lasso_coef <- coef(cv_lasso, s = "lambda.min")
  important_features <- which(as.matrix(lasso_coef)[-1, 1] != 0)
  
  if(length(important_features) > 0) {
    cat("LASSO选出了", length(important_features), "个重要菌群:\n")
    cat("   Cluster", important_features, "\n")
    
    # 使用重要菌群重新训练随机森林
    lasso_data <- balanced_abundance[, important_features, drop = FALSE]
    rf_important <- randomForest(x = lasso_data,
                                 y = balanced_group,
                                 ntree = 500,
                                 importance = TRUE)
    
    pred_important <- predict(rf_important, type = "prob")
    if(is.null(dim(pred_important))) {
      pred_important_disease <- as.numeric(pred_important == "disease")
    } else {
      pred_important_disease <- pred_important[, "disease"]
    }
    
    roc_important <- roc(balanced_group, pred_important_disease)
    auc_important <- auc(roc_important)
    
    cat("模型2 AUC (仅重要菌群):", round(auc_important, 3), "\n")
    
    # 保存重要菌群信息
    lasso_importance <- data.frame(
      Cluster = important_features,
      Coefficient = as.matrix(lasso_coef)[important_features + 1, 1]
    )
    write.csv(lasso_importance, "cluster/deep_analysis/data/lasso_selected_clusters.csv", row.names = FALSE)
    
  } else {
    cat("LASSO没有选出重要特征\n")
    roc_important <- NULL
    auc_important <- NULL
  }
  
  # ==================== 模型3：只使用差异显著的菌群 ====================
  cat("\n训练模型3：只使用差异显著的菌群...\n")
  
  # 运行差异分析
  diff_pvalues <- sapply(1:ncol(balanced_abundance), function(i) {
    healthy_vals <- balanced_abundance[balanced_group == "healthy", i]
    disease_vals <- balanced_abundance[balanced_group == "disease", i]
    wilcox.test(healthy_vals, disease_vals)$p.value
  })
  
  sig_clusters <- which(diff_pvalues < 0.05)
  
  if(length(sig_clusters) > 0) {
    cat("差异显著的菌群:", length(sig_clusters), "个\n")
    cat("   Cluster", sig_clusters, "\n")
    
    sig_data <- balanced_abundance[, sig_clusters, drop = FALSE]
    rf_sig <- randomForest(x = sig_data,
                           y = balanced_group,
                           ntree = 500,
                           importance = TRUE)
    
    pred_sig <- predict(rf_sig, type = "prob")
    if(is.null(dim(pred_sig))) {
      pred_sig_disease <- as.numeric(pred_sig == "disease")
    } else {
      pred_sig_disease <- pred_sig[, "disease"]
    }
    
    roc_sig <- roc(balanced_group, pred_sig_disease)
    auc_sig <- auc(roc_sig)
    
    cat("模型3 AUC (仅差异显著):", round(auc_sig, 3), "\n")
  } else {
    cat("没有发现差异显著的菌群\n")
    roc_sig <- NULL
    auc_sig <- NULL
  }
  
  # ==================== 模型4：仅使用多样性指标 ====================
  cat("\n训练模型4：仅使用多样性指标...\n")
  
  # 计算多样性指标
  calculate_shannon <- function(x) {
    x <- x[x > 0]
    if(length(x) <= 1) return(0)
    p <- x / sum(x)
    -sum(p * log(p))
  }
  
  calculate_simpson <- function(x) {
    x <- x[x > 0]
    if(length(x) <= 1) return(0)
    p <- x / sum(x)
    1 - sum(p^2)
  }
  
  diversity_features <- data.frame(
    Shannon = apply(balanced_abundance, 1, calculate_shannon),
    Simpson = apply(balanced_abundance, 1, calculate_simpson),
    Richness = rowSums(balanced_abundance > 0)
  )
  
  rf_diversity <- randomForest(x = diversity_features,
                               y = balanced_group,
                               ntree = 500,
                               importance = TRUE)
  
  pred_diversity <- predict(rf_diversity, type = "prob")
  if(is.null(dim(pred_diversity))) {
    pred_diversity_disease <- as.numeric(pred_diversity == "disease")
  } else {
    pred_diversity_disease <- pred_diversity[, "disease"]
  }
  
  roc_diversity <- roc(balanced_group, pred_diversity_disease)
  auc_diversity <- auc(roc_diversity)
  
  cat("模型4 AUC (仅多样性):", round(auc_diversity, 3), "\n")
  
  # ==================== 绘制ROC曲线比较图 ====================
  pdf("cluster/deep_analysis/figures/roc_comparison.pdf", width = 10, height = 8)
  
  # 设置颜色和线型
  plot(roc_all, col = "blue", lwd = 2, 
       main = "ROC Curve Comparison of Different Models",
       xlab = "1 - Specificity", ylab = "Sensitivity")
  
  if(exists("roc_important") && !is.null(roc_important)) {
    plot(roc_important, col = "red", lwd = 2, add = TRUE)
  }
  
  if(exists("roc_sig") && !is.null(roc_sig)) {
    plot(roc_sig, col = "green", lwd = 2, add = TRUE)
  }
  
  plot(roc_diversity, col = "purple", lwd = 2, add = TRUE)
  
  # 添加对角线
  abline(a = 0, b = 1, lty = 2, col = "gray")
  
  # 添加图例
  legend_text <- c(paste("All 10 clusters (AUC =", round(auc_all, 3), ")"))
  legend_colors <- c("blue")
  
  if(exists("auc_important") && !is.null(auc_important)) {
    legend_text <- c(legend_text, paste("LASSO-selected (AUC =", round(auc_important, 3), ")"))
    legend_colors <- c(legend_colors, "red")
  }
  
  if(exists("auc_sig") && !is.null(auc_sig)) {
    legend_text <- c(legend_text, paste("Significant only (AUC =", round(auc_sig, 3), ")"))
    legend_colors <- c(legend_colors, "green")
  }
  
  legend_text <- c(legend_text, paste("Diversity only (AUC =", round(auc_diversity, 3), ")"))
  legend_colors <- c(legend_colors, "purple")
  
  legend("bottomright", legend = legend_text, 
         col = legend_colors, lwd = 2, cex = 0.9)
  
  dev.off()
  
  cat("\nROC曲线比较图已保存到: cluster/deep_analysis/figures/roc_comparison.pdf\n")
  
  # ==================== 保存模型比较结果 ====================
  model_comparison <- data.frame(
    Model = c("All 10 clusters", 
              if(exists("auc_important") && !is.null(auc_important)) "LASSO-selected" else NULL,
              if(exists("auc_sig") && !is.null(auc_sig)) "Significant only" else NULL,
              "Diversity only"),
    AUC = c(auc_all,
            if(exists("auc_important") && !is.null(auc_important)) auc_important else NULL,
            if(exists("auc_sig") && !is.null(auc_sig)) auc_sig else NULL,
            auc_diversity)
  )
  
  write.csv(model_comparison, "cluster/deep_analysis/data/model_comparison.csv", row.names = FALSE)
  
  # ==================== 打印总结 ====================
  cat("\n========================================\n")
  cat("模型性能总结\n")
  cat("========================================\n")
  print(model_comparison)
  
  # 找出最佳模型
  best_model <- model_comparison[which.max(model_comparison$AUC), ]
  cat("\n最佳模型:", best_model$Model, " (AUC =", round(best_model$AUC, 3), ")\n")
  
  # ==================== 特征重要性汇总 ====================
  cat("\n========================================\n")
  cat("特征重要性 (所有10个菌群模型)\n")
  cat("========================================\n")
  importance_df <- data.frame(
    Cluster = 1:10,
    Importance = importance(rf_model_all)[, "MeanDecreaseGini"]
  )
  importance_df <- importance_df[order(-importance_df$Importance), ]
  print(importance_df)
  
  # 保存特征重要性
  write.csv(importance_df, "cluster/deep_analysis/data/feature_importance_all.csv", row.names = FALSE)
  
  # ==================== 保存所有模型 ====================
  save(rf_model_all, rf_important, rf_sig, rf_diversity,
       file = "cluster/deep_analysis/data/all_models.RData")
  
  cat("\n所有模型已保存到: cluster/deep_analysis/data/all_models.RData\n")
  
  # ==================== 生成报告 ====================
  sink("cluster/deep_analysis/model_comparison_report.txt")
  
  cat("========================================\n")
  cat("模型比较报告\n")
  cat("========================================\n\n")
  
  cat("分析日期:", date(), "\n\n")
  
  cat("1. 数据集信息\n")
  cat(sprintf("   - 平衡数据集大小: %d\n", nrow(balanced_abundance)))
  cat(sprintf("   - 健康组: %d\n", sum(balanced_group == "healthy")))
  cat(sprintf("   - 疾病组: %d\n\n", sum(balanced_group == "disease")))
  
  cat("2. 模型性能\n")
  for(i in 1:nrow(model_comparison)) {
    cat(sprintf("   - %s: AUC = %.3f\n", 
                model_comparison$Model[i], 
                model_comparison$AUC[i]))
  }
  cat("\n")
  
  cat("3. 最佳模型\n")
  cat(sprintf("   %s (AUC = %.3f)\n\n", best_model$Model, best_model$AUC))
  
  cat("4. 最重要特征（所有菌群模型）\n")
  for(i in 1:min(3, nrow(importance_df))) {
    cat(sprintf("   - Cluster %d: 重要性 = %.2f\n", 
                importance_df$Cluster[i], 
                importance_df$Importance[i]))
  }
  
  sink()
  
  cat("\n模型比较报告已保存到: cluster/deep_analysis/model_comparison_report.txt\n")
}

cat("\n========================================\n")
cat("ROC曲线比较完成！\n")
cat("========================================\n")

# ==================== 8. 生成综合分析报告 ====================
cat("\n=== 8. 生成综合分析报告 ===\n")

sink("cluster/deep_analysis/deep_analysis_report.txt")

cat("========================================\n")
cat("微生物组菌群深度分析报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 菌群功能分类\n")
for(i in 1:10) {
  cat(sprintf("   Cluster %d: %s (%s)\n", 
              cluster_functions$Cluster[i],
              cluster_functions$Function[i],
              cluster_functions$Clinical_Association[i]))
}
cat("\n")

cat("2. 疾病风险评分\n")
cat(sprintf("   Wilcoxon检验p值: %.2e\n", wilcox_result$p.value))
cat("   风险评分公式: (有害菌群) / (有益菌群)\n")
cat("   有害菌群: Cluster 1,2,3,8\n")
cat("   有益菌群: Cluster 4,5,6,9,10\n\n")

cat("3. 机器学习模型\n")
cat(sprintf("   LASSO选出的重要特征: %d个\n", nrow(lasso_importance)))
if(nrow(lasso_importance) > 0) {
  cat("   最重要的特征:\n")
  for(i in 1:min(3, nrow(lasso_importance))) {
    cat(sprintf("     - Cluster %d (系数: %.3f)\n", 
                lasso_importance$Cluster[i], 
                lasso_importance$Coefficient[i]))
  }
}
cat("\n")

cat("4. 菌群多样性\n")
cat(sprintf("   Shannon多样性: p = %.2e\n", shannon_test$p.value))
cat(sprintf("   Simpson多样性: p = %.2e\n", simpson_test$p.value))
cat("   结论: 疾病组与健康组的菌群多样性", 
    ifelse(shannon_test$p.value < 0.05, "有显著差异", "无显著差异"), "\n\n")


cat("5. 输出文件\n")
cat("   - 数据文件: cluster/deep_analysis/data/\n")
cat("   - 图表文件: cluster/deep_analysis/figures/\n")
cat("   - 报告文件: cluster/deep_analysis/deep_analysis_report.txt\n")

sink()

cat("综合分析报告已保存\n\n")

# ==================== 完成 ====================
cat("========================================\n")
cat("深入分析完成！\n")
cat("结果已保存到 'cluster/deep_analysis' 目录\n")
cat("========================================\n")