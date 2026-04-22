# ==================== 特征选择高级可视化 ====================
cat("\n=== 特征选择高级可视化 ===\n")

# 确保目录存在
if(!dir.exists("cluster/classifier/figures")) {
  dir.create("cluster/classifier/figures", recursive = TRUE)
}

library(ggplot2)
library(reshape2)
library(corrplot)
library(ggrepel)
library(ggthemes)

# ==================== 1. 特征重要性比较图 ====================
cat("\n生成特征重要性比较图...\n")

# 准备数据
importance_df <- data.frame(
  Feature = rownames(rf_importance_model$importance),
  Importance = rf_importance_model$importance[, "MeanDecreaseGini"],
  Type = ifelse(rownames(rf_importance_model$importance) %in% best_biomarkers, 
                "Selected", "Not Selected")
)

# 排序
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
importance_df$Feature <- factor(importance_df$Feature, levels = importance_df$Feature)

# 水平条形图 - 带颜色区分
p1 <- ggplot(importance_df, aes(x = Feature, y = Importance, fill = Type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Selected" = "#E41A1C", "Not Selected" = "#377EB8")) +
  labs(title = "Feature Importance for Disease Classification",
       subtitle = paste("Selected biomarkers:", paste(best_biomarkers, collapse = ", ")),
       x = "Feature", y = "Mean Decrease Gini") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "bottom")

ggsave("cluster/classifier/figures/feature_importance_selected.pdf", p1, width = 10, height = 6)

# ==================== 2. 特征选择方法比较图 ====================
cat("\n生成特征选择方法比较图...\n")

# 收集各方法选中的特征
method_features <- list(
  "RF Importance" = best_features,
  "Stepwise" = best_features_step,
  "LASSO" = best_features_lasso,
  "Combined" = best_combined_features
)

# 创建热图矩阵
all_features_unique <- unique(unlist(method_features))
method_matrix <- matrix(0, nrow = length(all_features_unique), ncol = length(method_features))
rownames(method_matrix) <- all_features_unique
colnames(method_matrix) <- names(method_features)

for(i in 1:length(method_features)) {
  method_matrix[method_features[[i]], i] <- 1
}

# 按行和排序
row_order <- order(rowSums(method_matrix), decreasing = TRUE)
method_matrix <- method_matrix[row_order, ]

# 热图
pdf("cluster/classifier/figures/feature_selection_comparison.pdf", width = 8, height = 10)
heatmap(method_matrix, 
        Rowv = NA, Colv = NA,
        col = c("white", "steelblue"),
        scale = "none",
        main = "Feature Selection Methods Comparison",
        xlab = "Method", ylab = "Feature",
        margins = c(8, 10))
dev.off()

# 使用ggplot2版本
method_long <- melt(method_matrix)
colnames(method_long) <- c("Feature", "Method", "Selected")

p2 <- ggplot(method_long, aes(x = Method, y = Feature, fill = as.factor(Selected))) +
  geom_tile() +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue"), 
                    labels = c("Not Selected", "Selected"),
                    name = "Selection Status") +
  labs(title = "Feature Selection Methods Comparison",
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

ggsave("cluster/classifier/figures/feature_selection_heatmap.pdf", p2, width = 8, height = 10)

# ==================== 3. 不同特征集性能雷达图 ====================
cat("\n=== 性能雷达图 ===\n")

# 确保目录存在
if(!dir.exists("cluster/classifier/figures")) {
  dir.create("cluster/classifier/figures", recursive = TRUE)
}

library(ggplot2)
library(reshape2)
library(pheatmap)

# ==================== 方法1：使用现有性能数据 ====================
cat("\n方法1：使用现有性能数据\n")

# 检查performance_results是否存在
if(!exists("performance_results") || nrow(performance_results) == 0) {
  cat("重新计算性能数据...\n")
  
  # 重新计算性能数据
  feature_sets <- list(
    "All Features" = selected_features,
    "Top 3" = head(importance_df$Feature, 3),
    "Top 5" = head(importance_df$Feature, 5)
  )
  
  performance_results <- data.frame()
  
  for(set_name in names(feature_sets)) {
    features <- feature_sets[[set_name]]
    
    if(length(features) >= 2) {
      X_subset <- X_selected[, features, drop = FALSE]
      
      set.seed(123)
      cv_control <- trainControl(method = "cv", number = 5, 
                                 summaryFunction = defaultSummary)
      
      model_cv <- tryCatch({
        train(x = X_subset, y = y,
              method = "rf", ntree = 200,
              trControl = cv_control)
      }, error = function(e) NULL)
      
      if(!is.null(model_cv)) {
        best_result <- model_cv$results[which.max(model_cv$results$Accuracy), ]
        
        performance_results <- rbind(performance_results, data.frame(
          Feature_Set = set_name,
          N_Features = length(features),
          Accuracy = best_result$Accuracy,
          Kappa = best_result$Kappa
        ))
      }
    }
  }
}

# 打印性能结果
print(performance_results)

# ==================== 创建雷达图数据 ====================
cat("\n创建雷达图数据...\n")

# 定义性能指标
metrics <- c("Accuracy", "Kappa", "Sensitivity", "Specificity", "Precision")

# 为每个特征集计算完整性能指标
radar_data <- data.frame(Metric = metrics)

for(i in 1:nrow(performance_results)) {
  set_name <- performance_results$Feature_Set[i]
  features <- feature_sets[[set_name]]
  X_subset <- X_selected[, features, drop = FALSE]
  
  # 计算每个类别的性能
  set.seed(123)
  train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
  
  train_X <- X_subset[train_idx, ]
  train_y <- y[train_idx]
  test_X <- X_subset[-train_idx, ]
  test_y <- y[-train_idx]
  
  # 训练模型
  rf_model_temp <- randomForest(x = train_X, y = train_y, ntree = 200)
  
  # 预测
  pred <- predict(rf_model_temp, test_X)
  conf_mat <- confusionMatrix(pred, test_y)
  
  # 计算宏平均性能
  sensitivity <- mean(conf_mat$byClass[, "Sensitivity"], na.rm = TRUE)
  specificity <- mean(conf_mat$byClass[, "Specificity"], na.rm = TRUE)
  precision <- mean(conf_mat$byClass[, "Precision"], na.rm = TRUE)
  
  radar_data[[set_name]] <- c(
    performance_results$Accuracy[i],
    performance_results$Kappa[i],
    sensitivity,
    specificity,
    precision
  )
}

# 处理缺失值
radar_data[is.na(radar_data)] <- 0

# ==================== 雷达图1：使用ggplot2 ====================
cat("\n生成ggplot2雷达图...\n")

# 转换为长格式
radar_long <- melt(radar_data, id.vars = "Metric", variable.name = "Feature_Set", value.name = "Score")

# 创建雷达图（极坐标）
p1 <- ggplot(radar_long, aes(x = Metric, y = Score, color = Feature_Set, group = Feature_Set)) +
  geom_polygon(aes(fill = Feature_Set), alpha = 0.2, size = 1) +
  geom_point(size = 2) +
  coord_polar() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Performance Comparison of Feature Sets",
       subtitle = "Radar Chart (Polar Coordinates)",
       x = "", y = "Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 10),
        legend.position = "bottom")

ggsave("cluster/classifier/figures/performance_radar_ggplot.pdf", p1, width = 10, height = 8)

# ==================== 雷达图2：使用平行坐标图 ====================
cat("\n生成平行坐标图...\n")

p2 <- ggplot(radar_long, aes(x = Metric, y = Score, color = Feature_Set, group = Feature_Set)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Performance Comparison of Feature Sets",
       subtitle = "Parallel Coordinates Plot",
       x = "", y = "Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("cluster/classifier/figures/performance_parallel.pdf", p2, width = 10, height = 6)

# ==================== 雷达图3：使用条形图矩阵 ====================
cat("\n生成条形图矩阵...\n")

p3 <- ggplot(radar_long, aes(x = Metric, y = Score, fill = Feature_Set)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Performance Comparison of Feature Sets",
       subtitle = "Grouped Bar Chart",
       x = "", y = "Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("cluster/classifier/figures/performance_bars.pdf", p3, width = 10, height = 6)

# ==================== 雷达图4：热图 ====================
cat("\n生成性能热图...\n")

# 创建矩阵
performance_matrix <- as.matrix(radar_data[, -1])
rownames(performance_matrix) <- radar_data$Metric

p4 <- pheatmap(performance_matrix,
               main = "Performance Comparison of Feature Sets",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               display_numbers = TRUE,
               number_format = "%.3f",
               color = colorRampPalette(c("white", "steelblue"))(50),
               fontsize_row = 10,
               fontsize_col = 10)

ggsave("cluster/classifier/figures/performance_heatmap.pdf", p4, width = 8, height = 6)

# ==================== 雷达图5：蜘蛛图（使用ggradar）====================
cat("\n尝试使用ggradar包...\n")

# 检查是否安装了ggradar
if(require(ggradar)) {
  # 准备数据
  radar_ggradar <- radar_data
  colnames(radar_ggradar)[1] <- "group"
  
  # 创建蜘蛛图
  p5 <- ggradar(radar_ggradar,
                values.radar = c(0, 0.5, 1),
                grid.min = 0, grid.mid = 0.5, grid.max = 1,
                group.line.width = 1.5,
                group.point.size = 3,
                legend.position = "bottom")
  
  ggsave("cluster/classifier/figures/performance_spider.pdf", p5, width = 10, height = 8)
} else {
  cat("ggradar包未安装，跳过\n")
}

# ==================== 雷达图6：简单的表格 + 条形图 ====================
cat("\n生成性能汇总表...\n")

# 创建性能汇总表
performance_table <- radar_data
performance_table[, -1] <- round(performance_table[, -1], 3)

# 保存为CSV
write.csv(performance_table, "cluster/classifier/data/performance_summary.csv", row.names = FALSE)

# 创建表格图
library(gridExtra)
table_grob <- tableGrob(performance_table, rows = NULL,
                        theme = ttheme_minimal(
                          core = list(fg_params = list(hjust = 0, x = 0.1)),
                          colhead = list(fg_params = list(fontface = "bold"))
                        ))

p6 <- ggplot() +
  annotation_custom(table_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  labs(title = "Performance Metrics by Feature Set") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave("cluster/classifier/figures/performance_table.pdf", p6, width = 10, height = 4)

# ==================== 生成性能比较报告 ====================
sink("cluster/classifier/performance_report.txt")

cat("========================================\n")
cat("特征集性能比较报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 特征集性能汇总\n\n")
print(performance_table)
cat("\n")

cat("2. 最佳特征集\n")
best_set <- performance_results[which.max(performance_results$Accuracy), ]
cat(sprintf("   - 特征集: %s\n", best_set$Feature_Set))
cat(sprintf("   - 特征数量: %d\n", best_set$N_Features))
cat(sprintf("   - 准确率: %.4f\n", best_set$Accuracy))
cat(sprintf("   - Kappa: %.4f\n", best_set$Kappa))
cat("\n")

cat("3. 各特征集性能排名\n")
performance_results <- performance_results[order(-performance_results$Accuracy), ]
for(i in 1:nrow(performance_results)) {
  cat(sprintf("   %d. %s: 准确率 = %.4f, Kappa = %.4f\n", 
              i, performance_results$Feature_Set[i], 
              performance_results$Accuracy[i],
              performance_results$Kappa[i]))
}

sink()

cat("\n报告已保存到: cluster/classifier/performance_report.txt\n")

# ==================== 打印摘要 ====================
cat("\n\n========================================\n")
cat("性能雷达图生成完成！\n")
cat("========================================\n")
cat("生成的文件:\n")
cat("  - performance_radar_ggplot.pdf (极坐标雷达图)\n")
cat("  - performance_parallel.pdf (平行坐标图)\n")
cat("  - performance_bars.pdf (分组条形图)\n")
cat("  - performance_heatmap.pdf (性能热图)\n")
cat("  - performance_table.pdf (性能表格)\n")
cat("  - performance_summary.csv (性能数据)\n")
cat("  - performance_report.txt (分析报告)\n")

# 打印最佳结果
cat("\n最佳特征集: ", best_set$Feature_Set, "\n")
cat("准确率: ", round(best_set$Accuracy, 4), "\n")

# ==================== 4. 生物标志物相关性网络 ====================
cat("\n生成生物标志物相关性网络...\n")

if(length(best_biomarkers) >= 3) {
  # 计算生物标志物之间的相关性
  biomarker_cor <- cor(X_selected[, best_biomarkers], method = "spearman")
  
  # 相关性热图
  pdf("cluster/classifier/figures/biomarker_correlation.pdf", width = 8, height = 7)
  corrplot(biomarker_cor, method = "color", type = "upper",
           addCoef.col = "black", tl.col = "black",
           title = "Biomarker Correlation Network",
           mar = c(0,0,1,0))
  dev.off()
  
  # 网络图
  library(igraph)
  threshold <- 0.3
  adj_matrix <- abs(biomarker_cor) > threshold
  diag(adj_matrix) <- FALSE
  
  if(sum(adj_matrix) > 0) {
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    
    # 添加节点属性
    V(g)$size <- importance_df$Importance[match(V(g)$name, importance_df$Feature)] * 2
    V(g)$color <- "lightblue"
    
    pdf("cluster/classifier/figures/biomarker_network.pdf", width = 8, height = 7)
    plot(g, 
         vertex.size = V(g)$size,
         vertex.label = V(g)$name,
         vertex.label.cex = 0.8,
         vertex.label.dist = 1,
         edge.width = abs(biomarker_cor[adj_matrix]) * 2,
         edge.color = ifelse(biomarker_cor[adj_matrix] > 0, "red", "blue"),
         main = "Biomarker Co-occurrence Network")
    dev.off()
  }
}

# ==================== 5. 生物标志物箱线图矩阵 ====================
cat("\n生成生物标志物箱线图矩阵...\n")

if(length(best_biomarkers) >= 2) {
  # 准备数据
  biomarker_box_data <- data.frame(
    Disease = y,
    X_selected[, best_biomarkers, drop = FALSE]
  )
  
  biomarker_long <- melt(biomarker_box_data, 
                         id.vars = "Disease", 
                         variable.name = "Biomarker", 
                         value.name = "Abundance")
  
  # 箱线图矩阵
  p3 <- ggplot(biomarker_long, aes(x = Disease, y = Abundance, fill = Disease)) +
    geom_boxplot() +
    facet_wrap(~ Biomarker, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = rainbow(length(levels(y)))) +
    labs(title = "Distribution of Selected Biomarkers by Disease",
         x = "Disease", y = "Relative Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  ggsave("cluster/classifier/figures/biomarker_boxplots.pdf", p3, width = 12, height = 8)
}

# ==================== 6. 特征选择稳定性图 ====================
cat("\n生成特征选择稳定性图...\n")

# Bootstrap特征选择稳定性
set.seed(123)
n_bootstrap <- 50
feature_stability <- matrix(0, nrow = ncol(X_selected), ncol = n_bootstrap)
colnames(feature_stability) <- paste0("Bootstrap", 1:n_bootstrap)
rownames(feature_stability) <- colnames(X_selected)

for(b in 1:n_bootstrap) {
  # Bootstrap采样
  boot_idx <- sample(1:nrow(X_selected), replace = TRUE)
  boot_X <- X_selected[boot_idx, ]
  boot_y <- y[boot_idx]
  
  # 随机森林特征重要性
  boot_rf <- randomForest(x = boot_X, y = boot_y, ntree = 200, importance = TRUE)
  boot_imp <- importance(boot_rf)[, "MeanDecreaseGini"]
  
  # 记录是否入选Top 3
  top_features <- names(sort(boot_imp, decreasing = TRUE)[1:3])
  feature_stability[top_features, b] <- 1
}

# 计算稳定性得分
stability_scores <- rowSums(feature_stability) / n_bootstrap
stability_df <- data.frame(
  Feature = names(stability_scores),
  Stability = stability_scores
)
stability_df <- stability_df[order(-stability_df$Stability), ]

# 稳定性图
p4 <- ggplot(stability_df, aes(x = reorder(Feature, Stability), y = Stability)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Selection Stability (Bootstrap)",
       subtitle = "Proportion of times feature was in Top 3",
       x = "Feature", y = "Stability Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("cluster/classifier/figures/feature_stability.pdf", p4, width = 8, height = 5)

# ==================== 7. 生物标志物与疾病的热图 ====================
cat("\n生成生物标志物-疾病热图...\n")

# 计算每个生物标志物在各类疾病中的平均丰度
biomarker_means <- aggregate(X_selected[, best_biomarkers], 
                             by = list(Disease = y), 
                             FUN = mean)
rownames(biomarker_means) <- biomarker_means$Disease
biomarker_means <- biomarker_means[, -1]

# 标准化
biomarker_zscore <- scale(biomarker_means)

# 热图
pdf("cluster/classifier/figures/biomarker_disease_heatmap.pdf", width = 8, height = 6)
pheatmap(biomarker_zscore,
         main = "Biomarker Expression by Disease",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize_row = 10,
         fontsize_col = 8)
dev.off()

# ==================== 8. 特征选择性能曲线 ====================
cat("\n生成特征选择性能曲线...\n")

# 测试不同特征数量下的性能
max_features <- min(10, ncol(X_selected))
feature_performance <- data.frame()

for(n in 1:max_features) {
  top_features <- importance_df$Feature[1:n]
  X_subset <- X_selected[, top_features, drop = FALSE]
  
  set.seed(123)
  cv_control <- trainControl(method = "cv", number = 5)
  model_cv <- train(x = X_subset, y = y,
                    method = "rf", ntree = 200,
                    trControl = cv_control)
  
  feature_performance <- rbind(feature_performance, data.frame(
    N_Features = n,
    Accuracy = max(model_cv$results$Accuracy)
  ))
}

# 性能曲线
p5 <- ggplot(feature_performance, aes(x = N_Features, y = Accuracy)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = length(best_biomarkers), linetype = "dashed", color = "red") +
  annotate("text", x = length(best_biomarkers) + 0.5, y = max(feature_performance$Accuracy), 
           label = "Selected Set", color = "red") +
  labs(title = "Classification Performance vs Number of Features",
       x = "Number of Features", y = "Accuracy") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("cluster/classifier/figures/performance_curve.pdf", p5, width = 8, height = 5)

# ==================== 9. 生物标志物重要性排名图 ====================
cat("\n生成生物标志物重要性排名图...\n")

# 创建排名图
importance_rank <- importance_df
importance_rank$Rank <- 1:nrow(importance_rank)
importance_rank$Selected <- importance_rank$Feature %in% best_biomarkers

p6 <- ggplot(importance_rank, aes(x = Rank, y = Importance, color = Selected, label = Feature)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1), color = "gray", alpha = 0.5) +
  geom_text_repel(data = subset(importance_rank, Selected == TRUE), 
                  aes(label = Feature), size = 3, nudge_x = 2) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  labs(title = "Feature Importance Ranking",
       subtitle = paste("Selected biomarkers:", paste(best_biomarkers, collapse = ", ")),
       x = "Rank", y = "Importance (Mean Decrease Gini)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

ggsave("cluster/classifier/figures/importance_ranking.pdf", p6, width = 10, height = 6)

# ==================== 10. 生成综合报告 ====================
sink("cluster/classifier/visualization_report.txt")

cat("========================================\n")
cat("特征选择与生物标志物可视化报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 生成的可视化文件\n")
cat("   - feature_importance_selected.pdf: 特征重要性（标记选中特征）\n")
cat("   - feature_selection_comparison.pdf: 各方法特征选择比较\n")
cat("   - feature_selection_heatmap.pdf: 特征选择热图\n")
cat("   - performance_radar.pdf: 不同特征集性能雷达图\n")
cat("   - biomarker_correlation.pdf: 生物标志物相关性热图\n")
cat("   - biomarker_network.pdf: 生物标志物网络图\n")
cat("   - biomarker_boxplots.pdf: 生物标志物箱线图矩阵\n")
cat("   - feature_stability.pdf: 特征选择稳定性图\n")
cat("   - biomarker_disease_heatmap.pdf: 生物标志物-疾病热图\n")
cat("   - performance_curve.pdf: 性能曲线\n")
cat("   - importance_ranking.pdf: 特征重要性排名图\n")
cat("\n")

cat("2. 最佳生物标志物组合\n")
cat(sprintf("   - 特征数量: %d\n", length(best_biomarkers)))
cat("   - 生物标志物:\n")
for(i in 1:length(best_biomarkers)) {
  cat(sprintf("     %d. %s\n", i, best_biomarkers[i]))
}
cat("\n")

cat("3. 性能指标\n")
best_perf <- performance_results[performance_results$Feature_Set == best_set$Feature_Set, ]
cat(sprintf("   - 准确率: %.4f\n", best_perf$Accuracy))
cat(sprintf("   - Kappa: %.4f\n", best_perf$Kappa))

sink()

cat("\n报告已保存到: cluster/classifier/visualization_report.txt\n")
cat("\n所有可视化已保存到: cluster/classifier/figures/\n")

# ==================== 打印摘要 ====================
cat("\n\n========================================\n")
cat("可视化生成完成！\n")
cat("========================================\n")
cat("生成的可视化文件:\n")
cat("  1. 特征重要性比较图\n")
cat("  2. 特征选择方法比较热图\n")
cat("  3. 性能雷达图\n")
cat("  4. 生物标志物相关性网络\n")
cat("  5. 生物标志物箱线图矩阵\n")
cat("  6. 特征选择稳定性图\n")
cat("  7. 生物标志物-疾病热图\n")
cat("  8. 特征选择性能曲线\n")
cat("  9. 特征重要性排名图\n")
cat("\n最佳生物标志物组合:", paste(best_biomarkers, collapse = ", "))