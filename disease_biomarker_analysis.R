# ==================== 疾病特异性生物标志物识别 ====================
cat("\n=== 疾病特异性生物标志物识别 ===\n")

# 创建分析目录
if(!dir.exists("cluster/biomarker_analysis")) {
  dir.create("cluster/biomarker_analysis")
}
if(!dir.exists("cluster/biomarker_analysis/figures")) {
  dir.create("cluster/biomarker_analysis/figures")
}
if(!dir.exists("cluster/biomarker_analysis/data")) {
  dir.create("cluster/biomarker_analysis/data")
}

library(ggplot2)

# 确保必要的对象存在
if(!exists("cluster_abundance_rel")) {
  cat("创建菌群组成数据...\n")
  counts <- assay(tse_2w, "counts")
  if(!exists("species_clusters")) {
    species_cor <- cor(t(counts), method = "spearman")
    dist_species <- as.dist(1 - abs(species_cor))
    hc_species <- hclust(dist_species, method = "ward.D2")
    species_clusters <- cutree(hc_species, k = 10)
  }
  
  cluster_abundance <- matrix(0, ncol(tse_2w), 10)
  colnames(cluster_abundance) <- paste0("Cluster", 1:10)
  rownames(cluster_abundance) <- colnames(tse_2w)
  
  for(i in 1:10) {
    cluster_species <- names(species_clusters[species_clusters == i])
    cluster_abundance[, i] <- colSums(counts[cluster_species, ])
  }
  cluster_abundance_rel <- cluster_abundance / rowSums(cluster_abundance)
}

# 获取简化的疾病分类
if(!"condition_simplified" %in% colnames(colData(tse_2w))) {
  # 创建简化分类
  condition <- colData(tse_2w)$condition
  simplified <- condition
  
  simplify_map <- list(
    "healthy" = "Healthy",
    "diabesity|diabetes|T2D" = "T2D",
    "gestational_diabetes" = "Gestational_DM",
    "MS|MS_obesity" = "Metabolic_Syndrome",
    "NAFLD|NAFLD_obesity" = "NAFLD",
    "obesity" = "Obesity",
    "overweight" = "Overweight",
    "bariatric_surgery" = "Post_Surgery",
    "asthma_obesity" = "Asthma_Obesity",
    "opioid_use" = "Opioid_Use"
  )
  
  for(pattern in names(simplify_map)) {
    matched <- grepl(pattern, simplified, ignore.case = TRUE)
    simplified[matched] <- simplify_map[[pattern]]
  }
  colData(tse_2w)$condition_simplified <- simplified
}

# ==================== 1. 疾病特异性生物标志物 ====================
cat("\n=== 1. 识别疾病特异性生物标志物 ===\n")

# 获取所有疾病类型（排除健康组）
disease_types <- unique(colData(tse_2w)$condition_simplified)
disease_types <- disease_types[disease_types != "Healthy" & !is.na(disease_types)]
disease_types <- disease_types[sapply(disease_types, function(d) sum(colData(tse_2w)$condition_simplified == d) >= 10)]

cat("分析的疾病类型:", paste(disease_types, collapse = ", "), "\n")

# 存储结果
biomarker_results <- list()

for(disease in disease_types) {
  cat("\n分析:", disease, "\n")
  
  # 该疾病组 vs 所有其他疾病组（不包括健康组）
  idx_this <- colData(tse_2w)$condition_simplified == disease
  idx_other <- colData(tse_2w)$condition_simplified %in% disease_types & 
              colData(tse_2w)$condition_simplified != disease
  
  if(sum(idx_this) < 5 || sum(idx_other) < 5) next
  
  disease_results <- data.frame()
  
  for(cluster in 1:10) {
    vals_this <- cluster_abundance_rel[idx_this, cluster]
    vals_other <- cluster_abundance_rel[idx_other, cluster]
    
    # Wilcoxon检验
    test <- wilcox.test(vals_this, vals_other)
    
    # 计算效应量（AUC）
    auc_value <- tryCatch({
      roc(rep(c(1,0), c(length(vals_this), length(vals_other))), 
          c(vals_this, vals_other))$auc[1]
    }, error = function(e) NA)
    
    # 计算fold change
    fc <- (mean(vals_this) + 0.0001) / (mean(vals_other) + 0.0001)
    log2fc <- log2(fc)
    
    disease_results <- rbind(disease_results, data.frame(
      Disease = disease,
      Cluster = cluster,
      Mean_This = mean(vals_this),
      Mean_Other = mean(vals_other),
      Log2FC = log2fc,
      FC = fc,
      AUC = auc_value,
      P_value = test$p.value,
      N_This = sum(idx_this),
      N_Other = sum(idx_other)
    ))
  }
  
  disease_results$FDR <- p.adjust(disease_results$P_value, method = "fdr")
  disease_results$Significant <- disease_results$FDR < 0.05
  disease_results$Biomarker_Type <- ifelse(disease_results$Log2FC > 0, "Enriched", "Depleted")
  
  biomarker_results[[disease]] <- disease_results
  
  # 显示显著标志物
  sig_markers <- disease_results[disease_results$Significant, ]
  if(nrow(sig_markers) > 0) {
    cat("  显著标志物:\n")
    for(i in 1:nrow(sig_markers)) {
      cat(sprintf("    Cluster %d: %s (Log2FC = %.2f, p = %.2e)\n",
                  sig_markers$Cluster[i],
                  sig_markers$Biomarker_Type[i],
                  sig_markers$Log2FC[i],
                  sig_markers$P_value[i]))
    }
  } else {
    cat("  未发现显著标志物\n")
  }
}

# 合并所有结果
all_biomarkers <- do.call(rbind, biomarker_results)
write.csv(all_biomarkers, "cluster/biomarker_analysis/data/all_biomarkers.csv", row.names = FALSE)

# ==================== 2. 疾病特异性标志物可视化 ====================

# 热图：疾病特异性标志物
if(nrow(all_biomarkers) > 0) {
  library(reshape2)
  
  # 创建Log2FC矩阵
  fc_matrix <- acast(all_biomarkers, Disease ~ Cluster, value.var = "Log2FC")
  
  # 创建显著性矩阵
  sig_matrix <- acast(all_biomarkers, Disease ~ Cluster, value.var = "Significant")
  
  # 热图
  pdf("cluster/biomarker_analysis/figures/disease_specific_biomarkers.pdf", 
      width = 10, height = 8)
  pheatmap(fc_matrix,
           main = "Disease-Specific Microbial Biomarkers",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = sig_matrix,
           number_format = function(x) ifelse(x, "*", ""),
           color = colorRampPalette(c("blue", "white", "red"))(50),
           fontsize_row = 10,
           fontsize_col = 8)
  dev.off()
  
  # 条形图：每个疾病最重要的标志物
  top_markers <- all_biomarkers[all_biomarkers$Significant, ]
  top_markers <- top_markers[order(abs(top_markers$Log2FC), decreasing = TRUE), ]
  top_markers <- top_markers[!duplicated(paste(top_markers$Disease, top_markers$Cluster)), ]
  
  if(nrow(top_markers) > 0) {
    pdf("cluster/biomarker_analysis/figures/top_biomarkers.pdf", width = 10, height = 6)
    ggplot(top_markers, aes(x = reorder(paste(Disease, Cluster), Log2FC), 
                            y = Log2FC, fill = Biomarker_Type)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("Enriched" = "red", "Depleted" = "blue")) +
      labs(title = "Top Disease-Specific Microbial Biomarkers",
           x = "Disease - Cluster", y = "Log2 Fold Change") +
      theme_minimal()
    dev.off()
  }
}

# ==================== 3. 疾病严重程度梯度分析 ====================
cat("\n\n=== 3. 疾病严重程度梯度分析 ===\n")

# 定义疾病严重程度顺序（基于临床认知）
severity_order <- c(
  "Healthy",
  "Overweight",
  "Gestational_DM",
  "Pre_Diabetes_Overweight",
  "Obesity",
  "T2D",
  "Metabolic_Syndrome",
  "NAFLD",
  "Asthma_Obesity",
  "Opioid_Use",
  "Post_Surgery"
)

# 只保留存在的疾病类型
existing_diseases <- unique(colData(tse_2w)$condition_simplified)
severity_order <- severity_order[severity_order %in% existing_diseases]

cat("严重程度梯度顺序:\n")
for(i in seq_along(severity_order)) {
  cat(sprintf("  %d. %s\n", i, severity_order[i]))
}

# ==================== 4. 计算疾病进展评分 ====================

# 使用PCA第一主成分作为综合评分
pca_result <- prcomp(cluster_abundance_rel, scale. = TRUE)
disease_score <- pca_result$x[, 1]

# 创建评分数据框
score_df <- data.frame(
  Sample = rownames(cluster_abundance_rel),
  Condition = colData(tse_2w)$condition_simplified,
  PC1_Score = disease_score
)

# 按严重程度排序
score_df$Severity_Level <- factor(score_df$Condition, levels = severity_order)
score_df <- score_df[!is.na(score_df$Severity_Level), ]

# 计算每个疾病类型的平均评分
severity_means <- aggregate(PC1_Score ~ Severity_Level, data = score_df, 
                            FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
severity_means <- do.call(data.frame, severity_means)
colnames(severity_means) <- c("Disease", "Mean_Score", "SD_Score", "N")

# 保存
write.csv(severity_means, "cluster/biomarker_analysis/data/severity_scores.csv", row.names = FALSE)

# ==================== 5. 梯度可视化 ====================

# 箱线图
pdf("cluster/biomarker_analysis/figures/disease_severity_gradient.pdf", 
    width = 12, height = 6)
boxplot(PC1_Score ~ Severity_Level, data = score_df,
        main = "Disease Severity Gradient Based on Microbial Composition",
        xlab = "Disease Severity (Increasing →)",
        ylab = "PC1 Score (Microbial Dysbiosis Index)",
        col = colorRampPalette(c("green", "yellow", "red"))(length(severity_order)),
        las = 2,
        cex.axis = 0.8)
abline(h = 0, lty = 2, col = "gray", lwd = 2)
dev.off()

# 点线图（均值和误差线）
pdf("cluster/biomarker_analysis/figures/severity_trend.pdf", width = 10, height = 6)
plot(1:length(severity_means$Mean_Score), severity_means$Mean_Score,
     type = "b", pch = 19, cex = 1.5, col = "blue",
     xlab = "Disease Severity (Increasing →)", 
     ylab = "Mean PC1 Score",
     main = "Disease Severity Trend",
     xaxt = "n", ylim = c(min(severity_means$Mean_Score - severity_means$SD_Score),
                          max(severity_means$Mean_Score + severity_means$SD_Score)))
axis(1, at = 1:length(severity_means$Disease), labels = severity_means$Disease, las = 2)
arrows(1:length(severity_means$Mean_Score), 
       severity_means$Mean_Score - severity_means$SD_Score,
       1:length(severity_means$Mean_Score), 
       severity_means$Mean_Score + severity_means$SD_Score,
       length = 0.05, angle = 90, code = 3, col = "gray")
abline(h = 0, lty = 2, col = "red")
dev.off()

# ==================== 6. 梯度相关性检验 ====================

# Spearman相关性检验
severity_numeric <- as.numeric(score_df$Severity_Level)
cor_test <- cor.test(severity_numeric, score_df$PC1_Score, method = "spearman")

cat("\n严重程度梯度检验结果:\n")
cat(sprintf("  Spearman相关系数: %.3f\n", cor_test$estimate))
cat(sprintf("  P值: %.2e\n", cor_test$p.value))
cat(sprintf("  结论: %s\n", ifelse(cor_test$p.value < 0.05, 
    "存在显著的疾病严重程度梯度", "未发现显著梯度")))

# ==================== 7. 各菌群与严重程度的相关性 ====================

cat("\n=== 各菌群与疾病严重程度的相关性 ===\n")

severity_correlations <- data.frame()

for(cluster in 1:10) {
  cor_test <- cor.test(severity_numeric, cluster_abundance_rel[, cluster], 
                       method = "spearman")
  
  severity_correlations <- rbind(severity_correlations, data.frame(
    Cluster = cluster,
    Spearman_rho = cor_test$estimate,
    P_value = cor_test$p.value
  ))
}

severity_correlations$FDR <- p.adjust(severity_correlations$P_value, method = "fdr")
severity_correlations$Significant <- severity_correlations$FDR < 0.05
severity_correlations <- severity_correlations[order(-abs(severity_correlations$Spearman_rho)), ]

print(severity_correlations)
write.csv(severity_correlations, "cluster/biomarker_analysis/data/severity_correlations.csv", row.names = FALSE)

# 可视化
pdf("cluster/biomarker_analysis/figures/cluster_severity_correlation.pdf", width = 8, height = 5)
barplot(severity_correlations$Spearman_rho,
        names.arg = paste("Cluster", severity_correlations$Cluster),
        main = "Correlation with Disease Severity",
        ylab = "Spearman's rho",
        col = ifelse(severity_correlations$Spearman_rho > 0, "red", "blue"),
        las = 2)
abline(h = 0, lty = 2)
dev.off()

# ==================== 8. 沿着梯度方向的菌群变化 ====================

cat("\n=== 沿着严重程度梯度的菌群变化 ===\n")

# 计算每个疾病类型的平均菌群组成
severity_composition <- aggregate(cluster_abundance_rel, 
                                  by = list(Severity = score_df$Severity_Level), 
                                  FUN = mean)
rownames(severity_composition) <- severity_composition$Severity
severity_composition <- severity_composition[, -1]

# 热图展示
pdf("cluster/biomarker_analysis/figures/severity_composition_heatmap.pdf", 
    width = 10, height = 6)
pheatmap(t(severity_composition),
         main = "Microbial Composition Along Disease Severity Gradient",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize_row = 8)
dev.off()

# ==================== 9. 生成报告 ====================
sink("cluster/biomarker_analysis/biomarker_report.txt")

cat("========================================\n")
cat("疾病特异性生物标志物与严重程度梯度分析报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 疾病特异性生物标志物\n")
cat(sprintf("   - 分析的疾病类型数: %d\n", length(biomarker_results)))
cat(sprintf("   - 发现的显著标志物数: %d\n", sum(all_biomarkers$Significant)))
cat("\n   各疾病的显著标志物:\n")
for(disease in names(biomarker_results)) {
  sig_count <- sum(biomarker_results[[disease]]$Significant)
  cat(sprintf("     %s: %d个\n", disease, sig_count))
}
cat("\n")

cat("2. 疾病严重程度梯度\n")
cat(sprintf("   - Spearman相关系数: %.3f\n", cor_test$estimate))
cat(sprintf("   - P值: %.2e\n", cor_test$p.value))
cat(sprintf("   - 结论: %s\n", ifelse(cor_test$p.value < 0.05, 
    "存在显著的疾病严重程度梯度", "未发现显著梯度")))
cat("\n")

cat("3. 与严重程度显著相关的菌群\n")
sig_severity <- severity_correlations[severity_correlations$Significant, ]
if(nrow(sig_severity) > 0) {
  for(i in 1:nrow(sig_severity)) {
    direction <- ifelse(sig_severity$Spearman_rho[i] > 0, "正相关（随严重程度增加）", 
                        "负相关（随严重程度减少）")
    cat(sprintf("   - Cluster %d: %s (rho = %.3f, p = %.2e)\n",
                sig_severity$Cluster[i], direction,
                sig_severity$Spearman_rho[i], sig_severity$P_value[i]))
  }
} else {
  cat("   未发现显著相关的菌群\n")
}

sink()

cat("\n报告已保存到: cluster/biomarker_analysis/biomarker_report.txt\n")
cat("\n=== 分析完成 ===\n")

# ==================== 10. 打印关键结果摘要 ====================
cat("\n\n========================================\n")
cat("关键结果摘要\n")
cat("========================================\n")

cat("\n【疾病特异性生物标志物】\n")
for(disease in names(biomarker_results)) {
  sig_markers <- biomarker_results[[disease]][biomarker_results[[disease]]$Significant, ]
  if(nrow(sig_markers) > 0) {
    cat(sprintf("\n%s:\n", disease))
    for(i in 1:nrow(sig_markers)) {
      cat(sprintf("  Cluster %d: %s (%.2f倍变化)\n", 
                  sig_markers$Cluster[i],
                  ifelse(sig_markers$Log2FC[i] > 0, "↑升高", "↓降低"),
                  abs(sig_markers$FC[i])))
    }
  }
}

cat("\n【严重程度梯度】\n")
cat(sprintf("相关系数: %.3f (p = %.2e)\n", cor_test$estimate, cor_test$p.value))
cat("\n与严重程度最相关的菌群:\n")
for(i in 1:min(3, nrow(severity_correlations))) {
  cat(sprintf("  Cluster %d: rho = %.3f (p = %.2e)\n",
              severity_correlations$Cluster[i],
              severity_correlations$Spearman_rho[i],
              severity_correlations$P_value[i]))
}