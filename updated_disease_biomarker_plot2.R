# ==================== 首先检查数据结构 ====================
cat("\n=== 检查数据结构 ===\n")

if(exists("all_biomarkers")) {
  cat("all_biomarkers 的列名:\n")
  print(colnames(all_biomarkers))
  
  # 检查是否有 P_value 列（可能是 P_value 或 p_value）
  if(!"P_value" %in% colnames(all_biomarkers) && "p_value" %in% colnames(all_biomarkers)) {
    all_biomarkers$P_value <- all_biomarkers$p_value
    cat("已将 p_value 重命名为 P_value\n")
  }
  
  # 检查是否有 Log2FC 列
  if(!"Log2FC" %in% colnames(all_biomarkers) && "log2FC" %in% colnames(all_biomarkers)) {
    all_biomarkers$Log2FC <- all_biomarkers$log2FC
  }
  
  # 检查是否有 Significant 列
  if(!"Significant" %in% colnames(all_biomarkers) && "significant" %in% colnames(all_biomarkers)) {
    all_biomarkers$Significant <- all_biomarkers$significant
  }
  
  cat("\n修复后的列名:\n")
  print(colnames(all_biomarkers))
  
} else {
  cat("错误：all_biomarkers 不存在，正在重新创建...\n")
  
  # 重新创建 biomarker 结果
  disease_types <- unique(colData(tse_2w)$condition_simplified)
  disease_types <- disease_types[disease_types != "Healthy" & !is.na(disease_types)]
  disease_types <- disease_types[sapply(disease_types, function(d) sum(colData(tse_2w)$condition_simplified == d) >= 10)]
  
  all_biomarkers <- data.frame()
  
  for(disease in disease_types) {
    idx_this <- colData(tse_2w)$condition_simplified == disease
    idx_other <- colData(tse_2w)$condition_simplified %in% disease_types & 
                colData(tse_2w)$condition_simplified != disease
    
    if(sum(idx_this) < 5 || sum(idx_other) < 5) next
    
    for(cluster in 1:10) {
      vals_this <- cluster_abundance_rel[idx_this, cluster]
      vals_other <- cluster_abundance_rel[idx_other, cluster]
      
      test <- wilcox.test(vals_this, vals_other)
      log2fc <- log2((mean(vals_this) + 0.0001) / (mean(vals_other) + 0.0001))
      
      all_biomarkers <- rbind(all_biomarkers, data.frame(
        Disease = disease,
        Cluster = cluster,
        Log2FC = log2fc,
        P_value = test$p.value,
        N_This = sum(idx_this),
        N_Other = sum(idx_other),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  all_biomarkers$FDR <- p.adjust(all_biomarkers$P_value, method = "fdr")
  all_biomarkers$Significant <- all_biomarkers$FDR < 0.05
  
  cat("重新创建完成，共", nrow(all_biomarkers), "行\n")
}

# ==================== 修复后的可视化 ====================
cat("\n=== 生成修复后的可视化 ===\n")

# 确保目录存在
if(!dir.exists("cluster/biomarker_analysis/figures")) {
  dir.create("cluster/biomarker_analysis/figures", recursive = TRUE)
}

library(ggplot2)
library(reshape2)

# ==================== 图1：热图（简化版，不使用P_value）====================
cat("生成热图...\n")

# 准备数据
heatmap_data <- all_biomarkers[, c("Disease", "Cluster", "Log2FC", "Significant")]
heatmap_data$Cluster <- factor(heatmap_data$Cluster)
heatmap_data$Disease <- factor(heatmap_data$Disease)

# 使用PNG格式（更稳定）
png("cluster/biomarker_analysis/figures/disease_specific_biomarkers.png", 
    width = 800, height = 600)

p1 <- ggplot(heatmap_data, aes(x = Cluster, y = Disease, fill = Log2FC)) +
  geom_tile() +
  geom_text(aes(label = ifelse(Significant, "*", "")), size = 6, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Log2 FC") +
  labs(title = "Disease-Specific Microbial Biomarkers",
       x = "Microbial Cluster", y = "Disease") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14))

print(p1)
dev.off()

# 也保存为PDF
pdf("cluster/biomarker_analysis/figures/disease_specific_biomarkers.pdf", 
    width = 10, height = 8)
print(p1)
dev.off()

cat("热图已保存\n")

# ==================== 图2：条形图（每个疾病最重要的标志物）====================
cat("生成条形图...\n")

top_markers <- all_biomarkers[all_biomarkers$Significant, ]

if(nrow(top_markers) > 0) {
  # 为每个疾病选择最显著的标志物（按p值排序）
  top_per_disease <- top_markers %>%
    group_by(Disease) %>%
    arrange(P_value, desc(abs(Log2FC))) %>%
    slice(1)
  
  top_per_disease$Label <- paste(top_per_disease$Disease, 
                                 "Cluster", top_per_disease$Cluster)
  top_per_disease$Direction <- ifelse(top_per_disease$Log2FC > 0, 
                                      "Enriched in Disease", "Depleted in Disease")
  
  png("cluster/biomarker_analysis/figures/top_biomarkers.png", 
      width = 800, height = 600)
  
  p2 <- ggplot(top_per_disease, aes(x = reorder(Label, Log2FC), 
                                    y = Log2FC, fill = Direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Enriched in Disease" = "red", 
                                 "Depleted in Disease" = "blue")) +
    labs(title = "Top Disease-Specific Microbial Biomarkers",
         x = "Disease - Cluster", y = "Log2 Fold Change") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 9))
  
  print(p2)
  dev.off()
  
  pdf("cluster/biomarker_analysis/figures/top_biomarkers.pdf", 
      width = 10, height = 6)
  print(p2)
  dev.off()
  
  cat("条形图已保存\n")
} else {
  cat("没有发现显著的生物标志物，跳过条形图\n")
}

# ==================== 图3：点图（显示所有标志物）====================
cat("生成点图...\n")

# 计算 -log10(P_value)
all_biomarkers$neg_log10_p <- -log10(all_biomarkers$P_value)
all_biomarkers$Direction <- ifelse(all_biomarkers$Log2FC > 0, "Enriched", "Depleted")

png("cluster/biomarker_analysis/figures/biomarkers_dotplot.png", 
    width = 1000, height = 700)

p3 <- ggplot(all_biomarkers, aes(x = factor(Cluster), y = Log2FC)) +
  geom_point(aes(size = neg_log10_p, color = Direction), alpha = 0.7) +
  facet_wrap(~ Disease, ncol = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Enriched" = "red", "Depleted" = "blue")) +
  scale_size_continuous(name = "-log10(P-value)") +
  labs(title = "Disease-Specific Microbial Biomarkers",
       x = "Microbial Cluster", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

print(p3)
dev.off()

pdf("cluster/biomarker_analysis/figures/biomarkers_dotplot.pdf", 
    width = 12, height = 8)
print(p3)
dev.off()

cat("点图已保存\n")

# ==================== 图4：每个疾病的标志物热图（单独）====================
cat("生成各疾病单独的热图...\n")

for(disease in unique(all_biomarkers$Disease)) {
  disease_data <- all_biomarkers[all_biomarkers$Disease == disease, ]
  
  png(paste0("cluster/biomarker_analysis/figures/biomarkers_", 
             gsub(" ", "_", disease), ".png"), 
      width = 600, height = 400)
  
  p4 <- ggplot(disease_data, aes(x = factor(Cluster), y = "1", fill = Log2FC)) +
    geom_tile() +
    geom_text(aes(label = ifelse(Significant, 
                                 sprintf("%.2f", Log2FC), "")), 
              size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, name = "Log2 FC") +
    labs(title = paste(disease, "- Microbial Biomarkers"),
         x = "Cluster", y = "") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p4)
  dev.off()
}

cat("各疾病单独热图已保存\n")

# ==================== 生成文本摘要 ====================
cat("\n=== 生成文本摘要 ===\n")

sink("cluster/biomarker_analysis/biomarker_summary.txt")

cat("========================================\n")
cat("疾病特异性生物标志物摘要\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 总体统计\n")
cat(sprintf("   - 分析的疾病类型数: %d\n", length(unique(all_biomarkers$Disease))))
cat(sprintf("   - 总比较次数: %d\n", nrow(all_biomarkers)))
cat(sprintf("   - 显著标志物数 (FDR < 0.05): %d\n", sum(all_biomarkers$Significant)))
cat("\n")

cat("2. 各疾病显著标志物\n")
for(disease in unique(all_biomarkers$Disease)) {
  sig_markers <- all_biomarkers[all_biomarkers$Disease == disease & all_biomarkers$Significant, ]
  cat(sprintf("\n  %s:\n", disease))
  
  if(nrow(sig_markers) > 0) {
    for(i in 1:nrow(sig_markers)) {
      direction <- ifelse(sig_markers$Log2FC[i] > 0, "↑ 升高", "↓ 降低")
      cat(sprintf("    Cluster %d: %s (Log2FC = %.2f, p = %.2e)\n", 
                  sig_markers$Cluster[i], direction,
                  sig_markers$Log2FC[i], sig_markers$P_value[i]))
    }
  } else {
    cat("    未发现显著标志物\n")
  }
}

sink()

cat("文本摘要已保存: biomarker_summary.txt\n")

# ==================== 打印到控制台 ====================
cat("\n\n========================================\n")
cat("显著生物标志物列表:\n")
cat("========================================\n")

sig_all <- all_biomarkers[all_biomarkers$Significant, ]
if(nrow(sig_all) > 0) {
  for(i in 1:nrow(sig_all)) {
    cat(sprintf("%s - Cluster %d: Log2FC = %.2f (p = %.2e)\n",
                sig_all$Disease[i], sig_all$Cluster[i],
                sig_all$Log2FC[i], sig_all$P_value[i]))
  }
} else {
  cat("未发现显著生物标志物\n")
}

cat("\n=== 修复完成 ===\n")