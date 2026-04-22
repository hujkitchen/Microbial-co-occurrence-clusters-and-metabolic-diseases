# ==================== 完整的多疾病类型综合分析 ====================
cat("\n=== 完整的多疾病类型综合分析 ===\n")

# 创建分析目录
if(!dir.exists("cluster/multi_disease_analysis")) {
  dir.create("cluster/multi_disease_analysis")
}
if(!dir.exists("cluster/multi_disease_analysis/figures")) {
  dir.create("cluster/multi_disease_analysis/figures")
}
if(!dir.exists("cluster/multi_disease_analysis/data")) {
  dir.create("cluster/multi_disease_analysis/data")
}

library(mia)
library(pheatmap)

zero_samples <- which(colSums(assay(tse_2w, "counts")) == 0)
if(length(zero_samples) > 0) {
cat("移除以下零计数样本：", names(zero_samples), "\n")
tse_2w <- tse_2w[, -zero_samples]
}

# ==================== 0. 数据准备 ====================
cat("\n=== 0. 数据准备 ===\n")

# 检查并创建 cluster_abundance_rel
if(!exists("cluster_abundance_rel")) {
  cat("创建菌群组成数据...\n")
  
  # 获取计数数据
  counts <- assay(tse_2w, "counts")
  
  # 检查是否有物种聚类结果
  if(!exists("species_clusters")) {
    cat("进行物种聚类...\n")
    # 计算物种相关性（使用更快的近似方法）
    species_cor <- cor(t(counts), method = "spearman")
    dist_species <- as.dist(1 - abs(species_cor))
    hc_species <- hclust(dist_species, method = "ward.D2")
    
    # 使用10个聚类
    species_clusters <- cutree(hc_species, k = 10)
    cat("物种聚类完成\n")
  }
  
  # 计算菌群组成
  cluster_abundance <- matrix(0, ncol(tse_2w), 10)
  colnames(cluster_abundance) <- paste0("Cluster", 1:10)
  rownames(cluster_abundance) <- colnames(tse_2w)
  
  for(i in 1:10) {
    cluster_species <- names(species_clusters[species_clusters == i])
    cluster_abundance[, i] <- colSums(counts[cluster_species, ])
  }
  
  cluster_abundance_rel <- cluster_abundance / rowSums(cluster_abundance)
  cat("菌群组成数据创建完成，维度:", dim(cluster_abundance_rel), "\n")
}

# ==================== 1. 提取和处理condition信息 ====================
cat("\n=== 1. 处理疾病分类信息 ===\n")

# 获取condition信息
condition <- colData(tse_2w)$condition

# 查看唯一值
unique_conditions <- unique(condition)
cat("检测到的唯一条件数量:", length(unique_conditions), "\n")
print(table(condition))

# 创建简化的疾病分类
simplified_condition <- condition

# 定义简化规则（按顺序应用）
simplify_rules <- list(
  # 健康
  list(pattern = "healthy", replacement = "Healthy"),
  # 糖尿病相关
  list(pattern = "diabesity", replacement = "T2D"),
  list(pattern = "diabetes", replacement = "T2D"),
  list(pattern = "T2D", replacement = "T2D"),
  list(pattern = "gestational_diabetes", replacement = "Gestational_Diabetes"),
  # 代谢综合征
  list(pattern = "^MS$", replacement = "Metabolic_Syndrome"),
  list(pattern = "MS_obesity", replacement = "Metabolic_Syndrome_Obesity"),
  # NAFLD
  list(pattern = "NAFLD", replacement = "NAFLD"),
  # 肥胖相关
  list(pattern = "obesity$", replacement = "Obesity"),
  list(pattern = "overweight_obesity", replacement = "Overweight_Obesity"),
  list(pattern = "overweight", replacement = "Overweight"),
  list(pattern = "pre-diabetes_overweight", replacement = "Pre_Diabetes_Overweight"),
  # 其他
  list(pattern = "bariatric_surgery", replacement = "Post_Surgery"),
  list(pattern = "asthma_obesity", replacement = "Asthma_Obesity"),
  list(pattern = "opioid_use", replacement = "Opioid_Use"),
  list(pattern = "T2D_hypertension_obesity", replacement = "T2D_Metabolic")
)

# 应用简化规则
for(rule in simplify_rules) {
  matched <- grepl(rule$pattern, simplified_condition, ignore.case = TRUE)
  simplified_condition[matched] <- rule$replacement
}

# 更新colData
colData(tse_2w)$condition_simplified <- simplified_condition

# 查看简化后的分布
cat("\n简化后的疾病分类:\n")
print(table(colData(tse_2w)$condition_simplified))

# ==================== 2. 不同疾病类型的菌群组成 ====================
cat("\n=== 2. 不同疾病类型的菌群组成 ===\n")

# 计算各疾病类型的平均菌群组成
disease_composition <- data.frame(
  Condition = colData(tse_2w)$condition_simplified,
  cluster_abundance_rel
)

# 计算均值
disease_means <- aggregate(. ~ Condition, data = disease_composition, FUN = mean)

# 保存
write.csv(disease_means, "cluster/multi_disease_analysis/data/disease_composition_means.csv", row.names = FALSE)

# 创建热图（只显示样本数>=5的疾病类型）
sample_counts <- table(disease_composition$Condition)
diseases_to_show <- names(sample_counts[sample_counts >= 5])

if(length(diseases_to_show) > 1) {
  disease_means_filtered <- disease_means[disease_means$Condition %in% diseases_to_show, ]
  disease_means_matrix <- as.matrix(disease_means_filtered[, -1])
  rownames(disease_means_matrix) <- disease_means_filtered$Condition
  
  pdf("cluster/multi_disease_analysis/figures/disease_composition_heatmap.pdf", 
      width = 12, height = 8)
  pheatmap(disease_means_matrix,
           main = "Microbial Cluster Composition by Disease Type",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.3f",
           color = colorRampPalette(c("navy", "white", "red"))(50),
           fontsize_row = 10,
           fontsize_col = 8)
  dev.off()
  cat("热图已保存\n")
}

# ==================== 3. 健康组vs各疾病组 ====================
cat("\n=== 3. 健康组 vs 各疾病组 ===\n")

# 提取健康组
healthy_idx <- colData(tse_2w)$condition_simplified == "Healthy"
healthy_abundance <- cluster_abundance_rel[healthy_idx, ]

# 分析各疾病组
disease_types <- unique(colData(tse_2w)$condition_simplified)
disease_types <- disease_types[disease_types != "Healthy"]

results_list <- list()

for(disease in disease_types) {
  disease_idx <- colData(tse_2w)$condition_simplified == disease
  if(sum(disease_idx) < 3) next
  
  disease_abundance <- cluster_abundance_rel[disease_idx, ]
  
  for(cluster in 1:10) {
    test <- wilcox.test(healthy_abundance[, cluster], disease_abundance[, cluster])
    log2fc <- log2((mean(disease_abundance[, cluster]) + 0.0001) / 
                   (mean(healthy_abundance[, cluster]) + 0.0001))
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Disease = disease,
      Cluster = cluster,
      Log2FC = log2fc,
      P_value = test$p.value,
      Healthy_n = nrow(healthy_abundance),
      Disease_n = sum(disease_idx)
    )
  }
}

if(length(results_list) > 0) {
  all_results <- do.call(rbind, results_list)
  all_results$FDR <- p.adjust(all_results$P_value, method = "fdr")
  all_results$Significant <- all_results$FDR < 0.05
  
  write.csv(all_results, "cluster/multi_disease_analysis/data/disease_vs_healthy.csv", row.names = FALSE)
  
  # 创建热图
  library(reshape2)
  fc_matrix <- acast(all_results, Disease ~ Cluster, value.var = "Log2FC")
  
  pdf("cluster/multi_disease_analysis/figures/disease_vs_healthy_heatmap.pdf", 
      width = 10, height = 6)
  pheatmap(fc_matrix,
           main = "Log2 Fold Change (Disease vs Healthy)",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           fontsize_row = 10)
  dev.off()
  
  cat("差异分析完成，显著结果数:", sum(all_results$Significant), "\n")
}

# ==================== 4. 生成报告 ====================
sink("cluster/multi_disease_analysis/report.txt")

cat("========================================\n")
cat("多疾病类型综合分析报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("1. 数据概况\n")
cat(sprintf("   - 总样本数: %d\n", ncol(tse_2w)))
cat(sprintf("   - 健康组: %d\n", sum(colData(tse_2w)$condition_simplified == "Healthy")))
cat("   - 疾病类型分布:\n")
for(d in names(table(colData(tse_2w)$condition_simplified))) {
  cat(sprintf("     %s: %d\n", d, sum(colData(tse_2w)$condition_simplified == d)))
}

if(exists("all_results")) {
  cat("\n2. 差异分析摘要\n")
  cat(sprintf("   - 分析的疾病类型数: %d\n", length(unique(all_results$Disease))))
  cat(sprintf("   - 总比较次数: %d\n", nrow(all_results)))
  cat(sprintf("   - 显著差异数 (FDR<0.05): %d\n", sum(all_results$Significant)))
}

sink()

cat("\n报告已保存到: cluster/multi_disease_analysis/report.txt\n")
cat("\n=== 分析完成 ===\n")