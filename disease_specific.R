# ==================== 按疾病类型分层分析（不使用合并的Disease组）====================
cat("\n========================================\n")
cat("按疾病类型分层分析\n")
cat("========================================\n")

library(lme4)

# 获取所有疾病类型（排除Healthy）
disease_types <- unique(simplified_condition)
disease_types <- disease_types[disease_types != "Healthy"]
cat("疾病类型:", paste(disease_types, collapse = ", "), "\n")

# 简化协变量
geo_info <- colData(tse_2w)$geo_loc_name
geo_counts <- sort(table(geo_info), decreasing = TRUE)
top10_geo <- names(geo_counts)[1:min(10, length(geo_counts))]
geo_simplified <- factor(ifelse(geo_info %in% top10_geo, geo_info, "Other"))

platform_info <- colData(tse_2w)$instrument
plat_counts <- sort(table(platform_info), decreasing = TRUE)
top5_plat <- names(plat_counts)[1:min(5, length(plat_counts))]
plat_simplified <- factor(ifelse(platform_info %in% top5_plat, platform_info, "Other"))

# ==================== 模型: 每种疾病 vs Healthy + 协变量 + Study随机效应 ====================
cat("\n--- 每种疾病 vs Healthy（调整Study + Geography + Platform）---\n")

all_disease_results <- data.frame()

for(disease in disease_types) {
  cat(sprintf("\n分析: %s\n", disease))
  
  # 选择该疾病和Healthy样本
  disease_idx <- simplified_condition == disease
  subset_idx <- disease_idx | healthy_idx
  
  if(sum(disease_idx) < 30) {
    cat(sprintf("  样本数不足 (n=%d)，跳过\n", sum(disease_idx)))
    next
  }
  
  cat(sprintf("  Healthy: %d, %s: %d\n", sum(healthy_idx[subset_idx]), disease, sum(disease_idx)))
  
  # 子集数据
  X_sub <- cluster_abundance_rel[subset_idx, ]
  y_sub <- factor(ifelse(simplified_condition[subset_idx] == disease, disease, "Healthy"),
                  levels = c("Healthy", disease))
  study_sub <- factor(study_info[subset_idx])
  geo_sub <- geo_simplified[subset_idx]
  plat_sub <- plat_simplified[subset_idx]
  
  # 对每个Cluster拟合模型
  for(cl in 1:k) {
    model_data <- data.frame(
      Abundance = X_sub[, cl],
      Group = y_sub,
      Study = study_sub,
      Geography = geo_sub,
      Platform = plat_sub
    )
    
    tryCatch({
      model <- lmer(Abundance ~ Group + Geography + Platform + (1|Study), 
                    data = model_data)
      coefs <- coef(summary(model))
      
      # 查找该疾病的系数行
      disease_coef_name <- paste0("Group", disease)
      disease_row <- which(rownames(coefs) == disease_coef_name)
      
      if(length(disease_row) > 0) {
        all_disease_results <- rbind(all_disease_results, data.frame(
          Disease = disease,
          Cluster = cl,
          Coefficient = coefs[disease_row, "Estimate"],
          SE = coefs[disease_row, "Std. Error"],
          t_value = coefs[disease_row, "t value"],
          P_value = 2 * pnorm(-abs(coefs[disease_row, "t value"])),
          N_Disease = sum(disease_idx),
          N_Total = sum(subset_idx)
        ))
      }
    }, error = function(e) {
      # 静默处理
    })
  }
}

# FDR校正（按疾病类型分别校正）
if(nrow(all_disease_results) > 0) {
  all_disease_results$FDR <- NA
  
  for(disease in unique(all_disease_results$Disease)) {
    disease_rows <- which(all_disease_results$Disease == disease)
    all_disease_results$FDR[disease_rows] <- p.adjust(
      all_disease_results$P_value[disease_rows], method = "fdr"
    )
  }
  
  all_disease_results$Significant <- all_disease_results$FDR < 0.05
  
  cat("\n========================================\n")
  cat("按疾病类型分层分析结果\n")
  cat("========================================\n")
  
  # 显示显著结果
  sig_results <- all_disease_results[all_disease_results$Significant, ]
  if(nrow(sig_results) > 0) {
    cat("\n显著结果 (FDR < 0.05):\n")
    print(sig_results[, c("Disease","Cluster","Coefficient","P_value","FDR","N_Disease")])
  } else {
    cat("\n未发现任何显著的疾病特异性Cluster差异\n")
  }
  
  # 汇总
  cat("\n========================================\n")
  cat("显著Cluster数汇总（按疾病）\n")
  cat("========================================\n")
  for(disease in unique(all_disease_results$Disease)) {
    n_sig <- sum(all_disease_results$Significant[all_disease_results$Disease == disease])
    n_total <- sum(all_disease_results$Disease == disease)
    cat(sprintf("  %s: %d / %d\n", disease, n_sig, n_total))
  }
  
  # 与未调整Wilcoxon比较
  cat("\n========================================\n")
  cat("对比：未调整Wilcoxon vs LMM调整\n")
  cat("========================================\n")
  
  # 计算未调整的疾病特异性Wilcoxon检验
  unadjusted_results <- data.frame()
  for(disease in unique(all_disease_results$Disease)) {
    d_idx <- simplified_condition == disease
    for(cl in 1:k) {
      vals_d <- cluster_abundance_rel[d_idx, cl]
      vals_h <- cluster_abundance_rel[healthy_idx, cl]
      test <- wilcox.test(vals_d, vals_h)
      unadjusted_results <- rbind(unadjusted_results, data.frame(
        Disease = disease, Cluster = cl, P_value = test$p.value
      ))
    }
  }
  unadjusted_results$FDR <- NA
  for(disease in unique(unadjusted_results$Disease)) {
    rows <- which(unadjusted_results$Disease == disease)
    unadjusted_results$FDR[rows] <- p.adjust(unadjusted_results$P_value[rows], method = "fdr")
  }
  
  # 合并比较
  compare_disease <- data.frame(
    Disease = all_disease_results$Disease,
    Cluster = all_disease_results$Cluster,
    Wilcoxon_P = unadjusted_results$P_value[match(
      paste(all_disease_results$Disease, all_disease_results$Cluster),
      paste(unadjusted_results$Disease, unadjusted_results$Cluster)
    )],
    Wilcoxon_FDR = unadjusted_results$FDR[match(
      paste(all_disease_results$Disease, all_disease_results$Cluster),
      paste(unadjusted_results$Disease, unadjusted_results$Cluster)
    )],
    LMM_P = all_disease_results$P_value,
    LMM_FDR = all_disease_results$FDR,
    N_Disease = all_disease_results$N_Disease
  )
  
  cat("\n前20行对比:\n")
  print(head(compare_disease, 20))
  
  # 汇总显著数
  cat("\n显著Cluster数汇总 (FDR < 0.05):\n")
  cat(sprintf("  未调整Wilcoxon: %d\n", sum(compare_disease$Wilcoxon_FDR < 0.05, na.rm = TRUE)))
  cat(sprintf("  LMM调整:        %d\n", sum(compare_disease$LMM_FDR < 0.05, na.rm = TRUE)))
  
  # 保存
  write.csv(all_disease_results, 
            "revision_analysis/data/TableS7_disease_specific_LMM.csv", 
            row.names = FALSE)
  write.csv(compare_disease, 
            "revision_analysis/data/TableS7_disease_comparison.csv", 
            row.names = FALSE)
  
  cat("\n结果已保存\n")
}