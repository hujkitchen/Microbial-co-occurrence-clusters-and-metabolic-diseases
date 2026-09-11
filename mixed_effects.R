# ==================== 混合效应模型分析====================
cat("\n=== 混合效应模型分析===\n")

library(lme4)

healthy_idx <- simplified_condition == "Healthy"
disease_idx <- !healthy_idx

study_adjusted_results <- data.frame()

cat("对每个Cluster拟合混合效应模型...\n")

for(cl in 1:k) {
  cat(sprintf("Cluster %d...\n", cl))
  
  model_data <- data.frame(
    Abundance = cluster_abundance_rel[, cl],
    Group = factor(ifelse(disease_idx, "Disease", "Healthy"), 
                   levels = c("Healthy", "Disease")),
    Study = factor(study_info)
  )
  
  tryCatch({
    model <- lmer(Abundance ~ Group + (1 | Study), data = model_data)
    summary_model <- summary(model)
    coefs <- coef(summary_model)
    
    # 查找GroupDisease行
    disease_row <- which(rownames(coefs) == "GroupDisease")
    
    if(length(disease_row) > 0) {
      # 提取系数和标准误
      coefficient <- coefs[disease_row, "Estimate"]
      se <- coefs[disease_row, "Std. Error"]
      t_value <- coefs[disease_row, "t value"]
      
      # 计算p值（使用正态近似，因为样本量大）
      p_value <- 2 * pnorm(-abs(t_value))
      
      study_adjusted_results <- rbind(study_adjusted_results, data.frame(
        Cluster = cl,
        Coefficient = coefficient,
        SE = se,
        t_value = t_value,
        P_value = p_value
      ))
      cat(sprintf("  系数=%.6f, SE=%.6f, t=%.3f, p=%.2e\n", 
                  coefficient, se, t_value, p_value))
    }
  }, error = function(e) {
    cat(sprintf("  失败: %s\n", e$message))
  })
}

# 结果
cat(sprintf("\n成功分析的Cluster数: %d\n", nrow(study_adjusted_results)))

if(nrow(study_adjusted_results) > 0) {
  study_adjusted_results$FDR <- p.adjust(study_adjusted_results$P_value, method = "fdr")
  study_adjusted_results$Significant <- study_adjusted_results$FDR < 0.05
  
  cat("\n========================================\n")
  cat("混合效应模型结果（Study作为随机效应）\n")
  cat("========================================\n")
  print(round(study_adjusted_results, 6))
  
  cat(sprintf("\n显著Cluster数 (LMM, FDR<0.05): %d / %d\n", 
              sum(study_adjusted_results$Significant), nrow(study_adjusted_results)))
  
  # 与Wilcoxon比较
  if(exists("diff_results") && nrow(diff_results) > 0) {
    cat(sprintf("显著Cluster数 (Wilcoxon, FDR<0.05): %d / %d\n", 
                sum(diff_results$Significant), nrow(diff_results)))
  }
  
  # 保存
  write.csv(study_adjusted_results, 
            "revision_analysis/data/TableS6_mixed_effects_results.csv", 
            row.names = FALSE)
  
  cat("\n结果已保存\n")
}