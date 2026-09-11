# ==================== 置换检验结果深入分析 ====================
cat("\n========================================\n")
cat("置换检验深入分析\n")
cat("========================================\n")


cat(sprintf("原始AUC: %.3f\n", auc_original))
cat(sprintf("置换AUC: %.3f +/- %.3f [%.3f, %.3f]\n", 
            mean(permutation_results$AUC),
            sd(permutation_results$AUC),
            min(permutation_results$AUC),
            max(permutation_results$AUC)))
cat(sprintf("差值: %.3f\n", auc_original - mean(permutation_results$AUC)))
cat(sprintf("经验p值: %.3f\n", p_perm))

# 可视化（避免特殊字符）
pdf("revision_analysis/figures/FigS_permutation_test_detailed.pdf", width = 10, height = 6)

par(mfrow = c(1, 2))

# Panel A: 直方图
hist(permutation_results$AUC, breaks = 20, 
     main = "Permutation Test (Study-Stratified Shuffle)",
     xlab = "AUC", col = "lightblue", border = "darkblue",
     xlim = c(0.7, 0.8))
abline(v = auc_original, col = "red", lwd = 3)
abline(v = mean(permutation_results$AUC), col = "blue", lwd = 2, lty = 2)

legend("topleft", 
       legend = c(
         paste("Observed =", round(auc_original, 3)),
         paste("Permuted mean =", round(mean(permutation_results$AUC), 3)),
         paste("Diff =", round(auc_original - mean(permutation_results$AUC), 3)),
         "p < 0.001 (effect negligible)"),
       col = c("red", "blue", "black", "white"), 
       lty = c(1, 2, 0, 0), lwd = c(3, 2, 0, 0), bty = "n", cex = 0.9)

# Panel B: 箱线图
boxplot(list(Observed = auc_original, 
             Permuted = permutation_results$AUC),
        main = "AUC: Observed vs Permuted",
        ylab = "AUC", col = c("red", "lightblue"),
        ylim = c(0.7, 0.8))

# 添加差值标注（使用ASCII字符）
diff_val <- round(auc_original - mean(permutation_results$AUC), 3)
text(1.5, auc_original + 0.005, 
     paste("Diff =", diff_val), 
     cex = 1.2, font = 2)

dev.off()

cat("\n图表已保存\n")