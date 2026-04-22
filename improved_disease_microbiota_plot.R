# ==================== 改进的疾病特异性菌群可视化 ====================
cat("\n=== 改进的疾病特异性菌群可视化 ===\n")


# 确保目录存在
if(!dir.exists("cluster/biomarker_analysis/figures")) {
  dir.create("cluster/biomarker_analysis/figures", recursive = TRUE)
}

library(ggplot2)
library(reshape2)
library(gridExtra)

# ==================== 方法1：使用 fmsb 包的雷达图（更稳定）====================
cat("生成雷达图（使用fmsb包）...\n")

# 检查是否安装了fmsb包
if(!require(fmsb)) {
  install.packages("fmsb")
  library(fmsb)
}

# 为每个疾病创建雷达图数据
for(disease in unique(all_biomarkers$Disease)) {
  disease_data <- all_biomarkers[all_biomarkers$Disease == disease, ]
  
  # 创建雷达图数据框
  radar_df <- data.frame(
    C1 = disease_data$Log2FC[disease_data$Cluster == 1],
    C2 = disease_data$Log2FC[disease_data$Cluster == 2],
    C3 = disease_data$Log2FC[disease_data$Cluster == 3],
    C4 = disease_data$Log2FC[disease_data$Cluster == 4],
    C5 = disease_data$Log2FC[disease_data$Cluster == 5],
    C6 = disease_data$Log2FC[disease_data$Cluster == 6],
    C7 = disease_data$Log2FC[disease_data$Cluster == 7],
    C8 = disease_data$Log2FC[disease_data$Cluster == 8],
    C9 = disease_data$Log2FC[disease_data$Cluster == 9],
    C10 = disease_data$Log2FC[disease_data$Cluster == 10]
  )
  
  # 添加最大最小值行
  max_val <- max(abs(radar_df))
  radar_df <- rbind(rep(max_val, 10), rep(-max_val, 10), radar_df)
  rownames(radar_df) <- c("Max", "Min", disease)
  
  # 绘制雷达图
  png(paste0("cluster/biomarker_analysis/figures/radar_", 
             gsub(" ", "_", disease), ".png"), 
      width = 800, height = 800)
  
  radarchart(radar_df,
             axistype = 1,
             pcol = rgb(0.2, 0.5, 0.8, 0.9),
             pfcol = rgb(0.2, 0.5, 0.8, 0.3),
             plwd = 2,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "grey",
             caxislabels = seq(-max_val, max_val, length.out = 5),
             cglwd = 0.8,
             vlcex = 1.2,
             title = paste(disease, "- Microbial Profile"))
  
  dev.off()
}

cat("雷达图已保存\n")

# ==================== 方法2：水平条形图（所有疾病对比）- 最推荐 ====================
cat("生成水平条形图...\n")

# 选择显著标志物
sig_markers <- all_biomarkers[all_biomarkers$Significant, ]

if(nrow(sig_markers) > 0) {
  # 创建标签和排序
  sig_markers$Label <- paste(sig_markers$Disease, "- Cluster", sig_markers$Cluster)
  sig_markers$Direction <- ifelse(sig_markers$Log2FC > 0, "Enriched in Disease", "Depleted in Disease")
  sig_markers$neg_log10_p <- -log10(sig_markers$P_value)
  
  # 按Log2FC排序
  sig_markers <- sig_markers[order(sig_markers$Log2FC), ]
  sig_markers$Label <- factor(sig_markers$Label, levels = sig_markers$Label)
  
  # 水平条形图
  p1 <- ggplot(sig_markers, aes(x = Label, y = Log2FC, fill = Direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Enriched in Disease" = "#E41A1C", 
                                 "Depleted in Disease" = "#377EB8")) +
    labs(title = "Significant Disease-Specific Microbial Biomarkers",
         x = "", y = "Log2 Fold Change") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.y = element_text(size = 10),
          legend.position = "bottom")
  
  ggsave("cluster/biomarker_analysis/figures/horizontal_bars.png", 
         p1, width = 10, height = 8, dpi = 300)
  ggsave("cluster/biomarker_analysis/figures/horizontal_bars.pdf", 
         p1, width = 10, height = 8)
  
  cat("水平条形图已保存\n")
}

# ==================== 方法3：点阵图（Lollipop chart）====================
cat("生成点阵图...\n")

if(nrow(sig_markers) > 0) {
  # 按疾病分组
  p2 <- ggplot(sig_markers, aes(x = Log2FC, y = reorder(Label, Log2FC))) +
    geom_segment(aes(x = 0, xend = Log2FC, y = reorder(Label, Log2FC), 
                     yend = reorder(Label, Log2FC)), color = "gray") +
    geom_point(aes(color = Direction, size = neg_log10_p)) +
    scale_color_manual(values = c("Enriched in Disease" = "#E41A1C", 
                                  "Depleted in Disease" = "#377EB8")) +
    scale_size_continuous(name = "-log10(P-value)", range = c(2, 8)) +
    labs(title = "Lollipop Chart of Disease-Specific Biomarkers",
         x = "Log2 Fold Change", y = "") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  ggsave("cluster/biomarker_analysis/figures/lollipop_chart.png", 
         p2, width = 10, height = 8, dpi = 300)
  ggsave("cluster/biomarker_analysis/figures/lollipop_chart.pdf", 
         p2, width = 10, height = 8)
  
  cat("点阵图已保存\n")
}

# ==================== 方法4：热图（简化版）====================
cat("生成热图...\n")

# 创建矩阵
fc_matrix <- acast(all_biomarkers, Disease ~ Cluster, value.var = "Log2FC")
sig_matrix <- acast(all_biomarkers, Disease ~ Cluster, value.var = "Significant")

# 转换为长格式用于ggplot
heatmap_long <- melt(fc_matrix)
colnames(heatmap_long) <- c("Disease", "Cluster", "Log2FC")
heatmap_long$Significant <- melt(sig_matrix)$value
heatmap_long$Cluster <- factor(heatmap_long$Cluster)

# 热图
p3 <- ggplot(heatmap_long, aes(x = Cluster, y = Disease, fill = Log2FC)) +
  geom_tile() +
  geom_text(aes(label = ifelse(Significant, sprintf("%.2f", Log2FC), "")), 
            size = 3, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, name = "Log2 FC") +
  labs(title = "Disease-Specific Microbial Biomarkers",
       x = "Microbial Cluster", y = "Disease") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave("cluster/biomarker_analysis/figures/heatmap.png", 
       p3, width = 10, height = 8, dpi = 300)
ggsave("cluster/biomarker_analysis/figures/heatmap.pdf", 
       p3, width = 10, height = 8)

cat("热图已保存\n")

# ==================== 方法5：每个疾病的簇状条形图 ====================
cat("生成各疾病簇状条形图...\n")

for(disease in unique(all_biomarkers$Disease)) {
  disease_data <- all_biomarkers[all_biomarkers$Disease == disease, ]
  
  # 标记显著性
  disease_data$Signif_Label <- ifelse(disease_data$Significant, 
                                      ifelse(disease_data$Log2FC > 0, "↑", "↓"), 
                                      "")
  disease_data$Color <- ifelse(disease_data$Log2FC > 0, "Enriched", "Depleted")
  
  p4 <- ggplot(disease_data, aes(x = factor(Cluster), y = Log2FC, fill = Color)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Signif_Label), vjust = ifelse(disease_data$Log2FC > 0, -0.5, 1.5),
              size = 5, fontface = "bold") +
    scale_fill_manual(values = c("Enriched" = "#E41A1C", "Depleted" = "#377EB8")) +
    labs(title = paste(disease, "- Microbial Biomarkers"),
         x = "Cluster", y = "Log2 Fold Change") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  ggsave(paste0("cluster/biomarker_analysis/figures/barplot_", 
                gsub(" ", "_", disease), ".png"), 
         p4, width = 8, height = 5, dpi = 300)
}

cat("各疾病条形图已保存\n")

# ==================== 方法6：汇总点阵图（所有疾病）====================
cat("生成汇总点阵图...\n")

# 计算每个疾病中显著变化的菌群数量
disease_summary <- aggregate(Significant ~ Disease, data = all_biomarkers, FUN = sum)
colnames(disease_summary)[2] <- "N_Significant"

# 计算每个疾病的平均Log2FC绝对值
disease_summary$Mean_abs_Log2FC <- aggregate(abs(Log2FC) ~ Disease, 
                                             data = all_biomarkers[all_biomarkers$Significant, ], 
                                             FUN = mean)$`abs(Log2FC)`

p5 <- ggplot(disease_summary, aes(x = N_Significant, y = reorder(Disease, N_Significant), 
                                   size = Mean_abs_Log2FC, color = N_Significant)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red", name = "# Significant\nBiomarkers") +
  scale_size_continuous(name = "Mean |Log2 FC|", range = c(3, 10)) +
  labs(title = "Summary of Disease-Specific Biomarkers",
       x = "Number of Significant Biomarkers", y = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("cluster/biomarker_analysis/figures/summary_dotplot.png", 
       p5, width = 8, height = 6, dpi = 300)

cat("汇总点阵图已保存\n")

# ==================== 生成索引文件 ====================
sink("cluster/biomarker_analysis/figures_index.html")

cat("<html>\n<head>\n<title>Disease-Specific Biomarkers Visualization</title>\n")
cat("<style>\n")
cat("body { font-family: Arial, sans-serif; margin: 20px; }\n")
cat("h1 { color: #333; }\n")
cat("h2 { color: #666; margin-top: 30px; }\n")
cat(".gallery { display: flex; flex-wrap: wrap; gap: 20px; }\n")
cat(".figure { border: 1px solid #ddd; padding: 10px; width: 300px; }\n")
cat(".figure img { width: 100%; height: auto; }\n")
cat(".figure p { margin: 5px 0; font-size: 12px; text-align: center; }\n")
cat("</style>\n</head>\n<body>\n")

cat("<h1>Disease-Specific Microbial Biomarkers</h1>\n")

cat("<h2>Main Figures</h2>\n")
cat("<div class='gallery'>\n")
cat("<div class='figure'><img src='horizontal_bars.png'><p>Horizontal Bar Chart</p></div>\n")
cat("<div class='figure'><img src='lollipop_chart.png'><p>Lollipop Chart</p></div>\n")
cat("<div class='figure'><img src='heatmap.png'><p>Heatmap</p></div>\n")
cat("<div class='figure'><img src='summary_dotplot.png'><p>Summary Dotplot</p></div>\n")
cat("</div>\n")

cat("<h2>Individual Disease Barplots</h2>\n")
cat("<div class='gallery'>\n")
for(disease in unique(all_biomarkers$Disease)) {
  filename <- paste0("barplot_", gsub(" ", "_", disease), ".png")
  cat(paste0("<div class='figure'><img src='", filename, "'><p>", disease, "</p></div>\n"))
}
cat("</div>\n")

cat("<h2>Radar Charts</h2>\n")
cat("<div class='gallery'>\n")
for(disease in unique(all_biomarkers$Disease)) {
  filename <- paste0("radar_", gsub(" ", "_", disease), ".png")
  cat(paste0("<div class='figure'><img src='", filename, "'><p>", disease, "</p></div>\n"))
}
cat("</div>\n")

cat("</body>\n</html>\n")

sink()

cat("\nHTML索引文件已保存: figures_index.html\n")
cat("\n=== 所有图形生成完成 ===\n")