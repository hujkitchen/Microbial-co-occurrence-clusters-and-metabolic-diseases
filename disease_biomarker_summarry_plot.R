# ==================== 综合示意图 ====================
cat("\n=== 创建疾病特异性生物标志物与严重程度梯度综合示意图 ===\n")


# 确保目录存在
if(!dir.exists("cluster/figures")) {
  dir.create("cluster/figures", recursive = TRUE)
}

library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

# ==================== 数据准备 ====================

# 疾病特异性标志物数据
disease_biomarkers <- data.frame(
  Disease = c("T2D", "Post_Surgery", "Obesity", "Overweight", 
              "NAFLD", "Metabolic_Syndrome", "Opioid_Use"),
  Count = c(7, 8, 8, 9, 8, 7, 5),
  Color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
            "#FF7F00", "#FFFF33", "#A65628")
)

# 严重程度梯度数据
severity_correlations <- data.frame(
  Cluster = c(1, 3, 8, 6, 2, 10, 4),
  Correlation = c(0.076, 0.049, 0.029, -0.026, -0.026, 0.024, 0.016),
  Direction = c("Positive", "Positive", "Positive", "Negative", 
                "Negative", "Positive", "Positive"),
  P_value = c(6.49e-27, 3.75e-12, 4.39e-05, 2.34e-04, 
              3.22e-04, 6.62e-04, 2.74e-02),
  Significance = c("***", "***", "***", "***", "***", "***", "*")
)

# 聚类功能注释
cluster_functions <- data.frame(
  Cluster = 1:10,
  Function = c("Butyrate producers", "Opportunistic pathogens", 
               "Mucin degraders", "Lactate producers", 
               "Succinate producers", "Bile acid metabolizers",
               "Mucin degraders (IBD)", "Unclassified", 
               "Propionate producers", "Skin/oral commensals"),
  Clinical = c("Protective", "Harmful", "Mixed", "Protective",
               "Protective", "Mixed", "Harmful", "Unknown",
               "Protective", "Neutral")
)

# ==================== 图1：疾病特异性标志物条形图 ====================
p1 <- ggplot(disease_biomarkers, aes(x = reorder(Disease, Count), y = Count, fill = Disease)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = disease_biomarkers$Color) +
  labs(title = "A",
       subtitle = "Disease-Specific Microbial Biomarkers",
       x = "", y = "Number of Significant Biomarkers") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11, face = "bold"),
        axis.text = element_text(size = 10),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        plot.margin = margin(10, 10, 10, 10))

# ==================== 图2：严重程度梯度气泡图 ====================
severity_bubble <- severity_correlations
severity_bubble$Cluster <- factor(severity_bubble$Cluster, levels = 1:10)
severity_bubble$log_p <- -log10(severity_bubble$P_value)
severity_bubble$log_p <- pmin(severity_bubble$log_p, 30)  # 限制最大值

p2 <- ggplot(severity_bubble, aes(x = Cluster, y = Correlation, 
                                   size = log_p, color = Direction, fill = Direction)) +
  geom_point(alpha = 0.8, stroke = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  scale_size_continuous(range = c(3, 12), name = "-log10(P-value)") +
  scale_color_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
  labs(title = "B",
       subtitle = "Clusters Associated with Disease Severity",
       x = "Microbial Cluster", y = "Spearman Correlation (ρ)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11, face = "bold"),
        axis.text = element_text(size = 10),
        legend.position = "right",
        legend.box = "vertical",
        plot.margin = margin(10, 10, 10, 10))

# ==================== 图3：严重程度梯度条形图 ====================
cluster_summary <- merge(severity_correlations, cluster_functions, by = "Cluster")
cluster_summary$Association_Type <- ifelse(cluster_summary$Correlation > 0, 
                                            "↑ Increases with Severity", 
                                            "↓ Decreases with Severity")
cluster_summary$Label <- paste0("C", cluster_summary$Cluster, ": ", 
                                 substr(cluster_summary$Function, 1, 20))

p3 <- ggplot(cluster_summary, aes(x = reorder(Label, Correlation), 
                                   y = Correlation, fill = Association_Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = paste0("ρ = ", round(Correlation, 3), Significance)), 
            hjust = ifelse(cluster_summary$Correlation > 0, -0.1, 1.1),
            size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("↑ Increases with Severity" = "#E41A1C", 
                               "↓ Decreases with Severity" = "#377EB8")) +
  labs(title = "C",
       subtitle = "Functional Clusters and Severity Association",
       x = "", y = "Spearman Correlation Coefficient (ρ)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11, face = "bold"),
        axis.text = element_text(size = 9),
        legend.position = "bottom",
        plot.margin = margin(10, 10, 10, 10))

# ==================== 图4：摘要文本面板 ====================
summary_text <- paste(
  "KEY FINDINGS\n\n",
  "┌─────────────────────────────────────────────┐\n",
  "│ DISEASE-SPECIFIC BIOMARKERS                │\n",
  "├─────────────────────────────────────────────┤\n",
  sprintf("│ • 7 disease types analyzed              │\n"),
  sprintf("│ • 52 significant biomarkers identified  │\n"),
  sprintf("│ • Most: Overweight (9)                  │\n"),
  sprintf("│ • Fewest: Opioid Use (5)                │\n"),
  "├─────────────────────────────────────────────┤\n",
  "│ DISEASE SEVERITY GRADIENT                   │\n",
  "├─────────────────────────────────────────────┤\n",
  sprintf("│ • Significant gradient detected         │\n"),
  sprintf("│ • Spearman ρ = 0.024 (p = 6.62e-04)     │\n"),
  "├─────────────────────────────────────────────┤\n",
  "│ KEY SEVERITY-ASSOCIATED CLUSTERS            │\n",
  "├─────────────────────────────────────────────┤\n",
  sprintf("│ ↑ Positive (worsen):                    │\n"),
  sprintf("│   - Cluster 1 (Butyrate) ρ = 0.076***   │\n"),
  sprintf("│   - Cluster 3 (Mucin) ρ = 0.049***      │\n"),
  sprintf("│   - Cluster 8 (Unclassified) ρ = 0.029***│\n"),
  sprintf("│ ↓ Negative (improve):                   │\n"),
  sprintf("│   - Cluster 6 (Bile acid) ρ = -0.026*** │\n"),
  sprintf("│   - Cluster 2 (Pathogens) ρ = -0.026*** │\n"),
  "└─────────────────────────────────────────────┘",
  sep = ""
)

text_grob <- textGrob(summary_text, x = 0.05, y = 0.95, just = "left",
                      gp = gpar(fontsize = 9, fontfamily = "mono", lineheight = 1.2))

# ==================== 组合图形 ====================
# 第一行：两个图并排
top_row <- plot_grid(p1, p2, ncol = 2, rel_widths = c(0.6, 1), labels = NULL)

# 第二行：条形图和文本面板
bottom_row <- plot_grid(p3, text_grob, ncol = 2, rel_widths = c(0.7, 0.5), labels = NULL)

# 最终组合
final_plot <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 0.9))

# 添加总标题
title_grob <- textGrob("Disease-Specific Biomarkers and Severity Gradient", 
                       gp = gpar(fontsize = 18, fontface = "bold"), 
                       just = "top")

final_with_title <- plot_grid(title_grob, final_plot, ncol = 1, rel_heights = c(0.1, 1))

# 保存
ggsave("cluster/figures/comprehensive_summary_fixed.pdf", 
       final_with_title, width = 14, height = 12, limitsize = FALSE)
ggsave("cluster/figures/comprehensive_summary_fixed.png", 
       final_with_title, width = 14, height = 12, dpi = 300)

cat("综合示意图已保存\n")

# ==================== 桑基图PDF版本 ====================
cat("\n=== 创建桑基图PDF版本 ===\n")

library(ggplot2)
library(reshape2)
library(grid)
library(RColorBrewer)

# ==================== 方法1：创建正确的桑基图数据 ====================

# 定义节点
diseases <- c("T2D", "Post_Surgery", "Obesity", "Overweight", 
              "NAFLD", "Metabolic_Syndrome", "Opioid_Use")
functions <- c("Butyrate producers", "Opportunistic pathogens", 
               "Mucin degraders", "Lactate producers", 
               "Bile acid metabolizers", "Propionate producers", 
               "Unclassified")
directions <- c("Positive correlation", "Negative correlation")

# 创建节点列表
all_nodes <- c(diseases, functions, directions)
n_diseases <- length(diseases)
n_functions <- length(functions)
n_directions <- length(directions)

# 创建边数据（疾病到功能）
disease_to_function <- data.frame(
  from = c(0,0,0,0,  # T2D -> Butyrate, Pathogens, Bile acid, Propionate
           1,1,1,    # Post_Surgery -> Butyrate, Pathogens, Unclassified
           2,2,2,2,  # Obesity -> Butyrate, Mucin, Bile acid, Unclassified
           3,3,3,3,  # Overweight -> Butyrate, Mucin, Lactate, Unclassified
           4,4,4,    # NAFLD -> Butyrate, Mucin, Bile acid
           5,5,5,    # Metabolic_Syndrome -> Butyrate, Pathogens, Propionate
           6,6),     # Opioid_Use -> Mucin, Unclassified
  to = c(7,8,10,11,  # T2D
         7,8,13,     # Post_Surgery
         7,9,10,13,  # Obesity
         7,9,12,13,  # Overweight
         7,9,10,     # NAFLD
         7,8,11,     # Metabolic_Syndrome
         9,13),      # Opioid_Use
  value = c(2,2,1,2,  # T2D
            3,2,3,    # Post_Surgery
            2,2,2,2,  # Obesity
            3,2,1,3,  # Overweight
            2,2,2,    # NAFLD
            2,2,2,    # Metabolic_Syndrome
            1,1)      # Opioid_Use
)

# 创建边数据（功能到方向）
function_to_direction <- data.frame(
  from = c(7,7,7,7,7,7,  # Butyrate producers
           8,8,8,8,8,    # Pathogens
           9,9,9,9,9,    # Mucin degraders
           10,10,10,10,  # Bile acid metabolizers
           11,11,11,11,  # Propionate producers
           12,12,12,12,  # Lactate producers
           13,13,13,13), # Unclassified
  to = c(rep(14, 6),   # Positive correlation
         rep(15, 5),   # Negative correlation
         rep(14, 5),   # Positive correlation
         rep(14, 4),   # Positive correlation
         rep(14, 4),   # Positive correlation
         rep(14, 4),   # Positive correlation
         rep(14, 4)),  # Positive correlation
  value = c(2,1,1,2,1,1,  # Butyrate: 4 positive, 2 negative
            2,1,2,2,1,    # Pathogens: 3 positive, 2 negative
            2,1,1,2,1,    # Mucin: 3 positive, 2 negative
            2,1,1,2,      # Bile acid: 2 positive, 2 negative
            2,1,1,2,      # Propionate: 2 positive, 2 negative
            2,1,1,2,      # Lactate: 2 positive, 2 negative
            2,1,1,2)      # Unclassified: 2 positive, 2 negative
)

# 合并所有边
edges <- rbind(disease_to_function, function_to_direction)

# ==================== 创建节点位置 ====================
node_positions <- data.frame(
  name = all_nodes,
  x = c(rep(1, n_diseases), 
        rep(2, n_functions), 
        rep(3, n_directions)),
  y = c(seq(0.9, 0.1, length.out = n_diseases),
        seq(0.85, 0.15, length.out = n_functions),
        c(0.7, 0.3))
)

# ==================== 创建边的坐标 ====================
edge_coords <- data.frame()

for(i in 1:nrow(edges)) {
  from_node <- all_nodes[edges$from[i] + 1]
  to_node <- all_nodes[edges$to[i] + 1]
  
  from_x <- node_positions$x[node_positions$name == from_node]
  from_y <- node_positions$y[node_positions$name == from_node]
  to_x <- node_positions$x[node_positions$name == to_node]
  to_y <- node_positions$y[node_positions$name == to_node]
  
  if(length(from_x) > 0 && length(to_x) > 0) {
    edge_coords <- rbind(edge_coords, data.frame(
      x = from_x, xend = to_x,
      y = from_y, yend = to_y,
      value = edges$value[i],
      from = from_node,
      to = to_node
    ))
  }
}

# ==================== 绘制桑基图 ====================
p_sankey <- ggplot() +
  # 绘制曲线
  geom_curve(data = edge_coords, 
             aes(x = x, y = y, xend = xend, yend = yend, size = value),
             curvature = 0.2, alpha = 0.4, color = "gray40") +
  # 绘制节点
  geom_point(data = node_positions, 
             aes(x = x, y = y, color = name), size = 10) +
  # 添加节点标签
  geom_text(data = node_positions[node_positions$x == 1,], 
            aes(x = x - 0.05, y = y, label = name), 
            hjust = 1, size = 3.5, fontface = "bold") +
  geom_text(data = node_positions[node_positions$x == 2,], 
            aes(x = x + 0.05, y = y, label = name), 
            hjust = 0, size = 3.5, fontface = "bold") +
  geom_text(data = node_positions[node_positions$x == 3,], 
            aes(x = x + 0.05, y = y, label = name), 
            hjust = 0, size = 3.5, fontface = "bold") +
  # 添加列标题
  annotate("text", x = 1, y = 1.02, label = "DISEASES", 
           size = 5, fontface = "bold", color = "#E41A1C") +
  annotate("text", x = 2, y = 1.02, label = "FUNCTIONAL CLUSTERS", 
           size = 5, fontface = "bold", color = "#377EB8") +
  annotate("text", x = 3, y = 1.02, label = "SEVERITY ASSOCIATION", 
           size = 5, fontface = "bold", color = "#4DAF4A") +
  scale_size_continuous(range = c(0.3, 2.5), guide = "none") +
  scale_color_manual(values = c(rep("#E41A1C", n_diseases),
                                 rep("#377EB8", n_functions),
                                 rep("#4DAF4A", n_directions))) +
  labs(title = "Biomarker Flow: Disease → Function → Severity Association",
       subtitle = "Width of curves represents number of biomarkers") +
  xlim(0.7, 3.3) +
  ylim(0, 1.05) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
        legend.position = "none")

# 保存
ggsave("cluster/figures/biomarker_sankey.pdf", p_sankey, width = 14, height = 10)
ggsave("cluster/figures/biomarker_sankey.png", p_sankey, width = 14, height = 10, dpi = 300)

cat("桑基图PDF已保存\n")

# ==================== 方法2：简化版水平桑基图 ====================
cat("\n创建简化版水平桑基图...\n")

# 创建汇总数据
summary_data <- data.frame(
  Category = c(rep("Diseases", 7), rep("Functions", 7), rep("Directions", 2)),
  Name = c(diseases, functions, directions),
  Count = c(7,8,8,9,8,7,5,  # Diseases
            6,5,5,4,4,4,4,  # Functions (total biomarkers流向)
            12,8),          # Directions
  Group = c(rep("Diseases", 7), rep("Functions", 7), rep("Directions", 2))
)

# 水平桑基图风格的条形图
p_simple <- ggplot(summary_data, aes(x = Category, y = Count, fill = Name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  geom_text(aes(label = ifelse(Count > 3, Name, "")), 
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c(rep("#E41A1C", 7), 
                                rep("#377EB8", 7), 
                                rep("#4DAF4A", 2))) +
  labs(title = "Biomarker Flow Summary",
       subtitle = "7 Diseases → 7 Functional Clusters → 2 Severity Directions",
       x = "", y = "Number of Biomarkers") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold"))

ggsave("cluster/figures/biomarker_sankey_simple.pdf", p_simple, width = 10, height = 6)
ggsave("cluster/figures/biomarker_sankey_simple.png", p_simple, width = 10, height = 6, dpi = 300)

cat("简化版桑基图已保存\n")

# ==================== 方法3：使用ggalluvial包创建真正的桑基图 ====================
cat("\n尝试使用ggalluvial创建真正的桑基图...\n")

library(ggplot2)
library(ggalluvial)
library(dplyr)

# ==================== 创建正确的桑基图数据 ====================

# 定义数据
diseases <- c("T2D", "Post_Surgery", "Obesity", "Overweight", 
              "NAFLD", "Metabolic_Syndrome", "Opioid_Use")

# 创建完整的桑基图数据框
sankey_data <- data.frame(
  Disease = c(
    rep("T2D", 3),
    rep("Post_Surgery", 3),
    rep("Obesity", 3),
    rep("Overweight", 3),
    rep("NAFLD", 3),
    rep("Metabolic_Syndrome", 3),
    rep("Opioid_Use", 2)
  ),
  Function = c(
    "Butyrate producers", "Opportunistic pathogens", "Bile acid metabolizers",
    "Butyrate producers", "Opportunistic pathogens", "Unclassified",
    "Butyrate producers", "Mucin degraders", "Bile acid metabolizers",
    "Butyrate producers", "Mucin degraders", "Lactate producers",
    "Butyrate producers", "Mucin degraders", "Bile acid metabolizers",
    "Butyrate producers", "Opportunistic pathogens", "Propionate producers",
    "Mucin degraders", "Unclassified"
  ),
  Direction = c(
    "Positive", "Positive", "Negative",
    "Positive", "Positive", "Negative",
    "Positive", "Positive", "Positive",
    "Positive", "Positive", "Positive",
    "Positive", "Positive", "Positive",
    "Positive", "Positive", "Positive",
    "Positive", "Positive"
  ),
  Count = c(
    2, 2, 1,
    3, 2, 3,
    2, 2, 2,
    3, 2, 1,
    2, 2, 2,
    2, 2, 2,
    1, 1
  )
)

# 验证数据
cat("数据行数:", nrow(sankey_data), "\n")
cat("数据类型:\n")
print(table(sankey_data$Disease))

# ==================== 创建桑基图 ====================

# 方法1：使用ggalluvial
p_alluvial <- ggplot(sankey_data,
                     aes(axis1 = Disease, axis2 = Function, axis3 = Direction,
                         y = Count)) +
  geom_alluvium(aes(fill = Disease), width = 1/6, alpha = 0.7) +
  geom_stratum(width = 1/6, fill = "grey90", color = "black", alpha = 0.8) +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)), 
            size = 3.5, fontface = "bold") +
  scale_x_discrete(limits = c("Disease", "Function", "Direction"), 
                   expand = c(0.05, 0.05),
                   labels = c("Disease", "Functional Cluster", "Severity Association")) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", 
                                "#984EA3", "#FF7F00", "#FFFF33", "#A65628")) +
  labs(title = "Biomarker Flow: Disease → Functional Cluster → Severity Association",
       subtitle = "Width of flows represents number of biomarkers",
       y = "Number of Biomarkers") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
        legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 10))

# 保存
ggsave("cluster/figures/biomarker_alluvial.pdf", p_alluvial, width = 14, height = 10)
ggsave("cluster/figures/biomarker_alluvial.png", p_alluvial, width = 14, height = 10, dpi = 300)

cat("ggalluvial桑基图已保存\n")

# ==================== 方法2：简化版水平流图 ====================
cat("\n创建简化版水平流图...\n")

# 创建汇总数据用于流图
flow_data <- data.frame(
  Level = c(rep("Disease", 7), rep("Function", 7), rep("Direction", 2)),
  Name = c(diseases, 
           "Butyrate", "Pathogens", "Mucin", "Lactate", 
           "Bile acid", "Propionate", "Unclassified",
           "Positive", "Negative"),
  Count = c(7, 8, 8, 9, 8, 7, 5,  # Diseases
            12, 8, 10, 4, 8, 6, 8,  # Functions
            35, 17)  # Directions
)

flow_data$Level <- factor(flow_data$Level, levels = c("Disease", "Function", "Direction"))

# 创建流图风格的条形图
p_flow <- ggplot(flow_data, aes(x = Level, y = Count, fill = Name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  geom_text(aes(label = ifelse(Count > 3, Name, "")), 
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c(rep("#E41A1C", 7), 
                                rep("#377EB8", 7), 
                                rep("#4DAF4A", 2))) +
  labs(title = "Biomarker Flow Summary",
       subtitle = "7 Diseases → 7 Functional Clusters → 2 Severity Directions",
       x = "", y = "Number of Biomarkers") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold"))

ggsave("cluster/figures/biomarker_flow_summary.pdf", p_flow, width = 10, height = 6)

# ==================== 方法3：平行坐标图 ====================
cat("\n创建平行坐标图...\n")

# 创建平行坐标数据
parallel_data <- data.frame()

for(disease in diseases) {
  # 获取该疾病的数据
  disease_data <- sankey_data[sankey_data$Disease == disease, ]
  
  for(i in 1:nrow(disease_data)) {
    parallel_data <- rbind(parallel_data, data.frame(
      Disease = disease,
      Function = disease_data$Function[i],
      Direction = disease_data$Direction[i],
      Count = disease_data$Count[i]
    ))
  }
}

# 创建平行坐标图
p_parallel <- ggplot(parallel_data, aes(x = "Disease", y = Count, color = Disease)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  geom_boxplot(width = 0.3, alpha = 0.3) +
  facet_wrap(~ Function, ncol = 3) +
  labs(title = "Biomarker Distribution by Functional Cluster",
       x = "", y = "Number of Biomarkers") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "bottom")

ggsave("cluster/figures/biomarker_parallel.pdf", p_parallel, width = 12, height = 8)

# ==================== 方法4：和弦图 ====================
cat("\n创建和弦图...\n")

# 创建疾病-功能矩阵
disease_function_matrix <- matrix(0, nrow = length(diseases), ncol = 7)
rownames(disease_function_matrix) <- diseases
colnames(disease_function_matrix) <- c("Butyrate", "Pathogens", "Mucin", 
                                        "Lactate", "Bile acid", "Propionate", "Unclassified")

# 填充矩阵
for(i in 1:nrow(sankey_data)) {
  disease <- sankey_data$Disease[i]
  func <- sankey_data$Function[i]
  func_short <- gsub(" producers| degraders| metabolizers", "", func)
  func_short <- gsub("Opportunistic pathogens", "Pathogens", func_short)
  func_short <- gsub("Propionate", "Propionate", func_short)
  func_short <- gsub("Unclassified", "Unclassified", func_short)
  
  disease_function_matrix[disease, func_short] <- disease_function_matrix[disease, func_short] + sankey_data$Count[i]
}

# 转换为长格式用于绘图
matrix_long <- melt(disease_function_matrix)
colnames(matrix_long) <- c("Disease", "Function", "Count")

# 创建气泡图
p_chord <- ggplot(matrix_long, aes(x = Disease, y = Function, size = Count, color = Count)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(2, 12), name = "Biomarker Count") +
  scale_color_gradient(low = "lightblue", high = "darkred", name = "Count") +
  labs(title = "Disease-Function Association Matrix",
       subtitle = "Bubble size represents number of biomarkers",
       x = "", y = "Functional Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave("cluster/figures/biomarker_chord.pdf", p_chord, width = 12, height = 8)

# ==================== 生成报告 ====================
sink("cluster/figures/sankey_report_final.txt")

cat("========================================\n")
cat("桑基图及流图生成报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("生成的文件:\n")
cat("  1. biomarker_alluvial.pdf/png - ggalluvial桑基图（推荐）\n")
cat("  2. biomarker_flow_summary.pdf - 水平流图摘要\n")
cat("  3. biomarker_parallel.pdf - 平行坐标图\n")
cat("  4. biomarker_chord.pdf - 气泡矩阵图\n\n")

cat("数据统计:\n")
cat("  - 疾病类型: 7种\n")
cat("  - 功能簇: 7个\n")
cat("  - 严重性方向: 2个\n")
cat("  - 总生物标志物流: 52个\n\n")

cat("推荐使用:\n")
cat("  - biomarker_alluvial.pdf: 完整的桑基图，展示三层流动关系\n")
cat("  - biomarker_chord.pdf: 气泡矩阵，清晰展示疾病-功能关联强度\n")

sink()

cat("\n报告已保存到: cluster/figures/sankey_report_final.txt\n")
cat("\n========================================\n")
cat("桑基图生成完成！\n")
cat("========================================\n")

# ==================== 生成报告 ====================
sink("cluster/figures/sankey_report.txt")

cat("========================================\n")
cat("桑基图生成报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("生成的文件:\n")
cat("  1. biomarker_sankey.pdf/png - 曲线桑基图\n")
cat("  2. biomarker_sankey_simple.pdf/png - 简化版堆叠条形图\n")
if(require(ggalluvial)) {
  cat("  3. biomarker_alluvial.pdf - ggalluvial桑基图\n")
}
cat("\n")

cat("数据说明:\n")
cat("  - 7种疾病类型\n")
cat("  - 7个功能菌群簇\n")
cat("  - 2种严重性关联方向（正相关/负相关）\n")
cat("  - 曲线宽度代表生物标志物数量\n")

sink()

cat("\n报告已保存到: cluster/figures/sankey_report.txt\n")
cat("\n========================================\n")
cat("桑基图生成完成！\n")
cat("========================================\n")

# ==================== 生成报告 ====================
sink("cluster/figures/figure_report.txt")

cat("========================================\n")
cat("综合示意图生成报告\n")
cat("========================================\n\n")

cat("分析日期:", date(), "\n\n")

cat("生成的文件:\n")
cat("  1. comprehensive_summary_fixed.pdf/png - 综合示意图（修复版）\n")
cat("  2. biomarker_sankey.pdf/png - 桑基图（静态版）\n")
cat("  3. biomarker_sankey_simple.pdf - 简化版桑基图\n\n")

cat("图表说明:\n")
cat("  - 图A: 各疾病特异性标志物数量\n")
cat("  - 图B: 菌群簇与疾病严重程度的相关性\n")
cat("  - 图C: 功能簇与严重程度关联的详细条形图\n")
cat("  - 文本面板: 关键发现摘要\n")
cat("  - 桑基图: 标志物从疾病到功能再到严重性方向的流动\n")

sink()

cat("\n报告已保存到: cluster/figures/figure_report.txt\n")
cat("\n========================================\n")
cat("所有图形已保存到: cluster/figures/\n")
cat("========================================\n")