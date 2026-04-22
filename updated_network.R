# ==================== 修复网络图 ====================
cat("\n=== 修复网络图 ===\n")

# 确保目录存在
if(!dir.exists("cluster/cluster_species/figures")) {
  dir.create("cluster/cluster_species/figures", recursive = TRUE)
}

library(igraph)
library(ggplot2)

# 获取计数数据
counts <- assay(tse_2w, "counts")

# 合并所有cluster的数据
all_species_df <- do.call(rbind, cluster_species_list)

# 选择Top30菌种（按丰度）
all_species <- all_species_df$Species
mean_abundance_all <- all_species_df$Mean_Abundance
top30_species <- head(all_species[order(mean_abundance_all, decreasing = TRUE)], 30)

cat("使用Top30菌种构建网络...\n")

# 计算相关性
cor_matrix <- cor(t(counts[top30_species, ]), method = "spearman")

# 使用阈值
threshold <- 0.4
adj_matrix <- abs(cor_matrix) > threshold & cor_matrix > 0
diag(adj_matrix) <- FALSE

# 创建图
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

# 移除孤立节点
g <- delete.vertices(g, degree(g) == 0)

cat("网络节点数:", vcount(g), "\n")
cat("网络边数:", ecount(g), "\n")

if(vcount(g) > 0) {
  
  # 添加节点属性
  V(g)$cluster <- sapply(V(g)$name, function(sp) {
    for(i in 1:10) {
      if(sp %in% cluster_species_list[[i]]$Species) return(i)
    }
    return(0)
  })
  
  # 计算节点大小（基于丰度）
  sp_abundance <- rowMeans(counts[V(g)$name, ])
  V(g)$size <- scales::rescale(log1p(sp_abundance), to = c(3, 8))
  
  # ==================== 方法1：使用igraph内置布局（最稳定）====================
  
  # 布局1：Fruchterman-Reingold
  set.seed(123)
  layout_fr <- layout_with_fr(g, niter = 1000)
  
  pdf("cluster/cluster_species/figures/species_network_fr.pdf", 
      width = 12, height = 10)
  
  plot(g, 
       layout = layout_fr,
       vertex.color = rainbow(10)[V(g)$cluster],
       vertex.size = V(g)$size * 2,
       vertex.label = V(g)$name,
       vertex.label.cex = 0.7,
       vertex.label.dist = 1,
       vertex.label.color = "black",
       edge.width = 1,
       edge.color = "gray50",
       main = "Species Co-occurrence Network (Fruchterman-Reingold)")
  
  legend("topright", legend = paste("Cluster", 1:10), 
         col = rainbow(10), pch = 19, cex = 0.8)
  
  dev.off()
  
  # 布局2：Kamada-Kawai
  pdf("cluster/cluster_species/figures/species_network_kk.pdf", 
      width = 12, height = 10)
  
  plot(g, 
       layout = layout_with_kk(g),
       vertex.color = rainbow(10)[V(g)$cluster],
       vertex.size = V(g)$size * 2,
       vertex.label = V(g)$name,
       vertex.label.cex = 0.7,
       vertex.label.dist = 1,
       vertex.label.color = "black",
       edge.width = 1,
       edge.color = "gray50",
       main = "Species Co-occurrence Network (Kamada-Kawai)")
  
  legend("topright", legend = paste("Cluster", 1:10), 
         col = rainbow(10), pch = 19, cex = 0.8)
  
  dev.off()
  
  # ==================== 方法2：手动提取坐标用ggplot2绘制 ====================
  cat("生成ggplot2版本网络图...\n")
  
  # 获取布局坐标
  set.seed(123)
  layout_coords <- layout_with_fr(g, niter = 1000)
  colnames(layout_coords) <- c("x", "y")
  
  # 创建节点数据框
  nodes_df <- data.frame(
    name = V(g)$name,
    x = layout_coords[, 1],
    y = layout_coords[, 2],
    cluster = V(g)$cluster,
    size = V(g)$size,
    abundance = sp_abundance[V(g)$name]
  )
  
  # 创建边数据框
  edges_df <- as_data_frame(g, what = "edges")
  edges_df$from_name <- edges_df$from
  edges_df$to_name <- edges_df$to
  
  # 添加坐标
  edges_df$x_from <- nodes_df$x[match(edges_df$from, nodes_df$name)]
  edges_df$y_from <- nodes_df$y[match(edges_df$from, nodes_df$name)]
  edges_df$x_to <- nodes_df$x[match(edges_df$to, nodes_df$name)]
  edges_df$y_to <- nodes_df$y[match(edges_df$to, nodes_df$name)]
  
  # 使用ggplot2绘制
  p <- ggplot() +
    # 绘制边
    geom_segment(data = edges_df, 
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
                 color = "gray70", alpha = 0.5, size = 0.5) +
    # 绘制节点
    geom_point(data = nodes_df, 
               aes(x = x, y = y, color = factor(cluster), size = size),
               alpha = 0.8) +
    # 绘制标签
    geom_text(data = nodes_df,
              aes(x = x, y = y, label = name),
              size = 3, vjust = -0.5, hjust = 0.5,
              check_overlap = TRUE) +
    scale_color_manual(values = rainbow(10), name = "Cluster") +
    scale_size_continuous(range = c(3, 8), guide = "none") +
    labs(title = "Species Co-occurrence Network (Top 30 species, r > 0.4)") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom")
  
  ggsave("cluster/cluster_species/figures/species_network_ggplot.pdf",
         p, width = 14, height = 12)
  
  cat("ggplot2网络图已保存\n")
  
  # ==================== 方法3：圆形布局（避免重叠）====================
  pdf("cluster/cluster_species/figures/species_network_circle.pdf", 
      width = 12, height = 10)
  
  plot(g, 
       layout = layout_in_circle(g),
       vertex.color = rainbow(10)[V(g)$cluster],
       vertex.size = V(g)$size * 1.5,
       vertex.label = V(g)$name,
       vertex.label.cex = 0.6,
       vertex.label.dist = 0.5,
       vertex.label.color = "black",
       edge.width = 0.8,
       edge.color = "gray50",
       main = "Species Co-occurrence Network (Circular layout)")
  
  legend("topleft", legend = paste("Cluster", 1:10), 
         col = rainbow(10), pch = 19, cex = 0.7, pt.cex = 1)
  
  dev.off()
  
  # ==================== 方法4：按Cluster分组的网络图 ====================
  cat("生成按Cluster分组的网络图...\n")
  
  # 只保留有cluster信息的节点
  g_filtered <- delete.vertices(g, V(g)$cluster == 0)
  
  if(vcount(g_filtered) > 0) {
    set.seed(123)
    layout_filtered <- layout_with_fr(g_filtered, niter = 1000)
    
    pdf("cluster/cluster_species/figures/species_network_by_cluster.pdf", 
        width = 14, height = 12)
    
    plot(g_filtered, 
         layout = layout_filtered,
         vertex.color = rainbow(10)[V(g_filtered)$cluster],
         vertex.size = 8,
         vertex.label = V(g_filtered)$name,
         vertex.label.cex = 0.7,
         vertex.label.dist = 1,
         vertex.label.color = "black",
         edge.width = 1,
         edge.color = "gray50",
         main = "Species Co-occurrence Network (Grouped by Cluster)")
    
    legend("topright", legend = paste("Cluster", 1:10), 
           col = rainbow(10), pch = 19, cex = 0.8)
    
    dev.off()
  }
  
  # ==================== 图5：简化版（节点大小统一）====================
  pdf("cluster/cluster_species/figures/species_network_simple.pdf", 
      width = 10, height = 8)
  
  par(mar = c(1, 1, 3, 1))
  
  plot(g, 
       layout = layout_with_fr(g, niter = 2000),
       vertex.color = rainbow(10)[V(g)$cluster],
       vertex.size = 6,  # 统一大小
       vertex.label = V(g)$name,
       vertex.label.cex = 0.6,
       vertex.label.dist = 0.5,
       vertex.label.degree = pi/4,
       vertex.frame.color = NA,
       edge.width = 0.5,
       edge.color = rgb(0.5, 0.5, 0.5, 0.5),
       main = "Species Co-occurrence Network")
  
  legend("bottomright", legend = paste("Cluster", 1:10), 
         col = rainbow(10), pch = 19, cex = 0.6, 
         bty = "n", pt.cex = 1)
  
  dev.off()
  
  # ==================== 计算网络统计 ====================
  cat("\n计算网络统计...\n")
  
  network_stats <- data.frame(
    Metric = c("Nodes", "Edges", "Density", "Transitivity", 
               "Average Degree", "Diameter", "Average Path Length"),
    Value = c(vcount(g), ecount(g), 
              round(edge_density(g), 4),
              round(transitivity(g), 4),
              round(mean(degree(g)), 2),
              diameter(g),
              round(mean_distance(g), 2))
  )
  
  print(network_stats)
  write.csv(network_stats, "cluster/cluster_species/figures/network_stats.csv", 
            row.names = FALSE)
  
  # 计算节点中心性
  centrality_df <- data.frame(
    Species = V(g)$name,
    Cluster = V(g)$cluster,
    Degree = degree(g),
    Betweenness = round(betweenness(g), 2),
    Closeness = round(closeness(g), 4),
    Eigenvector = round(eigen_centrality(g)$vector, 4)
  )
  
  # 按Degree排序
  centrality_df <- centrality_df[order(-centrality_df$Degree), ]
  write.csv(centrality_df, "cluster/cluster_species/figures/node_centrality.csv", 
            row.names = FALSE)
  
  cat("\nTop 5 中心节点:\n")
  print(head(centrality_df, 5))
  
} else {
  cat("网络中没有节点，请调整相关性阈值\n")
}

# ==================== 完成 ====================
cat("\n\n========================================\n")
cat("网络图生成完成！\n")
cat("生成的文件:\n")
cat("  - species_network_fr.pdf (力导向布局)\n")
cat("  - species_network_kk.pdf (Kamada-Kawai布局)\n")
cat("  - species_network_circle.pdf (圆形布局)\n")
cat("  - species_network_ggplot.pdf (ggplot2版本)\n")
cat("  - species_network_by_cluster.pdf (按Cluster分组)\n")
cat("  - species_network_simple.pdf (简化版)\n")
cat("  - network_stats.csv (网络统计)\n")
cat("  - node_centrality.csv (节点中心性)\n")
cat("========================================\n")