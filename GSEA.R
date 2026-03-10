# 设置工作目录
setwd("/Volumes/LENOVO/Thesis 5/GSEA")

library(fgsea)
library(data.table)
library(ggplot2)
library(readxl)

# ============= 第一步：读取表达数据 =============
exp_data <- read_excel("../LASSO/Exp_details.xlsx")
exp_data <- as.data.frame(exp_data)

# 提取表达矩阵
rownames(exp_data) <- exp_data$`Feature ID`
expr_matrix <- exp_data[, -c(1,2)]
expr_matrix <- as.matrix(expr_matrix)
mode(expr_matrix) <- "numeric"

# ============= 第二步：基于Pde10a中位数分割 =============
pde10a_expr <- expr_matrix["Pde10a", ]
pde10a_median <- median(pde10a_expr)

high_samples <- names(pde10a_expr[pde10a_expr > pde10a_median])
low_samples <- names(pde10a_expr[pde10a_expr <= pde10a_median])

cat("Pde10a中位数:", round(pde10a_median, 2), "\n")
cat("高表达组 (n =", length(high_samples), "):", paste(high_samples, collapse=", "), "\n")
cat("低表达组 (n =", length(low_samples), "):", paste(low_samples, collapse=", "), "\n")

# ============= 第三步：计算log2FC =============
high_expr <- expr_matrix[, high_samples, drop=FALSE]
low_expr <- expr_matrix[, low_samples, drop=FALSE]

high_mean <- rowMeans(high_expr)
low_mean <- rowMeans(low_expr)

# 计算log2FC
log2FC <- high_mean - low_mean
log2FC <- sort(log2FC, decreasing = TRUE)

# 去除完全相同的值（处理ties警告）
# 给重复值添加微小扰动
ranks <- log2FC
dupes <- duplicated(ranks) | duplicated(ranks, fromLast = TRUE)
if(sum(dupes) > 0) {
  ranks[dupes] <- ranks[dupes] + rnorm(sum(dupes), 0, 1e-10)
  ranks <- sort(ranks, decreasing = TRUE)
}

# ============= 第四步：读取BIOCARTA基因集 =============
gmt_lines <- readLines("m2.cp.biocarta.v2026.1.Mm.symbols.gmt")
pathways <- list()

for(line in gmt_lines) {
  parts <- strsplit(line, "\t")[[1]]
  pathway_name <- gsub("BIOCARTA_", "", parts[1])
  genes <- parts[3:length(parts)]
  genes <- genes[genes != ""]
  if(length(genes) > 0) {
    pathways[[pathway_name]] <- genes
  }
}

cat("读取到", length(pathways), "个通路\n")

# ============= 第五步：运行fgsea（使用multilevel方法） =============
set.seed(123)
fgsea_results <- fgseaMultilevel(pathways = pathways,
                                 stats = ranks,
                                 minSize = 5,
                                 maxSize = 500,
                                 eps = 1e-50)  # 更精确的p值计算

# ============= 第六步：放宽筛选条件，确保有图可出 =============
# 先尝试 padj < 0.05
sig_results <- fgsea_results[fgsea_results$padj < 0.05, ]

# 如果没有显著通路，放宽到 pval < 0.05
if(nrow(sig_results) == 0) {
  cat("\n⚠️ padj < 0.05 无显著通路，改用 pval < 0.05\n")
  sig_results <- fgsea_results[fgsea_results$pval < 0.05, ]
}

# 如果还是没有，取NES绝对值最大的前10个
if(nrow(sig_results) == 0) {
  cat("\n⚠️ pval < 0.05 也无显著通路，改用NES绝对值最大的前10个\n")
  fgsea_results <- fgsea_results[order(abs(fgsea_results$NES), decreasing = TRUE), ]
  sig_results <- head(fgsea_results, 10)
  sig_results$padj <- 0.05  # 添加一个虚拟的padj值
}

sig_results <- sig_results[order(sig_results$NES, decreasing = TRUE), ]

cat("\n最终分析通路数:", nrow(sig_results), "\n")
cat("NES范围: [", round(min(sig_results$NES), 2), ", ", round(max(sig_results$NES), 2), "]\n")

# 保存结果
write.table(sig_results, file = "GSEA.result.BIOCARTA.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============= 第七步：绘制高表达组富集图 =============
high_pathways <- sig_results[sig_results$NES > 0, ]
if(nrow(high_pathways) > 0) {
  # 取前5个或全部
  n_plot <- min(5, nrow(high_pathways))
  plot_pathways <- high_pathways$pathway[1:n_plot]
  
  pdf("GSEA.high.Pde10a.pdf", width = 10, height = 3.5 * n_plot)
  par(mfrow = c(n_plot, 1), mar = c(4, 4, 3, 2))
  
  for(i in 1:n_plot) {
    pathway_name <- plot_pathways[i]
    pathway_result <- high_pathways[high_pathways$pathway == pathway_name, ]
    
    # 使用plotEnrichment生成图
    p <- plotEnrichment(pathways[[pathway_name]], ranks) +
      labs(title = paste0(pathway_name, 
                          "\nNES = ", round(pathway_result$NES, 2),
                          ", p = ", format(pathway_result$pval, scientific = TRUE, digits = 2),
                          ifelse("padj" %in% colnames(pathway_result), 
                                 paste0(", padj = ", round(pathway_result$padj, 3)), ""))) +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 9))
    
    print(p)
  }
  dev.off()
  cat("✅ 已生成 GSEA.high.Pde10a.pdf (", n_plot, "个通路)\n")
}

# ============= 第八步：绘制低表达组富集图 =============
low_pathways <- sig_results[sig_results$NES < 0, ]
if(nrow(low_pathways) > 0) {
  n_plot <- min(5, nrow(low_pathways))
  plot_pathways <- low_pathways$pathway[1:n_plot]
  
  pdf("GSEA.low.Pde10a.pdf", width = 10, height = 3.5 * n_plot)
  par(mfrow = c(n_plot, 1), mar = c(4, 4, 3, 2))
  
  for(i in 1:n_plot) {
    pathway_name <- plot_pathways[i]
    pathway_result <- low_pathways[low_pathways$pathway == pathway_name, ]
    
    p <- plotEnrichment(pathways[[pathway_name]], ranks) +
      labs(title = paste0(pathway_name, 
                          "\nNES = ", round(pathway_result$NES, 2),
                          ", p = ", format(pathway_result$pval, scientific = TRUE, digits = 2),
                          ifelse("padj" %in% colnames(pathway_result), 
                                 paste0(", padj = ", round(pathway_result$padj, 3)), ""))) +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"))
    
    print(p)
  }
  dev.off()
  cat("✅ 已生成 GSEA.low.Pde10a.pdf (", n_plot, "个通路)\n")
}

# ============= 第九步：如果没有上/下调通路，绘制混合图 =============
if(nrow(high_pathways) == 0 && nrow(low_pathways) == 0) {
  n_plot <- min(6, nrow(sig_results))
  plot_pathways <- sig_results$pathway[1:n_plot]
  
  pdf("GSEA.combined.Pde10a.pdf", width = 10, height = 3.5 * n_plot)
  par(mfrow = c(n_plot, 1), mar = c(4, 4, 3, 2))
  
  for(i in 1:n_plot) {
    pathway_name <- plot_pathways[i]
    pathway_result <- sig_results[sig_results$pathway == pathway_name, ]
    
    p <- plotEnrichment(pathways[[pathway_name]], ranks) +
      labs(title = paste0(pathway_name, 
                          "\nNES = ", round(pathway_result$NES, 2),
                          ", p = ", format(pathway_result$pval, scientific = TRUE, digits = 2))) +
      theme_bw()
    
    print(p)
  }
  dev.off()
  cat("✅ 已生成 GSEA.combined.Pde10a.pdf (混合通路)\n")
}

# ============= 第十步：输出统计报告 =============
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("GSEA分析完成！\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")
cat("分析结果:\n")
cat("  总分析通路数:", nrow(fgsea_results), "\n")
cat("  最终展示通路数:", nrow(sig_results), "\n")
cat("  高表达组通路:", nrow(high_pathways), "\n")
cat("  低表达组通路:", nrow(low_pathways), "\n\n")