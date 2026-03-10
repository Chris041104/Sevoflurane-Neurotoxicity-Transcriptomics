# 设置工作目录
setwd("/Volumes/LENOVO/Thesis 5/SVM")

# 安装必要的包（如果未安装）
# install.packages("e1071")
# install.packages("kernlab")
# install.packages("caret")
# install.packages("readxl")

library(e1071)
library(kernlab)
library(caret)
library(limma)
library(readxl)

# 设置随机数种子，固定结果
set.seed(123)

# ============= 第一步：读取基因列表和表达数据 =============
# 读取Genes.txt中的基因列表
genes_to_analyze <- read.table("../LASSO/Genes.txt", header=FALSE, stringsAsFactors=FALSE)
genes_to_analyze <- genes_to_analyze$V1
cat("从Genes.txt中读取到", length(genes_to_analyze), "个基因\n")
print(head(genes_to_analyze))

# 读取Exp_details.xlsx文件
exp_data <- read_excel("../LASSO/Exp_details.xlsx")
cat("\n原始表达数据维度:", dim(exp_data), "\n")

# 将exp_data转换为数据框
exp_data <- as.data.frame(exp_data)

# 检查Feature ID列中是否有匹配的基因
matching_genes <- intersect(genes_to_analyze, exp_data$`Feature ID`)
cat("\n在表达数据中找到", length(matching_genes), "个匹配的基因\n")

if(length(matching_genes) == 0) {
  stop("错误：在表达数据中没有找到任何匹配的基因！")
}

# 提取匹配基因的表达数据
data_subset <- exp_data[exp_data$`Feature ID` %in% matching_genes, ]

# 将Feature ID设为行名，移除Gene ID列
rownames(data_subset) <- data_subset$`Feature ID`
data_subset <- data_subset[, -c(1,2)]

# 转化为matrix
dimnames=list(rownames(data_subset), colnames(data_subset))
data=matrix(as.numeric(as.matrix(data_subset)), nrow=nrow(data_subset), dimnames=dimnames)

# ============= 第二步：数据预处理 =============
cat("\n=== 数据预处理 ===\n")
cat("原始数据维度:", dim(data), "\n")

# 去除所有表达值都为0的基因
data = data[rowSums(data) > 0, ]
cat("去除全0基因后剩余", nrow(data), "个基因\n")

# 标准化
if(nrow(data) > 0) {
  data=normalizeBetweenArrays(data)
  cat("数据标准化完成\n")
} else {
  stop("错误：所有基因的表达值都为0，无法进行分析")
}

# 导出标准化后的数据
write.table(data.frame(ID=rownames(data), data), 
            file="normalize.txt", sep="\t", quote=F, row.names = F)

# ============= 第三步：定义分组（合并0.7MAC和2MAC）=============
# 获取所有样本名
all_samples <- colnames(data)

# 提取各组样本
CTRL_samples <- all_samples[grep("CTRL", all_samples)]
TREAT_samples <- all_samples[grep("MAC", all_samples)]

cat("\n=== 样本分组 ===\n")
cat("CTRL样本:", paste(CTRL_samples, collapse=", "), "\n")
cat("TREAT样本 (0.7MAC+2MAC):", paste(TREAT_samples, collapse=", "), "\n")

# 保存分组文件（兼容原代码）
write.table(data.frame(CTRL_samples), file="Control.txt", 
            sep="\t", quote=F, row.names=F, col.names=F)
write.table(data.frame(TREAT_samples), file="Treat.txt", 
            sep="\t", quote=F, row.names=F, col.names=F)

# 设置分组信息
conNum = length(CTRL_samples)
treatNum = length(TREAT_samples)
Type = c(rep(1, conNum), rep(2, treatNum))

# 按照CTRL、TREAT顺序排列数据
data1 = data[, CTRL_samples, drop=FALSE]
data2 = data[, TREAT_samples, drop=FALSE]
data = cbind(data1, data2)

# ============= 第四步：差异表达分析（生成diff.Wilcoxon.txt）=============
cat("\n=== 差异表达分析 ===\n")

# 计算差异表达
p_values <- apply(data, 1, function(x) {
  tryCatch({
    t.test(x[1:conNum], x[(conNum+1):ncol(data)])$p.value
  }, error = function(e) NA)
})

logFC <- apply(data, 1, function(x) {
  mean(x[(conNum+1):ncol(data)]) - mean(x[1:conNum])
})

# 创建差异表达结果表
deg_results <- data.frame(
  logFC = logFC,
  p_value = p_values,
  fdr = p.adjust(p_values, method="BH"),
  stringsAsFactors = FALSE
)
rownames(deg_results) <- rownames(data)

# 按logFC排序（从小到大）
deg_results <- deg_results[order(deg_results$logFC, decreasing=FALSE), ]

cat("\n差异表达基因统计（fdr<0.05）:", sum(deg_results$fdr < 0.05, na.rm=TRUE), "\n")
cat("\nlogFC最小的前10个基因（下调最显著）:\n")
print(head(deg_results[, c("logFC", "p_value", "fdr")], 10))
cat("\nlogFC最大的前10个基因（上调最显著）:\n")
print(tail(deg_results[, c("logFC", "p_value", "fdr")], 10))

# 保存差异表达结果（兼容原代码）
write.table(data.frame(ID=rownames(deg_results), deg_results), 
            file="diff.Wilcoxon.txt", sep="\t", quote=F, row.names=F)

# ============= 第五步：选择用于SVM的基因（原代码的筛选逻辑）=============
# 按照原代码逻辑：选取logFC最小和最大的各15个基因
genes = deg_results

# 选取前后15个（共30个基因）
gene_num <- min(15, floor(nrow(genes)/2))
selected_genes <- c(rownames(genes)[1:gene_num], 
                    rownames(genes)[(nrow(genes)-gene_num+1):nrow(genes)])

cat("\n=== 选择用于SVM-RFE的基因 ===\n")
cat("选择了logFC最小和最大的各", gene_num, "个基因，共", length(selected_genes), "个基因\n")
print(selected_genes)

# 提取这些基因的表达数据
data_svm <- data[selected_genes, , drop=FALSE]

# ============= 第六步：SVM-RFE分析（完全按照原代码）=============
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("SVM-RFE分析\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# 准备数据
x = t(data_svm)
y = as.numeric(as.factor(Type))

cat("\nSVM输入数据维度:", dim(x), "\n")
cat("样本数:", nrow(x), ", 基因数:", ncol(x), "\n")
cat("分组分布:\n")
print(table(Type))

# SVM-RFE分析（完全按照原代码）
feature_sizes <- c(seq(1, min(15, ncol(x)), by=2))
cat("\n测试的特征数量:", paste(feature_sizes, collapse=", "), "\n")

Profile = rfe(x = x,
              y = y,
              sizes = feature_sizes,
              rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
              method = "svmRadial")

# ============= 第七步：画图（完全按照原代码）=============
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x_vals = Profile$results$Variables
y_vals = Profile$results$RMSE
plot(x_vals, y_vals, 
     xlab="Variables", 
     ylab="RMSE (Cross-Validation)", 
     col="darkgreen",
     main="SVM-RFE Feature Selection")
lines(x_vals, y_vals, col="darkgreen", lwd=2)

# 标记最优特征数
wmin = which.min(y_vals)
wmin.x = x_vals[wmin]
wmin.y = y_vals[wmin]
points(wmin.x, wmin.y, col="blue", pch=16, cex=1.5)
text(wmin.x, wmin.y, paste0('N=', wmin.x), pos=2, col="red", font=2)

# 添加网格
grid()
dev.off()

cat("\n已生成 SVM-RFE.pdf\n")

# ============= 第八步：输出筛选的基因（完全按照原代码）=============
featureGenes = Profile$optVariables

cat("\n=== SVM-RFE筛选结果 ===\n")
cat("筛选出的基因数:", length(featureGenes), "\n")
cat("筛选出的基因:\n")
print(featureGenes)

# 输出基因列表
write.table(file="SVM-RFE.gene.txt", featureGenes, 
            sep="\t", quote=F, row.names=F, col.names=F)

# 表达矩阵
if(length(featureGenes) > 0) {
  # 确保基因名在数据中存在
  valid_genes <- intersect(featureGenes, rownames(data_svm))
  
  if(length(valid_genes) > 0) {
    sigExp = data_svm[valid_genes, , drop=FALSE]
    write.table(data.frame(ID=rownames(sigExp), sigExp), 
                file="SVM-RFE.Gene.Exp.txt", sep="\t", quote=F, col.names=TRUE, row.names=FALSE)
    cat("\n已保存表达谱到 SVM-RFE.Gene.Exp.txt\n")
  }
}

# ============= 第九步：检查目标基因 =============
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("检查目标基因\n")
cat(paste(rep("=", 60), collapse=""), "\n")

target_genes <- c("Nr5a2", "Nr4a2", "Pde10a")

for(gene in target_genes) {
  if(gene %in% rownames(data)) {
    in_selected <- gene %in% selected_genes
    in_svm <- gene %in% featureGenes
    
    cat("\n", gene, ":\n", sep="")
    cat("  在表达数据中: ✓\n")
    cat("  在差异表达筛选（前30）中:", if(in_selected) "✓" else "✗", "\n")
    cat("  被SVM-RFE选中:", if(in_svm) "✓" else "✗", "\n")
    
    if(gene %in% rownames(deg_results)) {
      idx <- which(rownames(deg_results) == gene)
      cat("  logFC:", round(deg_results[idx, "logFC"], 4), "\n")
      cat("  p-value:", format(deg_results[idx, "p_value"], scientific=TRUE, digits=3), "\n")
    }
  } else {
    cat("\n", gene, "不在表达数据中\n", sep="")
  }
}

# 如果三个基因都被选中
if(all(c("Nr5a2", "Nr4a2", "Pde10a") %in% featureGenes)) {
  cat("\n🎉 恭喜！Nr5a2、Nr4a2、Pde10a 同时被SVM-RFE选中！\n")
} else if(sum(c("Nr5a2", "Nr4a2", "Pde10a") %in% featureGenes) >= 2) {
  cat("\n👍 有", sum(c("Nr5a2", "Nr4a2", "Pde10a") %in% featureGenes), "个目标基因被选中\n")
}

# ============= 第十步：生成热图 =============
if(length(featureGenes) >= 2 && require(pheatmap, quietly=TRUE)) {
  # 准备热图数据
  valid_genes <- intersect(featureGenes, rownames(data_svm))
  
  if(length(valid_genes) >= 2) {
    heatmap_data <- data_svm[valid_genes, , drop=FALSE]
    
    # 创建注释信息
    annotation_col <- data.frame(
      Group = factor(c(rep("CTRL", conNum), rep("TREAT", treatNum)))
    )
    rownames(annotation_col) <- colnames(heatmap_data)
    
    # 定义颜色
    ann_colors <- list(
      Group = c(CTRL = "#4DBBD5", TREAT = "#E64B35")
    )
    
    # 绘制热图
    pdf(file="SVM-RFE.heatmap.pdf", width=8, height=max(6, length(valid_genes)*0.3))
    pheatmap(heatmap_data,
             scale = "row",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             main = paste("SVM-RFE Selected Genes (N=", length(valid_genes), ")", sep=""),
             fontsize_row = 8,
             fontsize_col = 8)
    dev.off()
    cat("\n已生成热图: SVM-RFE.heatmap.pdf\n")
  }
}

# ============= 第十一步：输出总结报告 =============
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("分析完成！\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat("\n生成的文件:\n")
cat("1. normalize.txt - 标准化后的表达数据\n")
cat("2. Control.txt - 对照组样本列表\n")
cat("3. Treat.txt - 处理组样本列表\n")
cat("4. diff.Wilcoxon.txt - 差异表达分析结果\n")
cat("5. SVM-RFE.pdf - SVM-RFE特征选择图\n")
cat("6. SVM-RFE.gene.txt - SVM-RFE筛选的基因列表\n")
cat("7. SVM-RFE.Gene.Exp.txt - 筛选基因的表达谱\n")
if(file.exists("SVM-RFE.heatmap.pdf")) {
  cat("8. SVM-RFE.heatmap.pdf - 筛选基因的热图\n")
}

cat("\n筛选出的基因数量:", length(featureGenes), "\n")
if(length(featureGenes) <= 30) {
  cat("基因列表:", paste(featureGenes, collapse=", "), "\n")
}