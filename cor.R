# 设置工作目录
setwd("/Volumes/LENOVO/Thesis 5/cor")

# 加载包
library(limma)
library(reshape2)
library(ggpubr)
library(ggExtra)
library(readxl)
library(ggplot2)

# ============= 第一步：读取表达数据 =============
exp_data <- read_excel("../LASSO/Exp_details.xlsx")
exp_data <- as.data.frame(exp_data)

# 提取基因表达矩阵
gene_exp <- exp_data[, -c(1,2)]
rownames(gene_exp) <- exp_data$`Feature ID`
mat <- as.matrix(gene_exp)
mode(mat) <- "numeric"

# 查看样本名
cat("样本名:\n")
print(colnames(mat))

# ============= 第二步：提取0.7MAC和2MAC样本 =============
all_samples <- colnames(mat)
X0.7MAC_samples <- all_samples[grepl("0.7MAC", all_samples)]
X2MAC_samples <- all_samples[grepl("2MAC", all_samples)]
treat_samples <- c(X0.7MAC_samples, X2MAC_samples)

# 只保留处理组样本
mat <- mat[, treat_samples, drop=FALSE]

cat("\n处理组样本 (n=6):\n")
print(colnames(mat))

# ============= 第三步：获取Pde10a基因表达量 =============
gene <- "Pde10a"
if(!gene %in% rownames(mat)) {
  stop("找不到Pde10a基因！")
}

gene_exp_data <- t(mat[gene, , drop=FALSE])
gene_exp_data <- as.data.frame(gene_exp_data)
colnames(gene_exp_data) <- gene
cat("\nPde10a表达量:\n")
print(gene_exp_data)

# ============= 第四步：读取免疫细胞浸润分数 =============
# 使用之前ssGSEA计算的结果
immune <- read.table("immScore_simple.txt", header=T, sep="\t", check.names=F, row.names=1)
immune <- t(immune)  # 转置，使样本为行
immune <- as.data.frame(immune)

cat("\n免疫细胞数据维度:", dim(immune), "\n")
print(head(immune))

# ============= 第五步：数据合并 =============
# 确保样本名匹配
sameSample <- intersect(rownames(immune), rownames(gene_exp_data))
cat("\n共同样本数:", length(sameSample), "\n")

if(length(sameSample) == 0) {
  # 如果样本名不匹配，尝试清理空格
  rownames(immune) <- gsub(" ", "", rownames(immune))
  rownames(gene_exp_data) <- gsub(" ", "", rownames(gene_exp_data))
  sameSample <- intersect(rownames(immune), rownames(gene_exp_data))
  cat("清理空格后共同样本数:", length(sameSample), "\n")
}

# 合并数据
rt <- cbind(immune[sameSample, , drop=FALSE], gene_exp_data[sameSample, , drop=FALSE])
cat("\n合并后数据维度:", dim(rt), "\n")
print(head(rt))

# ============= 第六步：相关性分析 =============
outTab <- data.frame()

for(i in colnames(rt)[1:(ncol(rt)-1)]) {
  x <- as.numeric(rt[, gene])
  y <- as.numeric(rt[, i])
  
  # 处理标准差为0的情况
  if(sd(y) == 0) { y[1] <- 0.00001 }
  
  # 计算相关性（使用Spearman）
  cor_test <- cor.test(x, y, method="spearman")
  
  outVector <- cbind(Gene=gene, Cell=i, 
                     cor=cor_test$estimate, 
                     pvalue=cor_test$p.value)
  outTab <- rbind(outTab, outVector)
  
  # 为p<0.05的细胞绘制散点图
  if(cor_test$p.value < 0.05) {
    outFile <- paste0("cor_", gsub("[ /]", "_", i), ".pdf")
    df1 <- as.data.frame(cbind(x, y))
    colnames(df1) <- c("Gene", "Cell")
    
    p1 <- ggplot(df1, aes(x=Gene, y=Cell)) + 
      xlab(paste0(gene, " expression")) + 
      ylab(i) +
      geom_point(size=3, alpha=0.7, color="#2166AC") + 
      geom_smooth(method="lm", formula=y ~ x, color="#B2182B", fill="#FDB366", alpha=0.3) + 
      theme_bw(base_size=14) +
      theme(plot.title=element_text(hjust=0.5)) +
      stat_cor(method='spearman', aes(x=Gene, y=Cell), 
               label.x.npc="center", label.y.npc="top")
    
    p2 <- ggMarginal(p1, type="density", 
                     xparams=list(fill="#FDB366", alpha=0.6), 
                     yparams=list(fill="#2166AC", alpha=0.6))
    
    pdf(file=outFile, width=5.5, height=5)
    print(p2)
    dev.off()
    cat("已生成:", outFile, "\n")
  }
}

# 保存结果
outTab$cor <- as.numeric(outTab$cor)
outTab$pvalue <- as.numeric(outTab$pvalue)
write.table(outTab, file="cor_result.txt", sep="\t", row.names=F, quote=F)
cat("\n相关性结果已保存到 cor_result.txt\n")
print(outTab)

# ============= 第七步：绘制棒棒糖图（Lollipop）=============
# 定义颜色函数
p.col <- c('gold', 'pink', 'orange', 'LimeGreen', 'darkgreen')

fcolor <- function(x, p.col) {
  color <- ifelse(x > 0.8, p.col[1], 
                  ifelse(x > 0.6, p.col[2], 
                         ifelse(x > 0.4, p.col[3],
                                ifelse(x > 0.2, p.col[4], p.col[5]))))
  return(color)
}

# 定义点大小函数
p.cex <- seq(2.5, 5.5, length=5)

fcex <- function(x) {
  x <- abs(x)
  cex <- ifelse(x < 0.1, p.cex[1],
                ifelse(x < 0.2, p.cex[2],
                       ifelse(x < 0.3, p.cex[3],
                              ifelse(x < 0.4, p.cex[4], p.cex[5]))))
  return(cex)
}

# 计算颜色和大小
points.color <- fcolor(x=outTab$pvalue, p.col=p.col)
outTab$points.color <- points.color
points.cex <- fcex(x=outTab$cor)
outTab$points.cex <- points.cex

# 按相关系数排序
outTab <- outTab[order(outTab$cor), ]
xlim <- ceiling(max(abs(outTab$cor)) * 10) / 10

# 绘制棒棒糖图
pdf(file="Lollipop_Pde10a.pdf", width=10, height=8)

layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0), nc=2), 
       widths=c(8, 2.2), heights=c(1, 2, 1, 2, 1))

par(bg="white", las=1, mar=c(5, 20, 2, 4), cex.axis=1.3, cex.lab=1.5)

# 主图
plot(1, type="n", 
     xlim=c(-xlim, xlim), 
     ylim=c(0.5, nrow(outTab)+0.5),
     xlab="Correlation Coefficient", 
     ylab="", 
     yaxt="n", 
     yaxs="i", 
     axes=FALSE)

# 添加背景色
rect(par('usr')[1], par('usr')[3], 
     par('usr')[2], par('usr')[4], 
     col="#F5F5F5", border="#F5F5F5")

# 添加网格线
grid(ny=nrow(outTab), col="white", lty=1, lwd=2)

# 添加连接线
segments(x0=outTab$cor, y0=1:nrow(outTab), 
         x1=0, y1=1:nrow(outTab), 
         lwd=3, col="gray50")

# 添加点
points(x=outTab$cor, y=1:nrow(outTab), 
       col=outTab$points.color, 
       pch=16, 
       cex=outTab$points.cex)

# 添加细胞标签
text(par('usr')[1], 1:nrow(outTab), 
     outTab$Cell, 
     adj=1, xpd=TRUE, cex=1.2)

# 添加p值标签
pvalue.text <- ifelse(outTab$pvalue < 0.001, 
                      '<0.001', 
                      sprintf("%.03f", outTab$pvalue))

redcutoff_cor <- 0
redcutoff_pvalue <- 0.05

text(par('usr')[2], 1:nrow(outTab), 
     pvalue.text, 
     adj=0, xpd=TRUE, 
     col=ifelse(abs(outTab$cor) > redcutoff_cor & outTab$pvalue < redcutoff_pvalue, 
                "red", "black"), 
     cex=1.2)

axis(1, tick=FALSE)

# 图例1：点大小代表相关系数
par(mar=c(0, 4, 3, 4))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("left", 
       legend=c(0.1, 0.2, 0.3, 0.4, 0.5), 
       col="black", 
       pt.cex=p.cex, 
       pch=16, 
       bty="n", 
       cex=1.5, 
       title="abs(cor)")

# 图例2：颜色代表p值
par(mar=c(0, 6, 4, 6), cex.axis=1.5, cex.main=2)
barplot(rep(1, 5), 
        horiz=TRUE, 
        space=0, 
        border=NA, 
        col=p.col, 
        xaxt="n", 
        yaxt="n", 
        xlab="", 
        ylab="", 
        main="pvalue")

axis(4, at=seq(0.5, 4.5, by=1), 
     labels=c(1, 0.8, 0.6, 0.4, 0.2), 
     tick=FALSE)

dev.off()
cat("\n已生成 Lollipop_Pde10a.pdf\n")

# ============= 第八步：输出统计报告 =============
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("Pde10a与免疫细胞相关性分析报告\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("样本信息:\n")
cat("  0.7MAC组: 3个样本\n")
cat("  2MAC组: 3个样本\n")
cat("  总样本数:", nrow(rt), "\n\n")

cat("免疫细胞类型数:", ncol(immune), "\n\n")

cat("相关性统计:\n")
sig_positive <- sum(outTab$cor > 0 & outTab$pvalue < 0.05)
sig_negative <- sum(outTab$cor < 0 & outTab$pvalue < 0.05)
cat("  显著正相关 (p<0.05):", sig_positive, "\n")
cat("  显著负相关 (p<0.05):", sig_negative, "\n")
cat("  总显著相关:", sig_positive + sig_negative, "\n\n")

if(sig_positive + sig_negative > 0) {
  cat("显著相关的免疫细胞:\n")
  print(outTab[outTab$pvalue < 0.05, c("Cell", "cor", "pvalue")])
}

cat("\n生成的文件:\n")
cat("  1. cor_result.txt - 相关性结果表格\n")
cat("  2. cor_*.pdf - 显著相关细胞的散点图\n")
cat("  3. Lollipop_Pde10a.pdf - 棒棒糖图\n")
cat(paste(rep("=", 70), collapse=""), "\n")