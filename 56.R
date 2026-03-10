
#install.packages("glmnet")

library(glmnet)
library(limma)

#设置随机数种子，固定结果
set.seed(123)
#目录
setwd("F:/3SP/56LASSO")

#读取输入文件
data=read.table("GSE30219.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#去除低表达的基因
#data=data[rowMeans(data)>1,]
#boxplot(data.frame(data),col="#4DBBD5")
#标准化
data=normalizeBetweenArrays(data)
#boxplot(data.frame(data),col="#4DBBD5")
#dev.off()

#导出
write.table(data.frame(ID=rownames(data),data),file="normalize.txt", sep="\t", quote=F, row.names = F)

#读入Control样本
Control=read.table("Control.txt", header=F, sep="\t", check.names=F)

#读入Treat样本
Treat=read.table("Treat.txt", header=F, sep="\t", check.names=F)

#设置分组信息
conNum=length(rownames(Control))
treatNum=length(rownames(Treat))
Type=c(rep(1,conNum), rep(2,treatNum))

#按照正常，治疗排序
data1 = data[,Control[,1]]
data2 = data[,Treat[,1]]
data = cbind(data1,data2)

#此处，可以提取自己感兴趣的基因
#比如，某些特定基因集，差异表达分析，单因素cox分析的结果
#以这个为例子，大家需要改成自己的基因
data=data[rowMeans(data)>12,]
data=data[c(1:20),]

#构建模型
x=as.matrix(t(data))
y=Type
fit=glmnet(x, y, family = "binomial")
cvfit=cv.glmnet(x, y, family="binomial",nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#输出筛选的特征基因
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

#导出LASSO中基因的表达谱
rt123 = data[lassoGene,]
write.table(data.frame(ID=rownames(data[lassoGene,]),data[lassoGene,]), file="LASSO.gene.exp.txt", sep="\t", quote=F, row.names=F, col.names=T)

