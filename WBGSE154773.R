#读取数据----
GSE154773 <- read.table("WBGSE154773",sep = "\t",row.names = 1,header = T)

#构建一个metadata含ID Group Batch
meta_dataGSE154773 <- data.frame(
  ID = colnames(GSE154773),
  Group = ifelse(grepl("CHS", colnames(GSE154773)), "Normal", "Lesional"),
  Batch = "GSE154773")

write.table(meta_dataGSE154773,file = "GSE154773metadata.txt",sep = "\t")

#过滤低质量基因
#至少75%的样本不表达
#GSE154773
keep <- rowSums(GSE154773>0) >= floor(0.75*ncol(GSE154773))
table(keep)
filterGSE154773 <- GSE154773[keep,]
meta_dataGSE154773 <- column_to_rownames(meta_dataGSE154773,var = "ID")
identical(colnames(filterGSE154773),rownames(meta_dataGSE154773))
save(list = c("filterGSE154773" ,"meta_dataGSE154773"),file="filterGSE154773.Rda")

#TPM和FPKM的计算
load(file = "filterGSE154773.Rda")
#矩阵
exp <- filterGSE154773
#分组信息
metadata <- meta_dataGSE154773
metadata <- metadata[order(factor(metadata$Group,levels = c("Lesional", "Normal"))),]
rownames(metadata) <- NULL
metadata <- column_to_rownames(metadata,var = "ID")
exp <- exp[,match(rownames(metadata),colnames(exp))]
identical(colnames(exp),rownames(metadata))
##读取基因长度
eff_length <- read.csv("gene_length_2.csv", row.names = 1, header = T)
rownames(eff_length)<-eff_length$gene_id 
eff_length$e <- do.call(rbind,strsplit(as.character(eff_length$gene_id),'\\.'))[,1]
GS <- AnnotationDbi::select(org.Hs.eg.db, keys =eff_length$e, 
                            keytype = "ENSEMBL", column = "SYMBOL")
GS <- GS[match(eff_length$e,GS$ENSEMBL),]
eff_length <- merge(eff_length,GS,by.x ="e",by.y =  "ENSEMBL",all.x = TRUE)
##选取两者的交集基因
gen <- intersect(rownames(exp), eff_length$SYMBOL)
exp <- exp[gen,]
eff_length <- eff_length[match(gen,eff_length$SYMBOL), ]

#计算TPM和FPKM
##countToFpkm函数
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
##count转换为FPKM值
exp_fpkms <- as.data.frame(apply(exp, 2, countToFpkm, effLen = eff_length$eff_length))
##FPKM值转换成TPM值
fpkmToTpm <- function(exp_fpkms)
{
  exp(log(exp_fpkms) - log(sum(exp_fpkms)) + log(1e6))
}

##计算TPM值
exp_tpms <- as.data.frame(apply(exp_fpkms,2,fpkmToTpm))

#保存四个分析基本数据----

save(list =c( "exp","filterGSE154773", "exp_fpkms", "exp_tpms", "metadata"),file = "GSE154773.final" )
library(stringr)
Group <- str_replace_all(metadata$Group, "Lesional", "HS")


  
#差异分析----
load(file = "GSE154773.final")
exp <- filterGSE154773
#DESeq2----
metadata <- metadata[match(colnames(exp),rownames(metadata)),]
Group <- metadata$Group
#正式分析
dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = metadata,
  design = ~ Group)
dds <- DESeq(dds)
contrast <- c("Group", "Lesional", "Normal")
#前面的为实验组
res <- results(dds, contrast)
#查看每一项结果的具体含义
mcols(res, use.names = TRUE)
#查看描述性结果
summary(res)
resOrdered <- res[order(res$pvalue),]   #按照padj值排序
DEG<- as.data.frame(resOrdered)
#添加change列标记基因上调下调
logFC_cutoff <- 1
type1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NORMAL"))
table(DEG$change)
DESeq2_DEG <- DEG
read.table(file='GSE154773 DESeq2_results_diff.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
DESeq2_DEG <- read.table(file = "GSE154773 DESeq2_results_diff.txt",header = T,sep = "\t")
#火山图
DEG$logP <- -log10(DEG$pvalue)
p1 <- ggplot(data = DEG, 
             aes(x = log2FoldChange, y = logP)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p1
###添加标签
x1 = DEG %>% 
  filter(change == "UP") %>% 
  head(5)
x2 = DEG %>% 
  filter(change == "DOWN") %>% 
  head(5)
label = rbind(x1,x2)
P1 <- p1 +
  geom_point(size = 3, shape = 1, data = label) +
  ggrepel::geom_label_repel(data = label, aes(label = rownames(label)), color="black")
P1


# 热图
cg = rownames(DEG)[DEG$change !="NORMAL"]
exp_diff <- exp[cg,]
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(exp_diff)
P2 <- pheatmap(exp_diff,
               annotation_col=annotation_col,
               scale = "row",
               show_rownames = F,
               show_colnames =F,
               color = colorRampPalette(c("navy", "white", "red"))(50),
               cluster_cols =F,
               fontsize = 10,
               fontsize_row=3,
               fontsize_col=3)
#limma----
##差异分析
design <- model.matrix(~0+factor(metadata$Group))
colnames(design)=levels(factor(metadata$Group))
rownames(design)=colnames(exp)
#设计对比 前面比后面
comp <- "Lesional-Normal"
cont.matrix <- makeContrasts(contrasts=c(comp),levels = design) 
dge <- DGEList(counts=exp)
dge <- calcNormFactors(dge)
#精确权重法（precision weights）VOOM
v <- voom(dge,design, normalize="quantile")
#limma差异分析
fit <- lmFit(v, design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
fit2
#提取差异结果
DEG <- topTable(fit2, coef=comp, n=Inf,adjust.method = "BH")
DEG <- na.omit(DEG)
head(DEG)
#筛选上下调，设定阈值
#fc_cutoff <- 2
#pvalue <- 0.05
#DEG$change <- "NOT"
#up <- intersect(which(DEG$logFC > fc_cutoff),
#which(DEG$adj.P.Val < pvalue))
#down <- intersect(which(DEG$logFC < -fc_cutoff),
#which(DEG$adj.P.Val < pvalue)) 
#DEG$change[up] <- "UP"
#DEG$change[down] <- "DOWN"
#添加change列标记基因上调下调
logFC_cutoff <- 1
type1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NORMAL"))
table(DEG$change)
limma_DEG <- DEG
write.table(limma_DEG, file='GSE154773 limma_results_diff.txt',sep = "\t",row.names = T,col.names = NA,quote = F)

#检查上下调是否出错
head(DEG)
exp1 <- c(t(exp_tpms[match("ASIP",rownames(exp_tpms)),]))
test <- data.frame(value=exp1,group=metadata$Group)
test
ggplot(data=test,aes(x=group,y=value,fill=group))+geom_boxplot()

#火山图
DEG$logP <- -log10(DEG$P.Value)
p3 <- ggplot(data = DEG, 
             aes(x = logFC, y = logP)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p3
###添加标签
x1 = DEG %>% 
  filter(change == "UP") %>% 
  head(5)
x2 = DEG %>% 
  filter(change == "DOWN") %>% 
  head(5)
label = rbind(x1,x2)
P3 <- p3 +
  geom_point(size = 3, shape = 1, data = label) +
  ggrepel::geom_label_repel(data = label, aes(label = rownames(label)), color="black")
P3
# 热图
cg = rownames(DEG)[DEG$change !="NORMAL"]
exp_diff <- exp[cg,]
#这里的问题 应该用TPM值
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(exp_diff)
P4 <- pheatmap(exp_diff,
               annotation_col=annotation_col,
               scale = "row",
               show_rownames = F,
               show_colnames =F,
               color = colorRampPalette(c("navy", "white", "red"))(50),
               cluster_cols =F,
               fontsize = 10,
               fontsize_row=3,
               fontsize_col=3)
#edgeR----
#实验设计矩阵(Design matrix)
design <- model.matrix(~0+factor(metadata$Group))
colnames(design)=levels(factor(metadata$Group))
rownames(design)=colnames(exp)
#构建DGE对象
dge <- DGEList(counts=exp, group=factor(metadata$Group))
dge
#标准化
dge <- calcNormFactors(dge) 
dge$samples$norm.factors

#估计离散值（Dispersion）
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
#差异表达检验
design
#1:-1 表示前面比后面 -1：1表示后面比前面
fit2 <- glmLRT(fit, contrast=c(1,-1))   
#-1对应normal，1对应处理组/疾病组

#提取结果
DEG <- topTags(fit2, n=nrow(exp))
DEG <- as.data.frame(DEG)

#筛选基因保留全部
logFC_cutoff <- 1
type1 = (DEG$PValue< 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$PValue < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NORMAL"))
table(DEG$change)
edgeR_DEG <- DEG
write.table(edgeR_DEG, file='GSE154773 edgeR_results_diff.txt',sep = "\t",row.names = T,col.names = NA,quote = F)

##使用intersect()函数取出两者的共有基因
comgene <- intersect(rownames(DEG),rownames(Ginfo))
##提取注释文件中共有基因的行
DEG <- DEG[comgene,]
Ginfo <- Ginfo[comgene,]
DEG$Gene <- as.character(Ginfo$genename)
#火山图
#利用矫正P值化
DEG$logP <- -log10(DEG$PValue)
p5 <- ggplot(data = DEG, 
             aes(x = logFC, y = logP)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
p5
###添加标签
x1 = DEG %>% 
  filter(change == "UP") %>% 
  head(5)
x2 = DEG %>% 
  filter(change == "DOWN") %>% 
  head(5)
label = rbind(x1,x2)
P5 <- p5 +
  geom_point(size = 3, shape = 1, data = label) +
  ggrepel::geom_label_repel(data = label, aes(label = rownames(label)), color="black")
P5
# 热图
cg = rownames(DEG)[DEG$change !="NORMAL"]
exp_diff <- exp[cg,]
#exp_diff <- exp_diff[,order(factor(metadata$Group,levels = c("Lesional", "Normal")))]
#这里的问题
annotation_col=data.frame(group=Group)
rownames(annotation_col)=colnames(exp_diff)
P6 <- pheatmap(exp_diff,
               annotation_col=annotation_col,
               scale = "row",
               show_rownames = F,
               show_colnames =F,
               color = colorRampPalette(c("navy", "white", "red"))(50),
               cluster_cols =F,
               fontsize = 10,
               fontsize_row=3,
               fontsize_col=3)
P6
head(edgeR_DEG)

#交集的基因与提取----
#上调基因
head(DESeq2_DEG)
#可以用filter函数
filtered_data <- subset(DESeq2_DEG, change == "UP")
DESeq2 <- rownames(filtered_data)
filtered_data <- subset(limma_DEG, change == "UP")
limma <- rownames(filtered_data)
filtered_data <- subset(edgeR_DEG, change == "UP")
edgeR <- rownames(filtered_data)
GSE154773UP <- Reduce(intersect, list(DESeq2,limma,edgeR))
#三元#
venn.diagram(x=list(DESeq2,limma,edgeR),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), # 填充色 配色https://www.58pic.com/
             
             category.names = c("DESeq2", "limma","edgeR") , #标签名
             
             cat.dist = 0.1, # 标签距离圆圈的远近
             
             cat.pos = c(0, 0, 0), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='GSE154773UP.tiff',# 文件保存
             
             imagetype="tiff",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
)

#下调基因
filtered_data <- subset(DESeq2_DEG, change == "DOWN")
DESeq2 <- rownames(filtered_data)
filtered_data <- subset(limma_DEG, change == "DOWN")
limma <- rownames(filtered_data)
filtered_data <- subset(edgeR_DEG, change == "DOWN")
edgeR <- rownames(filtered_data)
GSE154773DOWN <- Reduce(intersect, list(DESeq2,limma,edgeR))
#三元#
venn.diagram(x=list(DESeq2,limma,edgeR),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), # 填充色 配色https://www.58pic.com/
             
             category.names = c("DESeq2", "limma","edgeR") , #标签名
             
             cat.dist = 0.1, # 标签距离圆圈的远近
             
             cat.pos = c(0, 0, 0), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='GSE154773DOWN.tiff',# 文件保存
             
             imagetype="tiff",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
)
#上下调基因
filtered_data <- subset(DESeq2_DEG, change != "NORMAL")
DESeq2 <- rownames(filtered_data)
filtered_data <- subset(limma_DEG, change != "NORMAL")
limma <- rownames(filtered_data)
filtered_data <- subset(edgeR_DEG, change != "NORMAL")
edgeR <- rownames(filtered_data)
GSE154773 <- Reduce(intersect, list(DESeq2,limma,edgeR))
venn.diagram(x=list(DESeq2,limma,edgeR),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), # 填充色 配色https://www.58pic.com/
             
             category.names = c("DESeq2", "limma","edgeR") , #标签名
             
             cat.dist = 0.1, # 标签距离圆圈的远近
             
             cat.pos = c(0, 0, 0), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='GSE154773COM.tiff',# 文件保存
             
             imagetype="tiff",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
)

save(list=c("DESeq2_DEG","limma_DEG","edgeR_DEG","GSE154773","GSE154773DOWN","GSE154773UP"),
     file = "DEG GSE154773")

#富集分析----
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
genelist <- bitr(GSE154773, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
#GO分析
ego <- enrichGO(gene = genelist$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.8, 
                qvalueCutoff =0.8,
                readable = TRUE)
write.table(ego,file="GO.txt",sep="\t",quote=F,row.names = F)
#可视化
## 柱状图
barplot(ego, showCategory = 20,label_format = 100,color = "pvalue")
##气泡图
dotplot(ego, showCategory = 20,label_format = 100)
## 分类展示
barplot(ego, drop = TRUE, showCategory =10,label_format = 100,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(ego,showCategory = 10,label_format =100,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')

#KEGG
ekegg <- enrichKEGG(gene =genelist$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.8,
                    qvalueCutoff =0.8)
write.table(ekegg,file="KEGG.txt",sep="\t",quote=F,row.names = F)
## 柱状图
barplot(ekegg, showCategory = 15,color = "pvalue")
## 气泡图
dotplot(ekegg, showCategory = 10,color = "pvalue")

#机器学习筛选关键变量----
#数据清洗
load(file = "GSE154773.final")
load(file = "DEG GSE154773")
mydata <- exp_tpms[GSE154773,]
mydata <- as.data.frame(t(mydata))
library(stringr)
metadata$Group <- str_replace_all(metadata$Group, "Lesional", "HS")
metadata$Group <- str_replace_all(metadata$Group, "Normal", "C")
mydata$group <- metadata[match(rownames(mydata),rownames(metadata)),"Group"]
mydata <- select(mydata,"group", everything())
table(mydata$group)
#将二分类变量中的字符转变成数字
mydata$group <- ifelse(mydata$group == "HS",1,0) #皮损为1，正常为0

#随机森林筛选变量
library(randomForest)
# 定义结局变量和自变量
y <- mydata[, 1]#提取第一列作为二分类结局变量
x <- mydata[, -1]#其余变量为要筛选的自变量
# 随机森林筛选变量 ----------------------------------------------------------------
rf_model <- randomForest(x, # 是一个包含预测变量的矩阵或数据框，用来训练随机森林模型。
                         y, # 是一个向量，包含因变量的观测值，用来训练随机森林模型。
                         ntree = 5000, #是一个整数，指定生成的决策树的数量，通常设置为一个较大的值（常用1000），以获得更好的预测性能。
                         importance = TRUE)#是否计算变量的重要性。
# 可视化变量重要性排序结果
#MSE
varImpPlot(rf_model, 
           type=1, #选择1或者2
           main="Variable Importance",#标题名称
           n.var = 15,#绘制的top变量个数
           scale = T,#是否显示横坐标
           cex = 1.2#字体大小
           )
#Node
varImpPlot(rf_model, 
           type=2, #选择1或者2
           main="Variable Importance",#标题名称
           n.var = 20,#绘制的top变量个数
           scale = T,#是否显示横坐标
           cex = 1.2#字体大小
)

varImpPlot(rf_model,n.var = 25,cex = 1.2)

# 提取变量重要性
importance <- importance(rf_model)

# ①提取%IncMSE前10的变量
importance <- as.data.frame(importance)
top10_IncMSE <- head(rownames(importance[order(-importance$`%IncMSE`),]), 15)
top10_IncMSE#可以用这个

# 构建筛选获得变量的表达矩阵
top10_IncMSE_data <- mydata[,top10_IncMSE]



# Lasso回归筛选变量 -------------------------------------------------------------
library(glmnet)
set.seed(12345)
lasso_model <- glmnet(x, # 包含自变量的矩阵
                      y, # 因变量向量
                      family = "binomial",# 表示因变量为二元分类变量
                      alpha = 1) # 表示采用L1正则化，即Lasso回归。
print(lasso_model) 
# 结果详解
# Df: 模型的自由度，也就是选择的自变量个数。
# %Dev: 模型的拟合优度，用对数似然比来表示，越接近1说明模型拟合效果越好。
# Lambda: 模型的正则化参数，它是通过交叉验证来选取的，用来控制模型的复杂度，Lambda越大表示模型的正则化强度越大，选择的自变量也越少。

#可视化
plot(lasso_model,# 拟合的Lasso模型
     xvar = "lambda",# 表示将正则化参数lambda作为横坐标
     label = F)# 不在图中显示变量名标签

coef_lasso <- coef(lasso_model,# 拟合的Lasso模型
                   s = 0.071930) # s为Lambda大小，Lambda越大表示模型的正则化强度越大，选择的自变量也越少。
coef_lasso


# 交叉验证选择合适的Lambda
x <- model.matrix(~., data = x) 
cv_model <- cv.glmnet(x, y, family = "binomial",alpha = 1,nfolds = 10)

#可视化
plot(cv_model)

# 根据交叉验证结果，选择lambda值
#我们可以选择平均误差最小的那个λ，即lambda.min
#也可以选择平均误差在一个标准差以内的最大的λ，即lambda.1se。
lambda_min <- cv_model$lambda.min
lambda_min
lambda_1se <- cv_model$lambda.1se
lambda_1se

# 根据lambda值，确定哪些变量应该被保留
coef_cv <- coef(lasso_model, s = lambda_min)
coef_cv 

# 根据回归系数计算OR值
exp(coef_cv)#自然对数的指数函数，它的定义为 exp(x) = e^x，其中 e 是自然对数的基数（约等于 2.718）。


# 把既定Lambda下有意义的变量提取出来
coef_cv <- as.matrix(coef_cv)# 先转为矩阵
coef_cv <- data.frame(coef_cv)# 再转为数据框

coef_cv$OR <- exp(coef_cv$s1)# 计算每一个变量的OR值
nonzero_vars <- rownames(coef_cv[coef_cv$OR != 1, ])# OR值不为1的就是我们要筛选的变量
nonzero_vars <- nonzero_vars [-1]#去除截距项，根据你自己获得的变量数量修改

# 构建筛选获得变量的表达矩阵
lasso_data <- mydata[,nonzero_vars]

RF1 <- intersect(nonzero_vars,top10_IncMSE)

# SVM算法变量筛选 ---------------------------------------------------------------
library(e1071)#支持向量机
library(caret)#集成了上百种分类和回归算法
library(pROC)#ROC曲线可视化
library(ggplot2)#绘图所需


X <- mydata[,-1] # 定义要选择的变量范围
Y <- as.numeric(as.factor(mydata$group)) # 定义预测变量
# ①设置重复交叉验证
control <- trainControl(method = "repeatedcv",   
                        # 交叉验证方法,这里选择重复的k折交叉验证
                        number = 10,     
                        # k值,将数据分成5份,可以根据样本量选择5-10
                        repeats = 5,     
                        # 重复次数,可选2-5次,这里选择1次
                        search = "random"   
                        # 选择重复交叉验证时,每轮选择 folds 的方式,可选"systematic"连续选择与"random"随机选择,这里选择随机     
)


# ②使用 SVM-RFE 特征选择方法
set.seed(12345)   # 为了保持结果的一致性
svm_rfe <- rfe(X,#输入待筛选变量的表达矩阵
               Y,#输入因变量
               #,sizes参数用于指定我们希望选择的特征数量。
               #它决定了RFE算法会考虑选择哪些特征数量的方案。
               #比如sizes = c(1:5, 10)，则RFE算法会考虑选择1个特征,2个特征,3个特征,4个特征,5个特征和10个特征的特征选择方案。
               sizes = 1:20,
               #选择1-10个特征
               rfeControl = rfeControl(functions = caretFuncs,
                                       #使用caret包中的评价指标
                                       method = "repeatedcv",
                                       #使用重复的留一交叉验证
                                       number = 10,
                                       #推荐为10，10折交叉验证
                                       repeats = 5,
                                       #推荐为5，重复5次
                                       verbose = FALSE),
               #不打印详细输出
               method = "svmLinear",
               #使用SVM分类器
               trControl = control,
               #trainControl设置
               preProc = c("center", "scale")
               #特征预处理,标准化
)

# ③提取前10个最相关的变量
svm_rfe_ranking <- svm_rfe$variables
head(svm_rfe_ranking)

# Overall:特征的总体重要性评分,越高表示该特征对模型性能贡献越大,重要性越高。
# var:特征的名称。
# Variables:特征数量,显示在选取该数量特征时,各特征的重要性评分。这里显示在选择21个特征的情况下,各特征的重要性。
# Resample:如果进行了交叉验证,显示交叉验证的fold信息。这里显示为Fold1.Rep1,表示第一折交叉验证的结果。


varImp(svm_rfe)#查看变量重要性评分

varImp_dataframe <- data.frame(Gene = row.names(varImp(svm_rfe))[1:20],
                               importance = varImp(svm_rfe)[1:20, 1])


# 删除NA行
varImp_dataframe <- na.omit(varImp_dataframe)

# 绘制柱状图
mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')

ggplot(varImp_dataframe, aes(x = reorder(Gene, -importance), y = importance , fill = Gene)) + 
  # 将Gene映射到x轴,按importance降序排序  
  geom_col() +
  # 指定geom_col()为柱状图
  ggtitle("Hub Genes") +
  # 设置标题  
  theme(panel.border = element_blank(),# 删除边框
        axis.text.x = element_text(size = 12, color = "black"),# 修改x轴和y轴文本大小和颜色
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),# 移动x轴标题,增加距离
        plot.title = element_text(margin = margin(b = 20)),# 移动图标题,增加距离  
        panel.grid.major = element_line(color = "grey", size = 0.2)) +# 增加灰色网格线
  # 修改x轴名称为"Gene" ，修改y轴名称为"Importance"
  xlab("Gene") + ylab("Importance") +
  scale_fill_manual(values = mycolors)# 使用自定义的mycolors向量填充颜色  


# 查看重要性前10的基因
top_10_vars <- svm_rfe_ranking$var[1:30]
top_10_vars
RS <- intersect(nonzero_vars,top_10_vars)

# 提取前10变量的表达矩阵
top10_SVM_data <- mydata[,top_10_vars]

# ④提取最优子集中的变量
X_plot = svm_rfe$results$Variables
Y_plot = svm_rfe$results$RMSE
plot(X_plot, Y_plot,  
     xlab="Variable Number",  
     ylab="RMSE (Cross-Validation)",  
     col="#7DEE44",    # 点的颜色
     pch=16,           # 点的形状, ici选择实心圆点       
     cex=1.5,          # 点的大小
     lwd=2,            # 线的宽度
     type="b",         # 同时绘制点与线
      )# y轴范围
lines(X_plot, Y_plot, col="#DF294C", lwd=2)   # 额外绘制一条红色粗线

abline(h=min(Y_plot), col="skyblue")   # 添加一条水平参考线
grid(col="grey",lwd=1,lty=3)   # 添加网格线 

legend("topright",c("Training RMSE","Cross-Validation RMSE"),
       col=c("#7DEE44","#DF294C"),pch=c(16,NA),lwd=2,bg="white")   # 添加图例

# 找到RMSE最小的点
# RMSE是Root Mean Square Error(均方根误差)的缩写,是评价回归模型效果的一个很重要指标。
# RMSE通过测量预测值与实际值的离差来评价模型的准确性,值越小表示模型越准确。
wmin <- which.min(svm_rfe$results$RMSE)
wmin

# 在图上标记RMSE最小的点
points(wmin, svm_rfe$results$RMSE[wmin], col = "orange", pch = 16, cex=2)  
text(wmin, svm_rfe$results$RMSE[wmin],  
     paste0("N=", wmin), pos = 2, col = "orange", cex=1)  

# 提取最优子集中变量的名称
Target_Genes <- svm_rfe$optVariables
Target_Genes

# 提取最优子集的表达矩阵
Best_SVM_data <- mydata[,Target_Genes]
#韦恩图-----
#SVM
COM <- Reduce(intersect, list(Target_Genes,nonzero_vars,top10_IncMSE))
# 查看重要性前10的基因
top_10_vars <- svm_rfe_ranking$var[1:20]
top_10_vars
#lasso
nonzero_vars
#RF
top10_IncMSE <- head(rownames(importance[order(-importance$`%IncMSE`),]), 15)
top10_IncMSE#可以用这个


#二元#
library(RColorBrewer)
library (VennDiagram) 
venn.diagram(x=list(nonzero_vars,top10_IncMSE),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFFCC','#CCFFFF'), # 填充色 配色https://www.58pic.com/
             
             category.names = c("LASSO", "RF") , #标签名
             
             cat.dist = 0.02, # 标签距离圆圈的远近
             
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='E:/桌面/HS 数据集/seq/pre deal/whole blood/RF LASSO1.tiff',# 文件保存
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
             compression = "lzw"# 压缩算法
             
)

#三元
venn.diagram(x=list(Target_Genes,nonzero_vars,top10_IncMSE),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFFCC','#CCFFFF',"#FFCCCC"), # 填充色 配色https://www.58pic.com/
             
             category.names = c("SVM", "LASSO","RF") , #标签名
             
             cat.dist = 0.1, # 标签距离圆圈的远近
             
             cat.pos = c(0, 0, 0), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='SVM LASSO RF PRR16 CCR7.tiff',# 文件保存
             
             imagetype="tiff",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
)



#表达量探究-----
exp <- exp_tpms
table(metadata$Group)
metadata$Group <- str_replace_all(metadata$Group, "Lesional", "HS")
metadata$Group <- str_replace_all(metadata$Group, "Normal", "HC")

exp1 <- as.data.frame(t(exp))
metadata <- metadata[match(rownames(exp1),rownames(metadata)),]
identical(rownames(exp1),rownames(metadata))
exp1$group <- metadata$Group
exp1 <- select(exp1,group,everything())
WBexp1 <- exp1
skinexp1 <- exp1

skinexp2 <- skinexp1[!(rownames(skinexp1) == "HS001_Skin"), ]
skinexp2 <- skinexp2[!(rownames(skinexp2) == "HS035_Skin"), ]

WBexp1<- WBexp1[!(rownames(WBexp1) == "HS030_WB"), ]

  

gene <- c("CCR7","PRR16")
table(exp1$group)
exp1$group <- factor(exp1$group,levels = c("HC","HS")) 


#col <-c("#5CB85C","#337AB7","#F0AD4E","#D9534F")
col <-c("#337AB7","#D9534F")
plist2<-list()
for (i in 1:length(gene)){
  bar_tmp<-exp1[,c(gene[i],"group")]#循环提取每个基因表达信息
  colnames(bar_tmp)<-c(gene[i],"group")#统一命名
  my_comparisons1 <- list(c("HC","HS")) #设置比较组
  pb1<-ggboxplot(bar_tmp,#ggboxplot画箱线图
                 x="group",#x轴为组别
                 y=gene[i],#y轴为表达量
                 color="group",#用样本分组填充
                 fill=NULL,
                 add = "jitter",#添加散点
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size=0.01,
                 font.label = list(size=30), 
                 palette = col)+ theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())#坐标轴修饰
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))#横坐标文字设置
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))#标题设置
  pb1<-pb1+theme(legend.position = "NA")+#（因为有组图，横坐标分组了，所以不需要设置legend）
    geom_signif(comparisons = list(c("HC","HS")),
                size = 0.1,
                textsize = 3,
                test = "wilcox.test")
  
  plist2[[i]]<-pb1 #将画好的图储存于plist2列表，并不断赋值循环直到结束
}

plot_grid(plist2[[1]],plist2[[2]],,ncol=2)#ncol=5表示图片排为几列

library(cowplot)
plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
          plist2[[4]],plist2[[5]],plist2[[6]],
          plist2[[7]],plist2[[8]],plist2[[9]],
          plist2[[10]],ncol=5)#ncol=5表示图片排为几列

plist2[[1]]
plist2[[2]]


#ROC曲线----
library(readxl)
library(pROC) 

#数据清洗
Target_Genes1 <- c("group",COM)
Best_data <- mydata[,Target_Genes1]
Best_data$group <- factor(Best_data$group,levels = c(0,1),labels = c("C","HS"))

#计算ROC
ROC_var_SHC2 <- roc(Best_data$group,Best_data$SHC2)
auc(ROC_var_SHC2)#查看曲线下面积AUC
ci(auc(ROC_var_SHC2))#查看95%置信区间

#绘制单个ROC
plot(1-ROC_var_CCR7$specificities,ROC_var_CCR7$sensitivities,type = "l",col="red",lty=1,
     xlab = "1-Specificity",ylab = "Sensitivity",lwd=2)
abline(0,1)
legend(0.38,0.08,#增加图例（坐标位置）
       c("AUC of CCR7:0.9429(95%CI 0.8684-1)"),
       lty = 1,#线段类型
       lwd = 2,#线段宽度
       col = "red",
       bty = "0")

#循环计算各变量AUC
#存储每个变量的 AUC 和 CI
auc_values <- vector("list", length = ncol(Best_data) - 1)
ci_values <- vector("list", length = ncol(Best_data) - 1)

#进行循环
for (i in 2:ncol(Best_data)) {  # 从第二列开始，因为第一列是 group
  variable <- colnames(Best_data)[i]
  
  # 计算 ROC
  roc_result <- roc(Best_data$group, Best_data[, i], levels = c("C", "HS"))
  
  # 存储 AUC 和 CI
  auc_values[[i - 1]] <- auc(roc_result)
  ci_values[[i - 1]] <- ci(auc(roc_result))
  # 绘制 ROC 曲线
  plot(1 - roc_result$specificities, roc_result$sensitivities, type = "l", col = "red", lty = 1,
       xlab = "1-Specificity", ylab = "Sensitivity", lwd = 2,
       main = paste("ROC Curve for", variable))  # 添加变量名到标题
  abline(0, 1)
  
  # 添加图例
  legend(0.3, 0.08, c(paste("AUC of", variable, ":", round(auc(roc_result), 4),
                             " (95%CI", round(ci(roc_result)[1], 4), "-", round(ci(roc_result)[3], 4), ")")),
         lty = 1, lwd = 2, col = "red", bty = "0")
}

#多个ROC汇总
#计算ROC
colnames(Best_data)
#PRR16
ROC_var_PRR16 <- roc(Best_data$group,Best_data$PRR16)
auc(ROC_var_PRR16)#查看曲线下面积AUC
ci(auc(ROC_var_PRR16))#查看95%置信区间
#CCR7
ROC_var_CCR7 <- roc(Best_data$group,Best_data$CCR7)
auc(ROC_var_CCR7)#查看曲线下面积AUC
ci(auc(ROC_var_CCR7))#查看95%置信区间

#画图
plot(1-ROC_var_PRR16$specificities,ROC_var_PRR16$sensitivities,type = "l",col="#FFCCCC",lty=1,
     xlab = "1-Specificity",ylab = "Sensitivity",lwd=2)
lines(1-ROC_var_CCR7$specificities,ROC_var_CCR7$sensitivities,type = "l",col="#54EDE6",lty=1,lwd=2)
abline(0,1)
legend(0.23,0.13,
       c("AUC of ROC_var_PRR16:0.9429(95%CI 0.8684-1)",
         "AUC of ROC_var_CCR7:0.9095(95%CI 0.7991-1)"),
       lty = c(1,1,1),
       lwd = c(2,2,2),
       col = c('#CCFFFF',"#FFCCCC"),
       bty = "0")

#Nomogram-logistic----
library(VIM)#缺失值可视化
library(rms)#拟合模型
library(nomogramFormula)#计算列线图得分
library(pROC)#绘制ROC曲线，计算AUC和95%置信区间
library(rmda)#临床决策曲线和临床影响曲线

#检查你的输入数据
aggr(Best_data, col=c('skyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(mydata), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Percentage"))
summary(Best_data) 
table(Best_data$group)
Best_data$group <- factor(Best_data$group,levels = c(0,1),labels = c("C","HS"))

#开始构建logistic回归模型
#数据打包 
dd = datadist(Best_data)
option <- options(datadist = "dd")

#拟合模型
colnames(Best_data)#查看列名，选择你要构建模型的变量
paste(colnames(Best_data), collapse = "+")
head(Best_data)
formula <- as.formula(group ~ PRR16+CCR7+SHC2 )
#在R中，变量名称中包含特殊字符（如破折号）时，需要使用反引号（`）将变量名括起来，以确保R正确识别变量名。
model <- lrm(formula, # 回归模型的公式,指定自变量和因变量
             data = Best_data, # 包含所有变量的数据框
             x=TRUE, # logistic回归也称为"广义线性模型",这个参数指定响应变量的二分类
             y=TRUE) # 参数y也是指定因变量是二分类的
model#查看模型具体情况
#计算Logistic回归模型中每个自变量的比值比(Odds Ratio, OR)
OR <- exp(model$coefficients)
OR

#可视化
Nomogram_1 <- nomogram(model,
                       fun = function(x)1/(1+exp(-x)),
                       lp=F,
                       fun.at = c(0.1,0.3,0.5,0.7,0.9),
                       funlabel = "Risk")
plot(Nomogram_1)


# 模型评估
#AUC计算
model_ROC <- glm(formula,data = mydata,family = binomial())
#lrm()函数属于Design包,专门用于拟合Logistic回归模型。
#glm()函数属于stats包,是用于拟合广义线性模型(Generalized Linear Models)的泛用函数,可以拟合Logistic回归、Poisson回归等多种模型。

#type=“response”给出具体的预测概率，而type=“class”按规定的阈值给出分类
Best_data$predvalue <- predict(model_ROC,type="response")
# 计算AUC和95%CI
ROC <- roc(Best_data$group,Best_data$predvalue) 
auc(ROC)
ci(auc(ROC))

