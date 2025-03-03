#读取文件
GSE154773 <- read.table("GSE154773",sep = "\t",row.names = 1,header = T)
meta_dataGSE154773 <- data.frame(
  ID = colnames(GSE154773),
  Group = ifelse(grepl("CHS", colnames(GSE154773)), "Normal", "Lesional"),
  Batch = "GSE154773")
meta_dataGSE154773 <- column_to_rownames(meta_dataGSE154773,var = "ID")
#过滤低质量基因
#至少75%的样本不表达
#GSE154773
keep <- rowSums(GSE154773>0) >= floor(0.75*ncol(GSE154773))
table(keep)
filterGSE154773 <- GSE154773[keep,]
identical(colnames(filterGSE154773),rownames(meta_dataGSE154773))
save(list = c("filterGSE154773" ,"meta_dataGSE154773"),file="filterGSE154773.Rda")
#矩阵
exp <- filterGSE154773
#分组信息
metadata <- meta_dataGSE154773
metadata <- metadata[order(factor(metadata$Group,levels = c("Lesional", "Normal"))),]
exp <- exp[,match(rownames(metadata),colnames(exp))]
identical(colnames(exp),rownames(metadata))
head(exp)
table(metadata$Group)

##读取基因长度
eff_length <- read.csv("gene_length_2.csv", row.names = 1, header = T)
rownames(eff_length)<-eff_length$gene_id 
rownames(eff_length) <- do.call(rbind,strsplit(as.character(eff_length$gene_id),'\\.'))[,1]
eff_length[1:3,]
#转化呈ensembleID
eff_length <- rownames_to_column(eff_length,var = "id")
library(org.Hs.eg.db)
GS <- select(org.Hs.eg.db, keys = eff_length$id, 
             keytype = "ENSEMBL", column = "SYMBOL")
colnames(GS)[1] <- "id"
eff_length<- merge(eff_length,GS,by="id")
##选取两者的交集基因
gen <- intersect(rownames(exp), eff_length$SYMBOL)
rt <- exp[gen,]
eff_length <- eff_length[match(gen,eff_length$SYMBOL), ]

##countToFpkm函数
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
##count转换为FPKM值
exp_fpkms <- as.data.frame(apply(rt, 2, countToFpkm, effLen = eff_length$eff_length))
write.table(exp_fpkms, "GSE154773exp_fpkms.txt", sep="\t", quote=F, row.names=T)

#3. FPKM值转换成TPM值
##FPKM转TPM
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
##计算TPM值
exp_tpms <- as.data.frame(apply(exp_fpkms,2,fpkmToTpm))
write.table(exp_tpms, "GSE154773exp_tpms.txt", sep="\t", quote=F, row.names=T)


head(DEG)
exp1 <- c(t(exp_tpms[match("ADAM10",rownames(exp_tpms)),]))
test <- data.frame(value=exp1,group=metadata$Group)
test
ggplot(data = test, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#85B22E", "#E29827")) +  # 使用 scale_fill_manual 设置箱线图的填充颜色
  
  geom_signif(comparisons = list(c("Normal", "Lesional")),
              size = 1,
              textsize = 4,
              test = "t.test")+
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )

 ##画图

gene <- c("ENSG00000137845","ENSG00000151694")
gene <- as.vector(gene)
#加载样本信息
Exp_plot<- metadata
Exp_plot$sam=metadata$Group
Exp_plot$sam <- factor(Exp_plot$sam,levels=c("Normal","Normal"))
col <-c("#5CB85C","#337AB7")
#保存四个分析基本数据----

save(list =c( "exp", "exp_fpkms", "exp_tpms", "metadata"),file = "GSE154773.final" )
load(file = "GSE154773.final")

#DESeq2----
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
resOrdered <- res[order(res$padj),]   #按照padj值排序
DEG<- as.data.frame(resOrdered)
#添加change列标记基因上调下调
logFC_cutoff <- 2
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NORMAL"))
table(DEG$change)
DESeq2_DEG <- DEG
write.table(DESeq2_DEG, file='GSE154773 DESeq2_results_diff.txt',sep = "\t",row.names = T,col.names = NA,quote = F)

#火山图
DEG$logP <- -log10(DEG$padj)
p1 <- ggplot(data = DEG, 
             aes(x = log2FoldChange, y = logP)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
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
annotation_col=data.frame(group=metadata$Group)
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
logFC_cutoff <- 2
type1 = (DEG$adj.P.Val < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$adj.P.Val < 0.05)&(DEG$logFC > logFC_cutoff)
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
DEG$logP <- -log10(DEG$adj.P.Val)
p3 <- ggplot(data = DEG, 
             aes(x = logFC, y = logP)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
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
exp_diff <- exp_tpms[cg,]
exp_diff <- na.omit(exp_diff)
#这里的问题 应该用TPM值
annotation_col=data.frame(group=metadata$Group)
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
logFC_cutoff <- 2
type1 = (DEG$FDR < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$FDR < 0.05)&(DEG$logFC > logFC_cutoff)
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
DEG$logP <- -log10(DEG$FDR)
p5 <- ggplot(data = DEG, 
             aes(x = logFC, y = logP)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
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
annotation_col=data.frame(group=metadata$Group)
rownames(annotation_col)=colnames(exp_diff)
P6 <- pheatmap(exp_diff,
               annotation_col=annotation_col,
               scale = "row",
               show_rownames = F,
               show_colnames =F,
               color = colorRampPalette(c("navy", "white", "red"))(50),
               cluster_cols =T,
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

GS1 <- select(org.Hs.eg.db, keys = GSE155176, 
              keytype = "ENSEMBL", column = "SYMBOL")
GS1 <- na.omit(GS1)
table(duplicated(GS1$SYMBOL))
GS2 <- GS1[!duplicated(GS1$SYMBOL),]
#二元#

venn.diagram(x=list(GSE154773,GS2$SYMBOL),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.5, #透明度
             
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#FFFFCC','#CCFFFF'), # 填充色 配色https://www.58pic.com/
             
             category.names = c("GSE154773", "GSE155176") , #标签名
             
             cat.dist = 0.02, # 标签距离圆圈的远近
             
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename='E:/桌面/HS 数据集/seq/pre deal/双GSE.tiff',# 文件保存
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
             compression = "lzw"# 压缩算法
             
)

load(file = "DEG GSE155176")

load(file = "DEG GSE154773")

COM <- Reduce(intersect, list(GS2$SYMBOL,GSE154773))
up <-intersect(COM,GSE154773UP)
down <-intersect(COM,GSE154773DOWN)
