library(tidyverse)
library(patchwork)
library(limma)
library(edgeR)
library(DESeq2)
library(pheatmap)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library (VennDiagram) 
#GSE151243_20HS L VS PL----
GSE151243 <- read.table("GSE151243_20HS L VS PL .txt",sep = "\t",row.names = 1,header = T)
#观察行名
colnames(GSE151243)
#提取 'HidradenitisSuppurativa_xL'的行
result_df <- GSE151243[, grep("HidradenitisSuppurativa_\\d+L", 
                              colnames(GSE151243), 
                              value = TRUE)]
colnames(result_df)
# 列名替换
# 使用sub函数进行替换
new_col_names <- sub("HidradenitisSuppurativa",
                     "GSE151243", 
                     colnames(result_df))

# 将替换后的列名应用到数据框
colnames(result_df) <- new_col_names
#将GSE151243皮损写出
write.table(result_df,file = "GSE151243",sep = "\t")

#GSE189266_10HS L VS PL V NL.txt----

GSE189266 <- read.table("GSE189266_10HS L VS PL V NL.txt",sep = "\t",row.names = 1,header = T)
#观察行名
colnames(GSE189266)
#保留带sam的行名
result_df <- GSE189266[, grep("sam", 
                              colnames(GSE189266), 
                              value = TRUE)]
colnames(result_df)
CF <- read.table("GSE189266_series.txt",sep = "\t")
CF <- as.data.frame(t(CF))
colnames(CF) <- CF[1, ]
CF <- CF[-1, ]
data <- CF
data$Time <- sub("timepoint: ", "", data$Time)
data$Skincondition <- sub("subytpe: ", "", data$Skincondition)

# 合并两列内容
data$sample <- paste(data$Time, data$Skincondition, sep = " ")
#data为临床信息
#提取只包含Baseline Lesional的样本


#GSE155176_19 L VS 13N-L(10cm) VS 16 H.csv----
GSE155176 <- read.csv("GSE155176_19 L VS 13N-L(10cm) VS 16 H.csv.gz",sep = ",",row.names = 1,header = T)
colnames(GSE155176)
#保留Lesional_TxNaive|Normal_Skin列
df <- GSE155176[,grep("Lesional_PRE|Lesional_TxNaive|Normal_Skin", colnames(GSE155176), value = TRUE)]
#计算多少个26LS（N7 P19及H（16
grep("Lesional_PRE|Lesional_TxNaive", colnames(GSE155176), value = TRUE)
grep("Normal_Skin", colnames(GSE155176), value = TRUE)
#将GSE155176皮损以及健康写出
write.table(df,file = "GSE155176",sep = "\t")
#GSE154773_22HS L VS 10 H----
GSE154773 <- read.table("GSE154773_22HS L VS 10 H .txt.gz",sep = "\t",row.names = 1,header = T)
GSE154773metadata <- read.table("GSE154773_series_matrix.txt",sep = "\t",row.names = 1,header = F)
GSE154773metadata <- as.data.frame(t(GSE154773metadata))
colnames(GSE154773metadata) <- c("sample","GSM")
GSE154773metadata$sample <- gsub(":", "_", GSE154773metadata$sample)

GSE154773metadata$sample <- gsub("whole blood", "WB", GSE154773metadata$sample)
GSE154773metadata$sample <- gsub("skin", "Skin", GSE154773metadata$sample)
GSE154773 <- GSE154773[,match(GSE154773metadata$sample,colnames(GSE154773))]

identical(GSE154773metadata$sample,colnames(GSE154773))

head(GSE154773metadata)
GSE154773metadata$type <- sapply(strsplit(as.character(GSE154773metadata$sample), "_"), function(x) x[2])
GSE154773metadata$group <- ifelse(grepl("CH", GSE154773metadata$sample), "HC", "HS")
write.table(GSE154773metadata,file = "GSE154773metadata.txt",sep = "\t")
df1 <- GSE154773metadata %>% filter(type == "Skin")
df2 <- GSE154773metadata %>% filter(type == "WB")
df1 <- df1 %>% filter(sample != "HS035_Skin")
df1$sample2 <- df2$sample
write.table(df1,file = "GSE154773metadata2.txt",sep = "\t")


df <- GSE154773[,grep("Skin", colnames(GSE154773), value = TRUE)]
grep("C", colnames(df), value = TRUE)
write.table(df,file = "GSE154773",sep = "\t")

#构建一个metadata含ID Group Batch----
#GSE155176
meta_data <- data.frame(
  ID = colnames(GSE155176),
  Group = ifelse(grepl("Normal_Skin", colnames(GSE155176)), "Normal", "Lesional"),
  Batch = "GSE155176")
#GSE154773
meta_data2 <- data.frame(
  ID = colnames(GSE154773),
  Group = ifelse(grepl("CHS", colnames(GSE154773)), "Normal", "Lesional"),
  Batch = "GSE154773")
#GSE151243
meta_data3 <- data.frame(
  ID = colnames(GSE151243),
  Group = "Lesional",
  Batch = "GSE151243")
metadata <- rbind(meta_data,meta_data2,meta_data3)
write.table(metadata,file = "metadata",sep = "\t")
#GSE154773和GSE155176差异分析----
#读取信息
GSE155176 <- read.table("GSE155176",sep = "\t",row.names = 1,header = T)
GSE154773 <- read.table("GSE154773",sep = "\t",row.names = 1,header = T)
#GSE155176
meta_data <- data.frame(
  ID = colnames(GSE155176),
  Group = ifelse(grepl("Normal_Skin", colnames(GSE155176)), "Normal", "Lesional"),
  Batch = "GSE155176")
meta_dataGSE155176 <- column_to_rownames(meta_dataGSE155176,var = "ID")

#GSE154773
meta_data2 <- data.frame(
  ID = colnames(GSE154773),
  Group = ifelse(grepl("CHS", colnames(GSE154773)), "Normal", "Lesional"),
  Batch = "GSE154773")
meta_dataGSE154773 <- column_to_rownames(meta_data2,var = "ID")

#读取GSE151243、GSE155176、GSE154773进行合并----
GSE151243 <- read.table("GSE151243",sep = "\t",row.names = 1,header = T)
GSE155176 <- read.table("GSE155176",sep = "\t",row.names = 1,header = T)
GSE154773 <- read.table("GSE154773",sep = "\t",row.names = 1,header = T)
metadata <- read.table("metadata",sep = "\t",row.names = 1,header = T)
#分割ENSG ID GSE155176
rt <- rownames_to_column(GSE155176,"id")
id <- apply(rt, 1, function(x){
  strsplit(x, "[.]")[[1]][1]
})
id <- as.data.frame(id)
rt$id <- id$id
#合并GSE155176和GSE151243
p <- rownames_to_column(GSE151243,"id")
merged_df <- merge(rt, p, by = "id")
table(table(merged_df$id)>1)
table(is.na(merged_df$id))
#将ENSG转化为基因symbol
library(org.Hs.eg.db)
GS <- AnnotationDbi::select(org.Hs.eg.db, keys = merged_df$id, 
                      keytype = "ENSEMBL", column = "SYMBOL")
colnames(GS)[1] <- "id"

# 方法1
# 使用 merge 合并数据框1和数据框2
merged_data <- merge(merged_df, GS, by = "id", all = FALSE)
# 删除包含缺失值的行
table(is.na(merged_data$SYMBOL))
merged_data <- na.omit(merged_data)
# 删除重复行
table(table(merged_data$SYMBOL)>1)
merged_data <- merged_data[!duplicated(merged_data$SYMBOL), ]
merged_data <- merged_data[,-1]
p1 <- dplyr::select(merged_data, SYMBOL, everything())

#方法2
#查找 Ensemble ID 使用索引位置替换 Ensemble ID 为 Gene Symbol
merged_df1 <- merged_df
merged_df1$id <- GS$SYMBOL[match(merged_df$id, GS$id)]
table(is.na(merged_df1$id))
table(table(merged_df1$id)>1)
merged_df1 <- na.omit(merged_df1)
merged_df1 <- merged_df1[!duplicated(merged_df1$id), ]

#检验两者为什么不一样
setdiff(merged_data$SYMBOL,merged_df1$id)
identical(merged_df1$id,merged_data$SYMBOL)
p1 <- merged_data[merged_data$SYMBOL == "ST7-OT3", ]

#合并GSE154773
colnames(GSE154773)
p2 <- rownames_to_column(GSE154773,"SYMBOL")
merged_p <- merge(p1, p2, by = "SYMBOL")
merge1 <- column_to_rownames(merged_p,var="SYMBOL")
write.table(merge1,file = "merge",sep = "\t")

#读取合并矩阵----
#进行无效基因的过滤
merge <- read.table("merge",sep = "\t",row.names = 1,header = T)
metadata <- read.table("metadata",sep = "\t",row.names = 1,header = T)
dim(merge)
#表达矩阵预处理
#过滤低表达基因
#1.至少75%的样本不表达
#2.平均值count＜10的基因
#3.评价cpm＜10的基因
keep <- rowSums(merge>0) >= floor(0.75*ncol(merge))
table(keep)
filter <- merge[keep,]


#数据探索----
#PCA分析----
library(DESeq2)
metadata$Group <- factor(metadata$Group)
metadata$Batch <- factor(metadata$Batch)
metadata <- data.frame(metadata, row.names = NULL)
metadata <- column_to_rownames(metadata,var = "ID")
metadata <- metadata[match(colnames(filter),rownames(metadata)),]
identical(colnames(filter),rownames(metadata))
dds <- DESeqDataSetFromMatrix(countData = filter, 
                              colData = metadata, 
                              design = ~ Group + Batch)
# 使用 VST 转换数据
vsd <- vst(dds, blind=FALSE)

# 进行 PCA 分析并返回 PCA 数据
pcaData <- plotPCA(filter, intgroup=c("Group", "Batch"), returnData=TRUE)

# 计算主成分分析的方差解释百分比
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 创建 PCA 图形
P1 <- ggplot(pcaData, aes(PC1, PC2, color = Batch, shape = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

mod <- model.matrix(~factor(metadata$Group))
assay(vsd) <- ComBat(dat = assay(vsd), bath= metadata$Batch,mod = mod)

library(limma)
assay(vsd)<- limma::removeBatchEffect(assay(vsd), vsd$type)

pcaData <- plotPCA(vsd, intgroup=c("Group", "Batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
P2 <- ggplot(pcaData, aes(PC1, PC2, color=Batch, shape = Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


P1+P2

pca_result <- prcomp(pcaData, scale = TRUE)
pca_plot <- fviz_pca_ind(pca_result,
                         geom = "point",
                         pointsize = 3,  # 设置点的大小
                         col.ind = pcaData$Batch,  # 颜色代表 condition
                         shape.ind = pcaData$Group,  # 形状代表 type
                         addEllipses = TRUE,  # 添加椭圆
                         legend.title = "Legend Title")  # 设置图例标题
  
library("FactoMineR")
library("factoextra")
# 准备数据
pca <- PCA(t(filter), graph = FALSE)
# 作图
fviz_pca_ind(pca,
             label = "none", 
             habillage = batch,
             palette = c("#00AFBB", "#E7B800","#FB9A99"),
             addEllipses = TRUE)
pca_plot <- fviz_pca_ind(pca,
                         geom.ind = "point",  # 使用点来表示个体
                         habillage = group_factor,  # 形状代表组信息
                         palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 指定颜色
                         addEllipses = TRUE,  # 添加椭圆
                         ellipse.type = "convex",  # 椭圆类型
                         label = "none"  # 不显示标签
)
pca_plot     
library(tinyarray)
P1 <- draw_pca(exp = filter, group_list = as.factor(metadata$Batch))
P2 <- draw_pca(exp = filter, group_list = as.factor(metadata$Group))
#SVA包 count去批次----
library(sva)

expr_count_combat <- ComBat_seq(counts = as.matrix(filter), 
                                batch = metadata$Batch,
                                group = metadata$Group)
P3 <- draw_pca(exp = expr_count_combat, group_list = as.factor(metadata$Group))
P4 <- draw_pca(exp = expr_count_combat, group_list = as.factor(metadata$Batch))
P1+P2+P4+P3

#三种差异分析----
#GSE155176
#读取信息
GSE155176 <- read.table("GSE155176",sep = "\t",row.names = 1,header = T)
meta_dataGSE155176 <- data.frame(
  ID = colnames(GSE155176),
  Group = ifelse(grepl("Normal_Skin", colnames(GSE155176)), "Normal", "Lesional"),
  Batch = "GSE155176")
meta_dataGSE155176 <- column_to_rownames(meta_dataGSE155176,var = "ID")
#过滤低质量基因
#至少75%的样本不表达
#GSE155176
keep <- rowSums(GSE155176>0) >= floor(0.75*ncol(GSE155176))
table(keep)
filterGSE155176 <- GSE155176[keep,]
identical(colnames(filterGSE155176),rownames(meta_dataGSE155176))
save(list = c("filterGSE155176" ,"meta_dataGSE155176"),file="filterGSE155176.Rda")

#转化信息----
load(file = "filterGSE155176.Rda")
#矩阵
exp <- filterGSE155176
#分组信息
metadata <- meta_dataGSE155176
metadata <- metadata[order(factor(metadata$Group,levels = c("Lesional", "Normal"))),]
exp <- exp[,match(rownames(metadata),colnames(exp))]
identical(colnames(exp),rownames(metadata))
#计算TPM和FPKM
#分割
rt <- rownames_to_column(exp,"id")
id <- apply(rt, 1, function(x){
  strsplit(x, "[.]")[[1]][1]
})
id <- as.data.frame(id)
rt$id <- id$id
rt <- as.data.frame(rt)
rt <- column_to_rownames(rt,"id")
##读取基因长度
eff_length <- read.csv("gene_length_2.csv", row.names = 1, header = T)
rownames(eff_length)<-eff_length$gene_id 
rownames(eff_length) <- do.call(rbind,strsplit(as.character(eff_length$gene_id),'\\.'))[,1]
eff_length[1:3,]
##选取两者的交集基因
gen <- intersect(rownames(rt), rownames(eff_length))
rt <- rt[gen,]
eff_length <- eff_length[gen, ]

##countToFpkm函数
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
##count转换为FPKM值
exp_fpkms <- as.data.frame(apply(rt, 2, countToFpkm, effLen = eff_length$eff_length))
#3. FPKM值转换成TPM值
##FPKM转TPM
fpkmToTpm <- function(exp_fpkms)
{
  exp(log(exp_fpkms) - log(sum(exp_fpkms)) + log(1e6))
}

##计算TPM值
exp_tpms <- as.data.frame(apply(exp_fpkms,2,fpkmToTpm))

head(DEG)
exp1 <- c(t(exp_tpms[match("ENSG00000151694",rownames(exp_tpms)),]))
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

save(list =c( "exp", "exp_fpkms", "exp_tpms", "metadata"),file = "GSE155176.final" )
load(file = "GSE155176.final")

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
write.table(DESeq2_DEG, file='GSE155176 DESeq2_results_diff.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#基因id的转换
##加载基因注释文件
GinfoFile <- "gene_length_Table.txt"
Ginfo_0 <- read.table(GinfoFile,sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] # 只取mRNA信息
rownames(Ginfo) <- Ginfo[,1]

##使用intersect()函数取出两者的共有基因
comgene <- intersect(rownames(DESeq2_DEG),rownames(Ginfo))
##提取注释文件中共有基因的行
DEG <- DESeq2_DEG[comgene,]
Ginfo <- Ginfo[comgene,]
DEG$Gene <- as.character(Ginfo$genename)
#save(DEG, file = paste0(cancer_type,"DEG_symbol.Rdata"))
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
###添加标签
x1 = DEG %>% 
  filter(change == "UP") %>% 
  head(5)
x2 = DEG %>% 
  filter(change == "DOWN") %>% 
  head(5)
label = rbind(x1,x2)
p1
P1 <- p1 +
  geom_point(size = 3, shape = 1, data = label) +
  ggrepel::geom_label_repel(data = label, aes(label = Gene), color="black")
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
         cluster_cols =T,
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
write.table(limma_DEG, file='GSE155176 limma_results_diff.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
##使用intersect()函数取出两者的共有基因
comgene <- intersect(rownames(DEG),rownames(Ginfo))
##提取注释文件中共有基因的行
DEG <- DEG[comgene,]
Ginfo <- Ginfo[comgene,]
DEG$Gene <- as.character(Ginfo$genename)
#检查上下调是否出错
head(DEG)
exp1 <- c(t(exp_tpms[match("ENSG00000163220",rownames(exp_tpms)),]))
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
  ggrepel::geom_label_repel(data = label, aes(label = Gene), color="black")
P3
# 热图
cg = rownames(DEG)[DEG$change !="NORMAL"]
exp_diff <- exp[cg,]
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
write.table(edgeR_DEG, file='GSE155176 edgeR_results_diff.txt',sep = "\t",row.names = T,col.names = NA,quote = F)

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
  ggrepel::geom_label_repel(data = label, aes(label = Gene), color="black")
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
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
P6
head(edgeR_DEG)

#交集的基因与提取----
#上调基因
load(file = "DEG GSE155176")
head(DESeq2_DEG)
filtered_data <- subset(DESeq2_DEG, change == "UP")
DESeq2 <- rownames(filtered_data)
filtered_data <- subset(limma_DEG, change == "UP")
limma <- rownames(filtered_data)
filtered_data <- subset(edgeR_DEG, change == "UP")
edgeR <- rownames(filtered_data)
GSE155176UP <- Reduce(intersect, list(DESeq2,limma,edgeR))
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
             
             filename='GSE155176UP.tiff',# 文件保存
             
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
GSE155176DOWN <- Reduce(intersect, list(DESeq2,limma,edgeR))
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
             
             filename='GSE155176DOWN.tiff',# 文件保存
             
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
GSE155176 <- Reduce(intersect, list(DESeq2,limma,edgeR))
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
             
             filename='GSE155176COM.tiff',# 文件保存
             
             imagetype="tiff",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
)
save(list=c("DESeq2_DEG","limma_DEG","edgeR_DEG","GSE155176","GSE155176DOWN","GSE155176UP"),
     file = "DEG GSE155176")
load(file = "DEG GSE155176")
