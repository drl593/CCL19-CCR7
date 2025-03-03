load(file = "merge.final")




#免疫浸润CIBERSORT----
##导出文件之前需要把行名变成一列,不然后面就会有报错（本质是因为CIBERSORT函数对于文件读取比较单一）
write.table(cbind(rownames(exp_tpms), exp_tpms),"exp_tpms.txt",quote = F, sep = "\t", row.names=FALSE)
source("source.R")   #注释文件
library(e1071)
library(parallel)
library(preprocessCore)
#CIBERSORT计算
sig_matrix <- "LM22.txt"   #注释文件名
mixture_file = 'exp_tpms normalize.txt'   #表达数据文件名，需要修改
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=1000, QN=FALSE)

save(list = c("CIBERSORT", "CoreAlg", "doPerm", "exp", "exp_fpkms", "exp_tpms", 
              "metadata", "mixture_file", "res_cibersort", "sig_matrix"),file = "123")


boxplot(exp_tpms,outline=F, notch=F , las=2)
qx <- as.numeric(quantile(exp_tpms, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))  #数据的分布，样本分位数
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)#判断是否进行log的标准
LogC


# 均数不一致需标准化 
library(limma)
DEG_expr <- normalizeBetweenArrays(exp_tpms)
boxplot(DEG_expr,outline=FALSE, notch=F , las=2)
write.table(cbind(rownames(DEG_expr), DEG_expr),"exp_tpms normalize.txt",quote = F, sep = "\t", row.names=FALSE)

res_cibersort2 <- CIBERSORT(sig_matrix, mixture_file, perm=1000, QN=FALSE)
result <- as.data.frame(res_cibersort)
result2 <- as.data.frame(res_cibersort2)

#metadata预先处理
colnames(metadata) <- c("group","Batch")
table(metadata$group)
metadata$group <- gsub("HC", "HC", metadata$group)
metadata$group <- gsub("Lesional", "LS", metadata$group)
metadata <- metadata[order(factor(metadata$group,levels = c("HC","LS"))),]
result <- result[match(rownames(metadata),rownames(result)),]
result2 <- result2[match(rownames(metadata),rownames(result)),]

#开始绘图
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(reshape2)
library(tidyverse)

mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
                       '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
                       '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
                       '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
                       '#A7DCE2','#AFDE9C')

result1 <- result[,1:22]
pheatmap(result1,
         color = colorRampPalette(c("#4CB43C", "#FEFCFB", "#ED5467"))(100),
         border="skyblue",#边框颜色
         main = "Heatmap",#指定图表的标题
         show_rownames = T,#是否展示行名
         show_colnames = T,#是否展示列名
         cexCol = 1,#指定列标签的缩放比例。
         scale = 'row',#指定是否应按行方向或列方向居中和缩放，或不居中和缩放。对应的值为row, column和none。
         cluster_col=T,#分别指定是否按列和行聚类。
         cluster_row=F,
         angle_col = "45",#指定列标签的角度。
         legend = F,#指定是否显示图例。
         legend_breaks=c(-3,0,3),#指定图例中显示的数据范围为-3到3。
         fontsize_row = 10,#分别指定行标签和列标签的字体大小。
         fontsize_col = 10)

data <- cbind(rownames(result1),result1)
colnames(data)[1] <- "Samples"
data <- melt(data,id.vars = c("Samples"))
colnames(data) <- c('Samples','celltype','proportion')
ggplot(data,
       aes(Samples,proportion,fill=celltype))+geom_bar(stat="identity",position="fill")+#x 轴是变量 Samples，y 轴是变量 proportion，条形的填充颜色由变量 celltype 决定
       scale_fill_manual(values=mycolors)+#填入需要填充的颜色，22种免疫细胞
       ggtitle("Proportion of immune cells")+theme_gray()+theme(axis.ticks.length=unit(3,'mm'),axis.title.x=element_text(size=11))+
       theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+
       guides(fill=guide_legend(title="Types of immune cells"))

#箱式图
data2 <- cbind(result2,metadata)
colnames(data2)
data2 <- data2[,c(26,27,1:22)]
data2 <- rownames_to_column(data2)
data2 <- pivot_longer(data = data2,
                      cols = 4:25,
                      names_to = "celltype",
                      values_to = "proportion")
ggboxplot(data = data2,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",#绘制箱形图的响应变量。
          combine = TRUE,#是否将数据合并为一个箱形图。
          merge = FALSE,#是否将相同值的分组合并。
          color = "black",#箱形图边框的颜色。
          fill = "group",#箱形图填充色。
          palette = NULL,#颜色调色板。
          title = "TME Cell composition",#图形标题。
          xlab = NULL,#x 轴标签。
          ylab = "Cell composition",#y 轴标签
          bxp.errorbar = FALSE,#是否在箱形图中绘制误差条。
          bxp.errorbar.width = 0.2,#误差条宽度。
          facet.by = NULL,#基于哪些变量进行分面
          panel.labs = NULL,#分面的标签
          short.panel.labs = TRUE,#是否将分面标签缩短
          linetype = "solid",#线条类型
          size = NULL,#图形大小。
          width = 0.8,#箱形图的宽度。
          notch = FALSE,#是否在箱形图中绘制刻度。
          outlier.shape = 20,#异常值标记的形状。
          select = NULL,#要绘制的变量
          remove = NULL,#不要绘制的变量。
          order = NULL,#箱形图的排序方式。
          error.plot = "pointrange",#如何绘制误差，可以是 "pointrange" 或 "errorbar"。
          label = NULL,#要添加的标签
          font.label = list(size = 12, color = "black"),#标签的字体属性
          label.select = NULL,#要添加标签的数据点
          repel = FALSE,#是否使用 repel 库的功能使标签互不重叠
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))  #这个函数的作用是将 x 轴文本旋转 90 度，并调整其对齐方式。



ggboxplot(data = data2,
          x = "celltype",#箱形图中的分组变量。
          y = "proportion",#绘制箱形图的响应变量。
          combine = TRUE,#是否将数据合并为一个箱形图。
          merge = FALSE,#是否将相同值的分组合并。
          color = "group",#箱形图边框的颜色。
          fill = NULL,#箱形图填充色。
          palette = c("#337AB7","#D9534F"),#颜色调色板。
          title = "IME Cell composition",#图形标题。
          xlab = NULL,#x 轴标签。
          ylab = "Cell composition",#y 轴标签
          bxp.errorbar = FALSE,#是否在箱形图中绘制误差条。
          bxp.errorbar.width = 0.2,#误差条宽度。
          facet.by = NULL,#基于哪些变量进行分面
          panel.labs = NULL,#分面的标签
          short.panel.labs = TRUE,#是否将分面标签缩短
          linetype = "solid",#线条类型
          size = NULL,#图形大小。
          width = 0.8,#箱形图的宽度。
          notch = FALSE,#是否在箱形图中绘制刻度。
          outlier.shape = 20,#异常值标记的形状。
          select = NULL,#要绘制的变量
          remove = NULL,#不要绘制的变量。
          order = NULL,#箱形图的排序方式。
          error.plot = "pointrange",#如何绘制误差，可以是 "pointrange" 或 "errorbar"。
          label = NULL,#要添加的标签
          font.label = list(size = 12, color = "black"),#标签的字体属性
          label.select = NULL,#要添加标签的数据点
          repel = TRUE,#是否使用 repel 库的功能使标签互不重叠
          label.rectangle = TRUE, ggtheme = theme_pubr())+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) +  #这个函数的作用是将 x 轴文本旋转 90 度，并调整其对齐方式。
  stat_compare_means(label = "p.signif",method = "t.test",ref.group = ".all.",hide.ns = T) 


#基因相关性热图的绘制
library(ggcorrplot)
library(tidyr)

tdata <- t(DEG_expr)#数据转置，行名是样本名，列名基因名
identical(rownames(metadata),rownames(result21))
identical(rownames(metadata),rownames(tdata))

result21 <- result2[metadata$group == "LS",1:22]
data <- tdata[metadata$group == "LS",c("AKR1B10", "IGFL2", "TMPRSS6", "SLAMF7", "PLEKHG7", "WNK2")]
corr_mat <- cor(data, result21, method = "spearman")
corr_mat <- t(corr_mat)

ggcorrplot(corr_mat, 
           # 显示颜色图例
           show.legend = T, 
           # 设置相关系数矩阵的颜色，其中第一个为最小值对应的颜色，
           # 中间的白色为0对应的颜色，最后一个为最大值对应的颜色
           colors = c("#2166AC", "white", "#B2182B"), 
           # 设置数字显示的位数
           digits = 2, 
           # 显示变量标签（默认为TRUE）
           lab = T) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#修改横轴标签方向

#免疫细胞相关性图
M <- round(cor(result21),2) # 计算相关性矩阵并保留两位小数

pheatmap::pheatmap(M)
ggcorrplot(M, 
           # 显示颜色图例
           show.legend = T, 
           # 设置相关系数矩阵的颜色，其中第一个为最小值对应的颜色，
           # 中间的白色为0对应的颜色，最后一个为最大值对应的颜色
           colors = c("#2166AC", "white", "#B2182B"), 
           # 设置数字显示的位数
           digits = 2, 
           # 显示变量标签（默认为TRUE）
           lab = T) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))#修改横轴标签方向

#患者barplot图
mycol <- ggplot2::alpha(rainbow(ncol(result21)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(result21)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(result21)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-8, # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(result21), 
       xpd = T,
       fill = mycol,
       cex = 0.4, 
       border = NA, 
       y.intersp = 0.6,
       x.intersp = 0.1,
       bty = "n")

#棒棒糖图的绘制
identical(rownames(data),rownames(result21))

#封装函数spearman相关性分析
calculate_correlation <- function(Gene_expr, cibersort_result) {
  cor_matrix <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
  
  for (i in 1:ncol(cibersort_result)) {
    result <- cor.test(Gene_expr, cibersort_result[, i], method = "spearman")
    new_row <- data.frame("Gene" = "Gene", "im_cell" = colnames(cibersort_result)[i], "Cor" = result$estimate, "p-value" = result$p.value)
    cor_matrix <- rbind(cor_matrix, new_row)
  }
  return(cor_matrix)
}

# 新建一个空的数据框保存结果
results <- data.frame("Gene" = character(),"im_cell" = character(),"Cor" = numeric(),"p-value" = numeric(), stringsAsFactors = FALSE)
# 使用for循环遍历数据框中的每一列，并计算相关性
for (i in 1:ncol(data)) {
  print(i)
  gene_expr <- data[, i]
  corr_result <- calculate_correlation(gene_expr, result21)
  corr_result$Gene<- colnames(data)[i]  # 获取当前列的列名
  # 将每次计算的结果添加到新的数据框中
  results <- rbind(results,corr_result)}

colnames(results)


#循环绘图
plist2<-list()

unique_genes <- unique(results$Gene)

for (i in 1:length(unique_genes)){
  
  # 选择目标基因的数据
  mydata <- results[results$Gene == unique_genes[i], ]
  
  # 数据处理
  mydata$im_cell <- reorder(mydata$im_cell, mydata$Cor)
  
  # 绘制相关性图
  plot <- ggplot(mydata, aes(x = Cor, y = im_cell)) +
    geom_segment(aes(x = 0, xend = Cor, y = im_cell, yend = im_cell), color = "black") +
    geom_point(aes(size = abs(Cor), colour = p.value), alpha = 0.2) +
    geom_text(aes(x = max(abs(Cor)) + 0.05, y = im_cell, label = sprintf("%.3f", p.value)), 
              size = 3, hjust = 0, color = ifelse(mydata$p.value < 0.05, "red", "black"))+
    scale_colour_gradient(low ="#ED063F", high ="#2166AC" ) +
    scale_size_continuous(range = c(2, 6)) +
    #scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5),
          plot.title = element_text(size = 14, face = "bold",hjust = 0.5),
          axis.text.y = element_text(hjust = 1)) +
    labs(title =unique_genes[i],x="Correlation Coefficient",y=NULL)
  
  
  plist2[[i]]<-plot
}

lapply(1:6, function(i) plist2[[i]])




 
