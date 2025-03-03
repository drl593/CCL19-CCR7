load(file = "GSE154773.final")

library(tidyverse)

identical(colnames(exp_tpms),rownames(metadata))
table(metadata$Group)
#保证两个数据框顺序
metadata$group <- metadata$Group
metadata$group <- gsub("Normal", "HC", metadata$group)
metadata$group <- gsub("Lesional", "LS", metadata$group)
  
  
#表达量图
exp1 <- as.data.frame(t(exp))
identical(rownames(exp1),rownames(metadata))
exp1$group <- metadata$group
exp1 <- exp1[order(factor(exp1$group,levels = c("HC","LS"))),]
exp1 <- select(exp1,group,everything())

genelist <-c("CCR7")
table(exp1$group)

gene <- intersect(genelist,colnames(exp1))

#画图
col <-c("#337AB7","#D9534F")
plist2<-list()
library(ggpubr)
library(ggsignif)
library(cowplot)
for (i in 1:length(gene)){
  bar_tmp<-exp1[,c(gene[i],"group")]#循环提取每个基因表达信息
  colnames(bar_tmp)<-c(gene[i],"group")#统一命名
  my_comparisons1 <- list(c("HC","LS")) #设置比较组
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
    geom_signif(comparisons = list(c("HC","LS")),
                size = 0.1,
                textsize = 3,
                test = "wilcox.test")
  
  plist2[[i]]<-pb1 #将画好的图储存于plist2列表，并不断赋值循环直到结束
}

#整合

plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
          plist2[[4]],plist2[[5]],plist2[[6]],
          ncol=3)#ncol=5表示图片排为几列 
