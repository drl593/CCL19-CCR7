#GSE154773.final
load(file = "GSE154773.final")
exp1 <- as.data.frame(t(exp_tpms))
exp1 <-rownames_to_column(exp1,var = "sample")
metadata <- rownames_to_column(metadata,var = "sample")
exp1 <- merge(exp1,metadata,by="sample")
gene <- c("CXCL17", "SLAMF7", "TCN1", "ADGRF1", "TMPRSS11D", "S100A9", 
          "S100A8", "AKR1B10")
exp1 <- exp1[,c("sample","Group","MYEOV", "TMPRSS6", "TAC1", "CXCL17", "CCER2", "SLAMF7", "TCN1", 
                "PI3", "WNK2", "ADGRF1", "TMPRSS11D", "IGFL2", "S100A7A", "S100A9", 
                "S100A8", "PLEKHG7", "AKR1B10")]
colnames(exp1)[2] <- "group"
exp1$group <- factor(exp1$group,levels = c("Normal","Lesional")) 


col <-c("#337AB7","#D9534F")
plist2<-list()
for (i in 1:length(gene)){
  bar_tmp<-exp1[,c(gene[i],"group")]#循环提取每个基因表达信息
  colnames(bar_tmp)<-c(gene[i],"group")#统一命名
  my_comparisons1 <- list(c("Normal","Lesional")) #设置比较组
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
    geom_signif(comparisons = list(c("Normal", "Lesional")),
                size = 0.1,
                textsize = 3,
                test = "wilcox.test")
  
  plist2[[i]]<-pb1 #将画好的图储存于plist2列表，并不断赋值循环直到结束
}
plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
          plist2[[4]],plist2[[5]],plist2[[6]],
          plist2[[7]],plist2[[8]],#plist2[[9]],
          #plist2[[10]],plist2[[11]],plist2[[12]],
          #plist2[[13]],plist2[[14]],plist2[[15]],
          #plist2[[16]],plist2[[17]],
          ncol=4)#ncol=5表示图片排为几列
