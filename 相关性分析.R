com <- skinexp2[,c("CCL19","CCL21")]

com$WB <- rownames(WBexp1)
WBexp1$group
WBexp1$group <- factor(WBexp1$group,levels = c("HS","HHC")) 
WBexp1 <- WBexp1[order(WBexp1$group),]
com$WB <- rownames(WBexp1)
com$group <- skinexp2$group
com$CCR7 <- WBexp1[,"CCR7"]
com$SKIN_CCR7 <- skinexp2[,"CCR7"]
com2 <- com[!(com$group=="HHC"),]


#基因相关性热图的绘制
library(ggcorrplot)
library(tidyr)



data <- com2[,c("CCL19", "CCR7")]


ggplot(com2,aes_string(x = "SKIN_CCR7",y = 'CCR7'))+
  geom_point(size = 5,color = '#2570A4',alpha = 3)+
  theme_bw()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 12),axis.ticks.length = unit(0.2,'cm'),axis.ticks = element_line(size = 1),panel.border = element_rect(size = 1.5),panel.grid = element_blank())+
  geom_smooth(method = 'lm',se = T,color = 'black',size = 2.0,fill = '#7A8991')+ stat_cor(method = "spearman",digits = 3,size=6)


result <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,header = T)
result <- result[!grepl("CH", rownames(result)), ]
result <- result[!(rownames(result) == "HS030_WB"), ]

identical(rownames(result),com2$WB)

com2 <- cbind(com2,result)
colnames(com2)
com2 <- com2[!(com2$WB == "HS034_WB"), ]
ggplot(com2,aes_string(x = "CCL21",y = 'Monocytes'))+
  geom_point(size = 5,color = '#2570A4',alpha = 3)+
  theme_bw()+
  theme(axis.title = element_text(size = 16),axis.text = element_text(size = 12),axis.ticks.length = unit(0.2,'cm'),axis.ticks = element_line(size = 1),panel.border = element_rect(size = 1.5),panel.grid = element_blank())+
  geom_smooth(method = 'lm',se = T,color = 'black',size = 2.0,fill = '#7A8991')+ stat_cor(method = "pearson",digits = 3,size=6)+
  ylab("Monocytes score")

com2 <- com2[,-c("P.value" ,"Correlation" ,"RMSE","group","WB"  )]
com2 <- select(com2, -one_of("P.value" ,"Correlation" ,"RMSE","group","WB" ))

target_gene <- "CCL19"
#设置你要进行相关性分析的范围
gene_list <-colnames(com2)

for (target_gene in gene_list) {
  print(paste("Calculating for target gene:", target_gene))
  
  # 循环读取
  for (gene in gene_list) {
    cor_R <- cor(x = com2[, target_gene], y = com2[, gene], method = 'spearman')
    cor_P <- round(cor.test(x = com2[, target_gene], y = com2[, gene])$p.value, 4)  # 四舍五入至4位小数
    temp_result <- data.frame(target_gene = target_gene,
                              gene_symbol = gene,
                              cor_R = cor_R,
                              cor_P = cor_P)
    result <- rbind(result, temp_result)
  }}

write.table(result,file = "相关性结果",sep = "\t")
