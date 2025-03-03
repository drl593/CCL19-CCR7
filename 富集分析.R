library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

load(file = "DEG GSE155176")
load(file = "DEG GSE154773")




GS <- AnnotationDbi::select(org.Hs.eg.db, keys = GSE154773, 
             keytype = "SYMBOL", column = "ENSEMBL")
GS <- na.omit(GS)
table(duplicated(GS))
GS1 <- select(org.Hs.eg.db, keys = GSE155176, 
              keytype = "ENSEMBL", column = "SYMBOL")
GS1 <- na.omit(GS1)
table(duplicated(GS1$SYMBOL))
COM <- Reduce(intersect, list(GS$ENSEMBL,GSE155176))
COM1 <- intersect(GS1$SYMBOL,GSE154773)
load(file = "DEG GSE155176")
load(file = "DEG GSE154773")




GS <- select(org.Hs.eg.db, keys = GSE154773, 
             keytype = "SYMBOL", column = "ENSEMBL")
GS <- na.omit(GS)
table(duplicated(GS))
GS1 <- select(org.Hs.eg.db, keys = GSE155176, 
              keytype = "ENSEMBL", column = "SYMBOL")
GS1 <- na.omit(GS1)
table(duplicated(GS1$SYMBOL))
COM <- Reduce(intersect, list(GS$ENSEMBL,GSE155176))
COM1 <- intersect(GS1$SYMBOL,GSE154773)

genelist <- bitr(COM1, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

#GO分析
ego <- enrichGO(gene = genelist$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
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
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.05)
write.table(ekegg,file="KEGG.txt",sep="\t",quote=F,row.names = F)
## 柱状图
barplot(ekegg, showCategory = 20,color = "pvalue")
## 气泡图
dotplot(ekegg, showCategory = 20)

###1: 提取logFC值，并储存在一个向量中
geneList = DEG[,2]
###2: 对geneList进行命名
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
###3: 根据logFC值降序排列
geneList = sort(geneList, decreasing = TRUE)



#GSEA转化  
genelist1 <- bitr(genelist$id, fromType="SYMBOL",
                  toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist1$change <- genelist[match(genelist1$SYMBOL,genelist$id),"log2FoldChange"]
genelist1 <- genelist1[order(genelist1$change,decreasing = T),]
genelist2 <- bitr(COM, fromType="ENSEMBL",
                  toType="ENTREZID", OrgDb='org.Hs.eg.db')
head(genelist1)


#2. GO的GSEA富集分析：gseGO
gsego <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               ont          = "CC",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
#将其中的基因名变成symbol ID
go_result <- setReadable(gsego, OrgDb = org.Hs.eg.db)
#转换成数据框
go_result_df <- as.data.frame(go_result)
write.table(go_result_df,file="GSEA_GO_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#可视化展示

#3. KEGG的GSEA富集分析：gseKEGG
gsekk <- gseKEGG(geneList     = geneList,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
#转换成数据框
gsekk_result_df <- as.data.frame(gsekk)
#可视化展示
gseaplot2(gsekk, geneSetID = c("hsa01100"))

#3. MSigDb的GSEA富集分析：GSEA
##MSigDb数据集下载网址
##https://www.gsea-msigdb.org/gsea/downloads.jsp
# 对应输入基因名，ENTREZID对应Entrez Gene ID
# msigdb.v7.0.symbols.gmt对应Gene Symbol名
##gmt文件储存的文件夹名
msigdb_GMTs <- "msigdb_v7.0_GMTs"
##指定用于GSEA富集分析的gmt文件
msigdb <- "c7.all.v7.0.entrez.gmt"

#提取文件
genelist$change<- DESeq2_DEG[match(genelist$SYMBOL,rownames(DESeq2_DEG)),"log2FoldChange"]
genelist <- genelist[order(genelist$change,decreasing = T),]
head(genelist)
G1 <- genelist1[,3]
names(G1) = as.character(genelist1[,'ENTREZID'])
head(G1)
##读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

GSEA3<-GSEA(G1,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
KEGG_result_df <- as.data.frame(GSEA2)

#单个图绘制
gseaplot2(GSEA3,1,color="red",pvalue_table = T)
gseaplot2(KEGG,3,color="red",pvalue_table = T)

#汇总结果
gseaplot2(GSEA2, geneSetID = 1:5)
gseaplot2(GSEA2, geneSetID = 1:3, subplots = 1)
gseaplot2(GSEA2, geneSetID = 1:3, subplots = 1:2)

write.table(COM1,file="200COM.txt",sep="\t",quote=F,row.names = F)

