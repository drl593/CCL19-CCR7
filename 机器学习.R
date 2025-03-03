library(glmnet)
library(org.Hs.eg.db)
library(tidyverse)
#GSE155176----
load(file = "GSE155176.final")
#数据读取与分析
#将ensemble ID 转化为gene symbol

exp_tpms <- rownames_to_column(exp_tpms,var = "id")
GS <- AnnotationDbi::select(org.Hs.eg.db, keys = exp_tpms$id, 
                            keytype = "ENSEMBL", column = "SYMBOL")
colnames(GS)[1] <- "id"
exp_tpms$id <- GS[match(exp_tpms$id,GS$id),2]
exp_tpms <- exp_tpms[!duplicated(exp_tpms$id),]
exp_tpms <- na.omit(exp_tpms)
table(is.na(exp_tpms))
rownames(exp_tpms) <- NULL
exp_tpms <- column_to_rownames(exp_tpms,var = "id")
#选择筛选的两百个基因
COM <- read.table(file = "200COM.txt",header = T)

#转置
texp_tpms <- as.data.frame(t(exp_tpms[c(COM$x),]))
texp_tpms$group <- metadata[match(rownames(texp_tpms),rownames(metadata)),"Group"]
texp_tpms <- select(texp_tpms,"group", everything())
table(texp_tpms$group)

#响应变量的更改
texp_tpms$group <- ifelse(texp_tpms$group == "Lesional",1,0) #皮损为1，正常为0


# 定义结局变量和自变量(注意区别，此处转化为矩阵)
y <- as.matrix(texp_tpms[,1])  # 提取第一列作为二分类结局变量
x <- as.matrix(texp_tpms[,-1])  # 第2至第51列变量为要筛选的自变量

#


#随机森林----
