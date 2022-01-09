#-----------------------数据整合
tumor_miExpr1<-read.csv(file= "../00rawdata/tumor_crc/CRC_125905_cellline_rawcount.csv")
tumor_miExpr2<-read.csv(file= "../00rawdata/tumor_crc/CRC_67004_cellline_rawcount.csv")
tumor_miExpr <- merge(tumor_miExpr1,tumor_miExpr2,by="gene_id",all=FALSE)

normal_miExpr1<- read.csv(file= "../00rawdata/normal/GBM_122488_normal_rawcount.csv")
normal_miExpr2<- read.csv(file= "../00rawdata/normal/CRC_71008_normal_rawcount.csv")
normal_miExpr3<- read.csv(file = '../00rawdata/normal/healthy_normal_130512_rawdata.csv')
normal_miExpr <- merge(normal_miExpr1,normal_miExpr2,by="gene_id",all=FALSE)
normal_miExpr <- merge(normal_miExpr,normal_miExpr3,by="gene_id",all=FALSE)

countData <- merge(normal_miExpr,tumor_miExpr,all=FALSE,by="gene_id")
rownames(countData)<- countData$gene_id
countData<- countData[,-1]
countData[is.na(countData)] <- 0

#----------------------差异分析
library("DESeq2")
condition<- factor(c(rep("normal",41),rep("cancer",21)))
colData <- data.frame(row.names =  colnames(countData), condition)
colData
dds <- DESeqDataSetFromMatrix(countData=countData , colData=colData,design = ~condition )
dds2 <- DESeq(dds)
resultsNames(dds2)
res1 <- results(dds2,contrast=c("condition", "normal", "cancer"))
table(res1$padj<0.01) #取P值小???0.05的结???
res1 <- res1[order(res1$pvalue),]
summary(res1)
diff_gene_deseq <-subset(res1, padj < 0.01 & abs(log2FoldChange) > 1  )#this
diff_gene_deseq2 <- data.frame(diff_gene_deseq)
diff_gene_deseq2$gene_id <- rownames(diff_gene_deseq2)
diff_gene2 <- rownames(diff_gene_deseq2)
diff_gene2 = data.frame(diff_gene2)
colnames(diff_gene2)<- "gene_id"

CRC_diff_1 <- diff_gene2


save(CRC_diff_1,file = "./test_data/CRC_0.01_200.Rdata")



d1<-CRC_diff_1
load('./diff_data/CRC_0.01_200.Rdata')
CRC_diff_1<-apply(CRC_diff_1,1,as.character)
intersect(CRC_diff_1,diff_gene2$gene_id)#有1个基因不同
#"hsa-miR-3975"
#"hsa-miR-764"
load('./CRC_var2_h20_51refgenes.Rdata')
